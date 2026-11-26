#!/usr/bin/env python3
"""V7D: identify the actual SOURCE of SRP074349 bridge runs vs selected HDF5.

READ ONLY. No expression matrices loaded, no accession/patient tokens printed.
This distinguishes an identifier-namespace problem from non-extracted runs.
"""
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py

ROOT = Path(sys.argv[1]).resolve()
REC = ROOT / "26_recount3_external_primary_v1"
S5 = REC / "stage5_output"
STUDY = "SRP074349"
RUN_RE = re.compile(r"(?:SRR|ERR|DRR)\d+", re.IGNORECASE)

def tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        cols = rdr.fieldnames or []
        return cols, list(rdr)

def dec(x):
    return bytes(x).decode("utf-8").strip() if isinstance(x, bytes) else str(x).strip()

def report(label, c):
    if not c:
        print(label, "EMPTY")
    else:
        print(label, " ".join(f"{repr(k)}={v}" for k,v in sorted(c.items(), key=lambda x: str(x[0]))))

def prefix(s):
    match = re.match(r"^([A-Za-z]+)", str(s).strip())
    return match.group(1).upper() if match else "OTHER"

print("V7D_SRP074349_BRIDGE_VS_EXPRESSION_SOURCE_V1")
bridgepath = S5 / "03_bridge_165_availability.tsv"
_, bridge = tsv(bridgepath)
pairs = [r for r in bridge if r["study"] == STUDY]
assert len(bridge) == 165 and len(pairs) == 143
b_libs = [r["library_id"].strip() for r in pairs]
b_samples = [r["biosample"].strip() for r in pairs]
b_runs = [RUN_RE.findall(r["recount3_run_ids"] or "") for r in pairs]
assert len(set(b_libs)) == 143
assert all(len(x)==int(r["n_runs"]) for x,r in zip(b_runs,pairs))
assert all(len(x)==len(set(x)) for x in b_runs)
all_runs = [run.upper() for x in b_runs for run in x]
b_runset = set(all_runs)
print("BRIDGE_EXPERIMENTS",len(pairs))
print("BRIDGE_RUN_REFERENCES",len(all_runs))
print("BRIDGE_DISTINCT_RUNS",len(b_runset))
report("BRIDGE_RUN_MULTIPLICITY",Counter(len(x) for x in b_runs))
report("BRIDGE_CLASS",Counter(r["bridge_class"] for r in pairs))
report("BRIDGE_LIBRARY_ID_PREFIX",Counter(prefix(x) for x in b_libs))
report("BRIDGE_RUN_ID_PREFIX",Counter(prefix(x) for x in b_runset))
print("BRIDGE_BIOSAMPLES_POPULATED",sum(bool(x) for x in b_samples))
print("BRIDGE_DISTINCT_NONEMPTY_BIOSAMPLES",len({x for x in b_samples if x}))

h5path = REC / "stage2_output/studies/SRP074349/selected_run_coverage_G026.h5"
with h5py.File(h5path, "r") as h:
    h_runs = [dec(x).upper() for x in h["run_id"][:]]
    h_libs = [dec(x) for x in h["library_id"][:]]
    assert len(h_runs)==len(h_libs)==144
hruns = set(h_runs)
hlibs = set(h_libs)
h_run_lib = dict(zip(h_runs,h_libs))
print("H5_RUN_COLUMNS",len(h_runs))
print("H5_DISTINCT_RUNS",len(hruns))
print("H5_DISTINCT_LIBRARIES",len(hlibs))
report("H5_RUN_PREFIX",Counter(prefix(x) for x in hruns))
report("H5_LIBRARY_PREFIX",Counter(prefix(x) for x in hlibs))

print("EXACT_BRIDGE_H5_RUN_OVERLAP",len(b_runset & hruns))
print("EXACT_BRIDGE_H5_LIBRARY_OVERLAP",len(set(b_libs) & hlibs))
print("EXACT_BRIDGE_RUN_TOKENS_MISSING_FROM_H5",len(b_runset - hruns))
print("BRIDGE_EXPERIMENTS_WITH_ANY_RUN_IN_H5",sum(bool(set(rs)&hruns) for rs in b_runs))
print("BRIDGE_EXPERIMENTS_WITH_ALL_RUNS_IN_H5",sum(bool(rs) and set(rs)<=hruns for rs in b_runs))
print("H5_RUNS_NOT_IN_BRIDGE",len(hruns-b_runset))
print("BRIDGE_LIBRARIES_MATCH_H5_RUN_TOKENS",len(set(b_libs)&hruns))
print("BRIDGE_RUN_TOKENS_MATCH_H5_LIBRARY_IDS",len(b_runset&hlibs))

# Historical selection and quality records. Exact joins only.
records = (
    ("STAGE2_LIBRARY_QC", REC / "stage2_output/library_qc_all.tsv"),
    ("STAGE2_RUN_QC", REC / "stage2_output/run_qc_all.tsv"),
    ("STAGE2_STUDY_RUN_QC", REC / "stage2_output/studies/SRP074349/run_qc.tsv"),
    ("STAGE5_PREANALYSIS", S5 / "04_library_preanalysis_decision_3639.tsv"),
    ("STAGE41_CANDIDATES", ROOT / "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/03_candidate_integrated_evidence_3602.tsv"),
)
for label,path in records:
    print("TABLE",label,"PRESENT",path.is_file())
    if not path.is_file():
        continue
    cols,rows=tsv(path)
    current=[r for r in rows if r.get("study")==STUDY]
    print("TABLE",label,"STUDY_ROWS",len(current))
    print("TABLE",label,"HAS_LIBRARY_ID","library_id" in cols,"HAS_RUN_ID","run_id" in cols,"HAS_BIOSAMPLE","biosample" in cols)
    for key,source in (("library_id",set(b_libs)),("run_id",b_runset),("biosample",{x for x in b_samples if x})):
        if key not in cols:
            continue
        vals={r.get(key,"").strip() for r in current}
        if key=="run_id":
            vals={v.upper() for v in vals}
        print("TABLE",label,"BRIDGE_"+key+"_OVERLAP",len(vals&source))
        matched=[r for r in current if (r.get(key,"").strip().upper() if key=="run_id" else r.get(key,"").strip()) in source]
        print("TABLE",label,"ROWS_MATCHING_"+key,len(matched))
    for k in ("preanalysis_disposition","stage41_partition_NOT_ADMISSION","technical_hold","run_count_status"):
        if k in cols:
            report("TABLE_"+label+"_"+k,Counter(r.get(k,"") for r in current))

# Native study metadata: locate run accessions without emitting identifier values.
native = REC / "raw_studies/SRP074349/metadata/sra.sra.SRP074349.MD.gz"
print("NATIVE_SRA_METADATA_PRESENT",native.is_file())
if native.is_file():
    with gzip.open(native,"rt",encoding="utf-8-sig",errors="replace",newline="") as f:
        rdr=csv.DictReader(f,delimiter="\t")
        cols=rdr.fieldnames or []
        rows=list(rdr)
    print("NATIVE_SRA_ROWS",len(rows))
    print("NATIVE_SRA_COLUMNS",";".join(cols[:60]))
    print("NATIVE_EXPERIMENT_ID_BRIDGE_MATCH",sum(r.get("experiment_acc","").strip() in set(b_libs) for r in rows))
    for col in cols:
        if not re.search("run|experiment|sample|biosample",col,re.I):
            continue
        matches_run=set()
        matches_lib=set()
        matched_rows_run=0
        matched_rows_lib=0
        for row in rows:
            val=(row.get(col) or "")
            toks={t.upper() for t in RUN_RE.findall(val)}
            exactparts={x.strip() for x in re.split(r"[;,|\s]+",val) if x.strip()}
            hits=toks&b_runset
            libs=exactparts&set(b_libs)
            if hits:
                matches_run.update(hits)
                matched_rows_run+=1
            if libs:
                matches_lib.update(libs)
                matched_rows_lib+=1
        if matches_run or matches_lib:
            print("NATIVE_ID_COLUMN",col,"MATCHED_BRIDGE_RUN_TOKENS",len(matches_run),
                  "ROWS_WITH_RUN_MATCH",matched_rows_run,
                  "MATCHED_BRIDGE_LIBRARIES",len(matches_lib),
                  "ROWS_WITH_LIBRARY_MATCH",matched_rows_lib)

# Inspect other *existing* data assets in this one study directory by filename only.
study_dir = REC/"stage2_output/studies/SRP074349"
if study_dir.is_dir():
    files=sorted(p.name for p in study_dir.iterdir() if p.is_file())
    print("STAGE2_STUDY_FILE_COUNT",len(files))
    print("STAGE2_STUDY_FILES",";".join(files[:40]))

print("NO_HDF5_EXPRESSION_VALUES_LOADED")
print("NO_INDIVIDUAL_SAMPLE_OR_RUN_IDENTIFIERS_PRINTED")
print("NO_SOURCE_OR_FAILED_RUN_MODIFICATIONS")
print("NO_MODEL_PROJECTION")
print("CLINICAL_ADMISSIONS 0")
print("V7D_SRP074349_BRIDGE_VS_EXPRESSION_SOURCE_V1_COMPLETE")
