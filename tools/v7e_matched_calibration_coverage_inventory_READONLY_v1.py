#!/usr/bin/env python3
"""V7E: fast READ-ONLY availability census for 454 matched calibration runs.

Scans run_id metadata axes (NOT coverage arrays) in existing Stage2 HDF5 files.
Records only aggregate counts, never individual run, sample or patient IDs.
Identifies existing extraction candidates without downloading or altering data.
"""
import csv
import gzip
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT=Path(sys.argv[1]).resolve()
REC=ROOT/"26_recount3_external_primary_v1"
STAGE2=REC/"stage2_output/studies"
S5=REC/"stage5_output"
RUN_RE=re.compile(r"(?:SRR|ERR|DRR)\d+",re.I)

def decode(x):
    return bytes(x).decode("utf-8").strip() if isinstance(x,(bytes,np.bytes_)) else str(x).strip()

def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

print("V7E_MATCHED_CALIBRATION_SOURCE_AVAILABILITY_V1")
bridge=rows(S5/"03_bridge_165_availability.tsv")
m=[r for r in bridge if r["study"]=="SRP074349"]
assert len(bridge)==165 and len(m)==143
wanted=set()
by_lib={}
for r in m:
    ids=[x.upper() for x in RUN_RE.findall(r["recount3_run_ids"] or "")]
    assert len(ids)==int(r["n_runs"]) and len(ids)==len(set(ids))
    k=r["library_id"].strip()
    assert k not in by_lib
    by_lib[k]=set(ids)
    wanted.update(ids)
assert len(wanted)==454
print("TARGET_STUDY SRP074349")
print("MATCHED_BRIDGE_LIBRARIES",len(by_lib))
print("MATCHED_DISTINCT_RUNS",len(wanted))

# Check the existing run-level frozen extraction inventory.
runqc=REC/"stage2_output/run_qc_all.tsv"
print("STAGE2_RUN_QC_PRESENT",runqc.is_file())
if runqc.is_file():
    a=rows(runqc)
    existing={r.get("run_id","").upper() for r in a}
    print("STAGE2_QC_TOTAL_ROWS",len(a))
    print("TARGET_RUNS_IN_STAGE2_QC",len(wanted&existing))
    print("TARGET_RUNS_MISSING_FROM_STAGE2_QC",len(wanted-existing))

h5files=sorted(STAGE2.glob("*/selected_run_coverage_G026.h5")) if STAGE2.is_dir() else []
print("STAGE2_SELECTED_H5_FILES",len(h5files))
matched=set()
h5_reads=0
bad=[]
file_hits=[]
for p in h5files:
    try:
        with h5py.File(p,"r") as h:
            if "run_id" not in h or "coverage_bp" not in h:
                bad.append("MISSING_EXPECTED_AXES")
                continue
            ids={decode(x).upper() for x in h["run_id"][:]}
            shape=h["coverage_bp"].shape
        h5_reads+=1
        hits=ids&wanted
        if hits:
            file_hits.append((p.parent.name,len(hits),shape))
            matched.update(hits)
    except Exception as e:
        bad.append(type(e).__name__)
print("STAGE2_H5_AXES_OPENED",h5_reads)
print("STAGE2_H5_READ_ERRORS",len(bad))
if bad:
    print("STAGE2_H5_ERROR_CLASSES",dict(Counter(bad)))
print("TARGET_RUNS_IN_ANY_STAGE2_H5",len(matched))
print("TARGET_RUNS_NOT_IN_STAGE2_H5",len(wanted-matched))
print("STAGE2_H5_FILES_WITH_TARGET_RUNS",len(file_hits))
for study,n,shape in file_hits:
    print("H5_MATCHING_STUDY",study,"RUNS",n,"SHAPE",str(shape))
print("BRIDGE_LIBRARIES_ALL_RUNS_IN_H5",sum(v<=matched for v in by_lib.values()))
print("BRIDGE_LIBRARIES_ANY_RUNS_IN_H5",sum(bool(v&matched) for v in by_lib.values()))

# Native sample membership (descriptive only): experiment IDs, no values.
native=REC/"raw_studies/SRP074349/metadata/sra.sra.SRP074349.MD.gz"
print("NATIVE_SRA_PRESENT",native.is_file())
if native.is_file():
    with gzip.open(native,"rt",encoding="utf-8-sig",errors="replace",newline="") as f:
        rdr=csv.DictReader(f,delimiter="\t")
        exp=[r.get("experiment_acc","").strip() for r in rdr]
    print("NATIVE_SRA_METADATA_ROWS",len(exp))
    print("BRIDGE_EXPERIMENTS_IN_NATIVE_SRA",len(set(by_lib)&set(exp)))

# Candidate provenance names only. Bound recursion to TWO folder levels.
# The user can then direct a targeted exact-run investigation if a cache exists.
roots=[
 ("STAGE2_STUDY",REC/"stage2_output/studies/SRP074349"),
 ("STAGE3",REC/"stage3_output"),
 ("STAGE5",REC/"stage5_output"),
 ("STUDY_RAW",REC/"raw_studies/SRP074349"),
]
tokens=("coverage","count","expression","gene","h5","parquet","bridge","calibr","quant")
for label,base in roots:
    print("CANDIDATE_ROOT",label,"PRESENT",base.is_dir())
    if not base.is_dir():continue
    observed=[]
    for current,folders,files in os.walk(base):
        depth=len(Path(current).relative_to(base).parts)
        if depth>=2:folders[:]=[]
        else: folders[:]=sorted(f for f in folders if not f.startswith("."))
        for filename in files:
            if any(token in filename.lower() for token in tokens):
                path=Path(current)/filename
                observed.append((str(path.relative_to(base)),int(path.stat().st_size)))
        if len(observed)>1200:
            break
    print("CANDIDATE_FILES",label,len(observed))
    for name,sz in sorted(observed)[:30]:
        print("CANDIDATE_FILE",label,name,"BYTES",sz)
    if len(observed)>30:print("CANDIDATE_FILE_LIST_TRUNCATED",label)

if len(wanted-matched)==454:
    print("NEXT_DECISION_NO_MATCHED_RUNS_IN_CURRENT_STAGE2_H5")
elif wanted<=matched:
    print("NEXT_DECISION_MATCHED_RUNS_FOUND_IN_STAGE2_H5")
else:
    print("NEXT_DECISION_PARTIAL_MATCHED_RUNS_FOUND_IN_STAGE2_H5")
print("NO_SOURCE_FILES_MODIFIED")
print("NO_EXPRESSION_VALUES_LOADED")
print("NO_DOWNLOADS_OR_MODEL_EXECUTION")
print("NO_CLINICAL_ADMISSIONS")
print("V7E_MATCHED_CALIBRATION_SOURCE_AVAILABILITY_V1_COMPLETE")
