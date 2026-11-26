#!/usr/bin/env python3
"""V7F: read-only inspect raw SRP074349 recount3 G026 gzip header vs 454 bridge runs.

Only reads the first 20 decompressed lines; never prints run IDs, gene IDs, 
individual data values, or patient-linked metadata. No data extraction.
"""
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT=Path(sys.argv[1]).resolve()
REC=ROOT/"26_recount3_external_primary_v1"
BRIDGE=REC/"stage5_output/03_bridge_165_availability.tsv"
RAW=REC/"raw_studies/SRP074349/gene_sums/sra.gene_sums.SRP074349.G026.gz"
REPORT=REC/"stage5_output/bridge_SRP074349/bridge_report.json"
RUN=re.compile(r"^(?:SRR|ERR|DRR)\d+$", re.I)
RUN_IN_TEXT=re.compile(r"(?:SRR|ERR|DRR)\d+", re.I)
ENSG=re.compile(r"^ENSG\d+",re.I)

def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def split_line(s,d):
    # Filename is G026.gz; delimiter is NOT assumed.
    return s.rstrip("\r\n").split(d)

print("V7F_RAW_SRP074349_G026_HEADER_MEMBERSHIP_V1")
with BRIDGE.open(encoding="utf-8-sig",newline="") as f:
    bridge=list(csv.DictReader(f,delimiter="\t"))
pairs=[r for r in bridge if r["study"]=="SRP074349"]
assert len(bridge)==165 and len(pairs)==143
wanted=set()
for r in pairs:
    ids=RUN_IN_TEXT.findall(r["recount3_run_ids"] or "")
    assert len(ids)==int(r["n_runs"]) and len(set(ids))==len(ids)
    wanted.update(x.upper() for x in ids)
assert len(wanted)==454
print("BRIDGE_TARGET_EXPERIMENTS",len(pairs))
print("BRIDGE_TARGET_UNIQUE_RUNS",len(wanted))
print("RAW_GZIP_EXISTS",RAW.is_file())
if not RAW.is_file():
    print("V7F_RAW_FILE_MISSING_NO_EXTRACTION")
    sys.exit(1)
print("RAW_COMPRESSED_BYTES",RAW.stat().st_size)
print("RAW_SHA256",sha256(RAW))
try:
    with gzip.open(RAW,"rt",encoding="utf-8-sig",errors="replace",newline="") as f:
        lines=[]
        for _ in range(20):
            line=f.readline()
            if line=="":break
            lines.append(line)
except Exception as exc:
    print("RAW_GZIP_READ_ERROR",type(exc).__name__)
    raise
print("RAW_FIRST_LINES_AVAILABLE",len(lines))
assert lines
for delim,label in (("\t","TAB"),(",","COMMA"),(" ","SPACE")):
    shapes=Counter()
    hits_per_line=[]
    max_run_per_line=0
    best_line=None
    best_runs=set()
    for i,line in enumerate(lines):
        cells=split_line(line,delim)
        shapes[len(cells)]+=1
        # exact accession values only, not substrings
        run_values=[cell.strip().upper().strip('"') for cell in cells if RUN.fullmatch(cell.strip().strip('"'))]
        runset=set(run_values)
        hit=len(runset & wanted)
        hits_per_line.append(hit)
        if hit>len(best_runs&wanted):
            best_line=i
            best_runs=runset
        max_run_per_line=max(max_run_per_line,len(runset))
    print("DELIMITER",label,
          "FIRST_LINE_FIELDS",len(split_line(lines[0],delim)),
          "FIELD_COUNTS",dict(sorted(shapes.items())),
          "MAX_RUN_FIELDS_ANY_LINE",max_run_per_line,
          "MAX_BRIDGE_RUN_MATCHES_ANY_LINE",max(hits_per_line))
    if label=="TAB":
        header=split_line(lines[0],delim)
        names=[x.strip().strip('"').upper() for x in header]
        hrun={x for x in names if RUN.fullmatch(x)}
        print("TAB_FIRST_LINE_UNIQUE_RUN_ACCESSIONS",len(hrun))
        print("TAB_FIRST_LINE_BRIDGE_RUN_OVERLAP",len(hrun&wanted))
        print("TAB_FIRST_LINE_BRIDGE_RUN_MISSING",len(wanted-hrun))
        print("TAB_FIRST_LINE_NON_BRIDGE_RUNS",len(hrun-wanted))
        print("TAB_FIRST_LINE_DUPLICATED_RUN_LABELS",sum(v-1 for k,v in Counter(names).items() if v>1 and RUN.fullmatch(k)))
        print("TAB_FIRST_COLUMN_GENE_PATTERN_ROWS_2_TO_20",sum(bool(ENSG.match(split_line(x,"\t")[0].strip().strip('"'))) for x in lines[1:]))
        print("TAB_DATA_ROW_WIDTH_MATCH_FIRST_LINE",sum(len(split_line(x,"\t"))==len(header) for x in lines[1:]))
        if len(hrun & wanted)==454:
            print("EXACT_ALL_454_RUNS_PRESENT_IN_TAB_HEADER_TRUE")
        else:
            print("EXACT_ALL_454_RUNS_PRESENT_IN_TAB_HEADER_FALSE")
    if best_line is not None:
        print("DELIMITER",label,"BEST_LINE_0_BASED",best_line,
              "RUN_VALUES",len(best_runs),
              "BRIDGE_OVERLAP",len(best_runs&wanted))

print("PRIOR_BRIDGE_REPORT_EXISTS",REPORT.is_file())
if REPORT.is_file():
    report=json.loads(REPORT.read_text(encoding="utf-8"))
    print("PRIOR_BRIDGE_REPORT_TOP_LEVEL_KEYS",";".join(sorted(map(str,report.keys()))[:35]))
    # Report only scalar numerical facts; never print sample IDs or paths.
    for k,v in sorted(report.items()):
        if isinstance(v,(int,float,bool)) and not isinstance(v,str):
            print("PRIOR_BRIDGE_NUMERIC",k,v)
print("NO_EXPRESSION_MATRIX_EXTRACTION")
print("NO_SAMPLE_OR_RUN_ACCESSIONS_PRINTED")
print("NO_FROZEN_OR_SOURCE_FILES_MODIFIED")
print("NO_MODEL_EXECUTION")
print("CLINICAL_ADMISSIONS 0")
print("V7F_RAW_SRP074349_G026_HEADER_MEMBERSHIP_V1_COMPLETE")
