#!/usr/bin/env python3
"""FUSCC remaining release gates: original run-QC membership, protocol and local evidence inventory."""
import csv,gzip,hashlib,json,sys
from pathlib import Path
from collections import Counter,defaultdict
root=Path(sys.argv[1]);repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1"
s41=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
study="SRP157974"
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def count(rs,f):
    assert all(f in r for r in rs),f
    return dict(sorted(Counter(r[f] for r in rs).items()))
def metric(rs,f):
    vals=[]
    for r in rs:
        try:vals.append(float(r[f]))
        except (KeyError,ValueError,TypeError):pass
    vals.sort()
    import math
    vals=[v for v in vals if math.isfinite(v)]
    return {"n_finite":len(vals),"n_other":len(rs)-len(vals),
            "min":vals[0] if vals else None,
            "median":(vals[(len(vals)-1)//2]+vals[len(vals)//2])/2 if vals else None,
            "max":vals[-1] if vals else None}
ledger=repo/"private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv"
assert sha(ledger)=="943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0"
lr=[r for r in read(ledger) if r["study"]==study]
keep={r["library_id_LOCAL_ONLY"] for r in lr}
assert len(lr)==len(keep)==334 and all(r["proposed_non_tumour_exclusion"]=="NO" and r["final_approval"]=="NO_UNSIGNED" for r in lr)
cand=s41/"03_candidate_integrated_evidence_3602.tsv"
full=s41/"01_full_source_experiment_identity_5221.tsv"
assert sha(cand)=="aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81"
assert sha(full)=="95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa"
cr=[r for r in read(cand) if r["study"]==study]
fr=[r for r in read(full) if r["study"]==study]
hold={r["library_id"] for r in cr if r["technical_hold"]=="YES"}
assert len(cr)==360 and len(fr)==448 and len(hold)==26
assert keep=={r["library_id"] for r in cr if r["technical_hold"]=="NO"}
assert not keep&hold
print("FUSCC_LOCKED_MEMBERSHIP",len(fr),len(cr),len(keep),len(hold))
native=rec/"raw_studies"/study/"metadata"/f"sra.sra.{study}.MD.gz"
assert sha(native)=="c2e1c7b5b1169e13059447e14bdc71ed3ffd4166dd90131b2791360ebab28155"
nr=read(native)
assert len(nr)==len({r["experiment_acc"] for r in nr})==448
assert {r["experiment_acc"] for r in nr}=={r["experiment_acc"] for r in fr}
print("NATIVE_SOURCE_HASH_PASS")
print("=== CLINICAL AND ASSAY EVIDENCE ===")
print("HISTORICAL_PATIENT_GRADE",count(lr,"historical_patient_grade"))
print("HISTORICAL_TREATMENT_GRADE",count(lr,"historical_treatment_grade"))
for label,ids in (("RETAINED334",keep),("HELD26",hold)):
    rs=[r for r in nr if r["experiment_acc"] in ids]
    print("SCOPE",label,"N",len(rs))
    for f in ("library_strategy","library_selection","library_layout","platform_model"):
        print(f,count(rs,f))
    print("TUMORTISSUE_IN_NATIVE_NAME",sum("tumortissue" in r["sample_name"].casefold() for r in rs))
    import re
    protocols=Counter(r["library_construction_protocol"].strip() for r in rs)
    print("PROTOCOL_NONEMPTY",sum(n for t,n in protocols.items() if t),"DISTINCT_NONEMPTY",sum(bool(t) for t in protocols))
    for t,n in protocols.items():
        if not t:continue
        print("PROTOCOL_GROUP",n,"SHA256",hashlib.sha256(t.encode()).hexdigest(),
              "RIBO_ZERO_OR_DEPLETION",bool(re.search(r"ribo.?zero|(?:ribosomal|rrna).{0,20}deplet|deplet.{0,30}(?:ribosomal|rrna)",t,re.I)),
              "POLYA_OR_OLIGODT",bool(re.search(r"poly\s*[\(-]?\s*a\b|polyadenylat|oligo.?dt",t,re.I)))
# Original authoritative QC inputs; do not rerun QC or release existing holds.
files={"library_qc":rec/"stage2_output/library_qc_all.tsv",
       "run_qc":rec/"stage2_output/run_qc_all.tsv",
       "hold_queue":rec/"cohort_readiness_v1/source_snapshot/05_run_hold_queue_409.tsv"}
for name,path in files.items():
    assert path.is_file()
    print("QC_SOURCE",name,"SHA256",sha(path))
lq=[r for r in read(files["library_qc"]) if r["study"]==study]
rq=[r for r in read(files["run_qc"]) if r["study"]==study]
hq=[r for r in read(files["hold_queue"]) if r["study"]==study]
print("QC_LIBRARY_ROWS",len(lq),"QC_RUN_ROWS",len(rq),"HOLD_QUEUE_ROWS",len(hq))
for label,ids in (("RETAINED334",keep),("HELD26",hold)):
    print("=== QC",label,"===")
    for name,data in (("library_qc",lq),("run_qc",rq),("hold_queue",hq)):
        rs=[r for r in data if r["library_id"] in ids]
        found={r["library_id"] for r in rs}
        print(name,"ROWS",len(rs),"COVERED_LIBRARIES",len(found),"ABSENT_LIBRARIES",len(ids-found))
        if name=="library_qc":
            for f in ("run_count_status","not_collapsed"):print(f,count(rs,f))
            for f in ("n_runs_expected","n_runs_expression_verified","min_run_n_genes_gt0","max_run_n_genes_gt0"):print(f,metric(rs,f))
        elif name=="run_qc":
            for f in ("basic_expression_flag","pred_metadata_status"):print(f,count(rs,f))
            for f in ("fraction_genes_detected","gene_body_coverage_bp"):print(f,metric(rs,f))
        else:
            for f in ("decision","reason_codes","automatic_proposal","triage_priority","proposed_action_not_applied"):print(f,count(rs,f))
            for f in ("gene_detection_fraction","top10_coverage_fraction","average_mapped_length_bp"):print(f,metric(rs,f))
print("HOLD_QUEUE_LIBRARIES_OUTSIDE_360",len({r["library_id"] for r in hq}-(keep|hold)))
print("=== EXISTING LOCAL FUSCC EVIDENCE FILE INVENTORY ===")
# Search only the local repository evidence/private roots; no content or identifiers.
# Exclude installed environment trees and raw sequence/expression files.
import os
for parent in (repo/"private_data",repo/"diagnostics",repo/"evidence"):
    for dp,dirs,fs in os.walk(parent):
        dirs[:]=[d for d in dirs if d not in {"environments",".git","__pycache__"}]
        for name in sorted(fs):
            rel=(Path(dp)/name).relative_to(repo)
            key=str(rel).casefold()
            if ("fuscc" in key or "srp157974" in key or "19783498" in key) and (name.lower().endswith((".tsv",".csv",".xlsx",".xls",".json",".txt",".zip",".gz"))):
                print("LOCAL_EVIDENCE_FILE",str(rel),"BYTES",(repo/rel).stat().st_size)
print("V8_FUSCC_RELEASE_GATE_REVIEW_COMPLETE")
print("NO_FILES_MODIFIED_NO_MODEL_EXECUTION_CLINICAL_ADMISSIONS_0")
