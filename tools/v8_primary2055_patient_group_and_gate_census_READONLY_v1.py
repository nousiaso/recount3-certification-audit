#!/usr/bin/env python3
"""Read-only cohort-wide historical evidence and candidate grouping census.
No patient or specimen identifiers printed; no approval, collapse or filtering.
"""
import csv, hashlib, io, re, sys, zipfile
from collections import Counter, defaultdict
from pathlib import Path

root=Path(sys.argv[1])
repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1"
s41=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def val(v):return str(v or "").strip()
def present(v):return val(v).casefold() not in {"","na","n/a","none","null","nan"}
def codes(v):
    s=val(v)
    # Only categorical uppercase evidence grades are displayed.
    return s if re.fullmatch(r"[A-Z][A-Z0-9_]{0,100}",s) else "EMPTY" if not present(s) else "OTHER_NONCODE_VALUE"
def groups_summary(groups):
    return {"groups":len(groups),"reused_groups":sum(len(v)>1 for v in groups.values()),
            "records_in_reused_groups":sum(len(v) for v in groups.values() if len(v)>1),
            "max_group_size":max(map(len,groups.values()),default=0)}

v5=rec/"RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"
assert sha(v5)=="a76a68465629dafba98f16f688a69e25fa3cfcb0282ab30add66206a85c62b36"
with zipfile.ZipFile(v5) as z:
    paths=[p for p in z.namelist() if p.endswith("/results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv")]
    assert len(paths)==1
    raw=z.read(paths[0])
assert hashlib.sha256(raw).hexdigest()=="06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49"
rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t"))
primary={r["library_id"]:r for r in rows}
assert len(rows)==len(primary)==2055 and len({r["study"] for r in rows})==77
assert all(r["FINAL_CLINICAL_APPROVAL"]=="NO_UNSIGNED" for r in rows)
candpath=s41/"03_candidate_integrated_evidence_3602.tsv"
assert sha(candpath)=="aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81"
crows=read(candpath)
c={(r["study"],r["library_id"]):r for r in crows}
assert len(c)==len(crows)==3602
assert all((r["study"],r["library_id"]) in c and c[(r["study"],r["library_id"])]["technical_hold"]=="NO" for r in rows)

# Reuse existing V6 precedence decisions; do not rerun their construction.
evpath=repo/"diagnostics/20261123_v6_evidence_precedence_dryrun/proposed_evidence_2055_UNSIGNED.tsv"
assert sha(evpath)=="08efe2c8f469edcbdfb98231dc87f705d801434e78dcbea85f57c4b94da88eb0"
evrows=read(evpath);ev={r["library_id"]:r for r in evrows}
assert len(ev)==len(evrows)==2055 and ev.keys()==primary.keys()
assert all(ev[k]["study"]==r["study"] for k,r in primary.items())
assert all(r["FINAL_APPROVAL"]=="NO_UNSIGNED" for r in evrows)
exclude={k for k,r in ev.items() if r["unsigned_evidence_bucket"]=="EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE"}
assert len(exclude)==22
print("V8_PRIMARY2055_PATIENT_AND_GATE_CENSUS")
print("LOCKED_INPUTS_PASS",2055,"STUDIES",77)
print("HISTORICAL_PROPOSED_NON_TUMOUR_EXCLUSIONS",len(exclude),"ADMIN_REMAINDER_NOT_APPROVED",2055-len(exclude))
print("HISTORICAL_GRADES_NOT_CURRENT_ADJUDICATION")
by=defaultdict(list)
for r in rows:by[r["study"]].append(r)

BIO=re.compile(r"(?<![A-Za-z0-9])(?:SAMN|SAMEA|SAMD)\d+(?![A-Za-z0-9])",re.I)
bios=defaultdict(set)
for r in rows:
    for token in BIO.findall(r["biosample"]):
        bios[token.upper()].add((r["study"],r["library_id"]))
assert all(BIO.search(r["biosample"]) for r in rows)
print("EXTERNAL_BIOSAMPLE_GROUPS",groups_summary(bios))
print("EXTERNAL_BIOSAMPLES_SHARED_ACROSS_STUDIES",sum(len({s for s,_ in g})>1 for g in bios.values()))
print("EXTERNAL_RECORDS_IN_CROSS_STUDY_BIOSAMPLE_GROUPS",len(set().union(*(g for g in bios.values() if len({s for s,_ in g})>1))))
print("NOTE_BIOSAMPLE_IDENTITY_IS_NOT_PATIENT_IDENTITY")

grades=("patient_evidence_grade","treatment_evidence_grade","tissue_evidence_grade")
gates=("patient_independence_certified","primary_tumour_certified","pretreatment_certified","assay_chemistry_certified","atlas_independence_certified")
for f in grades:
    print("V5_HISTORICAL_GRADE",f,dict(Counter(codes(r[f]) for r in rows)))
for f in gates:
    print("STAGE41_HISTORICAL_GATE",f,dict(Counter(codes(c[(r["study"],r["library_id"])][f]) for r in rows)))

print("=== STUDY_SCOPED_PATIENT_KEY_GROUPS_AND_EVIDENCE ===")
header=["study","admin_n","proposed_non_tumour_n","remaining_unsigned_n","key_missing_n","distinct_provisional_keys","reused_key_groups","records_in_reused_groups","max_key_group","biosample_proxy_grade_n","curated_key_unverified_grade_n","pretreatment_not_established_grade_n"]
print("\t".join(header))
allkeys=defaultdict(set);totalmissing=0
for s,rs in sorted(by.items(),key=lambda q:(-len(q[1]),q[0])):
    keys=defaultdict(set)
    missing=0
    for r in rs:
        k=val(r["patient_key"])
        if present(k):
            keys[k].add(r["library_id"])
            allkeys[(s,k)].add(r["library_id"])
        else:missing+=1
    totalmissing+=missing
    summary=groups_summary(keys)
    n_ex=sum(r["library_id"] in exclude for r in rs)
    line=[s,len(rs),n_ex,len(rs)-n_ex,missing,summary["groups"],summary["reused_groups"],summary["records_in_reused_groups"],summary["max_group_size"],
          sum(r["patient_evidence_grade"]=="BIOSAMPLE_PROXY_NOT_PATIENT_ID" for r in rs),
          sum(r["patient_evidence_grade"]=="CURATED_STUDY_PATIENT_KEY_UNVERIFIED" for r in rs),
          sum(r["treatment_evidence_grade"]=="PRETREATMENT_NOT_ESTABLISHED" for r in rs)]
    print("\t".join(map(str,line)))
print("ALL_STUDY_SCOPED_PROVISIONAL_KEYS",groups_summary(allkeys),"MISSING_RECORDS",totalmissing)
print("V8_PRIMARY2055_PATIENT_AND_GATE_CENSUS_COMPLETE")
print("NO_FILES_MODIFIED_NO_PATIENT_CERTIFICATION_NO_CLINICAL_ADMISSIONS")
