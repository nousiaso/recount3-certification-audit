#!/usr/bin/env python3
"""Prioritize substantive clinical review using existing joint evidence, not approvals."""
import csv,hashlib,sys
from collections import Counter,defaultdict
from pathlib import Path
root=Path(sys.argv[1]);repo=root/"RECOUNT3_CERTIFICATION_GIT"
def read_locked(rel,digest):
    p=repo/rel
    assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,rel
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
rows=read_locked("private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv","943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0")
group=read_locked("private_data/v8_repeated_group_evidence_addendum_v1/124_group_evidence_addendum_LOCAL_ONLY.tsv","12fc1fd8bee7f2d400c4ab9a3841f3460e494511035738d05e55f10a95579a1c")
fuscc=read_locked("private_data/v8_fuscc_author_annotation_addendum_v1/360_fuscc_author_annotation_evidence_LOCAL_ONLY.tsv","1cae682fe3ee006df02e22654f558c0ee2e28ab66fd5a422c67374609e1d1dd5")
assert len(rows)==len({r["library_id_LOCAL_ONLY"] for r in rows})==2055
assert all(r["final_approval"]=="NO_UNSIGNED" for r in rows)
rs=[r for r in rows if r["proposed_non_tumour_exclusion"]=="NO"]
assert len(rs)==2033
G={r["library_id_LOCAL_ONLY"]:r for r in group};F={r["experiment_acc_LOCAL_ONLY"]:r for r in fuscc if r["membership"]=="RETAINED334"}
assert len(G)==124 and len(F)==334
assert G.keys()<= {r["library_id_LOCAL_ONLY"] for r in rs} and F.keys()=={r["library_id_LOCAL_ONLY"] for r in rs if r["study"]=="SRP157974"}
assert all(r["final_clinical_approval"]=="NO_UNSIGNED" for r in group+fuscc)
by=defaultdict(list)
for r in rs:by[r["study"]].append(r)
assert len(by)==76 # one study has zero remaining after proposed exclusions
# Historical positive labels are screening evidence, not final clinical attestations.
# Sarcoma reciprocal evidence adds joint primary/timing support conditionally on paired material.
def tissue(r):return r["historical_tissue_grade"]=="PRIMARY_TUMOUR_STUDY_SOURCE_SUPPORTED" or r["sarcoma_untreated_primary_evidence"]=="RECIPROCAL_FROZEN_MICROARRAY_SOURCE_SUPPORTED"
def timing(r):return r["historical_treatment_grade"] in {"COHORT_PRETREATMENT_SOURCE_SUPPORTED","GEO_BASELINE_TEXT_SOURCE_SUPPORTED"} or r["sarcoma_untreated_primary_evidence"]=="RECIPROCAL_FROZEN_MICROARRAY_SOURCE_SUPPORTED"
def group_supported(r):return G.get(r["library_id_LOCAL_ONLY"],{}).get("group_evidence_status","").startswith("SOURCE_SUPPORTED_")
def group_hold(r):return G.get(r["library_id_LOCAL_ONLY"],{}).get("group_evidence_status","").startswith("HOLD_")
# Explicit contextual blockers from executed investigations and original methods.
context={"SRP158272":"HOLD_NO_EXACT_CLINICAL_BRIDGE","SRP157974":"CLINICAL_TABLE_S1_AND_CHEMISTRY_OPEN","SRP057793":"PATIENT_CROSSWALK_OPEN","SRP163173":"OUTCOME_SELECTED_PRETREATMENT_OPEN_ID_CONFLICT2","DRP005276":"PRETREATMENT_AND_ASSAY_CONTEXT_HOLD","SRP185703":"CROSSWALK_HOLD"}
print("V8_JOINT_EVIDENCE_PRIORITY_NOT_CERTIFICATION")
print("LOCKED_INPUTS_PASS",len(rows),"REMAINING",len(rs),"STUDIES",len(by))
print("HISTORICAL_POSITIVE_TREATMENT_LABELS_INCLUDE_BASELINE_NOT_NECESSARILY_ALL_PRIOR_THERAPY")
print("study\tremaining\tprimary_evidence\tpretreatment_or_baseline_evidence\tjoint_evidence\tjoint_without_group_hold\tgroup_supported\tgroup_hold\tbiosample_proxy\tcontext")
rank=[]
for s,a in by.items():
    joint=[r for r in a if tissue(r) and timing(r)]
    row=(s,len(a),sum(map(tissue,a)),sum(map(timing,a)),len(joint),sum(not group_hold(r) for r in joint),sum(map(group_supported,a)),sum(map(group_hold,a)),sum(r["historical_patient_grade"]=="BIOSAMPLE_PROXY_NOT_PATIENT_ID" for r in a),context.get(s,"REVIEW_REQUIRED"))
    rank.append(row)
for row in sorted(rank,key=lambda x:(-x[5],-x[4],-x[1],x[0])):print("\t".join(map(str,row)))
print("TOTAL",len(rs),"JOINT_EVIDENCE",sum(tissue(r) and timing(r) for r in rs),"JOINT_WITHOUT_GROUP_HOLD",sum(tissue(r) and timing(r) and not group_hold(r) for r in rs))
print("JOINT_HISTORICAL_GRADES")
for k,n in sorted(Counter((r["historical_tissue_grade"],r["historical_treatment_grade"],r["historical_patient_grade"]) for r in rs if tissue(r) and timing(r)).items()):
    print("\t".join(k),n)
print("V8_JOINT_EVIDENCE_PRIORITY_COMPLETE")
print("NO_FILES_MODIFIED_NO_CLINICAL_ADMISSIONS")
