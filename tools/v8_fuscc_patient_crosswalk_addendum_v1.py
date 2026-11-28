#!/usr/bin/env python3
"""Persist exact clinical patient/project identity evidence; no final admission."""
import csv,hashlib,sys,re,json
from pathlib import Path
from collections import Counter
import pandas as pd
assert len(sys.argv)==3
root=Path(sys.argv[1]).resolve();repo=root/"RECOUNT3_CERTIFICATION_GIT"
out=Path(sys.argv[2]).resolve()
assert out.parent==(repo/"private_data").resolve() and out.name=="v8_fuscc_patient_crosswalk_addendum_v1" and not out.exists()
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def norm(x):return str(x).strip().upper() if not pd.isna(x) else ""
def index(rows,k):
    d={r[k]:r for r in rows};assert len(d)==len(rows)
    return d
def verify_audit(d,required,marker,stdout="stdout.txt"):
    entries={}
    for l in (d/"SHA256SUMS.txt").read_text().splitlines():
        m=re.fullmatch(r"([0-9a-f]{64}) [ *](.+)",l);assert m
        p=Path(m[2]);p=(d/p).resolve() if not p.is_absolute() else p.resolve()
        assert p not in entries
        entries[p]=m[1]
    assert {p.resolve() for p in required}<=entries.keys()
    for p,h in entries.items():assert sha(p)==h,str(p)
    assert marker in (d/stdout).read_text().splitlines()
workbook=repo/"private_data/v8_fuscc_original_clinical_source_v1/mmc2.xlsx"
ann=repo/"private_data/v8_fuscc_author_annotation_addendum_v1/360_fuscc_author_annotation_evidence_LOCAL_ONLY.tsv"
ledger=repo/"private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv"
prep=repo/"private_data/v8_joint461_full_G026_QC_normalization_v1"
projection=repo/"private_data/v8_joint461_fixed500_compatibility_v1_0_1"
audit=repo/"diagnostics/20261128_v8_fuscc_original_clinical_table_v1"
verify_audit(audit,[workbook,audit/"stdout.txt",repo/"tools/v8_fuscc_original_clinical_table_READONLY_v1.py"],"INDEPENDENT_CANONICAL_CONTENT_REFERENCE_PASS")
verify_audit(repo/"diagnostics/20261128_v8_joint461_full_G026_QC_normalization_v1",[prep/"461_amended_technical_input_manifest_LOCAL_ONLY.tsv",prep/"input_provenance.json"],"V8_JOINT461_FULL_G026_PREPARATION_COMPLETE_PENDING_REVIEW")
pa=repo/"diagnostics/20261128_v8_joint461_fixed500_compatibility_v1_0_1"
verify_audit(pa,[projection/"SHA256SUMS.txt"],"V8_JOINT461_FIXED500_COMPATIBILITY_V1_0_1_COMPLETE_PENDING_REVIEW","03_projection_stdout.txt")
# Validate the three specific FUSCC projection outputs against the nested audit.
entries={}
for l in (projection/"SHA256SUMS.txt").read_text().splitlines():
    h,n=l.split(maxsplit=1);entries[Path(n.lstrip("* ")).name]=h
assert sha(workbook)=="1581533099f7a7ce15d5ab39551fa93855c25b30fc624160b3940700151becb7"
assert sha(ann)=="1cae682fe3ee006df02e22654f558c0ee2e28ab66fd5a422c67374609e1d1dd5"
assert sha(ledger)=="943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0"
t=pd.read_excel(workbook,sheet_name="Table S1",header=1,engine="openpyxl")
assert t.shape==(465,50)
cols=["Project_ID","RNA_Seqencing","mRNA_Subtype","Intrinsic_Subtype","Age_at_surgery",
      "Histology","Histology_extended","Grade","T","N","Size_cm","LN_examed","LN_positive",
      "Ki67","sTILs","iTILs","Followup_Month","RFS_Status","RFS_time_Days","RFS_time_Months",
      "Date_of_surgery","Chemotherapy","Paclitaxel","Radiotherapy","No_Chemotherapy",
      "BRCA1/2 MUT 20160718","ER_IHC_score","PR_IHC_score","ERBB2_IHC_score","ERBB2_FISH","Sex"]
assert hashlib.sha256(t[cols].to_csv(index=False,lineterminator="\n",float_format="%.10g").encode()).hexdigest()=="96ab7091b1bbdcc7ae10d9926ff4f2f4c998a434e3372e2b459493a30ef19646"
keys=[norm(x) for x in t.Project_ID]
assert all(keys) and len(set(keys))==465
t=t.assign(_key=keys)
rna=t[t.RNA_Seqencing.map(norm)=="YES"].set_index("_key",verify_integrity=True)
assert len(rna)==360
a=read(ann)
assert len(a)==len({r["experiment_acc_LOCAL_ONLY"] for r in a})==360
assert len({r["author_RNA_specimen_base_LOCAL_ONLY"] for r in a})==360
assert set(rna.index)=={norm(r["author_RNA_specimen_base_LOCAL_ONLY"]) for r in a}
assert all(r["final_clinical_approval"]=="NO_UNSIGNED" for r in a)
assert Counter(r["membership"] for r in a)=={"RETAINED334":334,"HISTORICAL_HOLD26":26}
kept={r["experiment_acc_LOCAL_ONLY"] for r in a if r["membership"]=="RETAINED334"}
l=index([r for r in read(ledger) if r["study"]=="SRP157974"],"library_id_LOCAL_ONLY")
m=index([r for r in read(prep/"461_amended_technical_input_manifest_LOCAL_ONLY.tsv") if r["study"]=="SRP157974"],"library_id_LOCAL_ONLY")
assert set(l)==set(m)==kept and len(kept)==334
assert all(r["final_approval"]=="NO_UNSIGNED" for r in l.values())
assert all(r["input_status"]=="INPUT_READY_FOR_FULL_G026_QC_NOT_RELEASED" and r["final_clinical_approval"]=="NO_UNSIGNED" for r in m.values())
for branch in ("CPM","TMM_CPM","UQ_CPM"):
    f=projection/f"SRP157974_{branch}_projection_QC_LOCAL_ONLY.tsv"
    assert sha(f)==entries[f.name]
    assert set(index(read(f),"sample_id"))==kept
rows=[]
# Preserve clinical annotations verbatim as parsed, no inferred treatment timing,
# no survival date or free-text export. Workbook remains the authoritative source.
clinical_cols=[c for c in cols if c not in {"Project_ID","Date_of_surgery"}]
clinical_cols+=["Exome_Seqencing","OncoScan_Array","HRD","LSTm","Telomeric AI","LOH","Ascat_ACF","Ascat_Ploidy","Surgery"]
def scalar(x):return "" if pd.isna(x) else x
for r in a:
    e=r["experiment_acc_LOCAL_ONLY"];k=norm(r["author_RNA_specimen_base_LOCAL_ONLY"]);c=rna.loc[k]
    assert norm(c["mRNA_Subtype"])==norm(r["author_subtype_old"])
    is_retained=e in kept
    suffix=r["rep_suffix_status"]
    assert suffix in {"AGREES","DISAGREES_UNRESOLVED"}
    if not is_retained:
        state="HISTORICAL_QC_HOLD_PRESERVED"
    elif suffix!="AGREES":
        state="HOLD_LIBRARY_SUFFIX_PROVENANCE"
    else:
        state="PROPOSED_PATIENT_LINKED_REVIEW_CANDIDATE"
    row={**r,"clinical_project_id_LOCAL_ONLY":k,
         "clinical_patient_crosswalk_status":"EXACT_ORIGINAL_TABLE_S1_PROJECT_ID_RNA_YES",
         "distinct_patient_project_record_within_cohort":"YES",
         "patient_independence_certified":"NO_UNSIGNED",
         "clinical_primary_tumour_evidence":"PUBLISHED_PRIMARY_TNBC_COHORT_AND_EXACT_TUMOUR_MATERIAL",
         "clinical_pretreatment_evidence":"PUBLISHED_COHORT_LEVEL_SUPPORT_NOT_INDIVIDUAL_TIMING",
         "technical_projection_status":"FIXED500_THREE_BRANCHES_EXECUTED" if is_retained else "HISTORICAL_QC_HOLD_NOT_PROJECTED",
         "proposed_review_disposition":state,"final_clinical_approval":"NO_UNSIGNED"}
    row.update({"clinical_"+cname:scalar(c[cname]) for cname in clinical_cols})
    rows.append(row)
counts=Counter(r["proposed_review_disposition"] for r in rows)
assert counts=={"PROPOSED_PATIENT_LINKED_REVIEW_CANDIDATE":331,"HOLD_LIBRARY_SUFFIX_PROVENANCE":3,"HISTORICAL_QC_HOLD_PRESERVED":26}
out.mkdir(mode=0o700)
target=out/"360_fuscc_patient_clinical_crosswalk_LOCAL_ONLY.tsv"
with target.open("x",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
target.chmod(0o600)
print("EXACT_PATIENT_PROJECT_LINKS",len(rows),"DISTINCT",len({r["clinical_project_id_LOCAL_ONLY"] for r in rows}))
print("PROPOSED_DISPOSITIONS",dict(counts))
print("PRIVATE_CROSSWALK_SHA256",sha(target))
print("V8_FUSCC_PATIENT_CROSSWALK_ADDENDUM_COMPLETE")
print("NO_FINAL_PATIENT_INDEPENDENCE_SIGNATURE_NO_CLINICAL_ADMISSIONS_NO_REFIT")
