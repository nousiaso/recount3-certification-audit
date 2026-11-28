#!/usr/bin/env python3
"""Original Jiang Table S1 intake; aggregate-only, no identity or clinical release."""
import csv,hashlib,sys,re,math,zipfile
from pathlib import Path
from collections import Counter
from importlib.metadata import version
import pandas as pd
import openpyxl
assert len(sys.argv)==3
root=Path(sys.argv[1]).resolve();repo=root/"RECOUNT3_CERTIFICATION_GIT";workbook=Path(sys.argv[2]).resolve()
assert workbook.parent==(repo/"private_data/v8_fuscc_original_clinical_source_v1").resolve()
assert workbook.is_file() and zipfile.is_zipfile(workbook),"NOT_XLSX_CONTAINER"
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""):h.update(block)
    return h.hexdigest()
def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def check(ok,msg):
    if not ok:raise ValueError(msg)
print("WORKBOOK_SHA256",sha(workbook),"BYTES",workbook.stat().st_size)
print("PARSER_VERSIONS","pandas",pd.__version__,"openpyxl",version("openpyxl"))
with pd.ExcelFile(workbook,engine="openpyxl") as book:
    print("WORKSHEETS",book.sheet_names)
    check("Table S1" in book.sheet_names,"ORIGINAL_TABLE_S1_NOT_PRESENT")
    t=pd.read_excel(book,sheet_name="Table S1",header=1)
print("TABLE_S1_SHAPE",t.shape)
print("TABLE_S1_COLUMNS",list(t.columns))
check(t.shape==(465,50),"ORIGINAL_EXPECTED_TABLE_SHAPE_MISMATCH")
cols=["Project_ID","RNA_Seqencing","mRNA_Subtype","Intrinsic_Subtype","Age_at_surgery",
      "Histology","Histology_extended","Grade","T","N","Size_cm","LN_examed","LN_positive",
      "Ki67","sTILs","iTILs","Followup_Month","RFS_Status","RFS_time_Days","RFS_time_Months",
      "Date_of_surgery","Chemotherapy","Paclitaxel","Radiotherapy","No_Chemotherapy",
      "BRCA1/2 MUT 20160718","ER_IHC_score","PR_IHC_score","ERBB2_IHC_score","ERBB2_FISH","Sex"]
check(all(c in t.columns for c in cols),"ORIGINAL_EXPECTED_COLUMNS_MISSING")
canonical=hashlib.sha256(t[cols].to_csv(index=False,lineterminator="\n",float_format="%.10g").encode("utf-8")).hexdigest()
expected="96ab7091b1bbdcc7ae10d9926ff4f2f4c998a434e3372e2b459493a30ef19646"
print("CANONICAL_31_COLUMN_SHA256",canonical,"INDEPENDENT_REFERENCE_MATCH",canonical==expected)
print("CANONICAL_REFERENCE_PARSER_VERSIONS_MATCH",pd.__version__=="2.2.3" and version("openpyxl")=="3.1.5")
def norm(x):return str(x).strip().upper() if not pd.isna(x) else ""
keys=[norm(x) for x in t.Project_ID]
check(all(keys) and len(set(keys))==465,"PROJECT_KEYS_MISSING_OR_DUPLICATE_AFTER_CASE_TRIM")
flags=[norm(x) for x in t.RNA_Seqencing]
print("RNA_AVAILABILITY_FLAGS",dict(Counter(flags)))
check(sum(x=="YES" for x in flags)==360,"RNA_AVAILABILITY_COUNT_MISMATCH")
t=t.copy();t["_key"]=keys
RNA=t[[x=="YES" for x in flags]].copy()
rn=set(RNA["_key"])
ann=repo/"private_data/v8_fuscc_author_annotation_addendum_v1/360_fuscc_author_annotation_evidence_LOCAL_ONLY.tsv"
check(sha(ann)=="1cae682fe3ee006df02e22654f558c0ee2e28ab66fd5a422c67374609e1d1dd5","AUTHOR_ADDENDUM_HASH_MISMATCH")
a=read(ann)
check(len(a)==360 and len({r["experiment_acc_LOCAL_ONLY"] for r in a})==360,"AUTHOR_EXPERIMENT_AXIS")
base=[norm(r["author_RNA_specimen_base_LOCAL_ONLY"]) for r in a]
check(all(re.fullmatch(r"FUSCCTNBC\d+",x) for x in base) and len(set(base))==360,"AUTHOR_BASE_AXIS")
print("EXACT_PROJECT_ID_CASE_TRIM_JOIN","MATCHED",len(set(base)&rn),"AUTHOR_ONLY",len(set(base)-rn),"CLINICAL_RNA_ONLY",len(rn-set(base)))
check(set(base)==rn,"CLINICAL_RNA_SET_NOT_EXACT")
print("MEMBERSHIP",dict(Counter(r["membership"] for r in a)))
check(Counter(r["membership"] for r in a)=={"RETAINED334":334,"HISTORICAL_HOLD26":26},"HISTORICAL_MEMBERSHIP_MISMATCH")
check(all(r["final_clinical_approval"]=="NO_UNSIGNED" for r in a),"UNEXPECTED_APPROVAL_CHANGE")
lookup=RNA.set_index("_key",verify_integrity=True)
ret=[norm(r["author_RNA_specimen_base_LOCAL_ONLY"]) for r in a if r["membership"]=="RETAINED334"]
held=[norm(r["author_RNA_specimen_base_LOCAL_ONLY"]) for r in a if r["membership"]=="HISTORICAL_HOLD26"]
for scope,ss in (("ALL360",RNA),("RETAINED334",lookup.loc[ret]),("HELD26",lookup.loc[held])):
    print("SCOPE",scope,"N",len(ss),"SEX",dict(Counter(norm(x) for x in ss.Sex)))
    # Only counts and distributions; never print a patient key, date, or row.
    for f in ("mRNA_Subtype","Intrinsic_Subtype","Grade","T","N","RFS_Status","Chemotherapy","Paclitaxel","Radiotherapy","No_Chemotherapy"):
        v=[norm(x) for x in ss[f]]
        # Do not emit arbitrary free text or rare categories from patient fields.
        print("FIELD",f,"NONMISSING",sum(bool(x) for x in v),"DISTINCT_NONMISSING",len({x for x in v if x}))
    for f in ("Age_at_surgery","Ki67","sTILs","iTILs","RFS_time_Days","RFS_time_Months"):
        v=pd.to_numeric(ss[f],errors="coerce")
        finite=v[v.map(lambda x:math.isfinite(x))]
        print("NUMERIC",f,"FINITE",len(finite),"MISSING_OR_NONNUMERIC",len(v)-len(finite),"MEDIAN",float(finite.median()) if len(finite) else None)
    rfs=pd.to_numeric(ss.RFS_Status,errors="coerce")
    print("RFS_STATUS_BINARY_COUNTS","ZERO",int((rfs==0).sum()),"ONE",int((rfs==1).sum()),"OTHER_OR_MISSING",int((~rfs.isin([0,1])).sum()))
print("OLD_SUBTYPE_CASE_TRIM_AGREEMENT",sum(norm(lookup.loc[k,"mRNA_Subtype"])==norm(r["author_subtype_old"]) for k,r in zip(base,a)))
print("RETAINED_REP_SUFFIX_DISAGREEMENTS",sum(r["rep_suffix_status"]!="AGREES" for r in a if r["membership"]=="RETAINED334"))
check(all(norm(x)=="FEMALE" for x in t.Sex),"EXPECTED_COHORT_SEX_CHECK_FAILED")
rfs=pd.to_numeric(RNA.RFS_Status,errors="coerce")
check(rfs.isin([0,1]).all() and int(rfs.sum())==49,"INDEPENDENT_RNA_RELAPSE_COUNT_MISMATCH")
if canonical!=expected:
    print("SOURCE_CANONICAL_HASH_REVIEW_REQUIRED_NO_ADMISSION")
else:
    print("INDEPENDENT_CANONICAL_CONTENT_REFERENCE_PASS")
print("V8_FUSCC_ORIGINAL_CLINICAL_TABLE_INTAKE_COMPLETE")
print("NO_PATIENT_RECORDS_EXPORTED_NO_FROZEN_INPUT_CHANGES_NO_CLINICAL_ADMISSIONS")
