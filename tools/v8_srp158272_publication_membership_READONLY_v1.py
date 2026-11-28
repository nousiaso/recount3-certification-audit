#!/usr/bin/env python3
"""SRP158272: reconcile native specimen evidence with frozen candidate partitions.
Read-only, no patient/specimen values printed; publication subset membership is not inferred.
"""
import csv,gzip,hashlib,io,json,re,sys,zipfile
from pathlib import Path
from collections import Counter,defaultdict

root=Path(sys.argv[1])
repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1"
study="SRP158272"
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""):h.update(block)
    return h.hexdigest()
def read(path):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))
def norm(x):return str(x or "").strip().casefold()
def present(x):return norm(x) not in {"","na","n/a","none","null","nan","unknown","not available","not applicable","--"}
def attrs(r):
    d=defaultdict(set);bad=0
    for a in r.get("sample_attributes","").split("|"):
        if not a.strip():continue
        if ";;" not in a:bad+=1;continue
        k,v=a.split(";;",1)
        k=re.sub(r"[^a-z0-9]+","_",norm(k)).strip("_")
        if k and present(v):d[k].add(norm(v))
    return d,bad

ledger=repo/"private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv"
assert sha(ledger)=="943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0"
lr=[r for r in read(ledger) if r["study"]==study]
keep={r["library_id_LOCAL_ONLY"] for r in lr}
assert len(lr)==len(keep)==129
assert all(r["proposed_non_tumour_exclusion"]=="NO" and r["final_approval"]=="NO_UNSIGNED" for r in lr)
s41=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
cand=s41/"03_candidate_integrated_evidence_3602.tsv"
full=s41/"01_full_source_experiment_identity_5221.tsv"
assert sha(cand)=="aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81"
assert sha(full)=="95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa"
cr=[r for r in read(cand) if r["study"]==study]
fr=[r for r in read(full) if r["study"]==study]
ci={r["library_id"]:r for r in cr}
fi={r["experiment_acc"]:r for r in fr}
assert len(ci)==len(cr) and len(fi)==len(fr)
assert keep<=ci.keys()<=fi.keys()
assert all(ci[k]["technical_hold"]=="NO" for k in keep)
hold={k for k,r in ci.items() if r["technical_hold"]=="YES"}
other=set(ci)-keep-hold
sourceonly=set(fi)-set(ci)
print("V8_SRP158272_PUBLICATION_MEMBERSHIP_REVIEW")
print("MEMBERSHIP",{"source":len(fi),"candidate":len(ci),"retained":len(keep),"technical_hold":len(hold),"other_candidates":len(other),"source_only":len(sourceonly)})
for f in ("historical_patient_grade","historical_treatment_grade","historical_tissue_grade","v6_unsigned_evidence_bucket"):
    print("HISTORICAL",f,dict(Counter(r[f] for r in lr)))
native=rec/"raw_studies"/study/"metadata"/f"sra.sra.{study}.MD.gz"
assert sha(native)=="bfbdf336decb8295c204d17f0268740a7c5f00eeaebd66fddc50ddbc7b03e874"
nr=read(native);ni=defaultdict(list)
for r in nr:ni[r["experiment_acc"]].append(r)
print("NATIVE_ROWS",len(nr),"DISTINCT_EXPERIMENTS",len(ni),"FULL_LEDGER_MISSING",len(fi.keys()-ni.keys()),"NATIVE_OUTSIDE_FULL_LEDGER",len(ni.keys()-fi.keys()))
assert fi.keys()<=ni.keys() and all(len(ni[k])==1 for k in keep)
print("NATIVE_SHA256",sha(native))
for label,ids in (("RETAINED",keep),("TECHNICAL_HOLD",hold),("OTHER_CANDIDATE",other),("SOURCE_ONLY",sourceonly)):
    rs=[r for k in sorted(ids) for r in ni[k]]
    print("\nSCOPE",label,"EXPERIMENTS",len(ids),"METADATA_ROWS",len(rs))
    if not rs:continue
    for f in ("library_strategy","library_selection","library_layout","platform_model"):
        print("ASSAY",f,dict(Counter(r[f] for r in rs)))
    av={k:attrs(ni[k][0])[0] for k in ids if len(ni[k])==1}
    keys=sorted({f for v in av.values() for f in v})
    print("ATTRIBUTE_UNPARSED_SEGMENTS",sum(attrs(r)[1] for r in rs))
    for f in keys:
        vals=[v[f] for v in av.values() if f in v]
        one=[next(iter(v)) for v in vals if len(v)==1]
        sizes=Counter(one)
        print("ATTRIBUTE",f,"POPULATED",len(vals),"MULTIVALUE",len(vals)-len(one),"DISTINCT_SINGLE_VALUES",len(sizes),"REUSED_SINGLE_VALUES",sum(n>1 for n in sizes.values()),"MAX_MULTIPLICITY",max(sizes.values(),default=0))
    for f in ("sample_name","sample_title","external_id","library_name","sample_acc","sample_description","experiment_title","library_construction_protocol"):
        values=[norm(r.get(f)) for r in rs if present(r.get(f))]
        c=Counter(values)
        print("NATIVE_FIELD",f,"POPULATED",len(values),"DISTINCT",len(c),"REUSED_VALUES",sum(n>1 for n in c.values()),"MAX_MULTIPLICITY",max(c.values(),default=0))
    # Only count text evidence by exact field/attribute scope; never print raw patient values.
    patterns={
      "FFPE":r"\bffpe\b|formalin.{0,25}paraffin",
      "FROZEN":r"\bfrozen\b",
      "PRIMARY":r"\bprimary\b",
      "METASTATIC":r"\bmetast(?:atic|asis|ases)\b",
      "NORMAL_OR_ADJACENT":r"\bnormal\b|\badjacent\b",
      "UNTREATED_OR_TREATMENT_NAIVE":r"\buntreated\b|treatment[- ]na[iï]ve|no prior treatment|no neoadjuvant",
      "PRETREATMENT_OR_BASELINE":r"pre[- ]?treatment|\bbaseline\b",
      "SURGERY":r"\bsurg(?:ery|ical)\b|\bresection\b",
      "RIBO_DEPLETION":r"ribo.?zero|ribosomal.{0,25}deplet|rrna.{0,25}deplet|deplet.{0,25}rrna",
      "POLYA":r"poly\s*[\(-]?\s*a\b|polyadenylat|oligo.?dt",
      "DECEASED":r"\bdeceased\b|\bdead\b",
    }
    text_fields=("sample_name","sample_title","sample_description","experiment_title","library_construction_protocol","design_description","study_description","study_abstract")
    for f in text_fields:
        counts={name:sum(bool(re.search(pat,norm(r.get(f)),re.I)) for r in rs) for name,pat in patterns.items()}
        print("TEXT_EVIDENCE",f,{k:v for k,v in counts.items() if v})
    for f in keys:
        counts={name:sum(any(re.search(pat,x,re.I) for x in a.get(f,set())) for a in av.values()) for name,pat in patterns.items()}
        if any(counts.values()):print("ATTRIBUTE_TEXT_EVIDENCE",f,{k:v for k,v in counts.items() if v})
    # Identify whether existing public clinical data could be joined via exact GEO aliases.
    for f in ("sample_name","sample_title","external_id","library_name","sample_description","experiment_title","sample_attributes"):
        gs=[set(re.findall(r"(?<![A-Za-z0-9])GSM\d+(?![A-Za-z0-9])",r.get(f,""),re.I)) for r in rs]
        if any(gs):print("GSM_ALIAS_FIELD",f,"COVERED",sum(bool(x) for x in gs),"MULTIVALUED",sum(len(x)>1 for x in gs),"DISTINCT",len(set().union(*gs)))
# Existing V3 public review, not a new study-level certification.
v3=rec/"RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"
with zipfile.ZipFile(v3) as z:
    paths=[p for p in z.namelist() if p.endswith("/results/02_ALL_176_STUDY_SOURCE_REVIEW_UNSIGNED.tsv")]
    assert len(paths)==1
    data=list(csv.DictReader(io.StringIO(z.read(paths[0]).decode("utf-8-sig")),delimiter="\t"))
assert len(data)==176
v=[r for r in data if r["study"]==study];assert len(v)==1
print("\nV3_ARCHIVE_SHA256",sha(v3))
print("V3_STUDY_REVIEW_COLUMNS",list(v[0]))
for f in ("public_source_review_status","overall_study_disposition","review_action_from_publication","FINAL_STUDY_CLINICAL_APPROVAL_AS_OF_V3","qa_missing_pretreatment"):
    x=v[0].get(f,"")
    print("V3_REVIEW",f,x if re.fullmatch(r"[A-Z0-9_./-]{1,100}",x) else "TEXT_SUPPRESSED")
print("V8_SRP158272_PUBLICATION_MEMBERSHIP_REVIEW_COMPLETE")
print("NO_PUBLICATION_SUBSET_MEMBERSHIP_INFERRED_NO_FILES_MODIFIED_NO_CLINICAL_ADMISSIONS")
