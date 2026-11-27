#!/usr/bin/env python3
"""Read-only exact SRA join and aggregate clinical/assay labels; no identifiers."""
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

root=Path(sys.argv[1])
study="SRP057793"
rec=root/"26_recount3_external_primary_v1"
s41=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

source=[r for r in read(s41/"01_full_source_experiment_identity_5221.tsv") if r["study"]==study]
c=[r for r in read(s41/"03_candidate_integrated_evidence_3602.tsv") if r["study"]==study]
all_ids={r["experiment_acc"] for r in source}
retained={r["library_id"] for r in c if r["technical_hold"]=="NO"}
held={r["library_id"] for r in c if r["technical_hold"]=="YES"}
other=all_ids-retained-held
assert len(source)==len(all_ids)==149 and len(c)==148
assert (len(retained),len(held),len(other))==(140,8,1)
assert not retained&held and (retained|held)<=all_ids
native=read(rec/"raw_studies"/study/"metadata"/f"sra.sra.{study}.MD.gz")
lookup={r["experiment_acc"]:r for r in native}
assert len(native)==len(lookup)==149 and set(lookup)==all_ids

def attrs(row):
    d=defaultdict(set)
    for segment in row["sample_attributes"].split("|"):
        if not segment.strip(): continue
        assert ";;" in segment
        k,v=segment.split(";;",1)
        if k.strip() and v.strip(): d[k.strip().lower()].add(v.strip())
    return d

attributes={k:attrs(v) for k,v in lookup.items()}
def category(k,field):
    v=attributes[k].get(field,set())
    if len(v)!=1:return "<MISSING_OR_MULTIVALUE>"
    return next(iter(v))

print("=== PUBLIC STUDY TITLE AND LIBRARY METADATA ===")
for field in ("study_title","library_strategy","library_source","library_selection","library_layout","platform_model"):
    print(field,dict(Counter(lookup[k][field] for k in retained)))
print("\n=== RETAINED AND HELD CLINICAL LABELS ===")
for name,ids in (("RETAINED",retained),("HELD",held)):
    print("SCOPE",name,"N",len(ids))
    for field in ("disease","tissue","sample_type","cell_type","isolate"):
        print(field,dict(Counter(category(k,field) for k in ids)))
    print("TISSUE_DISEASE_CROSSTAB")
    for (tissue,disease),n in sorted(Counter((category(k,"tissue"),category(k,"disease")) for k in ids).items()):
        print(repr(tissue),repr(disease),n)
print("\n=== IDENTIFIER RELATIONSHIPS WITHOUT IDENTIFIERS ===")
# String equality only. No inferred patient grouping from prefixes.
def norm(s):return str(s or "").strip().casefold()
for field in ("sample_name","sample_title","sample_description","experiment_title","external_id"):
    vals=[norm(lookup[k][field]) for k in retained]
    ids=[norm(category(k,"rna_seq_id")) for k in retained]
    print(field,"DISTINCT",len(set(vals)),
          "EXACT_RNA_SEQ_ID_EQUAL",sum(v==i and bool(i) for v,i in zip(vals,ids)),
          "CONTAINS_RNA_SEQ_ID",sum(bool(i) and bool(re.search(r"(?<![A-Za-z0-9])"+re.escape(i)+r"(?![A-Za-z0-9])",v)) for v,i in zip(vals,ids)))
print("\nV8_SRP057793_NATIVE_CLINICAL_CATEGORIES_COMPLETE")
print("NO_FILES_MODIFIED_NO_CLINICAL_ADMISSIONS")
