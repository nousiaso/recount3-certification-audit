#!/usr/bin/env python3
"""Exact GEO SRX clinical evidence census; no outputs or patient identifiers."""
import csv, gzip, hashlib, re, sys
from collections import Counter, defaultdict
from pathlib import Path

root=Path(sys.argv[1])
paths=dict(zip(("GSE71119","GSE75885"),map(Path,sys.argv[2:])))
assert len(sys.argv)==4
study="SRP057793"
s41=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
rec=root/"26_recount3_external_primary_v1"
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))
source=[r for r in read(s41/"01_full_source_experiment_identity_5221.tsv") if r["study"]==study]
c=[r for r in read(s41/"03_candidate_integrated_evidence_3602.tsv") if r["study"]==study]
all_ids={r["experiment_acc"] for r in source}
keep={r["library_id"] for r in c if r["technical_hold"]=="NO"}
hold={r["library_id"] for r in c if r["technical_hold"]=="YES"}
assert len(all_ids)==149 and len(keep)==140 and len(hold)==8 and not keep&hold
assert keep|hold<=all_ids
native=read(rec/"raw_studies"/study/"metadata"/f"sra.sra.{study}.MD.gz")
native={r["experiment_acc"]:r for r in native}
assert set(native)==all_ids
def seqid(row):
    ids={v.strip().casefold() for segment in row["sample_attributes"].split("|") if ";;" in segment for k,v in [segment.split(";;",1)] if k.strip().lower()=="rna_seq_id"}
    return next(iter(ids)) if len(ids)==1 else ""
def parse(path,series):
    samples={}; cur=None; seen_series=set()
    with gzip.open(path,"rt",encoding="utf-8",errors="strict") as f:
        for line in f:
            line=line.rstrip("\r\n")
            if line.startswith("^"):
                cur=None
                if line.startswith("^SERIES = "): seen_series.add(line.split(" = ",1)[1])
                if line.startswith("^SAMPLE = "):
                    gsm=line.split(" = ",1)[1]
                    assert re.fullmatch(r"GSM\d+",gsm) and gsm not in samples
                    cur=defaultdict(list);samples[gsm]=cur
            elif cur is not None and line.startswith("!Sample_") and " = " in line:
                k,v=line.split(" = ",1);cur[k].append(v.strip())
    assert series in seen_series and samples, "Incorrect or empty GEO family source"
    return samples
def relations(r):
    return set(re.findall(r"\bSRX\d+\b"," ".join(r.get("!Sample_relation",[]))))
def traits(r):
    d=defaultdict(list)
    for v in r.get("!Sample_characteristics_ch1",[]):
        if ":" in v:
            k,val=v.split(":",1);d[k.strip().lower()].append(val.strip())
    return d
def flags(r):
    # Narrow evidence discovery only, not clinical adjudication.
    text=" ".join(r.get("!Sample_treatment_protocol_ch1",[])).casefold()
    return {
        "TREATMENT_PROTOCOL_PRESENT":bool(text.strip()),
        "EXPLICIT_UNTREATED_PRIMARY_PHRASE":bool(re.search(r"\buntreated primary tumou?r\b",text)),
        "PRIMARY_TUMOUR_PHRASE":bool(re.search(r"\bprimary tumou?r\b",text)),
        "BEFORE_TREATMENT_PHRASE":bool(re.search(r"\b(?:before|prior to) (?:any )?treatment\b",text)),
        "METASTATIC_OR_RECURRENCE_WORD":bool(re.search(r"\b(?:metastatic|metastasis|recurrence|recurrent)\b",text)),
        "EXTRACT_PROTOCOL_PRESENT":bool(r.get("!Sample_extract_protocol_ch1")),
    }
combined=defaultdict(dict)
for series,path in paths.items():
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    samples=parse(path,series)
    exact=defaultdict(dict); ambiguous=0
    for gsm,r in samples.items():
        ids=relations(r)
        if len(ids)==1:
            srx=next(iter(ids))
            if srx in all_ids: exact[srx][gsm]=r
        elif ids&all_ids: ambiguous+=1
    print("\nSERIES",series,"SHA256",h.hexdigest(),"GSM_RECORDS",len(samples),
          "EXACT_SOURCE_SRX",len(exact),"AMBIGUOUS_RELATION_RECORDS",ambiguous,
          "MULTIPLE_GSM_PER_EXACT_SRX",sum(len(x)>1 for x in exact.values()))
    print("RETAINED_EXACT",len(keep&exact.keys()),"HELD_EXACT",len(hold&exact.keys()))
    for srx,rs in exact.items():
        for gsm,r in rs.items():
            combined[srx][(series,gsm)]=r
    selected=[(s,r) for s,rs in exact.items() if s in keep for r in rs.values()]
    print("RETAINED_MATCHED_GEO_RECORDS",len(selected))
    print("GEO_TRAIT_KEYS",dict(Counter(k for _,r in selected for k in traits(r))))
    print("EVIDENCE_FLAGS",dict(Counter(k for _,r in selected for k,v in flags(r).items() if v)))
    print("TITLE_EXACT_NATIVE_RNA_SEQ_ID",sum(bool(seqid(native[s])) and any(t.casefold()==seqid(native[s]) for t in r.get("!Sample_title",[])) for s,r in selected))
    print("METASTASIS_TRAIT",dict(Counter(v.casefold() if v.casefold() in {"yes","no","unknown","na","n/a","not available"} else "OTHER_LABEL" for _,r in selected for v in traits(r).get("metastasis",[]))))
print("\nUNION_RETAINED_EXACT",len(keep&combined.keys()),"UNMATCHED",len(keep-combined.keys()))
print("RETAINED_IN_BOTH_SERIES",sum(len({g for g,_ in combined[s]})==2 for s in keep&combined.keys()))
print("RETAINED_DISTINCT_GSM_MULTIPLICITY",dict(Counter(len({gsm for _,gsm in combined[s]}) for s in keep&combined.keys())))
for key in flags({}):
    counts=Counter()
    for s in keep:
        vals=[flags(r)[key] for r in combined.get(s,{}).values()]
        counts["NO_EXACT_GEO" if not vals else "ALL_RECORDS_POSITIVE" if all(vals) else "MIXED_RECORDS" if any(vals) else "NO_POSITIVE_RECORD"]+=1
    print("RETAINED_UNION_FLAG",key,dict(counts))
print("V8_SRP057793_GEO_EXACT_CLINICAL_CENSUS_COMPLETE")
print("NO_FILES_MODIFIED_NO_CLINICAL_ADMISSIONS")
