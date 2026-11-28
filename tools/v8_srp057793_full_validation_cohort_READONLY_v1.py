#!/usr/bin/env python3
"""Full published-cohort scope reconciliation, including held/source-only entries.
No patient-ID inference from specimen labels; aggregate-only output."""
import csv,gzip,hashlib,re,sys
from collections import Counter,defaultdict
from pathlib import Path
root=Path(sys.argv[1]);repo=root/"RECOUNT3_CERTIFICATION_GIT"
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
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

def norm(v): return str(v).strip().casefold()
def one(values):
    v={norm(x) for x in values if norm(x) not in {"","--","na","n/a","none","not available"}}
    return next(iter(v)) if len(v)==1 else ""
def trait(r,k): return one(traits(r).get(k,[]))
def title(r): return one(r.get("!Sample_title",[]))
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()
def field_text(r,k): return " ".join(r.get(k,[])).strip()

def locked(p,h):
    assert sha(p)==h,"SOURCE_HASH_MISMATCH "+str(p)
    return p
def exact_index(rows,key):
    d={r[key]:r for r in rows};assert len(d)==len(rows) and all(d)
    return d
study="SRP057793"
s41=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
full=read(locked(s41/"01_full_source_experiment_identity_5221.tsv","95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa"))
cand=read(locked(s41/"03_candidate_integrated_evidence_3602.tsv","aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81"))
all_ids={r["experiment_acc"] for r in full if r["study"]==study}
c=exact_index([r for r in cand if r["study"]==study],"library_id")
keep={e for e,r in c.items() if r["technical_hold"]=="NO"}
held={e for e,r in c.items() if r["technical_hold"]=="YES"}
assert (len(all_ids),len(keep),len(held))==(149,140,8) and set(c)==keep|held and set(c)<=all_ids
native=exact_index(read(locked(root/"26_recount3_external_primary_v1/raw_studies/SRP057793/metadata/sra.sra.SRP057793.MD.gz","9cdac726793a9d8be3275bd7d1ef4a0a8adb759dfb4e445e586b3077cd5596f6")),"experiment_acc")
assert set(native)==all_ids
def seqid(row):
    return one([v for segment in row["sample_attributes"].split("|") if ";;" in segment for k,v in [segment.split(";;",1)] if k.strip().lower()=="rna_seq_id"])
cache=repo/"private_data/v8_srp057793_geo_source_v1"
a=parse(locked(cache/"GSE71119_family.soft.gz","249d5e0a299bf35a93802ca7e3a61e2aa8a36f100930300ac9c87400eae87858"),"GSE71119")
b=parse(locked(cache/"GSE71118_family.soft.gz","ea11c2ff661b47da23ed416e3e188cb80a6c7171edbd5d5389d22f65f6a4cc8c"),"GSE71118")
assert (len(a),len(b))==(136,312)
bm=defaultdict(list)
for g,r in b.items():
    if title(r):bm[title(r)].append(g)
pairs=[];status=Counter();exp_seen=set()
for g,r in a.items():
    ex=relations(r)
    assert len(ex)==1 and next(iter(ex)) in all_ids
    e=next(iter(ex));assert e not in exp_seen;exp_seen.add(e)
    scope="RETAINED" if e in keep else "HISTORICAL_HOLD" if e in held else "SOURCE_ONLY"
    token=trait(r,"paired microarray")
    if not token:status[(scope,"NO_UNIQUE_PAIR_TOKEN")]+=1;continue
    if len(bm[token])!=1:status[(scope,"MICROARRAY_TITLE_ABSENT_OR_AMBIGUOUS")]+=1;continue
    mg=bm[token][0];m=b[mg]
    if title(r)!=seqid(native[e]):status[(scope,"RNA_NATIVE_TITLE_MISMATCH")]+=1;continue
    if trait(m,"paired rna-seq")!=title(r):status[(scope,"RECIPROCAL_MISMATCH")]+=1;continue
    untreated=bool(re.search(r"\buntreated primary tumou?rs?\b",field_text(m,"!Sample_treatment_protocol_ch1"),re.I))
    frozen=(trait(r,"material support")==trait(m,"material support")=="frozen")
    status[(scope,"RECIPROCAL_PAIR")]+=1
    pairs.append((e,g,mg,r,m,scope,untreated,frozen))
assert len(exp_seen)==136
assert len({x[0] for x in pairs})==len(pairs)
print("FULL_GSE71119_PAIR_STATUS",dict(status))
print("FULL_RECIPROCAL_PAIRS",len(pairs),"DISTINCT_MICROARRAY_GSM",len({x[2] for x in pairs}))
print("FULL_PAIRED_UNTREATED_FROZEN_BY_SCOPE",dict(Counter((x[5],x[6],x[7]) for x in pairs)))
print("MATCH_TO_PUBLISHED_95_CASE_COUNT",len(pairs)==95,"CARDINALITY_NOT_PATIENT_ID_PROOF")
ret=[x for x in pairs if x[5]=="RETAINED" and x[6] and x[7]]
assert len(ret)==90
prep=repo/"private_data/v8_joint461_full_G026_QC_normalization_v1"
import json
prov=json.loads((prep/"input_provenance.json").read_text())
m=read(locked(prep/"461_amended_technical_input_manifest_LOCAL_ONLY.tsv",prov["amended_manifest_sha256"]))
prepared={r["library_id_LOCAL_ONLY"] for r in m if r["study"]==study}
assert prepared=={x[0] for x in ret}
print("EXACT_PROJECTED90_MEMBERSHIP_PASS",len(prepared))
print("FULL_PAIR_METASTASIS_LABEL_AGREEMENT",dict(Counter((x[5],trait(x[3],"metastasis")==trait(x[4],"metastasis"),bool(trait(x[3],"metastasis") and trait(x[4],"metastasis"))) for x in pairs)))
# New paired microarray aliases were not included in the earlier RNA GEO alias screen.
atlas=read(locked(root/"14_projection_use_analysis_ready_atlas_v1_0_3/03_projection_analysis_ready_core.tsv.gz","7f6852e829f7aaa11a26f86b44247c9c8c0840958494b7766db1ad6abd180092"))
ids={r["sample_id"].strip().upper() for r in atlas}
assert len(atlas)==len(ids)==9465
print("PAIRED_MICROARRAY_GSM_ATLAS_HITS_FULL",sum(x[2] in ids for x in pairs))
print("PAIRED_MICROARRAY_GSM_ATLAS_HITS_RETAINED90",sum(x[2] in ids for x in ret))
# Field coverage without exposing values: determine whether the expanded scope has
# explicit patient keys and/or linked clinical endpoint annotations.
for label,rr in [("ALL_PAIRED_RNA",[x[3] for x in pairs]),("ALL_PAIRED_MICROARRAY",[x[4] for x in pairs]),("RETAINED90_RNA",[x[3] for x in ret]),("RETAINED90_MICROARRAY",[x[4] for x in ret])]:
    fields=sorted(set().union(*(traits(r).keys() for r in rr))) if rr else []
    for k in fields:
        vals=[trait(r,k) for r in rr]
        nonempty=[v for v in vals if v]
        print("TRAIT_COVERAGE",label,repr(k),"N",len(rr),"NONEMPTY",len(nonempty),"DISTINCT",len(set(nonempty)))
    print("PATIENT_FIELD_NAMES",label,[k for k in fields if re.search(r"patient|donor|subject|conticabase|individual",k,re.I)])
print("V8_SRP057793_FULL_VALIDATION_COHORT_REVIEW_COMPLETE")
print("NO_PATIENT_ID_INFERENCE_NO_CLINICAL_ADMISSIONS_NO_REFIT")
