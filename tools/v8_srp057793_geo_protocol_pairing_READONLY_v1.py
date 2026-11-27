#!/usr/bin/env python3
"""Read-only GEO protocol and reciprocal cross-assay evidence review."""
import csv, gzip, hashlib, re, sys
from collections import Counter, defaultdict
from pathlib import Path

root=Path(sys.argv[1])
paths=dict(zip(("GSE71119","GSE75885","GSE71118"),map(Path,sys.argv[2:])))
assert len(sys.argv)==5
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
patterns={
    "UNTREATED_PRIMARY":r"\buntreated primary tumou?r[s]?\b",
    "PRIMARY_TUMOUR":r"\bprimary tumou?r[s]?\b",
    "BEFORE_TREATMENT":r"\b(?:before|prior to) (?:any |systemic |neo.adjuvant )?treatment\b",
    "RIBO_DEPLETION":r"\bribo.?zero\b|\bribosomal rna.deplet|\brrna.deplet|\bdepletion of ribosomal rna",
    "POLYA":r"\bpoly\s*[\(-]?\s*a\b|\bpolyadenylat|\boligo.?dt\b",
    "TRUSEQ":r"\btruseq\b",
}
def flag(text,key): return bool(re.search(patterns[key],text,re.I))
expected={"GSE71119":"249d5e0a299bf35a93802ca7e3a61e2aa8a36f100930300ac9c87400eae87858","GSE75885":"1c59dee75da58ddcff2bc2761f395e3ed27dd5fde83055cc11835dadcb04e6f4"}
geo={}
for series,path in paths.items():
    digest=sha(path)
    if series in expected: assert digest==expected[series], "Previously locked GEO source changed"
    geo[series]=parse(path,series)
    print("SOURCE",series,"SHA256",digest,"GSM_RECORDS",len(geo[series]))
exact={}
for series in expected:
    out=defaultdict(list)
    for gsm,r in geo[series].items():
        ids=relations(r)
        if len(ids)==1:
            s=next(iter(ids))
            if s in all_ids: out[s].append(r)
    assert all(len(v)==1 for v in out.values())
    exact[series]={s:v[0] for s,v in out.items()}
a,b=exact["GSE71119"],exact["GSE75885"]
assert (len(a),len(b),len(keep&a.keys()),len(keep&b.keys()),len(keep&a.keys()&b.keys()),len(keep&(a.keys()|b.keys())))==(136,112,129,106,97,138)
print("\n=== DIRECT RNA-SEQ GEO PROTOCOL EVIDENCE ===")
for series,rs in exact.items():
    selected=[r for s,r in rs.items() if s in keep]
    print("SERIES",series,"RETAINED",len(selected))
    for f in ("!Sample_extract_protocol_ch1","!Sample_treatment_protocol_ch1","!Sample_growth_protocol_ch1","!Sample_description","!Sample_source_name_ch1"):
        counts=Counter(field_text(r,f) for r in selected)
        print("FIELD",f,"NONEMPTY",sum(n for t,n in counts.items() if t),"DISTINCT_NONEMPTY",sum(bool(t) for t in counts))
        for k in patterns:
            print("PHRASE",k,sum(n for t,n in counts.items() if flag(t,k)))
        # Repeated public study protocol text only. No per-record descriptions or IDs.
        if f.endswith("protocol_ch1"):
            for t,n in sorted(counts.items(),key=lambda x:-x[1]):
                if not t: continue
                print("PROTOCOL_GROUP",n,"SHA256",hashlib.sha256(t.encode()).hexdigest())
                if n>=5 and not re.search(r"\b(?:GSM|SRX|SRR|SAMN)\d+\b|\bpatient\s*(?:id|number|#)\s*[:=]?\s*[A-Za-z0-9_-]+",t,re.I):
                    print("SHARED_PUBLIC_PROTOCOL",repr(t[:4000]))
                else: print("PROTOCOL_TEXT_LOCAL_REVIEW_ONLY")
    for k in ("material support","paired material support","cinsarc"):
        if any(k in traits(r) for r in selected):
            print("TRAIT",k,dict(Counter(trait(r,k) or "MISSING_OR_MULTIVALUE" for r in selected)))
print("\n=== OVERLAPPING GEO CLINICAL ANNOTATIONS ===")
both=keep&a.keys()&b.keys()
counts=Counter()
for s in both:
    av,bv=trait(a[s],"metastasis"),trait(b[s],"metastasis")
    counts[(av if av in {"yes","no"} else "OTHER",bv if bv in {"yes","no"} else "OTHER")]+=1
print("METASTASIS_CROSSTAB_GSE71119_GSE75885",dict(counts))
print("\n=== PAIRED MICROARRAY EVIDENCE ===")
micro=geo["GSE71118"]
titlemap=defaultdict(list)
for r in micro.values():
    if title(r):titlemap[title(r)].append(r)
status=Counter(); linked=[]; used=Counter()
for s in sorted(keep):
    r=a.get(s)
    if r is None: status["NO_GSE71119_RECORD"]+=1;continue
    token=trait(r,"paired microarray")
    if not token:status["NO_UNIQUE_PAIRED_MICROARRAY_TOKEN"]+=1;continue
    found=titlemap.get(token,[])
    if len(found)!=1:status["MICROARRAY_TITLE_NOT_UNIQUE_OR_ABSENT"]+=1;continue
    m=found[0]
    # Explicit reciprocal cross-assay link and exact native RNA specimen title.
    if not title(r) or title(r)!=seqid(native[s]):status["RNA_NATIVE_TITLE_MISMATCH"]+=1;continue
    if trait(m,"paired rna-seq")!=title(r):status["RECIPROCAL_LINK_MISSING_OR_MISMATCH"]+=1;continue
    status["RECIPROCAL_PAIR"]+=1;linked.append((s,r,m));used[token]+=1
assert sum(status.values())==140
print("PAIR_STATUS",dict(status))
print("REUSED_MICROARRAY_TITLE_GROUPS",sum(n>1 for n in used.values()))
for f in ("!Sample_treatment_protocol_ch1","!Sample_extract_protocol_ch1"):
    cnt=Counter(field_text(m,f) for _,_,m in linked)
    print("PAIRED_MICROARRAY_FIELD",f,"NONEMPTY",sum(n for t,n in cnt.items() if t),"DISTINCT_NONEMPTY",sum(bool(t) for t in cnt))
    for k in ("UNTREATED_PRIMARY","PRIMARY_TUMOUR","BEFORE_TREATMENT"):
        print("PHRASE",k,sum(n for t,n in cnt.items() if flag(t,k)))
    for t,n in sorted(cnt.items(),key=lambda x:-x[1]):
        if t and n>=5 and not re.search(r"\b(?:GSM|SRX|SRR|SAMN)\d+\b|\bpatient\s*(?:id|number|#)\s*[:=]?\s*[A-Za-z0-9_-]+",t,re.I):
            print("SHARED_PUBLIC_PAIRED_PROTOCOL",n,repr(t[:4000]))
print("PAIRED_MATERIAL_LABEL_CROSSTAB",dict(Counter((trait(r,"material support") or "MISSING",trait(m,"material support") or "MISSING") for _,r,m in linked)))
print("PAIRED_METASTASIS_CROSSTAB",dict(Counter((trait(r,"metastasis") if trait(r,"metastasis") in {"yes","no"} else "OTHER",trait(m,"metastasis") if trait(m,"metastasis") in {"yes","no"} else "OTHER") for _,r,m in linked)))
print("PAIRED_UNTREATED_PRIMARY_AND_UNIQUE_MICROARRAY_TITLE",sum(used[title(m)]==1 and any(flag(field_text(m,f),"UNTREATED_PRIMARY") for f in ("!Sample_treatment_protocol_ch1","!Sample_extract_protocol_ch1")) for _,r,m in linked))
print("V8_SRP057793_GEO_PROTOCOL_PAIRING_COMPLETE")
print("NO_FILES_MODIFIED_NO_CLINICAL_ADMISSIONS")
