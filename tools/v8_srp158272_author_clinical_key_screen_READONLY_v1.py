#!/usr/bin/env python3
"""Test exact source RR identifier correspondence with author clinical table.
No demographic matching, numeric-only matching, padding changes, or identity certification.
"""
import csv,gzip,hashlib,re,sys
from pathlib import Path
from collections import Counter,defaultdict
root=Path(sys.argv[1]); repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1"
study="SRP158272"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p,sep="\t"):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter=sep))
def norm(x):return str(x or "").strip().upper()
source=Path(sys.argv[2])
raw=source.read_bytes()
assert hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()=="786115ac47e7f073638d1e7909e8189a114a27f8"
print("AUTHOR_CLINICAL_TABLE_SHA256",sha(source))
clinical=read(source,",")
assert len(clinical)==558
author={norm(r["RR-Number"]):r for r in clinical}
assert len(author)==558 and all(re.fullmatch(r"RR-\d+",k) for k in author)
ledger=repo/"private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv"
assert sha(ledger)=="943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0"
keep={r["library_id_LOCAL_ONLY"] for r in read(ledger) if r["study"]==study}
cand=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/03_candidate_integrated_evidence_3602.tsv"
assert sha(cand)=="aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81"
cr=[r for r in read(cand) if r["study"]==study]
hold={r["library_id"] for r in cr if r["technical_hold"]=="YES"}
assert len(keep)==129 and len(hold)==5 and len(cr)==134 and not keep&hold
assert {r["library_id"] for r in cr}==keep|hold
native=rec/"raw_studies"/study/"metadata"/f"sra.sra.{study}.MD.gz"
assert sha(native)=="bfbdf336decb8295c204d17f0268740a7c5f00eeaebd66fddc50ddbc7b03e874"
nr=read(native); ni={r["experiment_acc"]:r for r in nr}
assert len(ni)==len(nr)==139 and keep|hold<=ni.keys()
other=set(ni)-keep-hold
assert len(other)==5
token=re.compile(r"(?<![A-Z0-9])RR-\d+(?![A-Z0-9])",re.I)
def fields(r):
    d={k:r.get(k,"") for k in ("sample_name","sample_title","external_id","library_name","sample_description","experiment_title","sample_acc")}
    for a in r.get("sample_attributes","").split("|"):
        if ";;" not in a:continue
        k,v=a.split(";;",1)
        if norm(k) in {"ISOLATE","READ1","READ2"}:
            d.setdefault("attribute_"+norm(k),"")
            d["attribute_"+norm(k)]+=" "+v
    return d
allmap={}
for scope,ids in (("RETAINED",keep),("HELD",hold),("SOURCE_ONLY",other)):
    print("SCOPE",scope,"N",len(ids))
    vals={e:fields(ni[e]) for e in ids}
    union=defaultdict(set)
    for f in sorted({f for d in vals.values() for f in d}):
        hits={e:set(token.findall(norm(d.get(f,"")))) for e,d in vals.items()}
        exact={e:x&author.keys() for e,x in hits.items()}
        for e,x in hits.items():union[e].update(x)
        print("FIELD",f,"RR_TOKEN_PRESENT",sum(bool(x) for x in hits.values()),
              "LITERAL_WHOLE_FIELD_MATCH",sum(norm(d.get(f,"")) in author for d in vals.values()),
              "AUTHOR_TOKEN_MATCH",sum(bool(x) for x in exact.values()),
              "MULTIPLE_RR_TOKENS",sum(len(x)>1 for x in hits.values()),
              "UNMATCHED_RR_TOKENS",sum(len(x-author.keys()) for x in hits.values()))
    match={e:next(iter(x)) for e,x in union.items() if len(x)==1 and x<=author.keys()}
    allmap.update(match)
    print("CONSISTENT_EXACT_AUTHOR_TOKEN",len(match),"DISTINCT_AUTHOR_KEYS",len(set(match.values())),
          "AMBIGUOUS_SOURCE_TOKENS",sum(len(x)>1 for x in union.values()))
    # Raw categorical code counts only; semantics/timing require source documentation.
    for f in ("NAC","Metastatic","DCISONLY","Other_Invasive","IDC+","ILC+","AD","BC"):
        print("MATCHED_AUTHOR_RAW_CODE",f,dict(Counter(str(author[k][f]).strip() if str(author[k][f]).strip() in {"0","1","-9","0.0","1.0","-9.0",""} else "OTHER" for k in match.values())))
print("ALL_SCOPE_EXACT_TOKEN_MATCHES",len(allmap),"DISTINCT_AUTHOR_KEYS",len(set(allmap.values())),
      "REUSED_AUTHOR_KEYS",sum(n>1 for n in Counter(allmap.values()).values()))
print("V8_SRP158272_AUTHOR_CLINICAL_KEY_SCREEN_COMPLETE")
print("NO_DEMOGRAPHIC_LINKAGE_NO_FILES_MODIFIED_NO_CLINICAL_ADMISSIONS")
