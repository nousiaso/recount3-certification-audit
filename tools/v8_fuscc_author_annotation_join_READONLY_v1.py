#!/usr/bin/env python3
"""Join pinned author annotations to the frozen FUSCC membership; aggregate stdout only.
Author sample IDs corroborate RNA specimen membership, not independent person identity.
"""
import csv,gzip,hashlib,json,re,sys,math,statistics
from pathlib import Path
from collections import Counter,defaultdict
root=Path(sys.argv[1])
repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1"
cache=Path(sys.argv[2])
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()
def read(p,sep="\t"):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter=sep))
def norm(x):return str(x or "").strip().casefold()
def numeric(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None
def describe(v):
    a=[x for x in v if x is not None]
    return {"n_finite":len(a),"n_missing":len(v)-len(a),"min":min(a) if a else None,
            "median":statistics.median(a) if a else None,"max":max(a) if a else None}
assert len(sys.argv)==3
table=cache/"Supplementary_Table_5.csv"
batch=cache/"RNAseq_batch_info.csv"
assert blob(table)=="9d8a51caee63a7ac3fb82c3e4bfa9e82f645641a"
assert blob(batch)=="41e705a589bc6e04572dc32c89b1b905259b632a"
print("AUTHOR_REPOSITORY_COMMIT","3857ad51af9ce4f13ff445763cdbc8cd82ea355f")
print("AUTHOR_SOURCE_SHA256",table.name,sha(table))
print("AUTHOR_SOURCE_SHA256",batch.name,sha(batch))
t=read(table,","); b=read(batch,",")
assert len(t)==360 and len(b)==448
assert set(t[0])=={"sample_id","subtype_old","subtype_new","tumor_purity"}
assert set(b[0])=={"sample_id","batch"}
token=re.compile(r"(?<![A-Za-z0-9])FUSCCTNBC\d+(?![A-Za-z0-9])",re.I)
def base(v):
    hits={x.upper() for x in token.findall(str(v or ""))}
    return next(iter(hits)) if len(hits)==1 else ""
def author_label(v):
    m=re.fullmatch(r"(FUSCCTNBC\d+)([._]PT)?([._]REP)?",str(v).strip(),re.I)
    assert m,"UNRECOGNIZED_AUTHOR_SAMPLE_LABEL"
    return (m[1].upper(),"NORMAL" if m[2] else "TUMOUR",bool(m[3]))
T={}
for r in t:
    k,kind,rep=author_label(r["sample_id"])
    assert kind=="TUMOUR" and not rep and k not in T
    T[k]=r
B={}
for r in b:
    k,kind,rep=author_label(r["sample_id"])
    assert (k,kind) not in B
    assert r["batch"] in {"batch1","batch2","batch3"}
    B[(k,kind)]=(r,rep)
assert len(T)==360 and sum(k[1]=="TUMOUR" for k in B)==360
assert {k[0] for k in B if k[1]=="TUMOUR"}==set(T)
assert sum(k[1]=="NORMAL" for k in B)==88
assert {k[0] for k in B if k[1]=="NORMAL"}<=set(T)

ledger=repo/"private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv"
assert sha(ledger)=="943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0"
lr=[r for r in read(ledger) if r["study"]=="SRP157974"]
keep={r["library_id_LOCAL_ONLY"] for r in lr}
assert len(lr)==len(keep)==334 and all(r["final_approval"]=="NO_UNSIGNED" for r in lr)
s41=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
c=s41/"03_candidate_integrated_evidence_3602.tsv"
assert sha(c)=="aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81"
cr=[r for r in read(c) if r["study"]=="SRP157974"]
hold={r["library_id"] for r in cr if r["technical_hold"]=="YES"}
assert len(cr)==360 and len(hold)==26 and not keep&hold
assert keep=={r["library_id"] for r in cr if r["technical_hold"]=="NO"}
native=rec/"raw_studies/SRP157974/metadata/sra.sra.SRP157974.MD.gz"
assert sha(native)=="c2e1c7b5b1169e13059447e14bdc71ed3ffd4166dd90131b2791360ebab28155"
nr=read(native); assert len(nr)==len({r["experiment_acc"] for r in nr})==448
assert keep|hold<={r["experiment_acc"] for r in nr}
mapping={}
issues=Counter()
for r in nr:
    e=r["experiment_acc"]
    k=base(r["library_name"])
    if not k:issues["LIBRARY_BASE_MISSING_OR_AMBIGUOUS"]+=1;continue
    # Patient base is a source RNA specimen token; verify selected tumour isolate independently.
    if e in keep|hold:
        iso=[v for a in r["sample_attributes"].split("|") if ";;" in a for key,v in [a.split(";;",1)] if norm(key)=="isolate"]
        if len(iso)!=1 or base(iso[0])!=k:issues["TUMOUR_ISOLATE_BASE_DISAGREES"]+=1;continue
    name=norm(r["library_name"])
    # Validate PT marker against the existing exact 360/88 source partition.
    pt=bool(re.search(r"(?:[._])pt(?:[._]|$)",name))
    kind="TUMOUR" if e in keep|hold else "NORMAL"
    if pt!=(kind=="NORMAL"):issues["PT_MARKER_PARTITION_CONFLICT"]+=1;continue
    key=(k,kind)
    if key in mapping:issues["DUPLICATE_BASE_MATERIAL"]+=1;continue
    rep=bool(re.search(r"(?:[._])rep(?:[._]|$)",name))
    mapping[key]=(r,rep)
print("NATIVE_BASE_MATERIAL_ISSUES",dict(issues))
print("NATIVE_AUTHOR_BASE_MATERIAL_MATCHES",len(mapping.keys()&B.keys()),"SOURCE_ONLY",len(mapping.keys()-B.keys()),"AUTHOR_ONLY",len(B.keys()-mapping.keys()))
assert not issues and len(mapping)==448 and set(mapping)==set(B),"SOURCE_PARTITION_NOT_EXACT"
print("EXACT_AUTHOR_TUMOUR_BASE_MATCHES",len({k[0] for k in mapping if k[1]=="TUMOUR"}&T.keys()))
assert {k[0] for k in mapping if k[1]=="TUMOUR"}==set(T)
print("REP_SUFFIX_AGREEMENT",sum(rep==B[k][1] for k,(_,rep) in mapping.items()),"REP_SUFFIX_DISAGREEMENT",sum(rep!=B[k][1] for k,(_,rep) in mapping.items()))
# Base+material matches are explicit; mismatched REP suffix is not silently declared full library-name identity.
for label,ids in (("RETAINED334",keep),("HELD26",hold)):
    keys=[k for k,(r,_) in mapping.items() if r["experiment_acc"] in ids]
    assert len(keys)==len(ids) and all(k[1]=="TUMOUR" for k in keys)
    print("SCOPE",label,"N",len(keys))
    print("AUTHOR_TUMOUR_SUBTYPE_OLD",dict(Counter(T[k[0]]["subtype_old"] for k in keys)))
    print("AUTHOR_TUMOUR_SUBTYPE_NEW",dict(Counter(T[k[0]]["subtype_new"] for k in keys)))
    print("AUTHOR_SUBTYPE_CHANGED",sum(T[k[0]]["subtype_old"]!=T[k[0]]["subtype_new"] for k in keys))
    print("AUTHOR_TUMOUR_PURITY_REPORTED",describe([numeric(T[k[0]]["tumor_purity"]) for k in keys]))
    print("AUTHOR_BATCH",dict(Counter(B[k][0]["batch"] for k in keys)))
    print("AUTHOR_REP_SUFFIX_DISAGREEMENT",sum(mapping[k][1]!=B[k][1] for k in keys))
    print("NATIVE_SELECTION_BY_AUTHOR_BATCH",dict(sorted(Counter((mapping[k][0]["library_selection"],B[k][0]["batch"]) for k in keys).items())))
    print("AUTHOR_PURITY_BY_BATCH")
    for ba in ("batch1","batch2","batch3"):
        print(ba,describe([numeric(T[k[0]]["tumor_purity"]) for k in keys if B[k][0]["batch"]==ba]))
print("V8_FUSCC_AUTHOR_ANNOTATION_JOIN_COMPLETE")
print("NO_PATIENT_CERTIFICATION_NO_EXPRESSION_PROCESSING_NO_CLINICAL_ADMISSIONS")
