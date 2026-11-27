#!/usr/bin/env python3
"""Build private unsigned evidence ledger; immutable inputs and explicit provenance."""
import csv, gzip, hashlib, io, json, re, sys, zipfile
from collections import Counter, defaultdict
from pathlib import Path

root=Path(sys.argv[1])
repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1"
s41=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def val(v):return str(v or "").strip()
def present(v):return val(v).casefold() not in {"","na","n/a","none","null","nan"}
def codes(v):
    s=val(v)
    # Only categorical uppercase evidence grades are displayed.
    return s if re.fullmatch(r"[A-Z][A-Z0-9_]{0,100}",s) else "EMPTY" if not present(s) else "OTHER_NONCODE_VALUE"
def groups_summary(groups):
    return {"groups":len(groups),"reused_groups":sum(len(v)>1 for v in groups.values()),
            "records_in_reused_groups":sum(len(v) for v in groups.values() if len(v)>1),
            "max_group_size":max(map(len,groups.values()),default=0)}

v5=rec/"RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"
assert sha(v5)=="a76a68465629dafba98f16f688a69e25fa3cfcb0282ab30add66206a85c62b36"
with zipfile.ZipFile(v5) as z:
    paths=[p for p in z.namelist() if p.endswith("/results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv")]
    assert len(paths)==1
    raw=z.read(paths[0])
assert hashlib.sha256(raw).hexdigest()=="06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49"
rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t"))
primary={r["library_id"]:r for r in rows}
assert len(rows)==len(primary)==2055 and len({r["study"] for r in rows})==77
assert all(r["FINAL_CLINICAL_APPROVAL"]=="NO_UNSIGNED" for r in rows)
candpath=s41/"03_candidate_integrated_evidence_3602.tsv"
assert sha(candpath)=="aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81"
crows=read(candpath)
c={(r["study"],r["library_id"]):r for r in crows}
assert len(c)==len(crows)==3602
assert all((r["study"],r["library_id"]) in c and c[(r["study"],r["library_id"])]["technical_hold"]=="NO" for r in rows)

# Reuse existing V6 precedence decisions; do not rerun their construction.
evpath=repo/"diagnostics/20261123_v6_evidence_precedence_dryrun/proposed_evidence_2055_UNSIGNED.tsv"
assert sha(evpath)=="08efe2c8f469edcbdfb98231dc87f705d801434e78dcbea85f57c4b94da88eb0"
evrows=read(evpath);ev={r["library_id"]:r for r in evrows}
assert len(ev)==len(evrows)==2055 and ev.keys()==primary.keys()
assert all(ev[k]["study"]==r["study"] for k,r in primary.items())
assert all(r["FINAL_APPROVAL"]=="NO_UNSIGNED" for r in evrows)
exclude={k for k,r in ev.items() if r["unsigned_evidence_bucket"]=="EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE"}
assert len(exclude)==22

assert len(sys.argv)==3
out=Path(sys.argv[2]).resolve()
assert out.parent== (repo/"private_data").resolve()
assert out.name=="v8_primary2055_integrated_unsigned_evidence_ledger_v1"
assert not out.exists(), "OUTPUT_EXISTS"
# Attach original V6 source evidence without reinterpreting its clinical grades.
def load_evidence(rel,digest,expected_ids):
    path=repo/rel
    assert sha(path)==digest,"SOURCE_HASH_MISMATCH "+rel
    data=read(path)
    assert data
    # Recognize only an explicit accession column, and prove exact membership.
    fields=[k for k in ("experiment_acc","library_id") if k in data[0]]
    assert fields,"EXPLICIT_EXPERIMENT_KEY_ABSENT "+repr(list(data[0]))
    keys=[{r[k] for r in data} for k in fields]
    assert all(x==expected_ids for x in keys),"EVIDENCE_MEMBERSHIP_MISMATCH "+rel
    key=fields[0]
    mapping={r[key]:r for r in data}
    assert len(data)==len(mapping)==len(expected_ids)
    assert all(r["study"]==primary[k]["study"] for k,r in mapping.items())
    return mapping
hist={k for k,r in ev.items() if r["unsigned_evidence_bucket"]=="HISTORICAL_EVIDENCE_REVIEW_PENDING"}
unattested={k for k in hist if ev[k]["historical_tissue_grade"]=="PRIMARY_TUMOUR_NOT_ATTESTED"}
assert len(hist)==1696 and len(unattested)==599
hrel="diagnostics/20261124_v6_historical1696_integrated_v1/1696_clinical_source_evidence_LOCAL_ONLY.tsv"
qrel="diagnostics/20261124_v6_599_integrated_review/599_accession_integrated_evidence_LOCAL_ONLY.tsv"
hmap=load_evidence(hrel,"d2cecf185f15dcd8c72ad8a2c4edc7e6481b5e0767c75245401f7db1df7762dc",hist)
qmap=load_evidence(qrel,"c77599f38ecfcf16be5f8a4648f0386421f65de215a9a3a3c13bcb965f82291a",unattested)
print("V6_EXACT_LOCAL_EVIDENCE_ATTACHMENTS",len(hmap),len(qmap))
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

# Reconstruct exact links for the new row-level ledger from already verified sources.
# This is a new evidence integration, not a repeat of the aggregate diagnostic.
cache=repo/"private_data/v8_srp057793_geo_source_v1"
sources={"GSE71119":("249d5e0a299bf35a93802ca7e3a61e2aa8a36f100930300ac9c87400eae87858",136),
         "GSE75885":("1c59dee75da58ddcff2bc2761f395e3ed27dd5fde83055cc11835dadcb04e6f4",117),
         "GSE71118":("ea11c2ff661b47da23ed416e3e188cb80a6c7171edbd5d5389d22f65f6a4cc8c",312)}
g={}
for series,(digest,n) in sources.items():
    path=cache/f"{series}_family.soft.gz"
    assert sha(path)==digest
    g[series]=parse(path,series)
    assert len(g[series])==n
sar={k for k,r in primary.items() if r["study"]=="SRP057793"}
assert len(sar)==140
nativepath=rec/"raw_studies/SRP057793/metadata/sra.sra.SRP057793.MD.gz"
with gzip.open(nativepath,"rt",encoding="utf-8-sig",newline="") as f:
    native_rows=list(csv.DictReader(f,delimiter="\t"))
native={r["experiment_acc"]:r for r in native_rows}
assert len(native)==len(native_rows)==149 and sar<=native.keys()
def exact(series):
    d={}
    for gsm,r in g[series].items():
        ids=relations(r)
        if len(ids)==1 and next(iter(ids)) in sar:
            k=next(iter(ids))
            assert k not in d
            d[k]=(gsm,r)
    return d
a,b=exact("GSE71119"),exact("GSE75885")
assert (len(a),len(b),len(a.keys()&b.keys()),len(a.keys()|b.keys()))==(129,106,97,138)
micro=defaultdict(list)
for gsm,r in g["GSE71118"].items():
    if title(r):micro[title(r)].append((gsm,r))
paired={}
for k,(gsm,r) in a.items():
    assert title(r)==seqid(native[k]) and trait(r,"material support")=="frozen tumor"
    # Verify the exact shared protocol hash observed in the executed review.
    assert hashlib.sha256(field_text(r,"!Sample_extract_protocol_ch1").encode()).hexdigest()=="a709cfdbb2d7723f169eed9abdab367aa6e20a4b079833acbdaf96d9e651969e"
    token=trait(r,"paired microarray")
    if not token:continue
    found=micro.get(token,[])
    if len(found)!=1:continue
    mgsm,m=found[0]
    if trait(m,"paired rna-seq")!=title(r):continue
    assert trait(m,"material support")=="frozen tumor"
    assert norm(field_text(m,"!Sample_treatment_protocol_ch1"))=="rna extracted from frozen tissue of the untreated primary tumor"
    paired[k]=(mgsm,m)
for k,(gsm,r) in b.items():
    assert hashlib.sha256(field_text(r,"!Sample_extract_protocol_ch1").encode()).hexdigest()=="32641e7629bdde972cba20b04a81d7da2eaa2252b70b9471627c1705d52f6100"
assert len(paired)==90 and len({gsm for gsm,_ in paired.values()})==90
discord={k for k,(_,m) in paired.items() if trait(a[k][1],"metastasis")!=trait(m,"metastasis")}
assert len(discord)==1
assert all(trait(a[k][1],"metastasis")==trait(b[k][1],"metastasis") for k in a.keys()&b.keys())
print("SARCOMA_EVIDENCE_EXACT_ATTACHMENTS",len(a.keys()|b.keys()),len(paired),"METASTASIS_DISCORDANCE",len(discord))
# Preserve V5 key groups and calculate actual post-proposal grouping separately.
allgroups=defaultdict(set);remaining=defaultdict(set)
for k,r in primary.items():
    assert present(r["patient_key"])
    group=(r["study"],r["patient_key"].strip())
    allgroups[group].add(k)
    if k not in exclude:remaining[group].add(k)
assert groups_summary(allgroups)=={"groups":1969,"reused_groups":54,"records_in_reused_groups":140,"max_group_size":6}
print("ADMIN_GROUPS",groups_summary(allgroups))
print("AFTER_PROPOSED_EXCLUSIONS_GROUPS_NOT_PATIENT_COUNTS",groups_summary(remaining))
# No historical or new evidence is converted to an approval.
# Original record fields kept in private JSON to preserve provenance; no data printed.
ledger=[]
for k,r in sorted(primary.items(),key=lambda kv:(kv[1]["study"],kv[0])):
    s=r["study"]; group=(s,r["patient_key"].strip())
    c_row=c[(s,k)]
    assert all(c_row[f]=="NO" for f in ("patient_independence_certified","primary_tumour_certified","pretreatment_certified","assay_chemistry_certified","atlas_independence_certified"))
    record={
        "study":s,"library_id_LOCAL_ONLY":k,
        "historical_patient_key_LOCAL_ONLY":group[1],
        "historical_patient_grade":r["patient_evidence_grade"],
        "historical_treatment_grade":r["treatment_evidence_grade"],
        "historical_tissue_grade":r["tissue_evidence_grade"],
        "v6_unsigned_evidence_bucket":ev[k]["unsigned_evidence_bucket"],
        "proposed_non_tumour_exclusion":"YES" if k in exclude else "NO",
        "admin_provisional_group_size":len(allgroups[group]),
        "remaining_provisional_group_size":len(remaining.get(group,set())),
        "exact_atlas_accession_screen":"NO_DETECTED_EXACT_HIT_V8",
        "sarcoma_polyA_protocol_evidence":"SOURCE_SUPPORTED_FROZEN_POLYA" if k in a or k in b else "NATIVE_LABEL_ONLY" if k in sar else "NOT_REVIEWED_BY_THIS_ADDENDUM",
        "sarcoma_untreated_primary_evidence":"RECIPROCAL_FROZEN_MICROARRAY_SOURCE_SUPPORTED" if k in paired else "UNRESOLVED_BY_THIS_ROUTE" if k in sar else "NOT_REVIEWED_BY_THIS_ADDENDUM",
        "sarcoma_metastasis_exception":"YES" if k in discord else "NO_DISCORDANCE_DETECTED_IN_PAIR" if k in paired else "NOT_EVALUABLE",
        "sarcoma_GEO_GSM_LOCAL_ONLY":"|".join(x[k][0] for x in (a,b) if k in x),
        "sarcoma_paired_microarray_GSM_LOCAL_ONLY":paired[k][0] if k in paired else "",
        "patient_independence_review":"UNSIGNED",
        "cross_study_patient_independence_review":"UNSIGNED",
        "primary_sampled_lesion_review":"UNSIGNED",
        "pretreatment_review":"UNSIGNED",
        "assay_chemistry_review":"UNSIGNED",
        "technical_release_review":"UNSIGNED",
        "final_approval":"NO_UNSIGNED",
        "v5_original_JSON_LOCAL_ONLY":json.dumps(r,sort_keys=True,ensure_ascii=False),
        "stage41_original_JSON_LOCAL_ONLY":json.dumps(c_row,sort_keys=True,ensure_ascii=False),
        "v6_precedence_JSON_LOCAL_ONLY":json.dumps(ev[k],sort_keys=True,ensure_ascii=False),
        "v6_1696_source_JSON_LOCAL_ONLY":json.dumps(hmap[k],sort_keys=True,ensure_ascii=False) if k in hmap else "",
        "v6_599_source_JSON_LOCAL_ONLY":json.dumps(qmap[k],sort_keys=True,ensure_ascii=False) if k in qmap else "",
    }
    ledger.append(record)
assert len(ledger)==2055
assert sum(r["proposed_non_tumour_exclusion"]=="YES" for r in ledger)==22
assert all(r["final_approval"]=="NO_UNSIGNED" for r in ledger)
study_rows=[]
for s in sorted({r["study"] for r in ledger}):
    rs=[r for r in ledger if r["study"]==s]
    rem=[r for r in rs if r["proposed_non_tumour_exclusion"]=="NO"]
    sg={k:v for k,v in remaining.items() if k[0]==s}
    z=groups_summary(sg)
    study_rows.append({"study":s,"admin_n":len(rs),"proposed_non_tumour_n":len(rs)-len(rem),
        "remaining_unsigned_n":len(rem),"remaining_provisional_keys":z["groups"],"remaining_reused_groups":z["reused_groups"],"remaining_records_in_reused_groups":z["records_in_reused_groups"],
        "v6_1696_evidence_n":sum(bool(r["v6_1696_source_JSON_LOCAL_ONLY"]) for r in rs),"v6_599_evidence_n":sum(bool(r["v6_599_source_JSON_LOCAL_ONLY"]) for r in rs),
        "new_sarcoma_protocol_support_n":sum(r["sarcoma_polyA_protocol_evidence"]=="SOURCE_SUPPORTED_FROZEN_POLYA" for r in rs),
        "new_sarcoma_untreated_primary_support_n":sum(r["sarcoma_untreated_primary_evidence"]=="RECIPROCAL_FROZEN_MICROARRAY_SOURCE_SUPPORTED" for r in rs),
        "clinical_admissions":0})
assert len(study_rows)==77 and sum(r["remaining_unsigned_n"] for r in study_rows)==2033
# All validation above precedes new file creation. Private output only.
out.mkdir(mode=0o700)
def write(path,records):
    with path.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(records)
    path.chmod(0o600)
lf=out/"2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv"
sf=out/"77_study_aggregate_UNSIGNED.tsv"
write(lf,ledger);write(sf,study_rows)
print("PRIVATE_LEDGER_SHA256",sha(lf))
print("STUDY_AGGREGATE_SHA256",sha(sf))
print("\t".join(study_rows[0]))
for r in study_rows:print("\t".join(str(v) for v in r.values()))
print("V8_PRIMARY2055_INTEGRATED_UNSIGNED_LEDGER_COMPLETE")
print("NO_FROZEN_INPUT_MODIFICATIONS_NO_MODEL_EXECUTION_CLINICAL_ADMISSIONS_0")
