#!/usr/bin/env python3
"""New direct GEO source evidence for two prioritized studies; private addendum, no certification."""
import csv,gzip,hashlib,re,sys
from collections import Counter,defaultdict
from pathlib import Path
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

assert len(sys.argv)==4
root=Path(sys.argv[1]);repo=root/"RECOUNT3_CERTIFICATION_GIT"
cache=Path(sys.argv[2]);out=Path(sys.argv[3]).resolve()
assert out.parent==(repo/"private_data").resolve()
assert out.name=="v8_two_study_geo_clinical_evidence_v1" and not out.exists()
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))
ledger=repo/"private_data/v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv"
assert sha(ledger)=="943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0"
rows=read(ledger)
assert len(rows)==2055 and all(r["final_approval"]=="NO_UNSIGNED" for r in rows)
full=root/"41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/01_full_source_experiment_identity_5221.tsv"
assert sha(full)=="95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa"
fr=read(full)
atlaspath=root/"14_projection_use_analysis_ready_atlas_v1_0_3/03_projection_analysis_ready_core.tsv.gz"
assert sha(atlaspath)=="7f6852e829f7aaa11a26f86b44247c9c8c0840958494b7766db1ad6abd180092"
atlas={r["sample_id"].strip().upper() for r in read(atlaspath)}
assert len(atlas)==9465
config={"SRP002628":("GSE22260",14,"prostate cancer tissue","prostate adenocarcinoma","tissue was extracted from patients before treatment"),
        "SRP094795":("GSE91065",23,"uterine cervix tumor","uterine cervix tumor","all biopsies were collected before patient treatment and stored in rnalater (qiagen) at -80°c until mrna isolation.")}
records=[]
for study,(series,n,tissue_expected,source_expected,timing_expected) in config.items():
    keep={r["library_id_LOCAL_ONLY"] for r in rows if r["study"]==study and r["proposed_non_tumour_exclusion"]=="NO"}
    assert len(keep)==n
    source={r["experiment_acc"]:r for r in fr if r["study"]==study}
    assert keep<=source.keys()
    path=cache/f"{series}_family.soft.gz"
    digest=sha(path); samples=parse(path,series)
    print("SOURCE",study,series,"SHA256",digest,"GSM_RECORDS",len(samples),"FULL_SOURCE_EXPERIMENTS",len(source))
    by=defaultdict(list)
    for gsm,r in samples.items():
        sr=relations(r)
        for e in sr&keep:by[e].append((gsm,r,len(sr)))
    stats=Counter()
    for e in sorted(keep):
        hits=by.get(e,[])
        record={"study":study,"library_id_LOCAL_ONLY":e,"GEO_series":series,"GEO_source_SHA256":digest,
                "GEO_GSM_LOCAL_ONLY":"","source_link":"UNRESOLVED","tissue_exact_source_support":"NO","pretreatment_exact_source_support":"NO",
                "assay_protocol_evidence":"UNRESOLVED","atlas_GSM_hit":"NOT_EVALUABLE","biosample_relation_agreement":"NOT_EVALUABLE",
                "patient_independence_certified":"NO","final_clinical_approval":"NO_UNSIGNED"}
        if len(hits)==1 and hits[0][2]==1:
            gsm,r,_=hits[0]
            t=trait(r,"tissue")
            sn=one(r.get("!Sample_source_name_ch1",[]))
            treatment=norm(field_text(r,"!Sample_treatment_protocol_ch1"))
            extract=norm(field_text(r,"!Sample_extract_protocol_ch1"))
            bios=set(re.findall(r"\b(?:SAMN|SAMEA|SAMD)\d+\b",field_text(r,"!Sample_relation"),re.I))
            native_bios=set(re.findall(r"\b(?:SAMN|SAMEA|SAMD)\d+\b",source[e]["biosamples_all_sources"],re.I))
            agreement=bool(bios) and {x.upper() for x in bios}=={x.upper() for x in native_bios}
            tissue_ok=t==tissue_expected and sn==source_expected
            timing_ok=treatment.rstrip(".")==timing_expected.rstrip(".") 
            record.update({"GEO_GSM_LOCAL_ONLY":gsm,"source_link":"UNIQUE_EXACT_SRX_RELATION",
                           "tissue_exact_source_support":"YES" if tissue_ok else "NO",
                           "pretreatment_exact_source_support":"YES" if timing_ok else "NO",
                           "assay_protocol_evidence":"OLIGODT_OR_POLYA_EXPLICIT" if re.search(r"oligo\s*\(?d\s*t\)?|poly\s*a\s*selected|polyadenylat",extract,re.I) else "TRUSEQ_RNA_KIT_TEXT" if "truseq rna sample prep kit" in extract else "UNRESOLVED",
                           "atlas_GSM_hit":"YES" if gsm in atlas else "NO",
                           "biosample_relation_agreement":"YES" if agreement else "NO"})
            stats["unique_exact_SRX"]+=1
            stats["tissue_exact"]+=tissue_ok
            stats["pretreatment_exact"]+=timing_ok
            stats["biosample_agreement"]+=agreement
            stats["atlas_GSM_hits"]+=gsm in atlas
            stats["joint_tissue_timing_biosample_no_GSM_hit"]+=tissue_ok and timing_ok and agreement and gsm not in atlas
            stats["assay_"+record["assay_protocol_evidence"]]+=1
        elif hits:stats["ambiguous_or_multi_SRX_relation"]+=1
        else:stats["no_exact_SRX_relation"]+=1
        records.append(record)
    print("RETAINED_EVIDENCE",study,"N",n,dict(stats))
    # Do not print any patient or specimen values; count field coverage for later review.
    matched=[v[0][1] for v in by.values() if len(v)==1 and v[0][2]==1]
    print("SOURCE_TITLE_DISTINCT",len({title(r) for r in matched if title(r)}))
    print("EXPLICIT_PATIENT_FIELDS", {k:sum(bool(trait(r,k)) for r in matched) for k in ("patient","patient id","patient_id","subject","subject id","donor")})
    print("TREATMENT_PROTOCOL_HASH_COUNTS",dict(Counter(hashlib.sha256(field_text(r,"!Sample_treatment_protocol_ch1").encode()).hexdigest() for r in matched)))
assert len(records)==37 and len({r["library_id_LOCAL_ONLY"] for r in records})==37
out.mkdir(mode=0o700)
target=out/"37_direct_GEO_evidence_LOCAL_ONLY.tsv"
with target.open("x",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t",lineterminator="\n")
    w.writeheader();w.writerows(records)
target.chmod(0o600)
print("PRIVATE_ADDENDUM_SHA256",sha(target))
print("V8_TWO_STUDY_GEO_EVIDENCE_COMPLETE")
print("NO_FROZEN_INPUT_MODIFICATIONS_NO_MODEL_EXECUTION_NO_CLINICAL_ADMISSIONS")
