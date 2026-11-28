#!/usr/bin/env python3
"""Prepare a private technical-input manifest for 461 joint clinical evidence leads.
No projection, no expression normalization, no clinical admission. Raw G026 header only.
"""
import csv,gzip,hashlib,re,sys,statistics,math
from pathlib import Path
from collections import Counter,defaultdict
assert len(sys.argv)==3
root=Path(sys.argv[1]);repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1"
out=Path(sys.argv[2]).resolve()
assert out.parent==(repo/"private_data").resolve() and out.name=="v8_joint461_technical_input_manifest_v1" and not out.exists()
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def locked(p,h):
    assert sha(p)==h,"INPUT_HASH_MISMATCH "+str(p)
    return read(p)
def summary(xs):
    a=[]
    for x in xs:
        try:
            v=float(x)
            if math.isfinite(v):a.append(v)
        except (ValueError,TypeError):pass
    return {"finite":len(a),"missing":len(xs)-len(a),"min":min(a) if a else None,"median":statistics.median(a) if a else None,"max":max(a) if a else None}
private=repo/"private_data"
rows=locked(private/"v8_primary2055_integrated_unsigned_evidence_ledger_v1/2055_integrated_unsigned_evidence_LOCAL_ONLY.tsv","943542b682b9aaaf46d31576e18b3de673dbca9d35b5e99c8f11379b028b35c0")
direct=locked(private/"v8_two_study_geo_clinical_evidence_v1/37_direct_GEO_evidence_LOCAL_ONLY.tsv","1cd5a57589377fa91d2a0bc5a22bf47b742ce6387003fbfb0313ffdab45f209f")
fuscc=locked(private/"v8_fuscc_author_annotation_addendum_v1/360_fuscc_author_annotation_evidence_LOCAL_ONLY.tsv","1cae682fe3ee006df02e22654f558c0ee2e28ab66fd5a422c67374609e1d1dd5")
assert len(rows)==len({r["library_id_LOCAL_ONLY"] for r in rows})==2055
assert len(direct)==37 and all(r["final_clinical_approval"]=="NO_UNSIGNED" for r in direct)
assert all(r["source_link"]=="UNIQUE_EXACT_SRX_RELATION" and r["tissue_exact_source_support"]==r["pretreatment_exact_source_support"]==r["biosample_relation_agreement"]=="YES" and r["atlas_GSM_hit"]=="NO" for r in direct)
D={r["library_id_LOCAL_ONLY"]:r for r in direct};assert len(D)==37
F={r["experiment_acc_LOCAL_ONLY"]:r for r in fuscc if r["membership"]=="RETAINED334"}
assert len(F)==334
selected={}
for r in rows:
    k=r["library_id_LOCAL_ONLY"];s=r["study"]
    assert r["final_approval"]=="NO_UNSIGNED"
    if r["proposed_non_tumour_exclusion"]!="NO":continue
    if k in F:
        assert s=="SRP157974" and r["historical_tissue_grade"]=="PRIMARY_TUMOUR_STUDY_SOURCE_SUPPORTED" and r["historical_treatment_grade"]=="COHORT_PRETREATMENT_SOURCE_SUPPORTED"
        evidence="COHORT_PRIMARY_PRETREATMENT_SOURCE_SUPPORT"
    elif s=="SRP057793" and r["sarcoma_untreated_primary_evidence"]=="RECIPROCAL_FROZEN_MICROARRAY_SOURCE_SUPPORTED":
        evidence="RECIPROCAL_PAIRED_UNTREATED_PRIMARY_SUPPORT"
    elif k in D:
        assert s==D[k]["study"]
        evidence="DIRECT_GEO_TISSUE_AND_PRETREATMENT_SUPPORT"
    else:continue
    assert k not in selected
    selected[k]=(r,evidence)
expected={"SRP157974":334,"SRP057793":90,"SRP002628":14,"SRP094795":23}
assert Counter(r["study"] for r,_ in selected.values())==expected and len(selected)==461
qcl=locked(rec/"stage2_output/library_qc_all.tsv","3a91c52a2a8cf00a6de7e88fd0c1c9f255152c70bc7cd6d3171d92a8e7279b")
qcr=locked(rec/"stage2_output/run_qc_all.tsv","a1c8fbacd32dca7e7f429ec2dc9a21d5a9b866dd7aefa6f6e2570de8c4445159")
hq=locked(rec/"cohort_readiness_v1/source_snapshot/05_run_hold_queue_409.tsv","d275c77591e7b93582708746e45d69ae817c8e2d13b3bb5b7f63e2e56fa56f0c")
L=defaultdict(list);R=defaultdict(list);H=defaultdict(list);run_owner=defaultdict(set)
for r in qcl:L[(r["study"],r["library_id"])].append(r)
for r in qcr:
    key=(r["study"],r["library_id"]);R[key].append(r)
    run_owner[r["run_id"].strip().upper()].add(key)
for r in hq:H[(r["study"],r["library_id"])].append(r)
mappath=rec/"stage5_output/00_gene_row_mapping.tsv"
assert sha(mappath)=="68be266dd6bd30bfa79a2bd8ac79f0e49c6f5d0f19d2d083e2ca57850ef17391"
mp=read(mappath)
assert len(mp)==63856 and [int(r["row_zero_based"]) for r in mp]==list(range(63856))
print("LOCKED_461_MEMBERSHIP_PASS",dict(expected),"G026_LOCKED_ROW_MAP",len(mp))
RUN=re.compile(r"(?:SRR|ERR|DRR)\d+",re.I)
manifest=[];agg=[]
for s,n in expected.items():
    ids=sorted(k for k,(r,_) in selected.items() if r["study"]==s)
    nativepath=rec/"raw_studies"/s/"metadata"/f"sra.sra.{s}.MD.gz"
    native=defaultdict(list)
    if nativepath.is_file():
        print("NATIVE_SOURCE",s,"SHA256",sha(nativepath))
        for r in read(nativepath):native[r["experiment_acc"]].append(r)
    raw=rec/"raw_studies"/s/"gene_sums"/f"sra.gene_sums.{s}.G026.gz"
    headers=[];header_status="RAW_FILE_MISSING";header_sha=""
    if raw.is_file():
        header_status="RAW_HEADER_NOT_FOUND"
        try:
            with gzip.open(raw,"rt",encoding="utf-8-sig",errors="strict",newline="") as f:
                for i in range(30):
                    line=f.readline()
                    if not line:break
                    v=line.rstrip("\r\n").split("\t")
                    if len(v)>1 and all(RUN.fullmatch(x.strip().strip('"')) for x in v[1:]):
                        headers=[x.strip().strip('"').upper() for x in v[1:]]
                        header_sha=hashlib.sha256(line.encode()).hexdigest()
                        header_status="UNIQUE_RUN_HEADER" if len(headers)==len(set(headers)) else "DUPLICATE_RUN_HEADER"
                        break
        except (OSError,UnicodeError,EOFError):header_status="RAW_HEADER_READ_ERROR"
    print("RAW_G026",s,"STATUS",header_status,"RUN_COLUMNS",len(headers),"HEADER_SHA256",header_sha or "NA","BYTES",raw.stat().st_size if raw.is_file() else 0)
    headerset=set(headers);study_rows=[];all_run_rows=[]
    for k in ids:
        r,evidence=selected[k];key=(s,k)
        lq=L[key];rq=R[key];hh=H[key];nr=native[k]
        all_run_rows.extend(rq)
        runs=[q["run_id"].strip().upper() for q in rq]
        valid_runs=bool(runs) and len(runs)==len(set(runs)) and all(RUN.fullmatch(x) and run_owner[x]=={key} for x in runs)
        expected_n=verified_n=None
        if len(lq)==1:
            try:
                expected_n=int(lq[0]["n_runs_expected"]);verified_n=int(lq[0]["n_runs_expression_verified"])
            except (ValueError,TypeError):pass
        count_ok=valid_runs and len(lq)==1 and expected_n==verified_n==len(runs) and lq[0]["run_count_status"]=="PRESENT_ALL"
        holds_clear=not hh
        raw_ok=header_status=="UNIQUE_RUN_HEADER" and valid_runs and set(runs)<=headerset
        status="INPUT_READY_FOR_FULL_G026_QC_NOT_RELEASED" if count_ok and holds_clear and raw_ok and len(nr)==1 else "HOLD_INPUT_CONTRACT_UNRESOLVED"
        row={"study":s,"library_id_LOCAL_ONLY":k,"joint_clinical_evidence":evidence,
             "historical_patient_grade":r["historical_patient_grade"],
             "patient_independence_certified":"NO","final_clinical_approval":"NO_UNSIGNED",
             "stage2_expected_runs":expected_n,"stage2_verified_runs":verified_n,
             "stage2_run_ids_LOCAL_ONLY":"|".join(sorted(set(runs))),
             "run_membership_and_count_ok":str(bool(count_ok)).upper(),
             "historical_hold_queue_rows":len(hh),"native_exact_experiment_rows":len(nr),
             "native_library_selection":nr[0]["library_selection"] if len(nr)==1 else "",
             "native_platform":nr[0]["platform_model"] if len(nr)==1 else "",
             "raw_g026_relative_path":str(raw.relative_to(root)),
             "raw_header_status":header_status,"raw_header_sha256":header_sha,
             "raw_expected_runs_present":str(bool(raw_ok)).upper(),
             "input_status":status}
        study_rows.append(row);manifest.append(row)
    counts=Counter(r["input_status"] for r in study_rows)
    print("STUDY",s,"N",n,"INPUT_STATUS",dict(counts))
    print("LIBRARY_SELECTION",dict(Counter(r["native_library_selection"] for r in study_rows)))
    print("PLATFORM",dict(Counter(r["native_platform"] for r in study_rows)))
    print("RUN_MULTIPLICITY",dict(Counter(r["stage2_expected_runs"] for r in study_rows)))
    print("RUN_QC_ROWS",len(all_run_rows),"GENE_DETECTION",summary([r["fraction_genes_detected"] for r in all_run_rows]),"GENE_BODY_BP",summary([r["gene_body_coverage_bp"] for r in all_run_rows]))
    print("RUN_BASIC_EXPRESSION_FLAGS",dict(Counter(r["basic_expression_flag"] for r in all_run_rows)))
    agg.append({"study":s,"n_joint_leads":n,"input_ready_not_released":counts["INPUT_READY_FOR_FULL_G026_QC_NOT_RELEASED"],"input_hold":counts["HOLD_INPUT_CONTRACT_UNRESOLVED"],"raw_run_columns":len(headers),"raw_header_status":header_status,"clinical_admissions":0})
assert len(manifest)==461 and len({r["library_id_LOCAL_ONLY"] for r in manifest})==461
out.mkdir(mode=0o700)
def write(name,data):
    f=out/name
    with f.open("x",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(data[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(data)
    f.chmod(0o600)
    print("OUTPUT_SHA256",name,sha(f))
write("461_technical_input_manifest_LOCAL_ONLY.tsv",manifest)
write("4_study_technical_input_aggregate_UNSIGNED.tsv",agg)
print("V8_JOINT461_TECHNICAL_INPUT_MANIFEST_COMPLETE")
print("HEADER_ONLY_NOT_FULL_EXPRESSION_QC_NO_MODEL_EXECUTION_NO_CLINICAL_ADMISSIONS")
