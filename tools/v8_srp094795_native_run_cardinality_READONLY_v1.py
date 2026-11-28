#!/usr/bin/env python3
"""Inspect the 22 cervical input holds. Aggregate evidence only, no status changes."""
import csv,gzip,hashlib,re,sys
from pathlib import Path
from collections import Counter,defaultdict
root=Path(sys.argv[1]);repo=root/"RECOUNT3_CERTIFICATION_GIT"
rec=root/"26_recount3_external_primary_v1";study="SRP094795"
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def verify_log_hash(path,log,label):
    matches=re.findall(r"^"+re.escape(label)+r" ([0-9a-f]{64})$",log,re.M)
    assert len(matches)==1,label
    assert sha(path)==matches[0],str(path)
audit=repo/"diagnostics/20261128_v8_joint461_technical_input_manifest_v1_0_1"
log=(audit/"stdout.txt").read_text()
assert "V8_JOINT461_TECHNICAL_INPUT_MANIFEST_V1_0_1_COMPLETE" in log.splitlines()
manifest=repo/"private_data/v8_joint461_technical_input_manifest_v1_0_1/461_technical_input_manifest_LOCAL_ONLY.tsv"
verify_log_hash(manifest,log,"OUTPUT_SHA256 "+manifest.name)
rows=read(manifest);assert len(rows)==461
m=[r for r in rows if r["study"]==study];assert len(m)==23
assert Counter(r["input_status"] for r in m)=={"HOLD_INPUT_CONTRACT_UNRESOLVED":22,"INPUT_READY_FOR_FULL_G026_QC_NOT_RELEASED":1}
assert all(r["final_clinical_approval"]=="NO_UNSIGNED" for r in m)
native=rec/"raw_studies"/study/"metadata"/f"sra.sra.{study}.MD.gz"
verify_log_hash(native,log,"NATIVE_SOURCE "+study+" SHA256")
nr=read(native);assert nr
qcr=rec/"stage2_output/run_qc_all.tsv"
verify_log_hash(qcr,log,"INPUT_SHA256_PASS "+qcr.name)
rq=read(qcr)
N=defaultdict(list);Q=defaultdict(list)
for r in nr:N[r["experiment_acc"]].append(r)
for r in rq:
    if r["study"]==study:Q[r["library_id"]].append(r)
print("NATIVE_ROWS",len(nr),"DISTINCT_EXPERIMENTS",len(N),"NATIVE_ROWS_PER_EXPERIMENT",dict(Counter(map(len,N.values()))))
print("NATIVE_COLUMNS",list(nr[0]))
RUN=re.compile(r"(?<![A-Za-z0-9])(?:SRR|ERR|DRR)\d+(?![A-Za-z0-9])",re.I)
def tokens(x):return {v.upper() for v in RUN.findall(str(x or ""))}
def attrs(row):
    result=defaultdict(set)
    for seg in row.get("sample_attributes","").split("|"):
        if ";;" in seg:
            k,v=seg.split(";;",1)
            result[k.strip().casefold()].add(v.strip())
    return dict(result)
held=[r for r in m if r["input_status"]=="HOLD_INPUT_CONTRACT_UNRESOLVED"]
for scope,records in (("ALL23",m),("HELD22",held)):
    print("SCOPE",scope,"N",len(records))
    for f in ("run_membership_and_count_ok","historical_hold_queue_rows","native_exact_experiment_rows","raw_expected_runs_present","stage2_expected_runs","stage2_verified_runs"):
        print("MANIFEST",f,dict(Counter(r[f] for r in records)))
    pairs=[(r,N[r["library_id_LOCAL_ONLY"]]) for r in records]
    assert all(len(a)==int(r["native_exact_experiment_rows"]) for r,a in pairs)
    print("STAGE2_RUN_SET_AGREEMENT",sum(set(r["stage2_run_ids_LOCAL_ONLY"].split("|"))=={q["run_id"].strip().upper() for q in Q[r["library_id_LOCAL_ONLY"]]} for r,_ in pairs))
    print("NATIVE_ROW_COUNT_EQUALS_EXPECTED_RUN_COUNT",sum(len(a)==int(r["stage2_expected_runs"]) for r,a in pairs))
    # Same specimen and library facts must be complete and concordant.
    for f in nr[0]:
        v=[[str(x.get(f,"") or "").strip() for x in a] for _,a in pairs]
        print("FIELD_CARDINALITY",f,"ALL_ROWS_POPULATED",sum(bool(x) and all(x) for x in v),"CONSTANT_NONEMPTY",sum(bool(x) and all(x) and len(set(x))==1 for x in v),"VARIES_WITHIN_EXPERIMENT",sum(len(set(x))>1 for x in v))
    ap=[[attrs(x) for x in a] for _,a in pairs]
    keys=sorted({k for a in ap for x in a for k in x})
    for k in keys:
        v=[[frozenset(x.get(k,set())) for x in a] for a in ap]
        print("ATTRIBUTE_CARDINALITY",k,"ALL_ROWS_POPULATED",sum(bool(x) and all(x) for x in v),"CONSTANT_NONEMPTY",sum(bool(x) and all(x) and len(set(x))==1 for x in v),"VARIES_WITHIN_EXPERIMENT",sum(len(set(x))>1 for x in v))
    # Test actual accession tokens; do not assume a native field is a run accession.
    for f in nr[0]:
        observed=[[tokens(x.get(f,"")) for x in a] for _,a in pairs]
        if not any(any(x) for x in observed):continue
        expected=[set(r["stage2_run_ids_LOCAL_ONLY"].split("|")) for r,_ in pairs]
        unions=[set().union(*x) for x in observed]
        exact=sum(u==e for u,e in zip(unions,expected))
        per_row=sum(u==e and all(len(z)==1 for z in v) and len(v)==len(e) for u,e,v in zip(unions,expected,observed))
        print("RUN_TOKEN_FIELD",f,"ANY_TOKEN",sum(bool(u) for u in unions),"EXACT_EXPECTED_RUN_SET",exact,"ONE_DISTINCT_EXPECTED_RUN_PER_ROW",per_row,"UNEXPECTED_TOKEN_GROUPS",sum(bool(u-e) for u,e in zip(unions,expected)))
print("V8_SRP094795_NATIVE_RUN_CARDINALITY_COMPLETE")
print("NO_FILES_MODIFIED_NO_HOLDS_RELEASED_NO_CLINICAL_ADMISSIONS")
