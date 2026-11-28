#!/usr/bin/env python3
"""Create the approved, separate FUSCC331 research stratum; no fitting."""
import csv,hashlib,json,math,os,re,sys,tempfile,subprocess
from collections import Counter
from pathlib import Path
from datetime import datetime,timezone
assert len(sys.argv)==3
os.umask(0o077)
root=Path(sys.argv[1]).resolve();repo=root/"RECOUNT3_CERTIFICATION_GIT"
out=Path(sys.argv[2]).resolve()
assert out.parent==(repo/"private_data").resolve() and out.name=="v8_fuscc331_clinical_lock_v1" and not out.exists()
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1048576),b""):h.update(block)
    return h.hexdigest()
def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f,delimiter="\t")
        assert rd.fieldnames and len(rd.fieldnames)==len(set(rd.fieldnames))
        rows=list(rd)
        assert all(None not in r and all(v is not None for v in r.values()) for r in rows)
        return rows
def index(rows,k):
    d={r[k]:r for r in rows};assert len(d)==len(rows) and all(d)
    return d
inputs={}
def pin(p,h=None):
    p=p.resolve();d=sha(p)
    if h is not None:assert d==h,"HASH_MISMATCH "+str(p)
    if str(p) in inputs:assert inputs[str(p)]==d
    inputs[str(p)]=d
    return p
def audit(d,marker=None,log="stdout.txt",required=()):
    m=pin(d/"SHA256SUMS.txt")
    entries={}
    for line in m.read_text().splitlines():
        x=re.fullmatch(r"([0-9a-f]{64}) [ *](.+)",line);assert x
        f=Path(x[2]);f=(d/f).resolve() if not f.is_absolute() else f.resolve()
        assert f not in entries
        entries[f]=x[1];pin(f,x[1])
    assert {Path(f).resolve() for f in required}<=entries.keys()
    if marker:
        f=(d/log).resolve()
        assert f in entries and marker in f.read_text().splitlines(),"AUDIT_MARKER_MISSING "+str(d)
    print("VERIFIED_AUDIT",d.name,len(entries))
    return entries
def finite(v):
    x=float(v);assert math.isfinite(x)
    return x
def write(p,rows):
    assert rows
    with p.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rows)
    p.chmod(0o400)
private=repo/"private_data";diag=repo/"diagnostics"
policy=pin(repo/"docs/V8_FUSCC331_PROPOSED_RELEASE_POLICY.md")
raw=policy.read_bytes()
assert hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()=="99f3748a496df19bbf873a4dbcb2a7e4992507c4"
approval=pin(repo/"docs/V8_FUSCC331_POLICY_APPROVAL_20261128.md")
assert 'responded "proceed"' in approval.read_text()
pin(Path(__file__))
cross=pin(private/"v8_fuscc_patient_crosswalk_addendum_v1/360_fuscc_patient_clinical_crosswalk_LOCAL_ONLY.tsv","8d34f94ba5834e6dc477fb740547ab46801cc89185e1b724ec49356d50400cd5")
book=pin(private/"v8_fuscc_original_clinical_source_v1/mmc2.xlsx","1581533099f7a7ce15d5ab39551fa93855c25b30fc624160b3940700151becb7")
audit(diag/"20261128_v8_fuscc_original_clinical_table_v1","INDEPENDENT_CANONICAL_CONTENT_REFERENCE_PASS",required=[book])
audit(diag/"20261128_v8_fuscc_patient_crosswalk_addendum_v1","V8_FUSCC_PATIENT_CROSSWALK_ADDENDUM_COMPLETE",required=[cross])
audit(diag/"20261128_v8_primary2055_atlas_exact_screen_v1","V8_PRIMARY2055_ATLAS_EXACT_SCREEN_COMPLETE")
ov=diag/"20261128_v8_fuscc_project_key_nonoverlap_v1"
audit(ov,"V8_FUSCC_PROJECT_KEY_NONOVERLAP_REVIEW_COMPLETE")
assert {"ATLAS_METADATA_UNION_PROJECT_KEYS_HIT 0","OTHER_EXTERNAL_PROJECT_KEYS_HIT 0"}<=set((ov/"stdout.txt").read_text().splitlines())
audit(diag/"20261128_v8_joint461_existing_projection_robustness_v1","V8_EXISTING_PROJECTION_ROBUSTNESS_REVIEW_COMPLETE")
for rel,h in [
("14_projection_use_analysis_ready_atlas_v1_0_3/03_projection_analysis_ready_core.tsv.gz","7f6852e829f7aaa11a26f86b44247c9c8c0840958494b7766db1ad6abd180092"),
("04_candidate_evidence_v1_0_4/04_candidate_sample_full_metadata.parquet","2e74aabd85a62c3f80835691a2ffa0c11f374bf262142349f579fd4e6f4fc36a"),
("05_study_sample_adjudication_v1_0_2_1/18_complete_adjudication_manifest.parquet","736e89b2969693f504d298c8e7274a1bcf5773a0f18667b114953c575ffd4a32")]:pin(root/rel,h)
prep=private/"v8_joint461_full_G026_QC_normalization_v1"
mp=prep/"461_amended_technical_input_manifest_LOCAL_ONLY.tsv"
pp=prep/"input_provenance.json"
audit(diag/"20261128_v8_joint461_full_G026_QC_normalization_v1","V8_JOINT461_FULL_G026_PREPARATION_COMPLETE_PENDING_REVIEW",required=[mp,pp])
prov=json.loads(pp.read_text())
assert (prov["membership"],prov["technical_runs"],prov["clinical_admissions"])==(461,483,0)
pin(mp,prov["amended_manifest_sha256"])
study="SRP157974"
for name,h in prov["studies"][study]["outputs_sha256"].items():
    assert Path(name).name==name
    pin(prep/study/name,h)
assert len(prov["studies"][study]["outputs_sha256"])==5
proj=private/"v8_joint461_fixed500_compatibility_v1_0_1"
pa=diag/"20261128_v8_joint461_fixed500_compatibility_v1_0_1"
audit(pa,"V8_JOINT461_FIXED500_COMPATIBILITY_V1_0_1_COMPLETE_PENDING_REVIEW","03_projection_stdout.txt",required=[proj/"SHA256SUMS.txt",pa/"01_integrity.txt"])
pe=audit(proj)
for line in (pa/"01_integrity.txt").read_text().splitlines():
    if line.startswith("REFERENCE_INPUT_SHA256 "):
        _,h,f=line.split(maxsplit=2);pin(Path(f),h)
print("FROZEN_REFERENCE_HASHES_REVERIFIED")
rows=read(cross)
by=index(rows,"experiment_acc_LOCAL_ONLY")
assert len(by)==len({r["clinical_project_id_LOCAL_ONLY"] for r in rows})==360
assert all(r["study"]==study and r["final_clinical_approval"]=="NO_UNSIGNED" and r["clinical_patient_crosswalk_status"]=="EXACT_ORIGINAL_TABLE_S1_PROJECT_ID_RNA_YES" for r in rows)
counts=Counter(r["proposed_review_disposition"] for r in rows)
assert counts=={"PROPOSED_PATIENT_LINKED_REVIEW_CANDIDATE":331,"HOLD_LIBRARY_SUFFIX_PROVENANCE":3,"HISTORICAL_QC_HOLD_PRESERVED":26}
keep={e for e,r in by.items() if r["proposed_review_disposition"]=="PROPOSED_PATIENT_LINKED_REVIEW_CANDIDATE"}
retained={e for e,r in by.items() if r["membership"]=="RETAINED334"}
assert len(retained)==334 and keep<=retained
assert all(by[e]["rep_suffix_status"]=="AGREES" for e in keep)
assert all(by[e]["rep_suffix_status"]=="DISAGREES_UNRESOLVED" for e in retained-keep)
man=index([r for r in read(mp) if r["study"]==study],"library_id_LOCAL_ONLY")
qcfile=prep/study/"library_expression_QC_LOCAL_ONLY.tsv"
assert str(qcfile.resolve()) in inputs
qc=index(read(qcfile),"library_id_LOCAL_ONLY")
assert set(man)==set(qc)==retained
assert all(r["input_status"]=="INPUT_READY_FOR_FULL_G026_QC_NOT_RELEASED" and r["final_clinical_approval"]=="NO_UNSIGNED" for r in man.values())
assert all(int(r["n_technical_runs"])==1 and r["UQ_evaluable"].upper()=="TRUE" for r in qc.values())
branches=("CPM","TMM_CPM","UQ_CPM");components=("NonSE",*("SE"+str(i) for i in range(1,10)))
scores={};comparisons={}
for branch in branches:
    f=proj/f"{study}_{branch}_projection_QC_LOCAL_ONLY.tsv"
    assert f.resolve() in pe
    scores[branch]=index(read(f),"sample_id")
    assert set(scores[branch])==retained
    for r in scores[branch].values():
        assert r["study"]==study and r["branch"]==branch
        v=[finite(r["fraction_"+c]) for c in components]
        h=[finite(r["raw_H_"+c]) for c in components]
        assert min(v)>=-1e-12 and min(h)>=-1e-12 and abs(sum(v)-1)<1e-8 and sum(h)>0
        assert max(abs(x-y/sum(h)) for x,y in zip(v,h))<1e-8
        for metric in ("normalized_KL_per_input_mass","reconstruction_RMSE","reconstruction_R2","reconstruction_cosine","normalized_score_entropy"):finite(r[metric])
for alt in branches[1:]:
    f=proj/f"{study}_CPM_vs_{alt}_LOCAL_ONLY.tsv"
    assert f.resolve() in pe
    comparisons[alt]=index(read(f),"sample_id")
    assert set(comparisons[alt])==retained
    for e,r in comparisons[alt].items():
        d=sum(abs(float(scores["CPM"][e]["fraction_"+c])-float(scores[alt][e]["fraction_"+c])) for c in components)/2
        assert r["comparison_evaluable"].upper()=="TRUE" and abs(d-finite(r["total_variation_distance"]))<1e-8
def add_prefixed(row,src,prefix,skip):
    for k,v in src.items():
        if k in skip:continue
        key=prefix+k;assert key not in row
        row[key]=v
manifest=[];dispositions=[]
for e,r in sorted(by.items()):
    approved=e in keep
    state="ADMITTED_UNDER_APPROVED_METADATA_POLICY" if approved else r["proposed_review_disposition"]
    dispositions.append({"study":study,"library_id_LOCAL_ONLY":e,"patient_project_id_LOCAL_ONLY":r["clinical_project_id_LOCAL_ONLY"],"original_disposition":r["proposed_review_disposition"],"release_disposition":state,"clinical_admission":"YES" if approved else "NO","release_version":"FUSCC331_V1"})
    if not approved:continue
    row=dict(r)
    # Original unsigned columns retained explicitly as historical evidence.
    row["release_version"]="FUSCC331_V1"
    row["release_disposition"]=state
    row["release_clinical_admission"]="YES"
    row["release_pretreatment_evidence"]="ACCEPTED_PUBLISHED_COHORT_LEVEL_NOT_INDIVIDUAL_TIMING"
    row["release_patient_unit"]="ONE_LIBRARY_PER_DISTINCT_ORIGINAL_CLINICAL_PROJECT_ID"
    row["release_independence_evidence"]="ACCEPTED_METADATA_NONOVERLAP_NOT_ABSOLUTE_BIOLOGICAL_INDEPENDENCE"
    row["release_primary_normalization"]="CPM"
    row["release_stratum"]="EXTERNAL_FUSCC_TNBC"
    add_prefixed(row,man[e],"input_",{"study","library_id_LOCAL_ONLY"})
    add_prefixed(row,qc[e],"expression_QC_",{"study","library_id_LOCAL_ONLY"})
    for branch in branches:
        add_prefixed(row,scores[branch][e],branch+"_",{"study","branch","sample_id"})
    for alt in branches[1:]:
        add_prefixed(row,comparisons[alt][e],"CPM_vs_"+alt+"_",{"sample_id"})
    manifest.append(row)
assert len(manifest)==331 and len(dispositions)==360
assert len({r["clinical_project_id_LOCAL_ONLY"] for r in manifest})==331
assert Counter(r["clinical_admission"] for r in dispositions)=={"YES":331,"NO":29}
batch=Counter(r["author_sequencing_batch_at_base_level"] for r in manifest)
assert set(batch)=={"batch1","batch2","batch3"}
coverage={}
for k in ("clinical_sTILs","clinical_iTILs","clinical_Ki67","clinical_RFS_Status","clinical_RFS_time_Days","clinical_Grade","clinical_T","clinical_N","author_tumour_purity_reported_percent"):
    coverage[k]=sum(bool(r[k].strip()) for r in manifest)
events=sum(float(r["clinical_RFS_Status"])==1 for r in manifest)
assert all(float(r["clinical_RFS_Status"]) in (0,1) for r in manifest)
# Check every recorded input again before publishing; do not change source files.
for f,h in inputs.items():assert sha(Path(f))==h,"INPUT_CHANGED "+f
tmp=Path(tempfile.mkdtemp(prefix=".v8_fuscc331_lock_staging_",dir=out.parent))
try:
    write(tmp/"331_clinical_projection_manifest_LOCAL_ONLY.tsv",manifest)
    write(tmp/"360_disposition_ledger_LOCAL_ONLY.tsv",dispositions)
    summary={"version":"FUSCC331_V1","created_utc":datetime.now(timezone.utc).isoformat(),
             "git_head":subprocess.check_output(["git","-C",str(repo),"rev-parse","HEAD"],text=True).strip(),
             "approval_record":str(approval),"policy":str(policy),"clinical_admissions":331,"held":29,
             "historical_QC_holds":26,"library_provenance_holds":3,"original_atlas_unchanged":9465,
             "hypothetical_combined_population":9796,"primary_normalization":"CPM","sensitivity_branches":["TMM_CPM","UQ_CPM"],
             "batch_counts":dict(batch),"clinical_marginal_coverage":coverage,"RFS_events":events,
             "input_sha256":inputs,"limitations":["Cohort-level pretreatment evidence, not individual treatment dates","Metadata-based nonoverlap, not absolute biological independence","Research compatibility, not assay equivalence or direct spatial validation","Historical unsigned fields preserved; release-specific admission fields authoritative","No retrospective outcome or ecotype-based exclusion"]}
    summary["output_sha256"]={f.name:sha(f) for f in sorted(tmp.iterdir())}
    f=tmp/"release_provenance_LOCAL_ONLY.json";f.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n");f.chmod(0o400)
    f=tmp/"SHA256SUMS.txt";f.write_text("".join(sha(x)+"  "+x.name+"\n" for x in sorted(tmp.iterdir()) if x.is_file()));f.chmod(0o400)
    assert not out.exists()
    # Exclusive final directory creation; do not replace an existing lock.
    out.mkdir(mode=0o700)
    for f in tmp.iterdir():f.rename(out/f.name)
    tmp.rmdir()
    out.chmod(0o500)
except Exception:
    print("LOCK_CREATION_INCOMPLETE_PRESERVE_STAGING_AND_INVESTIGATE")
    raise
print("FUSCC331_ADMITTED",len(manifest),"HELD",29)
print("BATCH_COUNTS",dict(batch))
print("CLINICAL_MARGINAL_COVERAGE",coverage,"RFS_EVENTS",events)
print("MANIFEST_SHA256",sha(out/"331_clinical_projection_manifest_LOCAL_ONLY.tsv"))
print("V8_FUSCC331_CLINICAL_LOCK_V1_COMPLETE")
print("ORIGINAL_ARCHS4_UNCHANGED_NO_REFIT_NO_PUBLIC_PATIENT_DATA")
