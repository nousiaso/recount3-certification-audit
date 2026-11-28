#!/usr/bin/env bash
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
ENV=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human
SRC="$REPO/private_data/v8_joint461_full_G026_QC_normalization_v1"
PRE="$REPO/diagnostics/20261128_v8_joint461_full_G026_QC_normalization_v1"
V7="$REPO/diagnostics/20261127_v7r_SRP157215_calibration_v1"
P="$REPO/tools/v8_joint461_fixed500_compatibility_v1.R"
HELPER="$REPO/tools/v7r_project_configured_matched_full500_KL_v1.R"
OUT="$REPO/private_data/v8_joint461_fixed500_compatibility_v1"
AUDIT="$REPO/diagnostics/20261128_v8_joint461_fixed500_compatibility_v1"
cd "$REPO"
test ! -e "$OUT";test ! -e "$AUDIT"
(cd "$PRE";sha256sum -c SHA256SUMS.txt)
(cd "$V7";sha256sum -c SHA256SUMS.txt)
grep -Fxq V8_JOINT461_FULL_G026_PREPARATION_COMPLETE_PENDING_REVIEW "$PRE/stdout.txt"
mkdir "$AUDIT"
"$ENV/bin/python" - "$ROOT" "$SRC" "$V7" "$HELPER" <<'PY' | tee "$AUDIT/01_integrity.txt"
import hashlib,json,sys
from pathlib import Path
root,src,v7,helper=map(Path,sys.argv[1:])
repo=root/"RECOUNT3_CERTIFICATION_GIT"
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def check(p,h):
    assert sha(p)==h,str(p)
v7entries={}
for l in (v7/"SHA256SUMS.txt").read_text().splitlines():
    h,n=l.split(maxsplit=1);p=Path(n.lstrip("* "))
    p=p if p.is_absolute() else v7/p
    v7entries[p.resolve()]=h
check(helper,v7entries[helper.resolve()])
control=repo/"private_data/v7r_SRP157215_prepared_v1"
cp=control/"input_provenance.json";check(cp,v7entries[cp.resolve()])
c=json.loads(cp.read_text());check(control/"pair_manifest_LOCAL_ONLY.tsv",c["outputs_sha256"]["pair_manifest_LOCAL_ONLY.tsv"])
p=json.loads((src/"input_provenance.json").read_text())
assert p["version"]=="v8_joint461_full_G026_QC_normalization_v1" and (p["membership"],p["technical_runs"],p["clinical_admissions"])==(461,483,0)
assert set(p["studies"])=={"SRP157974","SRP057793","SRP002628","SRP094795"}
check(src/"461_amended_technical_input_manifest_LOCAL_ONLY.tsv",p["amended_manifest_sha256"])
check(src/"4_study_expression_QC_aggregate_UNSIGNED.tsv",p["aggregate_sha256"])
for s,v in p["studies"].items():
    assert len(v["outputs_sha256"])==5
    for name,h in v["outputs_sha256"].items():
        assert Path(name).name==name
        check(src/s/name,h)
    print("PREPARED_STUDY_HASHES_PASS",s,len(v["outputs_sha256"]))
rel=root/"RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
s9="10_preprojection_normalization_v1_0_1";s10="11_pan_cancer_fixed_anchor_projection_v1_0_4"
def find(stage,name):
    return next(x for x in (root/stage/name,rel/stage/name) if x.is_file())
gene=find(s9,"08_model_gene_availability.tsv")
check(gene,p["model_gene_ledger_sha256"])
assert p["model_gene_ledger_sha256"]==c["frozen_model_gene_audit_sha256"]
assets=[gene,find(s9,"01_canonical_library_manifest.tsv.gz"),find(s10,"16A_frozen_W_gene_symbol_harmonization.tsv"),find(s10,"19_primary_CPM_fixed_W_model_adequacy.tsv.gz"),Path("/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds")]
for name in ("10_model_gene_CPM_canonical.tsv.gz","11_model_gene_TMM_CPM_canonical.tsv.gz","12_model_gene_UQ_CPM_canonical.tsv.gz"):assets.append(find(s9,name))
for name in ("06_fixed_anchor_scaled_input_CPM.rds","07_fixed_anchor_scaled_input_TMM_CPM.rds","08_fixed_anchor_scaled_input_UQ_CPM.rds","09_raw_H_CPM.tsv.gz","11_raw_H_TMM_CPM.tsv.gz","13_raw_H_UQ_CPM.tsv.gz","10_fraction_scores_CPM.tsv.gz","12_fraction_scores_TMM_CPM.tsv.gz","14_fraction_scores_UQ_CPM.tsv.gz","03_pan_cancer_anchor_scaling_CPM.tsv.gz","04_pan_cancer_anchor_scaling_TMM_CPM.tsv.gz","05_pan_cancer_anchor_scaling_UQ_CPM.tsv.gz"):assets.append(find(s10,name))
for f in assets:print("REFERENCE_INPUT_SHA256",sha(f),str(f))
print("V8_PREPARED_INTEGRITY_PASS")
PY
grep -Fxq V8_PREPARED_INTEGRITY_PASS "$AUDIT/01_integrity.txt"
"$ENV/bin/Rscript" --vanilla -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("R_SYNTAX_PASS\n")' "$P" | tee "$AUDIT/02_R_syntax.txt"
mkdir "$OUT"
git rev-parse HEAD > "$AUDIT/git_head.txt"
"$ENV/bin/Rscript" --vanilla "$P" "$ROOT" "$SRC" "$OUT" 2>&1 | tee "$AUDIT/03_projection_stdout.txt"
grep -Fxq V8_JOINT461_FIXED500_COMPATIBILITY_COMPLETE_PENDING_REVIEW "$AUDIT/03_projection_stdout.txt"
# Verify frozen reference files did not change during execution.
"$ENV/bin/python" - "$AUDIT/01_integrity.txt" <<'PY'
import hashlib,sys
from pathlib import Path
for l in Path(sys.argv[1]).read_text().splitlines():
    if not l.startswith("REFERENCE_INPUT_SHA256 "):continue
    _,h,p=l.split(maxsplit=2)
    d=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""):d.update(b)
    assert d.hexdigest()==h,p
print("REFERENCE_INPUTS_UNCHANGED")
PY
for f in aggregate_reconstruction_adequacy.tsv aggregate_normalization_sensitivity.tsv aggregate_component_distribution_DESCRIPTIVE.tsv aggregate_determinism.tsv aggregate_ARCHS4_original_CPM_reference_DESCRIPTIVE.tsv;do test -s "$OUT/$f";cp "$OUT/$f" "$AUDIT/";done
(cd "$OUT";sha256sum ./*.tsv ./*.rds > SHA256SUMS.txt;sha256sum -c SHA256SUMS.txt)
(cd "$AUDIT";sha256sum "$P" "$HELPER" "$REPO/scripts/run_v8_joint461_fixed500_compatibility_v1.sh" "$OUT/SHA256SUMS.txt" ./*.txt aggregate_*.tsv > SHA256SUMS.txt;sha256sum -c SHA256SUMS.txt)
echo V8_JOINT461_FIXED500_AUDIT_CAPTURED
