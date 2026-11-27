#!/usr/bin/env bash
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
ENV=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human
PYTHON="$ENV/bin/python"
RSCRIPT="$ENV/bin/Rscript"
STUDY="${1:?Supply configured study accession}"
case "$STUDY" in SRP074349|SRP157215) ;; *) echo "Unknown study" >&2; exit 1;; esac
CONFIG="$REPO/configs/v7r_matched_calibration_studies_v1.json"
PRE="$REPO/diagnostics/20261127_v7r_${STUDY}_preflight_v1"
AUDIT="$REPO/diagnostics/20261127_v7r_${STUDY}_calibration_v1"
SOURCE="$REPO/private_data/v7r_${STUDY}_prepared_v1"
OUT="$REPO/private_data/v7r_${STUDY}_projection_v1"
P="$REPO/tools/v7r_prepare_configured_matched_G026_full499_v1.py"
R="$REPO/tools/v7r_project_configured_matched_full500_KL_v1.R"
RAW="$ROOT/26_recount3_external_primary_v1/raw_studies/$STUDY/gene_sums/sra.gene_sums.$STUDY.G026.gz"
cd "$REPO"

test -f "$PRE/SHA256SUMS.txt"
for d in "$AUDIT" "$SOURCE" "$OUT"; do test ! -e "$d"; done
(cd "$PRE" && sha256sum -c SHA256SUMS.txt)
grep -Fxq V7R_ASSAY_KEY_REGRESSION_PASS "$PRE/04_assay_key_regression_stdout.txt"
mkdir "$AUDIT"
"$PYTHON" "$P" "$ROOT" "$SOURCE" "$CONFIG" "$STUDY" 2>&1 | tee "$AUDIT/01_preparation_stdout.txt"
grep -Fxq V7R_RECOUNT3_ONLY_MATCHED_FULL499_PREPARATION_PASS "$AUDIT/01_preparation_stdout.txt"
"$PYTHON" - "$SOURCE" "$CONFIG" "$STUDY" <<'PY' | tee "$AUDIT/02_prepared_integrity.txt"
import hashlib,json,sys
from pathlib import Path
src,config=map(Path,sys.argv[1:3])
study=sys.argv[3]
p=json.loads((src/"input_provenance.json").read_text())
assert p["study"]==study and p["study_contract_sha256"]==hashlib.sha256(config.read_bytes()).hexdigest()
assert len(p["outputs_sha256"])==5
for name,digest in p["outputs_sha256"].items():
    assert Path(name).name==name
    assert hashlib.sha256((src/name).read_bytes()).hexdigest()==digest
print("V7R_PREPARED_INPUT_SHA256_PASS")
PY
mkdir "$OUT"
"$RSCRIPT" "$R" "$ROOT" "$SOURCE" "$OUT" "$CONFIG" "$STUDY" 2>&1 | tee "$AUDIT/03_projection_stdout.txt"
grep -Fxq V7R_FROZEN_ARCHS4_VS_RECOUNT3_MATCHED_499_GENES_COMPLETED_PENDING_REVIEW "$AUDIT/03_projection_stdout.txt"
test -s "$OUT/projection_checkpoint_v7r_v1_LOCAL_ONLY.rds"
for f in aggregate_transfer_metrics.tsv aggregate_normalization_robustness.tsv aggregate_nine_gene_sensitivity.tsv aggregate_determinism_sentinel.tsv; do
    test -s "$OUT/$f"
    cp "$OUT/$f" "$AUDIT/$f"
done
(cd "$AUDIT" && sha256sum "$P" "$R" "$CONFIG" "$REPO/scripts/run_v7r_configured_matched_calibration_v1.sh" "$SOURCE/input_provenance.json" ./*.txt aggregate_*.tsv > SHA256SUMS.txt && sha256sum -c SHA256SUMS.txt)
echo "V7R_CONFIGURED_MATCHED_CALIBRATION_COMPLETE_PENDING_REVIEW"
echo "NO_FROZEN_FILES_MODIFIED_CLINICAL_ADMISSIONS_0"
