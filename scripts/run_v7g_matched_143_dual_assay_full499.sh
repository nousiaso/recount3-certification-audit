#!/usr/bin/env bash
# V7G real dual-assay projection: child Bash process, not source.
# Prefer dispatch via sbatch scripts/slurm_v7g_matched_143_dual_assay.sbatch
set -euo pipefail
umask 077
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
PYTHON="${PYTHON:-python3}"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
PREP="$REPO/tools/v7g_prepare_143_matched_dual_assay_G026_ARCHS4_v1.py"
PROJ="$REPO/tools/v7g_project_143_dual_assay_frozenW_v1.R"
OUT="$REPO/private_data/v7g_143_matched_dual_assay_full499_v1"
AUDIT="$REPO/diagnostics/20261127_v7g_143_matched_dual_assay_full499_v1"

test -f "$PREP" && test -f "$PROJ" && test -x "$RSCRIPT"
test ! -e "$OUT"
test ! -e "$AUDIT"
"$PYTHON" - "$PREP" <<'PY'
import ast,importlib.util,sys
from pathlib import Path
ast.parse(Path(sys.argv[1]).read_text(encoding="utf-8"))
for name in ("h5py","numpy","pandas"):
    if importlib.util.find_spec(name) is None:
        raise SystemExit("MISSING_PYTHON_DEPENDENCY "+name)
print("V7G_PYTHON_SYNTAX_AND_DEPENDENCIES_PASS")
PY
"$RSCRIPT" -e 'parse(file=commandArgs(trailingOnly=TRUE)[[1]]);cat("V7G_R_SYNTAX_PASS\n")' "$PROJ"
"$RSCRIPT" -e 'stopifnot(requireNamespace("SpatialEcoTyper",quietly=TRUE), requireNamespace("NMF",quietly=TRUE), requireNamespace("data.table",quietly=TRUE));cat("V7G_R_LIBRARIES_PASS\n")'
mkdir -p "$AUDIT"

echo "V7G_RAW_G026_454_RUNS_AND_ARCHS4_143_GSM_PREPARATION_START"
"$PYTHON" -u "$PREP" "$ROOT" "$OUT" 2>&1 | tee "$AUDIT/01_python_preparation_stdout.txt"
echo "V7G_FROZEN_W_DUAL_ASSAY_PROJECTION_START"
"$RSCRIPT" "$PROJ" "$ROOT" "$OUT" 2>&1 | tee "$AUDIT/02_R_projection_stdout.txt"
cp "$OUT/aggregate_transfer_metrics.tsv" "$AUDIT/aggregate_transfer_metrics.tsv"
cp "$OUT/aggregate_normalization_sensitivity.tsv" "$AUDIT/aggregate_normalization_sensitivity.tsv"
cp "$OUT/aggregate_nine_gene_sensitivity.tsv" "$AUDIT/aggregate_nine_gene_sensitivity.tsv"
cp "$OUT/input_provenance.json" "$AUDIT/input_provenance_aggregate.json"

(
  cd "$AUDIT"
  sha256sum "$PREP" "$PROJ" *.txt aggregate_*.tsv input_provenance_aggregate.json > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)
echo "LOCAL_SPECIMEN_LINKED_RESULTS $OUT"
echo "AGGREGATE_METRICS $AUDIT/aggregate_transfer_metrics.tsv"
echo "V7G_143_MATCHED_DUAL_ASSAY_CALIBRATION_COMPLETE_PENDING_SCIENTIFIC_REVIEW"
