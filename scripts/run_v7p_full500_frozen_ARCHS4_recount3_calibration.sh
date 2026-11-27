#!/usr/bin/env bash
# V7P full same-specimen assay transfer, starts from certified V7J prepared inputs.
# Reuses 454->143 G026 aggregation WITHOUT re-streaming/re-extracting.
# Technical calibration only; no clinical admissions.
set -euo pipefail
umask 077
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
SCRIPT="$REPO/tools/v7p_frozen_ARCHS4_vs_recount3_143pairs_full500_KL_v1.R"
V7J="$REPO/private_data/v7j_143_frozen_archs4_vs_recount3_full499_v1"
PRE="$REPO/diagnostics/20261127_v7p_static_preflight_v1"
OUT="$REPO/private_data/v7p_143_frozen_ARCHS4_vs_recount3_full500_KL_v1"
AUDIT="$REPO/diagnostics/20261127_v7p_143_frozen_ARCHS4_vs_recount3_full500_KL_v1"
test -x "$PYTHON" && test -x "$RSCRIPT"
test -f "$SCRIPT" && test -f "$V7J/input_provenance.json"
test ! -e "$OUT" && test ! -e "$AUDIT"
test -f "$PRE/SHA256SUMS.txt"
(
 cd "$PRE"
 sha256sum -c SHA256SUMS.txt
)
grep -q '^V7P_R_SYNTAX_PASS$' "$PRE/01_R_syntax_stdout.txt"
grep -q '^V7P_R_PACKAGES_PASS$' "$PRE/02_R_packages_stdout.txt"
grep -q '^V7P_V7J_PREPARED_INPUT_REUSE_SHA256_PASS' "$PRE/03_V7J_reused_inputs_stdout.txt"
mkdir -p "$OUT" "$AUDIT"
echo "V7P_REUSE_CERTIFIED_V7J_PREPARED_EXPRESSION_NO_RAW_EXTRACTION"
"$RSCRIPT" "$SCRIPT" "$ROOT" "$V7J" "$OUT" 2>&1 | tee "$AUDIT/01_paired_full500_KL_stdout.txt"
grep -q '^V7P_FROZEN_ARCHS4_VS_RECOUNT3_143_PAIRS_499_GENES_COMPLETED_PENDING_REVIEW$' "$AUDIT/01_paired_full500_KL_stdout.txt"
for f in aggregate_transfer_metrics.tsv aggregate_normalization_robustness.tsv aggregate_nine_gene_sensitivity.tsv aggregate_determinism_sentinel.tsv; do
  test -s "$OUT/$f"
  cp "$OUT/$f" "$AUDIT/$f"
done
(
 cd "$AUDIT"
 sha256sum "$SCRIPT" "$V7J/input_provenance.json" 01_paired_full500_KL_stdout.txt aggregate_*.tsv > SHA256SUMS.txt
 sha256sum -c SHA256SUMS.txt
)
echo "V7P_PAIRED_ARCHS4_RECOUNT3_FULL500_KL_AGGREGATES_CAPTURED_PENDING_REVIEW"
echo "NO_FROZEN_FILES_MODIFIED_CLINICAL_ADMISSIONS_0"
