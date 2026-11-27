#!/usr/bin/env bash
# V7N one-shot read-only frozen-H four-feature diagnostic.
# No optimizer; no external recount3 values loaded. Immutable audit output.
set -euo pipefail
umask 077
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
SCRIPT="$REPO/tools/v7n_original_frozen_143_four_feature_KL_READONLY_v1.R"
AUDIT="$REPO/diagnostics/20261127_v7n_original_143_four_feature_KL_v1"
test -x "$RSCRIPT" && test -f "$SCRIPT" && test ! -e "$AUDIT"
"$RSCRIPT" -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("V7N_R_SYNTAX_PASS\n")' "$SCRIPT"
mkdir -p "$AUDIT"
"$RSCRIPT" "$SCRIPT" "$ROOT" 2>&1 | tee "$AUDIT/01_frozen_four_feature_KL_stdout.txt"
(
cd "$AUDIT"
sha256sum "$SCRIPT" 01_frozen_four_feature_KL_stdout.txt > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
)
echo "V7N_FOUR_FEATURE_FROZEN_H_KL_AUDIT_CAPTURED_PENDING_REVIEW"
