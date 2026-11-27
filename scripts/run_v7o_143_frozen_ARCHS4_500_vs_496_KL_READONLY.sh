#!/usr/bin/env bash
# V7O: immutable numerical diagnostic on frozen original ARCHS4 only.
# Executes official 143-context solver + controlled per-sample KL objectives
# with 500 vs 496 W rows, using two starting points. No recount3 model run.
set -euo pipefail
umask 077
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
CODE="$REPO/tools/v7o_143_frozen_ARCHS4_controlled_500_vs_496_KL_READONLY_v1.R"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
AUDIT="$REPO/diagnostics/20261127_v7o_143_frozen_ARCHS4_500_vs_496_KL_v1"

test -x "$RSCRIPT" && test -f "$CODE" && test ! -e "$AUDIT"
"$RSCRIPT" -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("V7O_R_SYNTAX_PASS\n")' "$CODE"
mkdir -p "$AUDIT"
"$RSCRIPT" "$CODE" "$ROOT" 2>&1 | tee "$AUDIT/01_500_vs_496_KL_stdout.txt"
(
cd "$AUDIT"
sha256sum "$CODE" 01_500_vs_496_KL_stdout.txt > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
)
echo "V7O_FROZEN_ARCHS4_KL_FEATURE_OPTIMIZER_DIAGNOSTIC_CAPTURED_PENDING_REVIEW"
