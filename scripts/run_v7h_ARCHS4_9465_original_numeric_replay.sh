#!/usr/bin/env bash
# Frozen ARCHS4 9,465-sample numerical replay preflight and three-branch replay.
# Run as a child bash under SLURM, do not source or rerun into same audit folder.
set -euo pipefail
umask 077
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
SRC="$ROOT/RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
SCRIPT="$REPO/tools/v7h_replay_frozen_ARCHS4_9465_stages09_10_READONLY_v1.R"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
AUDIT="$REPO/diagnostics/20261127_v7h_ARCHS4_9465_original_numeric_replay_v1"

test -f "$SCRIPT"
test -x "$RSCRIPT"
test -x "$PYTHON"
test ! -e "$AUDIT"

test "$(sha256sum "$SRC/09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py" | cut -d' ' -f1)" = "46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500"
test "$(sha256sum "$SRC/10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R" | cut -d' ' -f1)" = "16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c"

"$RSCRIPT" -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]])); cat("V7H_R_SYNTAX_PASS\n")' "$SCRIPT"
mkdir -p "$AUDIT"
"$PYTHON" "$SRC/09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py" --run-self-tests 2>&1 | tee "$AUDIT/01_original_stage09_selftests.txt"
"$RSCRIPT" "$SRC/10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R" --run-self-tests 2>&1 | tee "$AUDIT/02_original_stage10_selftests.txt"

echo "V7H_REPLAY_ALL_THREE_NORMALIZATION_AND_PROJECTION_BRANCHES_START"
"$RSCRIPT" "$SCRIPT" "$ROOT" 2>&1 | tee "$AUDIT/03_ARCHS4_9465_numerical_replay_stdout.txt"
(
cd "$AUDIT"
sha256sum "$SCRIPT" *.txt > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
)
echo "V7H_ARCHS4_ORIGINAL_NUMERICAL_REPLAY_CAPTURED_PENDING_REVIEW"
