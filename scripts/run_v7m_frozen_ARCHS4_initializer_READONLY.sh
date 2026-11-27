#!/usr/bin/env bash
# V7M: one-shot read-only comparison of original ARCHS4 cohort-dependent
# feature retention and initial H. No NMF optimization, no recount3 expression.
set -euo pipefail
umask 077
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
CODE="$REPO/tools/v7m_original_ARCHS4_9465_vs_143_initialization_READONLY_v1.R"
AUDIT="$REPO/diagnostics/20261127_v7m_ARCHS4_9465_vs_143_initializer_v1"

test -x "$RSCRIPT"
test -f "$CODE"
test ! -e "$AUDIT"
"$RSCRIPT" -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("V7M_R_SYNTAX_PASS\n")' "$CODE"
mkdir -p "$AUDIT"
"$RSCRIPT" "$CODE" "$ROOT" 2>&1 | tee "$AUDIT/01_original_frozen_context_stdout.txt"
(
  cd "$AUDIT"
  sha256sum "$CODE" 01_original_frozen_context_stdout.txt > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)
echo "V7M_FROZEN_REFERENCE_NMF_INITIALIZER_AUDIT_CAPTURED_PENDING_REVIEW"
