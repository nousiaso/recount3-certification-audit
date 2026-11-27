#!/usr/bin/env bash
# V7K source inspection only; creates immutable small diagnostic stdout.
# Does not run SpatialEcoTyper projection or load biomedical matrices.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
SCRIPT="$REPO/tools/v7k_inspect_installed_DeconvoluteSE_SOURCE_READONLY_v1.R"
AUDIT="$REPO/diagnostics/20261127_v7k_DeconvoluteSE_installed_source_v1"
test -x "$RSCRIPT"
test -f "$SCRIPT"
test ! -e "$AUDIT"
"$RSCRIPT" -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("V7K_R_SOURCE_SYNTAX_PASS\n")' "$SCRIPT"
mkdir -p "$AUDIT"
"$RSCRIPT" "$SCRIPT" 2>&1 | tee "$AUDIT/01_installed_DeconvoluteSE_source_stdout.txt"
(
  cd "$AUDIT"
  sha256sum "$SCRIPT" 01_installed_DeconvoluteSE_source_stdout.txt > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)
echo "V7K_INSTALLED_SOURCE_INSPECTION_CAPTURED_PENDING_REVIEW"
