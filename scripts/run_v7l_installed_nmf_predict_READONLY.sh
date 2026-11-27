#!/usr/bin/env bash
# V7L: read-only installed low-level SpatialEcoTyper NMF code inspection.
# No model run / expression read. One-shot immutable diagnostic.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
SCRIPT="$REPO/tools/v7l_inspect_installed_nmf_predict_READONLY_v1.R"
AUDIT="$REPO/diagnostics/20261127_v7l_installed_nmf_predict_source_v1"
test -x "$RSCRIPT"
test -f "$SCRIPT"
test ! -e "$AUDIT"
"$RSCRIPT" -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("V7L_R_SYNTAX_PASS\n")' "$SCRIPT"
mkdir -p "$AUDIT"
"$RSCRIPT" "$SCRIPT" 2>&1 | tee "$AUDIT/01_installed_nmf_predict_stdout.txt"
(
cd "$AUDIT"
sha256sum "$SCRIPT" 01_installed_nmf_predict_stdout.txt > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
)
echo "V7L_INSTALLED_NMF_PREDICT_SOURCE_CAPTURED_PENDING_REVIEW"
