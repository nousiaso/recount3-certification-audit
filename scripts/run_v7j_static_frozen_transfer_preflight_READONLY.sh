#!/usr/bin/env bash
# V7J preflight only: verified scripts and original function fidelity.
# No expression values loaded; no model run. Preserve one immutable audit.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
CHECK="$REPO/tools/v7j_static_frozen_code_contract_READONLY_v1.py"
PREP="$REPO/tools/v7j_prepare_recount3_454runs_143pairs_frozen_ARCHS4_v1.py"
PROJ="$REPO/tools/v7j_project_recount3_vs_original_frozen_ARCHS4_143pairs_v1.R"
AUDIT="$REPO/diagnostics/20261127_v7j_static_frozen_transfer_preflight_v1"

test -x "$PYTHON" && test -x "$RSCRIPT"
test -f "$CHECK" && test -f "$PREP" && test -f "$PROJ"
test ! -e "$AUDIT"
mkdir -p "$AUDIT"
"$PYTHON" -u "$CHECK" "$ROOT" 2>&1 | tee "$AUDIT/01_static_fidelity_stdout.txt"
"$RSCRIPT" -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("V7J_R_SOURCE_SYNTAX_PASS\n")' "$PROJ" 2>&1 | tee "$AUDIT/02_R_parse_stdout.txt"
"$RSCRIPT" -e 'stopifnot(requireNamespace("SpatialEcoTyper",quietly=TRUE),requireNamespace("data.table",quietly=TRUE),requireNamespace("NMF",quietly=TRUE),requireNamespace("jsonlite",quietly=TRUE)); cat("V7J_R_PACKAGE_PREFLIGHT_PASS\n")' 2>&1 | tee "$AUDIT/03_R_package_stdout.txt"
"$PYTHON" - <<'PY' 2>&1 | tee "$AUDIT/04_python_package_stdout.txt"
import h5py,numpy,pandas
print("V7J_PYTHON_PACKAGE_PREFLIGHT_PASS")
PY
(
  cd "$AUDIT"
  sha256sum "$CHECK" "$PREP" "$PROJ" ./*.txt > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)
echo "V7J_READONLY_STATIC_PREFLIGHT_CAPTURED_PENDING_REVIEW"
