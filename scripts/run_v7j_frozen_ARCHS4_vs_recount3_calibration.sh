#!/usr/bin/env bash
# V7J full matched-assay calibration, run ONLY after V7J static preflight PASS.
# No source overwrites; all specimen-linked output stays under private_data/.
set -euo pipefail
umask 077
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
CHECK="$REPO/tools/v7j_static_frozen_code_contract_READONLY_v1.py"
PREP="$REPO/tools/v7j_prepare_recount3_454runs_143pairs_frozen_ARCHS4_v1.py"
PROJ="$REPO/tools/v7j_project_recount3_vs_original_frozen_ARCHS4_143pairs_v1.R"
PREF="$REPO/diagnostics/20261127_v7j_static_frozen_transfer_preflight_v1"
OUT="$REPO/private_data/v7j_143_frozen_archs4_vs_recount3_full499_v1"
AUDIT="$REPO/diagnostics/20261127_v7j_143_frozen_archs4_vs_recount3_full499_v1"

test -x "$PYTHON" && test -x "$RSCRIPT"
test -f "$CHECK" && test -f "$PREP" && test -f "$PROJ"
test ! -e "$OUT"
test ! -e "$AUDIT"
test -f "$PREF/01_static_fidelity_stdout.txt"
test -f "$PREF/02_R_parse_stdout.txt"
test -f "$PREF/03_R_package_stdout.txt"
test -f "$PREF/04_python_package_stdout.txt"
grep -q '^V7J_STATIC_FIDELITY_CONTRACT_PASS$' "$PREF/01_static_fidelity_stdout.txt"
grep -q '^V7J_R_SOURCE_SYNTAX_PASS$' "$PREF/02_R_parse_stdout.txt"
grep -q '^V7J_R_PACKAGE_PREFLIGHT_PASS$' "$PREF/03_R_package_stdout.txt"
grep -q '^V7J_PYTHON_PACKAGE_PREFLIGHT_PASS$' "$PREF/04_python_package_stdout.txt"

mkdir -p "$AUDIT"
# Recheck exact source code and source-input hashes even if files changed
# between static preflight and SLURM execution.
"$PYTHON" -u "$CHECK" "$ROOT" 2>&1 | tee "$AUDIT/00_reverified_static_contract.txt"

echo "V7J_G026_PREPARATION_START"
"$PYTHON" -u "$PREP" "$ROOT" "$OUT" 2>&1 | tee "$AUDIT/01_recount3_preparation_stdout.txt"
echo "V7J_ORIGINAL_STAGE10_PROJECTION_START"
"$RSCRIPT" "$PROJ" "$ROOT" "$OUT" 2>&1 | tee "$AUDIT/02_frozen_archs4_vs_recount3_stdout.txt"

cp "$OUT/aggregate_transfer_metrics.tsv" "$AUDIT/aggregate_transfer_metrics.tsv"
cp "$OUT/aggregate_normalization_robustness.tsv" "$AUDIT/aggregate_normalization_robustness.tsv"
cp "$OUT/aggregate_nine_gene_sensitivity.tsv" "$AUDIT/aggregate_nine_gene_sensitivity.tsv"
cp "$OUT/aggregate_determinism_sentinel.tsv" "$AUDIT/aggregate_determinism_sentinel.tsv"
cp "$OUT/input_provenance.json" "$AUDIT/aggregate_input_provenance.json"
(
cd "$AUDIT"
sha256sum "$CHECK" "$PREP" "$PROJ" ./*.txt ./aggregate_*.tsv ./aggregate_input_provenance.json > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
)
echo "V7J_FROZEN_ARCHS4_RECOUNT3_TECHNICAL_TRANSFER_CAPTURED_PENDING_SCIENTIFIC_REVIEW"
echo "NO_FROZEN_FILES_MODIFIED_CLINICAL_ADMISSIONS_0"
