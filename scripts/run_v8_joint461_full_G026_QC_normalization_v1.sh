#!/usr/bin/env bash
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
P="$REPO/tools/v8_joint461_full_G026_QC_normalization_v1.py"
HELPER="$REPO/tools/v7r_prepare_configured_matched_G026_full499_v1.py"
OUT="$REPO/private_data/v8_joint461_full_G026_QC_normalization_v1"
AUDIT="$REPO/diagnostics/20261128_v8_joint461_full_G026_QC_normalization_v1"
cd "$REPO"
test ! -e "$OUT"; test ! -e "$AUDIT"
"$PYTHON" - "$P" "$HELPER" <<'PY'
import ast,sys
from pathlib import Path
import numpy,pandas,h5py
for p in sys.argv[1:]:ast.parse(Path(p).read_text())
print("SYNTAX_AND_DEPENDENCIES_PASS")
PY
mkdir "$AUDIT"
git rev-parse HEAD > "$AUDIT/git_head.txt"
"$PYTHON" -u "$P" "$ROOT" "$OUT" 2>&1 | tee "$AUDIT/stdout.txt"
grep -Fxq V8_JOINT461_FULL_G026_PREPARATION_COMPLETE_PENDING_REVIEW "$AUDIT/stdout.txt"
cp "$OUT/4_study_expression_QC_aggregate_UNSIGNED.tsv" "$AUDIT/"
(cd "$AUDIT";sha256sum "$P" "$HELPER" "$REPO/scripts/run_v8_joint461_full_G026_QC_normalization_v1.sh" "$OUT/461_amended_technical_input_manifest_LOCAL_ONLY.tsv" "$OUT/input_provenance.json" git_head.txt stdout.txt 4_study_expression_QC_aggregate_UNSIGNED.tsv > SHA256SUMS.txt;sha256sum -c SHA256SUMS.txt)
echo V8_JOINT461_FULL_G026_AUDIT_CAPTURED
