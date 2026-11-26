#!/usr/bin/env bash
# V7C full-499 matched cross-assay projection. Execute as child Bash, not source.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
PYTHON="${PYTHON:-python3}"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
PREP="$REPO/tools/v7c_prepare_matched_143_G026_499.py"
PROJ="$REPO/tools/v7c_project_matched_143_frozenW.R"
OUT="$REPO/private_data/v7c_srp074349_matched_transfer_143_v1"
AUDIT="$REPO/diagnostics/20261127_v7c_srp074349_matched_transfer_v1"

test -f "$PREP"
test -f "$PROJ"
test -x "$RSCRIPT"
test ! -e "$OUT"
test ! -e "$AUDIT"

"$PYTHON" - "$PREP" <<'PY'
import ast, sys, importlib.util
from pathlib import Path
ast.parse(Path(sys.argv[1]).read_text(encoding="utf-8"))
for name in ("h5py","numpy","pandas"):
    if importlib.util.find_spec(name) is None:
        raise SystemExit("REQUIRED_PYTHON_MODULE_MISSING "+name)
print("PYTHON_SYNTAX_AND_MODULES_PASS")
PY

"$RSCRIPT" -e 'stopifnot(requireNamespace("SpatialEcoTyper",quietly=TRUE),requireNamespace("NMF",quietly=TRUE),requireNamespace("data.table",quietly=TRUE));cat("R_PACKAGES_PASS\n")'

mkdir -p "$AUDIT"
echo "V7C_COMPUTATIONAL_CALIBRATION_STARTED"
"$PYTHON" -u "$PREP" "$ROOT" "$OUT" 2>&1 | tee "$AUDIT/01_python_preparation_stdout.txt"
"$RSCRIPT" "$PROJ" "$ROOT" "$OUT" 2>&1 | tee "$AUDIT/02_R_projection_stdout.txt"
cp "$OUT/aggregate_transfer_metrics.tsv" "$AUDIT/aggregate_transfer_metrics.tsv"

(
  cd "$AUDIT"
  sha256sum "$PREP" "$PROJ" *.txt aggregate_transfer_metrics.tsv > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)
echo "LOCAL_PAIR_RESULTS $OUT"
echo "AGGREGATE_TRANSFER_METRICS $AUDIT/aggregate_transfer_metrics.tsv"
echo "V7C_MATCHED_143_FROZEN_W_TRANSFER_CAPTURED_PENDING_INTERPRETATION"
