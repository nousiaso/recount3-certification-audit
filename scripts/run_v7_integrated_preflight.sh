#!/usr/bin/env bash
# Run inside a child Bash: bash scripts/run_v7_integrated_preflight.sh
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
SCRIPT="$REPO/tools/v7_integrated_candidate_gates_READONLY_v1.py"
OUT="$REPO/private_data/v7_integrated_candidate_gates_v1"
AUDIT="$REPO/diagnostics/20261127_v7_integrated_candidate_gates_v1"

EXPECTED=b93996199e15e1b7d2439ed2db2aa94808275610e943f6b9c3fc1f7d1e31f702
ACTUAL=$(sha256sum "$SCRIPT" | awk '{print $1}')
test "$ACTUAL" = "$EXPECTED" || { echo 'SCRIPT_HASH_MISMATCH'; exit 1; }
test ! -e "$OUT"
test ! -e "$AUDIT"
mkdir -p "$AUDIT"
python3 - "$SCRIPT" <<'PY'
import ast,sys
from pathlib import Path
ast.parse(Path(sys.argv[1]).read_text(encoding='utf-8'))
print('PYTHON_SYNTAX_PASS')
PY
python3 -u "$SCRIPT" "$ROOT" --output "$OUT" 2>&1 | tee "$AUDIT/stdout.txt"
(
  cd "$AUDIT"
  sha256sum stdout.txt > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)
echo "LOCAL_ONLY_RESULTS $OUT"
echo "AUDIT_CAPTURE $AUDIT"
echo 'V7_FIRST_COHORT_WIDE_PREFLIGHT_CAPTURED_PENDING_REVIEW'
