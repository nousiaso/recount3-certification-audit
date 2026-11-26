#!/usr/bin/env bash
# V7D: READ-ONLY mapping diagnostic; preserves failed V7C artifacts.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
SCRIPT="$REPO/tools/v7d_srp074349_bridge_expression_identity_READONLY_v1.py"
AUDIT="$REPO/diagnostics/20261127_v7d_srp074349_bridge_expression_identity_v1"
PYTHON="${PYTHON:-python3}"

test -f "$SCRIPT"
test ! -e "$AUDIT"
"$PYTHON" - "$SCRIPT" <<'PY'
import ast,sys
from pathlib import Path
ast.parse(Path(sys.argv[1]).read_text(encoding="utf-8"))
print("PYTHON_SYNTAX_PASS")
PY
mkdir -p "$AUDIT"
"$PYTHON" -u "$SCRIPT" "$ROOT" 2>&1 | tee "$AUDIT/stdout.txt"
(
cd "$AUDIT"
sha256sum "$SCRIPT" stdout.txt > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
)
echo "V7D_BRIDGE_EXPRESSION_IDENTITY_CAPTURED_PENDING_REVIEW"
