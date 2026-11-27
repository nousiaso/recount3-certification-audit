#!/usr/bin/env bash
# V7I: strictly read-only transfer contract before revising/running V7G.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
SCRIPT="$REPO/tools/v7i_post_v7h2_transfer_contract_READONLY_v1.py"
AUDIT="$REPO/diagnostics/20261127_v7i_post_v7h2_transfer_contract_v1"
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
test -f "$SCRIPT"
test -x "$PYTHON"
test ! -e "$AUDIT"
"$PYTHON" - "$SCRIPT" <<'PY'
import ast,sys
from pathlib import Path
ast.parse(Path(sys.argv[1]).read_text(encoding="utf-8"))
print("V7I_PYTHON_SYNTAX_PASS")
PY
mkdir -p "$AUDIT"
"$PYTHON" -u "$SCRIPT" "$ROOT" 2>&1 | tee "$AUDIT/01_readonly_preflight_stdout.txt"
(
cd "$AUDIT"
sha256sum "$SCRIPT" 01_readonly_preflight_stdout.txt > SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
)
echo "V7I_POST_V7H2_TRANSFER_CONTRACT_CAPTURED_PENDING_REVIEW"
