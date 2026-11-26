#!/usr/bin/env bash
# V7F: read-only raw recount3 gzip header membership vs 454 bridge runs.
# Execute as a child Bash process (not source); no extraction/download.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
SCRIPT="$REPO/tools/v7f_raw_srp074349_gene_sums_header_READONLY_v1.py"
AUDIT="$REPO/diagnostics/20261127_v7f_raw_gene_sums_membership_v1"
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
echo "V7F_RAW_G026_HEADER_MEMBERSHIP_CAPTURED_PENDING_REVIEW"
