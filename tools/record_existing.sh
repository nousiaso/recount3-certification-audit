#!/usr/bin/env bash
# Import already executed diagnostic code and its recorded output. Does not run the code.
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo 'Usage: tools/record_existing.sh YYYYMMDD_TOPIC /abs/or/relative/script.py /abs/or/relative/stdout.txt' >&2
    exit 2
fi
id=$1
src=$2
log=$3
[[ "$id" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo 'Unsafe diagnostic ID' >&2; exit 2; }
[[ -f "$src" && -f "$log" ]] || { echo 'Script or log missing' >&2; exit 2; }
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
out="$repo/diagnostics/$id"
if [[ -e "$out/script.py" || -e "$out/stdout.txt" ]]; then
    echo 'A diagnostic with this ID already exists; choose a new version.' >&2
    exit 3
fi
mkdir -p "$out"
cp -- "$src" "$out/script.py"
cp -- "$log" "$out/stdout.txt"
(
 cd "$out"
 sha256sum script.py stdout.txt > SHA256SUMS.txt
)
cat > "$out/interpretation.md" <<'EOF'
# Imported existing diagnostic

**Provenance:** already executed on the cluster; code and stdout copied as supplied.
Check source-code identity against execution logs where available.

**Findings:** TO REVIEW. Do not infer clinical approval from a printed PASS.

**Open questions:** TO REVIEW.
EOF
printf 'Imported %s without running or modifying either original.\n' "$id"
printf 'Review %s then update AUDIT_REGISTRY.tsv and commit.\n' "$out"
