#!/usr/bin/env bash
# Read-only Git snapshot validation; does not rerun cohort diagnostics.
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
for d in "$repo"/diagnostics/*; do
  [[ -d "$d" ]] || continue
  if [[ -f "$d/SHA256SUMS.txt" ]]; then
    (cd "$d" && sha256sum -c SHA256SUMS.txt)
  fi
done
printf 'Snapshot hashes verified; clinical approval NOT inferred.\n'
