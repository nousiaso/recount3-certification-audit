#!/usr/bin/env bash
# Acquire the original G026 source from the official recount3 public S3 mirror.
# Accept only after all 22 bridge runs and all 63856 locked gene rows verify.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
DIR="$ROOT/26_recount3_external_primary_v1/raw_studies/SRP157215/gene_sums"
RAW="$DIR/sra.gene_sums.SRP157215.G026.gz"
URL="https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/sra/gene_sums/15/SRP157215/sra.gene_sums.SRP157215.G026.gz"
VERIFY="$REPO/tools/v7q_verify_SRP157215_G026_source_READONLY_v1.py"
AUDIT="$REPO/diagnostics/20261127_v7q_srp157215_source_acquisition_v1"
test -x "$PYTHON"
test -f "$VERIFY"
if [[ -e "$RAW" || -e "$AUDIT" ]]; then
    echo "REFUSE_EXISTING_RAW_OR_AUDIT" >&2
    exit 1
fi
mkdir -p "$DIR"
mkdir "$AUDIT"
TMP=$(mktemp "$DIR/.SRP157215.G026.XXXXXXXX.part")
trap 'rm -f -- "$TMP"' EXIT
printf '%s\n' "$URL" > "$AUDIT/00_source_url.txt"
curl -fL --retry 3 --connect-timeout 30 --max-time 1800 --output "$TMP" "$URL"
gzip -t "$TMP"
"$PYTHON" "$VERIFY" "$ROOT" "$TMP" | tee "$AUDIT/01_source_verification_stdout.txt"
grep -Fxq V7Q_SRP157215_RAW_SOURCE_AXIS_PASS "$AUDIT/01_source_verification_stdout.txt"
# Atomic no-overwrite publication on the same filesystem.
ln -- "$TMP" "$RAW"
(
    cd "$AUDIT"
    sha256sum "$VERIFY" "$REPO/scripts/run_v7q_acquire_SRP157215_G026_v1.sh" "$RAW" 00_source_url.txt 01_source_verification_stdout.txt > SHA256SUMS.txt
    sha256sum -c SHA256SUMS.txt
)
echo "V7Q_SRP157215_RAW_SOURCE_ACQUIRED_PENDING_REVIEW"
echo "NO_MODEL_EXECUTION_CLINICAL_ADMISSIONS_0"
