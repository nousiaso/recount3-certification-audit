#!/usr/bin/env python3
"""Historical V6 local provenance verification; Python 3.6 compatible.

Reads only existing files. Does not expose source records or modify inputs.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

print("V6_HISTORICAL_LOCAL_PROVENANCE_V1")
print("PYTHON_VERSION", sys.version.split()[0])

FILES = {
    "MATRIX_V1_SUPERSEDED": (
        "diagnostics/20261123_v6_study_review_matrix/"
        "study_review_matrix_176_UNSIGNED.tsv",
        "9b51f2ddf892e6e2d0d5ad3c502931f68f95ca1feec049fac32c5428853d188e",
        13565
    ),
    "HISTORICAL1696_STDOUT": (
        "diagnostics/20261124_v6_historical1696_integrated_v1/stdout.txt",
        "48fbf38938c2070e85f4c214bfbdd5fcbcbee0f675902d57e66d64fb2c2f4f96",
        692
    ),
    "HISTORICAL1696_SUMMARY": (
        "diagnostics/20261124_v6_historical1696_integrated_v1/summary.json",
        "3b3ebcbf16584a112e807b83034532bae8a390374084fe63d0dfce2ecd399d0e",
        1064
    ),
    "PROSTATE_GEO_SUMMARY": (
        "diagnostics/20261124_v6_prostate_geo_exact_v1/summary.json",
        "c7f643a2989486aaadb1b9643972e5b7c13584ca0581561a49453b7657dd38f3",
        768
    ),
    "PROSTATE_EXCEPTION_STDOUT": (
        "diagnostics/20261124_v6_prostate_patient_exception_v1/stdout.txt",
        "3b22f5ff12bf827b2a802ac8191aaabe55a612859f45d256c593b1f97c5a9127",
        378
    ),
    "PROSTATE_EXCEPTION_SUMMARY": (
        "diagnostics/20261124_v6_prostate_patient_exception_v1/summary.json",
        "307296075e21f6be8ccc3e2b852d01dc186452acaedf7d534f29b7ec64895900",
        753
    ),
    "HISTORICAL1696_SCRIPT": (
        "tools/v6_historical1696_integrated_evidence_v1.py",
        "5a4c4cbd46f746521aa0be054fada04cc426667ee9114af2531d25770830bb3a",
        8885
    ),
    "PROSTATE_EXCEPTION_SCRIPT": (
        "tools/v6_prostate_patient_exception_v1.py",
        "a01efc1955e1ce5d2e89f5f2ed5f28011b10a97824ed3b85d050ec2fe5e13932",
        5926
    ),
    "STRAY_UNREVIEWED_FILE": (
        "tatus --short",
        "7b1c5059598aaeaada3d2e5a1c4da7c48207feb4521a79f7d3cb0167eac461b0",
        31841
    ),
    "SRP218935_LAUNCHER": (
        "SRP218935_certify_archive_and_push_20261125.sh",
        "1aefe91f6616fd2d96b99f199156ee6f6fcd1297e7bdd7cb6c9eab19f3bdeca3",
        14788
    )
}

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

print("EXPECTED_FILES", len(FILES))

for label, (filename, expected_sha, expected_bytes) in FILES.items():
    path = Path(filename)

    assert path.is_file(), "Missing file group: " + label
    assert path.stat().st_size == expected_bytes, (
        "Byte count changed: " + label
    )
    assert sha256(path) == expected_sha, (
        "SHA256 mismatch: " + label
    )

    if path.suffix == ".json":
        with path.open(encoding="utf-8-sig") as f:
            obj = json.load(f)
        assert isinstance(obj, dict)

    if path.suffix == ".py":
        with path.open(encoding="utf-8-sig") as f:
            source = f.read()
        compile(source, filename, "exec")

    if label == "MATRIX_V1_SUPERSEDED":
        with path.open(encoding="utf-8-sig", newline="") as f:
            records = list(csv.DictReader(f, delimiter="\t"))
        assert len(records) == 176

    print("VERIFIED", label, "SHA256_PASS", expected_sha)

print("VERIFIED_FILE_COUNT", len(FILES))
print("HISTORICAL_1696_HASHES PASS")
print("PROSTATE_EXCEPTION_HASHES PASS")
print("MATRIX_V1 SUPERSEDED_NOT_RELEASED")
print("PROSTATE_PATIENT_CONFLICT UNRESOLVED")
print("EXISTING_CLINICAL_APPROVALS 0")
print("NO_ORIGINAL_FILES_MODIFIED")
print("NO_PATIENT_OR_SAMPLE_VALUES_PRINTED")
print("V6_HISTORICAL_LOCAL_PROVENANCE_V1_COMPLETE")
