#!/usr/bin/env python3
"""Locate original DRP005276 study metadata without reading patient rows."""

import hashlib
import os
import re
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"

locations = [
    ("RECOUNT3", root / "26_recount3_external_primary_v1"),
    ("STAGE41", root / "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"),
    ("AUDIT_EVIDENCE", repo / "evidence"),
    ("AUDIT_PRIVATE", repo / "private_data"),
]

strong = re.compile(
    r"DRP005276|DRA008030|E[-_]?GEAD[-_]?306|Kawahara",
    re.I
)

clinical = re.compile(
    r"clinical|pretreat|treatment|therapy|biosample|"
    r"sample_attributes|sdrf|idf",
    re.I
)

metadata_extensions = (
    ".xml", ".json", ".csv", ".tsv",
    ".xlsx", ".txt", ".md", ".gz", ".zip"
)

def safe_name(name):
    return re.sub(
        r"(?:SRX|ERX|DRX|SRR|ERR|DRR|GSM|SAMN|SAMD|SAMEA)\d+",
        "[ACCESSION]", name, flags=re.I
    )[:160]

print("DRP005276_ORIGINAL_METADATA_LOCATOR_V1")
total = 0
strong_total = 0
clinical_total = 0

for label, base in locations:
    if not base.is_dir():
        print("SOURCE_ROOT_MISSING", label)
        continue

    files_scanned = 0
    matched = []
    type_counts = Counter()

    for parent, dirs, filenames in os.walk(str(base)):
        depth = len(Path(parent).relative_to(base).parts)
        dirs[:] = sorted(
            d for d in dirs
            if d not in (".git", "__pycache__")
            and depth < 5
        )

        for name in sorted(filenames):
            files_scanned += 1

            if not name.lower().endswith(metadata_extensions):
                continue

            type_counts["METADATA_FILE"] += 1

            if strong.search(name):
                strong_total += 1
                matched.append(("STUDY", name))
            elif clinical.search(name):
                clinical_total += 1
                matched.append(("CLINICAL_NAME", name))

            if files_scanned >= 15000:
                break

        if files_scanned >= 15000:
            break

    total += files_scanned

    print(
        "ROOT", label,
        "FILES_VISITED", files_scanned,
        "METADATA_FILES", type_counts["METADATA_FILE"],
        "NAME_MATCHES", len(matched),
        "CAPPED", files_scanned >= 15000
    )

    for kind, name in matched[:25]:
        print(
            "POSSIBLE_METADATA", label,
            kind, safe_name(name)
        )

    print(
        "MATCHES_NOT_PRINTED",
        label, max(0, len(matched) - 25)
    )

print("TOTAL_FILES_VISITED", total)
print("STUDY_FILENAME_MATCHES", strong_total)
print("CLINICAL_FILENAME_MATCHES", clinical_total)
print("NO_METADATA_CONTENTS_READ")
print("NO_PATIENT_RECORDS_PRINTED")
print("NO_SOURCE_FILES_MODIFIED")
print("CLINICAL_ADMISSIONS 0")
print("DRP005276_ORIGINAL_METADATA_LOCATOR_V1_COMPLETE")
