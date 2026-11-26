#!/usr/bin/env python3
"""SRP185703 membership in V5 primary screen and V3 study reviews."""

import csv
import hashlib
import io
import re
import sys
import zipfile
from pathlib import Path

root = Path(sys.argv[1])
base = root / "26_recount3_external_primary_v1"
target = "SRP185703"

source = root / (
    "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/"
    "01_full_source_experiment_identity_5221.tsv"
)

v5 = base / "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"
v3 = base / "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def zip_table(archive, suffix):
    with zipfile.ZipFile(str(archive)) as z:
        names = [
            n for n in z.namelist()
            if n.endswith("/" + suffix)
        ]
        assert len(names) == 1, (archive.name, suffix)
        raw = z.read(names[0])

    reader = csv.DictReader(
        io.StringIO(raw.decode("utf-8-sig")),
        delimiter="\t"
    )
    return reader.fieldnames or [], list(reader)

def safe(value):
    value = str(value or "").strip()
    if re.fullmatch(r"[A-Z0-9_./-]{1,90}", value):
        return value
    return "[FREE_TEXT_NOT_PRINTED]" if value else "[EMPTY]"

assert sha256(source) == (
    "95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa"
)

with source.open(encoding="utf-8-sig", newline="") as f:
    original = list(csv.DictReader(f, delimiter="\t"))

source_ids = {
    r["experiment_acc"] for r in original
    if r["study"] == target
}

assert len(source_ids) == 75
assert all(re.fullmatch(r"(?:SRX|ERX|DRX)\d+", x) for x in source_ids)

print("SRP185703_PRIMARY2055_REVIEW_RECONCILIATION_V1")
print("SOURCE_EXPERIMENTS_VERIFIED", len(source_ids))
print("V5_ARCHIVE_SHA256", sha256(v5))
print("V3_ARCHIVE_SHA256", sha256(v3))

# --------------------------------------------
# Actual V5 primary-screen membership
# --------------------------------------------

fields, rows = zip_table(
    v5,
    "results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
)

print("\nV5_PRIMARY_SCREEN")
print("ROWS", len(rows))
print("FIELD_COUNT", len(fields))
assert len(rows) == 2055

study_fields = [
    x for x in fields
    if x.lower() in (
        "study", "study_id", "study_accession",
        "sra_study", "study_acc"
    )
]
print("STUDY_FIELD_CANDIDATES", ";".join(study_fields) or "[NONE]")

for field in study_fields:
    print(
        "STUDY_ROWS", field,
        sum(str(r.get(field) or "").strip() == target for r in rows)
    )

experiment_pattern = re.compile(r"(?:SRX|ERX|DRX)\d+", re.I)

id_columns = []

for field in fields:
    values = [
        str(r.get(field) or "").strip().upper()
        for r in rows
    ]
    valid = [x for x in values if experiment_pattern.fullmatch(x)]

    if valid:
        hits = set(valid) & source_ids

        id_columns.append(field)
        print(
            "EXPERIMENT_IDENTIFIER_FIELD", field,
            "VALID_ROWS", len(valid),
            "DISTINCT", len(set(valid)),
            "SRP185703_EXACT_INTERSECTION", len(hits)
        )

print("EXPERIMENT_ID_COLUMNS_DETECTED", len(id_columns))

if not id_columns:
    print("EXACT_EXPERIMENT_MEMBERSHIP NOT_EVALUABLE_NO_ID_COLUMN")
else:
    all_v5_ids = set()

    for row in rows:
        for field in id_columns:
            value = str(row.get(field) or "").strip().upper()
            if experiment_pattern.fullmatch(value):
                all_v5_ids.add(value)

    print(
        "SRP185703_EXACT_MEMBERSHIP_IN_PRIMARY2055",
        len(source_ids & all_v5_ids)
    )
    print(
        "SRP185703_NOT_IN_PRIMARY2055",
        len(source_ids - all_v5_ids)
    )

# --------------------------------------------
# V3 reviewed-study inventory
# --------------------------------------------

print("\nV3_HISTORICAL_STUDY_REVIEW")

for suffix in (
    "results/02_ALL_176_STUDY_SOURCE_REVIEW_UNSIGNED.tsv",
    "results/04_PUBLICATION_REVIEWED_STUDIES.tsv"
):
    columns, records = zip_table(v3, suffix)

    assert "study" in columns

    matches = [
        r for r in records if r.get("study") == target
    ]

    print(
        "REVIEW_TABLE", Path(suffix).name,
        "TOTAL_STUDIES", len(records),
        "SRP185703_ROWS", len(matches)
    )

    assert len(matches) <= 1

    if not matches:
        continue

    row = matches[0]

    for field in (
        "libraries",
        "patients",
        "qa_stage41_primary_screen",
        "qa_missing_pretreatment",
        "qa_clinical_holds",
        "qa_strict_atlas_overlap",
        "study_in_archs4_gse_list",
        "study_level_source_verification",
        "STUDY_REVIEW_DECISION",
        "overall_study_disposition",
        "public_source_review_status",
        "review_action_from_publication",
        "FINAL_STUDY_CLINICAL_APPROVAL_AS_OF_V3"
    ):
        if field in row:
            print(
                "REVIEW_FIELD", field,
                safe(row.get(field))
            )

    for field in (
        "public_source_tissue_evidence",
        "public_source_treatment_evidence",
        "public_source_patient_evidence",
        "review_caveat",
        "required_evidence"
    ):
        if field in row:
            value = str(row.get(field) or "")
            print(
                "REVIEW_TEXT", field,
                "POPULATED", bool(value.strip()),
                "LENGTH", len(value),
                "SHA256_PREFIX",
                hashlib.sha256(
                    value.encode("utf-8")
                ).hexdigest()[:16]
            )

print("\nINTERPRETATION_LIMITS")
print("PRIMARY_SCREEN_IS_ADMINISTRATIVE_NOT_ADMISSION")
print("HISTORICAL_REVIEW_NOT_NEW_CLINICAL_ATTESTATION")
print("NO_PATIENT_OR_EXPERIMENT_VALUES_PRINTED")
print("NO_FROZEN_DATA_MODIFIED")
print("PATIENT_INDEPENDENCE_UNSIGNED")
print("PRETREATMENT_UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("SRP185703_PRIMARY2055_REVIEW_RECONCILIATION_V1_COMPLETE")
