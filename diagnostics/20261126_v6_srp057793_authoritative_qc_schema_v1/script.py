#!/usr/bin/env python3
"""Read-only exact accession join to original Stage2 QC resources."""

import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
rec = root / "26_recount3_external_primary_v1"
s41 = root / (
    "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
)

def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return reader.fieldnames or [], list(reader)

_, source = read_tsv(
    s41 / "01_full_source_experiment_identity_5221.tsv"
)
_, candidates = read_tsv(
    s41 / "03_candidate_integrated_evidence_3602.tsv"
)

source = [r for r in source if r.get("study") == "SRP057793"]
candidates = [
    r for r in candidates if r.get("study") == "SRP057793"
]

assert len(source) == 149
assert len(candidates) == 148

source_ids = {r["experiment_acc"] for r in source}
candidate_ids = {r["library_id"] for r in candidates}

held_ids = {
    r["library_id"] for r in candidates
    if r.get("technical_hold") == "YES"
}
retained_ids = {
    r["library_id"] for r in candidates
    if r.get("technical_hold") == "NO"
}
source_only = source_ids - candidate_ids

assert len(source_ids) == 149
assert len(candidate_ids) == 148
assert len(held_ids) == 8
assert len(retained_ids) == 140
assert len(source_only) == 1

files = {
    "LIBRARY_QC_ALL": rec / "stage2_output/library_qc_all.tsv",
    "RUN_QC_ALL": rec / "stage2_output/run_qc_all.tsv",
    "STUDY_RUN_QC": rec / (
        "stage2_output/studies/SRP057793/run_qc.tsv"
    ),
    "RUN_HOLD_QUEUE": rec / (
        "cohort_readiness_v1/source_snapshot/05_run_hold_queue_409.tsv"
    ),
}

reason_rx = re.compile(
    r"reason|hold|filter|outlier|technical|"
    r"qc_flag|qc_status|failure|excluded",
    re.I
)

print("SRP057793_AUTHORITATIVE_QC_SCHEMA_V1")
print("SOURCE_EXPERIMENTS", len(source_ids))
print("HELD_EXPERIMENTS", len(held_ids))
print("RETAINED_EXPERIMENTS", len(retained_ids))
print("SOURCE_ONLY_EXPERIMENTS", len(source_only))

for label, path in files.items():
    print("\nFILE", label)
    print("PRESENT", path.is_file())

    if not path.is_file():
        continue

    fields, rows = read_tsv(path)

    print("ROWS", len(rows))
    print("FIELDS", len(fields))
    print("COLUMN_NAMES", ";".join(fields[:100]))

    matched = {}
    all_found = set()

    for field in fields:
        values = {
            str(r.get(field) or "").strip()
            for r in rows
        }

        common = values & source_ids
        if common:
            matched[field] = common
            all_found.update(common)

            print(
                "EXACT_IDENTIFIER_COLUMN", field,
                "SOURCE", len(common),
                "HELD", len(common & held_ids),
                "RETAINED", len(common & retained_ids),
                "SOURCE_ONLY", len(common & source_only)
            )

    print("EXACT_ID_COLUMNS", len(matched))
    print("TOTAL_SOURCE_ID_COVERAGE", len(all_found))
    print("TOTAL_HELD_ID_COVERAGE", len(all_found & held_ids))
    print(
        "TOTAL_SOURCE_ONLY_ID_COVERAGE",
        len(all_found & source_only)
    )

    reason_fields = [
        field for field in fields
        if reason_rx.search(field)
    ]

    print("QC_REASON_FIELD_CANDIDATES", ";".join(reason_fields[:50]))

    for field in reason_fields:
        held_count = 0
        retained_count = 0

        for row in rows:
            value = str(row.get(field) or "").strip()
            if not value:
                continue

            row_ids = {
                str(row.get(col) or "").strip()
                for col in matched
            }

            if row_ids & held_ids:
                held_count += 1
            if row_ids & retained_ids:
                retained_count += 1

        print(
            "QC_FIELD_COVERAGE",
            field,
            "HELD_MATCHED_ROWS", held_count,
            "RETAINED_MATCHED_ROWS", retained_count
        )

# Read the study-specific JSON schema, not its sensitive values.
study_json = rec / (
    "stage2_output/studies/SRP057793/study_qc.json"
)

print("\nSTUDY_QC_JSON_PRESENT", study_json.is_file())

if study_json.is_file():
    with study_json.open(encoding="utf-8") as f:
        obj = json.load(f)

    print("JSON_ROOT_TYPE", type(obj).__name__)
    if isinstance(obj, dict):
        print("JSON_TOP_LEVEL_KEYS", ";".join(sorted(obj)))

print("\nNO_INDIVIDUAL_ACCESSIONS_PRINTED")
print("QC_REASONS_NOT_ASSIGNED_WITHOUT_SOURCE_JOIN")
print("NO_FROZEN_FILES_MODIFIED")
print("PATIENT_INDEPENDENCE_UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("SRP057793_AUTHORITATIVE_QC_SCHEMA_V1_COMPLETE")
