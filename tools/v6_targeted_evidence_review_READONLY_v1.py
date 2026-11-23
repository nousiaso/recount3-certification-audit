#!/usr/bin/env python3
"""Review existing V6 grades and SRP009386 source; no mutations."""

import csv
import gzip
import hashlib
import sys
from pathlib import Path
from collections import Counter

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

p = repo / (
    "diagnostics/20261123_v6_evidence_precedence_dryrun/"
    "proposed_evidence_2055_UNSIGNED.tsv"
)
expected = (
    "08efe2c8f469edcbdfb98231dc87f705d801434e78"
    "dcbea85f57c4b94da88eb0"
)
assert hashlib.sha256(p.read_bytes()).hexdigest() == expected

with p.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 2055

for bucket in (
    "HISTORICAL_EVIDENCE_REVIEW_PENDING",
    "V5_TUMOUR_OR_CONTEXT_SPECIMEN_GATES_OPEN",
):
    subset = [
        r for r in rows
        if r["unsigned_evidence_bucket"] == bucket
    ]
    print("\nBUCKET", bucket, "N", len(subset))

    for field in (
        "historical_tissue_grade",
        "treatment_grade_UNREVIEWED",
        "patient_grade_UNREVIEWED",
    ):
        print("FIELD", field)
        for label, n in Counter(
            r[field] for r in subset
        ).most_common():
            print(" ", n, repr(label))

    print("TOP_STUDIES_BY_RECORD_COUNT")
    for study, n in Counter(
        r["study"] for r in subset
    ).most_common(20):
        print(" ", study, n)

print("\nSRP009386_EXISTING_V6_RECORDS")
for r in rows:
    if r["study"] == "SRP009386":
        print(
            r["library_id"],
            r["unsigned_evidence_bucket"],
            r["historical_tissue_grade"],
            r["patient_grade_UNREVIEWED"],
            sep="\t",
        )

target = {"SRX106046", "SRX106048"}
metadata = (
    rec / "raw_studies/SRP009386/metadata"
)
files = (
    sorted(metadata.glob("*.MD.gz"))
    if metadata.is_dir() else []
)

print("\nSRP009386_RAW_SOURCE_FILES", len(files))

found = set()
for path in files:
    with gzip.open(
        path, "rt", errors="replace", newline=""
    ) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "experiment_acc" not in (reader.fieldnames or []):
            print("UNSUPPORTED_SOURCE_SCHEMA", path.name)
            continue

        for row in reader:
            acc = row.get("experiment_acc", "")
            if acc not in target:
                continue
            found.add(acc)
            print("\nSOURCE_EXPERIMENT", acc)
            print("SOURCE_FILE", path.name)
            for key in (
                "sample_title",
                "sample_attributes",
                "sample_acc",
                "biosample",
            ):
                if key in row:
                    print(key, repr(row[key][:2500]))

print("\nSOURCE_MATCHED", sorted(found))
print("SOURCE_NOT_FOUND", sorted(target - found))
print("READ_ONLY_COMPLETE")
