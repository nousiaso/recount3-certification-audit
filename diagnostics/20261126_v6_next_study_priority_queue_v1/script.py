#!/usr/bin/env python3
"""Read-only schema and study-level census of historical priority queue."""

import csv
import io
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])

archive = root / (
    "26_recount3_external_primary_v1/"
    "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"
)

assert archive.is_file()

print("RECOUNT3_NEXT_STUDY_PRIORITY_QUEUE_V1")

with zipfile.ZipFile(str(archive)) as z:
    matches = [
        x for x in z.namelist()
        if x.endswith(
            "/results/06_REMAINING_STUDIES_PRIORITY_QUEUE.tsv"
        )
    ]
    assert len(matches) == 1

    data = z.read(matches[0]).decode("utf-8-sig")

reader = csv.DictReader(io.StringIO(data), delimiter="\t")
fields = reader.fieldnames or []
rows = list(reader)

assert fields and rows

print("QUEUE_ROWS", len(rows))
print("QUEUE_FIELD_COUNT", len(fields))
print("QUEUE_COLUMNS", ";".join(fields))

study_candidates = [
    x for x in fields
    if x.lower() in (
        "study", "study_id", "study_acc",
        "study_accession", "sra_study"
    )
]

priority_fields = [
    x for x in fields
    if re.search(
        r"priorit|rank|review_status|triage|"
        r"primary.*count|candidate.*count|"
        r"library.*count",
        x, re.I
    )
]

print("EXACT_STUDY_COLUMNS", ";".join(study_candidates))
print("PRIORITY_RELEVANT_COLUMNS", ";".join(priority_fields))

study_col = study_candidates[0] if study_candidates else None

if study_col:
    codes = Counter()

    for r in rows:
        value = str(r.get(study_col) or "").strip().upper()
        if re.fullmatch(r"(?:SRP|ERP|DRP)\d+", value):
            codes[value] += 1

    print("VALID_STUDY_ACCESSIONS", sum(codes.values()))
    print("DISTINCT_STUDIES", len(codes))
    print("DRP005276_QUEUE_ROWS", codes["DRP005276"])

    # Study-level accessions only, original file order.
    first = []
    for r in rows:
        value = str(r.get(study_col) or "").strip().upper()
        if re.fullmatch(r"(?:SRP|ERP|DRP)\d+", value):
            if value not in first:
                first.append(value)
        if len(first) == 15:
            break

    print(
        "FIRST_STUDIES_IN_QUEUE_FILE",
        ";".join(first)
    )
    print("QUEUE_ORDER_NOT_YET_VALIDATED_AS_PRIORITY_RANKING")
else:
    print("STUDY_COLUMN_NOT_IDENTIFIED")
    print("STUDY_SELECTION_DEFERRED_PENDING_SCHEMA_REVIEW")

print("NO_INDIVIDUAL_PATIENT_OR_SAMPLE_ROWS_PRINTED")
print("NO_COHORT_FLAGS_MODIFIED")
print("CLINICAL_ADMISSIONS 0")
print("RECOUNT3_NEXT_STUDY_PRIORITY_QUEUE_V1_COMPLETE")
