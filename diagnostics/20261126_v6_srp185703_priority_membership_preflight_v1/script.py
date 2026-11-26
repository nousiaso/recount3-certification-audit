#!/usr/bin/env python3
"""Read-only SRP185703 priority and source membership preflight."""

import csv
import hashlib
import io
import re
import sys
import zipfile
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
target = "SRP185703"

archive = root / (
    "26_recount3_external_primary_v1/"
    "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"
)

files = {
    "CANDIDATE": root / (
        "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/"
        "03_candidate_integrated_evidence_3602.tsv"
    ),
    "FULL_SOURCE": root / (
        "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/"
        "01_full_source_experiment_identity_5221.tsv"
    ),
    "SELECTED": repo / (
        "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
        "source_evidence_599_LOCAL_ONLY.tsv"
    )
}

hashes = {
    "CANDIDATE": "aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81",
    "FULL_SOURCE": "95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa",
    "SELECTED": "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb"
}

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def load(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return reader.fieldnames or [], list(reader)

def numeric(value):
    try:
        v = Decimal(str(value or "").strip())
        return v if v.is_finite() else None
    except InvalidOperation:
        return None

def safe(value):
    value = str(value or "").strip()
    if not value:
        return "[EMPTY]"
    if re.fullmatch(r"[A-Z0-9_/-]{1,90}", value):
        return value
    return "[OTHER_CATEGORY]"

print("SRP185703_PRIORITY_MEMBERSHIP_PREFLIGHT_V1")

assert archive.is_file()
print("QUEUE_ARCHIVE_SHA256", sha256(archive))

with zipfile.ZipFile(str(archive)) as z:
    matching = [
        x for x in z.namelist()
        if x.endswith(
            "/results/06_REMAINING_STUDIES_PRIORITY_QUEUE.tsv"
        )
    ]
    assert len(matching) == 1
    text = z.read(matching[0]).decode("utf-8-sig")

reader = csv.DictReader(io.StringIO(text), delimiter="\t")
queue = list(reader)
assert len(queue) == 147

target_rows = [
    r for r in queue if r.get("study") == target
]
assert len(target_rows) <= 1

print("STUDY", target)
print("IN_V3_REMAINING_QUEUE", len(target_rows) == 1)

if target_rows:
    row = target_rows[0]

    count_fields = (
        "libraries",
        "runs",
        "patients",
        "review_priority_score",
        "qa_priority_score",
        "qa_stage41_primary_screen",
        "qa_scientific_priority_review",
        "qa_clinical_holds",
        "qa_missing_pretreatment",
        "core_explicit_hold_libraries",
        "separate_context_libraries",
        "provisional_potential_core_libraries",
        "exact_archs4_libraries",
        "on_or_after_treatment_libraries"
    )

    for name in count_fields:
        value = numeric(row.get(name))
        print(
            "QUEUE_NUMERIC", name,
            str(value) if value is not None else "NOT_NUMERIC"
        )

    for name in (
        "STUDY_REVIEW_DECISION",
        "overall_study_disposition",
        "qa_action",
        "public_source_review_status",
        "review_action_from_publication",
        "FINAL_STUDY_CLINICAL_APPROVAL_AS_OF_V3"
    ):
        print(
            "QUEUE_STATUS", name,
            safe(row.get(name))
        )

    for score_name in (
        "review_priority_score",
        "qa_priority_score"
    ):
        current = numeric(row.get(score_name))
        comparable = [
            numeric(r.get(score_name))
            for r in queue
        ]
        comparable = [
            x for x in comparable if x is not None
        ]

        if current is not None:
            descending_rank = 1 + sum(
                x > current for x in comparable
            )
            print(
                "SCORE_DIAGNOSTIC", score_name,
                "DESCENDING_RANK", descending_rank,
                "NUMERIC_STUDIES", len(comparable),
                "RANK_IS_NOT_CLINICAL_APPROVAL"
            )

for label, path in files.items():
    assert path.is_file()
    assert sha256(path) == hashes[label], label

cf, candidates = load(files["CANDIDATE"])
ff, full = load(files["FULL_SOURCE"])
sf, selected = load(files["SELECTED"])

assert (len(candidates), len(full), len(selected)) == (
    3602, 5221, 599
)

cr = [r for r in candidates if r["study"] == target]
fr = [r for r in full if r["study"] == target]
sr = [r for r in selected if r["study"] == target]

def unique(rows, field):
    return {
        r[field] for r in rows if r.get(field)
    }

c_ids = unique(cr, "library_id")
f_ids = unique(fr, "experiment_acc")
s_ids = unique(sr, "experiment_acc")

print("\nSTAGE41_MEMBERSHIP")
print("CANDIDATE_ROWS", len(cr))
print("FULL_SOURCE_ROWS", len(fr))
print("HISTORICAL_SELECTED_ROWS", len(sr))
print("CANDIDATE_DISTINCT_IDS", len(c_ids))
print("FULL_SOURCE_DISTINCT_IDS", len(f_ids))
print("SELECTED_DISTINCT_IDS", len(s_ids))
print("CANDIDATE_NOT_IN_FULL_SOURCE", len(c_ids - f_ids))
print("SELECTED_NOT_IN_CANDIDATE", len(s_ids - c_ids))
print("FULL_SOURCE_NOT_IN_CANDIDATE", len(f_ids - c_ids))

print("\nSTAGE41_CANDIDATE_STATUS")

for name in (
    "technical_hold",
    "stage41_partition_NOT_ADMISSION",
    "tissue_review_status",
    "pretreatment_source_text_status",
    "patient_independence_certified",
    "primary_tumour_certified",
    "pretreatment_certified",
    "atlas_independence_certified",
    "cohort_approved"
):
    if name not in cf:
        print("FIELD_MISSING", name)
        continue

    counts = Counter(
        safe(r.get(name)) for r in cr
    )

    print(
        "FIELD", name,
        "COUNTS", repr(dict(sorted(counts.items())))
    )

print("\nGATES")
print("NO_PATIENT_OR_EXPERIMENT_VALUES_PRINTED")
print("NO_COHORT_DECISIONS_MODIFIED")
print("NO_FROZEN_FILES_MODIFIED")
print("CLINICAL_ADMISSIONS 0")
print("SRP185703_PRIORITY_MEMBERSHIP_PREFLIGHT_V1_COMPLETE")
