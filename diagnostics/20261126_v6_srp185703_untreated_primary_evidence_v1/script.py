#!/usr/bin/env python3
"""Targeted SRP185703 study evidence and V5 clinical-field coverage."""

import csv
import io
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

base = Path(sys.argv[1]) / "26_recount3_external_primary_v1"
target = "SRP185703"

v3 = base / "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"
v5 = base / "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"

def read_archive(path, suffix):
    with zipfile.ZipFile(str(path)) as z:
        matched = [
            n for n in z.namelist()
            if n.endswith("/" + suffix)
        ]
        assert len(matched) == 1
        raw = z.read(matched[0]).decode("utf-8-sig")

    reader = csv.DictReader(io.StringIO(raw), delimiter="\t")
    return reader.fieldnames or [], list(reader)

identifiers = re.compile(
    r"\b(?:SRX|ERX|DRX|SRR|ERR|DRR|GSM|"
    r"SAMN|SAMD|SAMEA)\d+\b", re.I
)

def redact(s):
    s = " ".join(str(s or "").split())
    s = identifiers.sub("[ACCESSION]", s)
    s = re.sub(
        r"\b[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}\b",
        "[EMAIL]", s
    )
    return s[:500] + (" [TRUNCATED]" if len(s) > 500 else "")

print("SRP185703_UNTREATED_PRIMARY_EVIDENCE_V1")

# Historical reviewed publication — one study-level row.
cols, rows = read_archive(
    v3, "results/04_PUBLICATION_REVIEWED_STUDIES.tsv"
)

study_rows = [r for r in rows if r.get("study") == target]
assert len(study_rows) == 1
review = study_rows[0]

print("STUDY_REVIEW_ROWS", len(study_rows))

for field in (
    "cancer",
    "publication_doi",
    "STUDY_EVIDENCE_URL",
    "public_source_url",
    "STUDY_REVIEW_NOTES",
    "public_source_tissue_evidence",
    "public_source_treatment_evidence",
    "public_source_patient_evidence",
    "public_source_scope",
    "review_action_from_publication",
    "review_caveat",
    "collection_procedures",
    "patient_key_basis",
    "required_evidence"
):
    if field in review:
        print("STUDY_EVIDENCE", field, redact(review[field]))
    else:
        print("STUDY_FIELD_ABSENT", field)

# Original V5 administrative primary screen.
fields, rows = read_archive(
    v5, "results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
)

assert len(rows) == 2055
selected = [r for r in rows if r.get("study") == target]
assert len(selected) == 75
assert len({r["library_id"] for r in selected}) == 75

print("\nV5_EXACT_STUDY_ROWS", len(selected))
print("V5_FIELD_COUNT", len(fields))

keywords = re.compile(
    r"treat|therap|pretreat|prior|untreat|"
    r"primar|metasta|recurr|tumou?r|"
    r"specimen|patient|subject|"
    r"sample|biops|patholog|source|"
    r"clinical|tissue|disease|stage|"
    r"collection|isolate|"
    r"resect|surg|assay|protocol",
    re.I
)

relevant = [f for f in fields if keywords.search(f)]

print("RELEVANT_V5_FIELD_COUNT", len(relevant))

for field in relevant:
    values = [
        str(r.get(field) or "").strip()
        for r in selected
    ]
    populated = [x for x in values if x]

    print(
        "V5_FIELD",
        field,
        "POPULATED", len(populated),
        "DISTINCT_NONEMPTY", len(set(populated))
    )

# Only explicitly controlled certification/status fields.
for field in (
    "pretreatment_source_text_status",
    "pretreatment_certified",
    "primary_tumour_certified",
    "patient_independence_certified",
    "tissue_review_status",
    "cohort_approved",
    "technical_hold"
):
    if field not in fields:
        print("STATUS_FIELD_ABSENT", field)
        continue

    counts = Counter(
        str(r.get(field) or "").strip()
        for r in selected
    )

    if all(re.fullmatch(r"[A-Z0-9_/-]{0,100}", k)
           for k in counts):
        print("V5_STATUS", field, repr(dict(sorted(counts.items()))))
    else:
        print("V5_STATUS_UNSAFE_FREE_TEXT", field)

print("\nNO_INDIVIDUAL_EXPERIMENT_OR_PATIENT_VALUES_PRINTED")
print("HISTORICAL_STUDY_REVIEW_NOT_PER_SPECIMEN_PROOF")
print("UNTREATED_PRIMARY_SUBSET_NOT_CERTIFIED")
print("PATIENT_INDEPENDENCE_UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("NO_SOURCE_FILES_MODIFIED")
print("SRP185703_UNTREATED_PRIMARY_EVIDENCE_V1_COMPLETE")
