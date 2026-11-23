#!/usr/bin/env python3
"""Corrected evidence-only study matrix. No scientific execution."""

import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

source = repo / (
    "diagnostics/20261123_v6_evidence_precedence_dryrun/"
    "proposed_evidence_2055_UNSIGNED.tsv"
)
ledger = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "11_176_STUDY_SOURCE_LEDGER_UNCHANGED_FROM_V4.tsv"
)

out = repo / "diagnostics/20261123_v6_study_review_matrix_v1_1"

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def read(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

assert sha(source) == (
    "08efe2c8f469edcbdfb98231dc87f705d801434e78"
    "dcbea85f57c4b94da88eb0"
)
assert sha(ledger) == (
    "b9e1ae8d91307cc47c5411f1b8b201db2a6cfe92"
    "a8db0d0f7d79917abaef2db7"
)

data = read(source)
studies = read(ledger)

study_ids = [r["study"] for r in studies]
assert len(studies) == len(set(study_ids)) == 176
assert len(data) == 2055
assert len({r["library_id"] for r in data}) == 2055
assert all(r["FINAL_APPROVAL"] == "NO_UNSIGNED" for r in data)

by = defaultdict(list)
for r in data:
    by[r["study"]].append(r)

# Count contributors BEFORE the 176-study reporting loop.
screen_studies = set(by)
assert screen_studies <= set(study_ids)
assert len(screen_studies) == 77
assert len(set(study_ids) - screen_studies) == 99

H = "HISTORICAL_EVIDENCE_REVIEW_PENDING"
N = "EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE"
P = "PRIMARY_TUMOUR_SPECIMEN_ONLY_OTHER_GATES_OPEN"
S = "HOLD_SAMPLED_SITE_METASTASIS_DISEASE_FLAG"
V = "V5_TUMOUR_OR_CONTEXT_SPECIMEN_GATES_OPEN"

output = []

for study in study_ids:
    records = by.get(study, [])
    buckets = Counter(r["unsigned_evidence_bucket"] for r in records)

    historical = [
        r for r in records
        if r["unsigned_evidence_bucket"] == H
    ]
    grades = Counter(r["historical_tissue_grade"] for r in historical)

    row = {
        "study": study,
        "n_admin_screen": len(records),
        "screen_contributor": "YES" if records else "NO",
        "n_non_tumour": buckets[N],
        "n_exact_primary_other_gates_open": buckets[P],
        "n_sarcoma_site_hold": buckets[S],
        "n_v5_context_other_gates_open": buckets[V],
        "n_historical_only": buckets[H],
        "n_historical_primary_unattested":
            grades["PRIMARY_TUMOUR_NOT_ATTESTED"],
        "n_historical_surgery_unconfirmed":
            grades["PRIMARY_SURGERY_CONTEXT_NOT_PATHOLOGY_CERTIFIED"],
        "n_historical_study_source_supported":
            grades["PRIMARY_TUMOUR_STUDY_SOURCE_SUPPORTED"],
        "n_pretreatment_not_established": sum(
            r["treatment_grade_UNREVIEWED"] ==
            "PRETREATMENT_NOT_ESTABLISHED"
            for r in records
        ),
        "n_curated_patient_key_unverified": sum(
            r["patient_grade_UNREVIEWED"] ==
            "CURATED_STUDY_PATIENT_KEY_UNVERIFIED"
            for r in records
        ),
        "n_biosample_proxy": sum(
            r["patient_grade_UNREVIEWED"] ==
            "BIOSAMPLE_PROXY_NOT_PATIENT_ID"
            for r in records
        ),
        "n_approved": 0,
    }
    output.append(row)

def total(field):
    return sum(r[field] for r in output)

assert total("n_admin_screen") == 2055
assert total("n_non_tumour") == 22
assert total("n_exact_primary_other_gates_open") == 5
assert total("n_sarcoma_site_hold") == 52
assert total("n_v5_context_other_gates_open") == 280
assert total("n_historical_only") == 1696

assert total("n_historical_primary_unattested") == 599
assert total("n_historical_surgery_unconfirmed") == 603
assert total("n_historical_study_source_supported") == 494

output.sort(key=lambda r: (
    -r["n_historical_primary_unattested"],
    -r["n_historical_surgery_unconfirmed"],
    -r["n_sarcoma_site_hold"],
    -r["n_admin_screen"],
    r["study"],
))

out.mkdir(parents=True, exist_ok=True)
target = out / "study_review_matrix_176_UNSIGNED.tsv"

with target.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(output[0]), delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(output)

summary = {
    "status": "UNSIGNED_EVIDENCE_ONLY",
    "n_inventory_studies": 176,
    "n_screen_studies": len(screen_studies),
    "n_other_studies": 99,
    "n_screen_libraries": 2055,
    "n_historical_unattested": 599,
    "n_historical_surgery_unconfirmed": 603,
    "n_historical_study_source_supported": 494,
    "n_v5_context_other_gates_open": 280,
    "clinical_admissions": 0,
    "source_sha256": sha(source),
    "matrix_sha256": sha(target),
}
(out / "summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)

print("CORRECTED_STUDY_MATRIX_PASS")
print("Inventory studies:", 176)
print("Screen contributors:", len(screen_studies))
print("Outside primary screen:", 99)
print("Historical unattested without V5 class:", 599)
print("Historical surgery context without V5 class:", 603)
print("Historical study-level support without V5 class:", 494)

print("\nTOP 20 STUDIES BY HISTORICAL PRIMARY ATTESTATION GAP")
print("study\tscreen_N\tunattested_no_V5\tsurgery_no_V5\tV5_context")

for r in output[:20]:
    print(
        r["study"],
        r["n_admin_screen"],
        r["n_historical_primary_unattested"],
        r["n_historical_surgery_unconfirmed"],
        r["n_v5_context_other_gates_open"],
        sep="\t"
    )

print("\nMATRIX_SHA256:", sha(target))
print("NO CLINICAL ADMISSIONS; ORIGINAL FILES UNCHANGED")
