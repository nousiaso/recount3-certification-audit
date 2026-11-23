#!/usr/bin/env python3
"""Build study-level evidence review matrix; no clinical decisions."""
import csv
import hashlib
import sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"
out = repo / "diagnostics/20261123_v6_study_review_matrix"

source = repo / (
    "diagnostics/20261123_v6_evidence_precedence_dryrun/"
    "proposed_evidence_2055_UNSIGNED.tsv"
)
ledger = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "11_176_STUDY_SOURCE_LEDGER_UNCHANGED_FROM_V4.tsv"
)

def digest(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()

assert digest(source) == (
    "08efe2c8f469edcbdfb98231dc87f705d801434e78"
    "dcbea85f57c4b94da88eb0"
), "V6 evidence file changed"

assert digest(ledger) == (
    "b9e1ae8d91307cc47c5411f1b8b201db2a6cfe92"
    "a8db0d0f7d79917abaef2db7"
), "176-study ledger changed"

def read(p):
    with p.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)

data = read(source)
studies = read(ledger)

assert len(data) == 2055
assert len(studies) == 176

study_ids = [r["study"] for r in studies]
assert len(set(study_ids)) == 176

by = defaultdict(list)
for r in data:
    assert r["study"] in study_ids
    assert r["FINAL_APPROVAL"] == "NO_UNSIGNED"
    by[r["study"]].append(r)

assert len(by) == 77
assert len({r["library_id"] for r in data}) == 2055

buckets = {
    "non_tumour": "EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE",
    "primary_exact": "PRIMARY_TUMOUR_SPECIMEN_ONLY_OTHER_GATES_OPEN",
    "site_hold": "HOLD_SAMPLED_SITE_METASTASIS_DISEASE_FLAG",
    "v5_context": "V5_TUMOUR_OR_CONTEXT_SPECIMEN_GATES_OPEN",
    "historical": "HISTORICAL_EVIDENCE_REVIEW_PENDING",
}
grades = {
    "primary_not_attested": "PRIMARY_TUMOUR_NOT_ATTESTED",
    "surgery_not_certified":
        "PRIMARY_SURGERY_CONTEXT_NOT_PATHOLOGY_CERTIFIED",
    "study_primary_supported":
        "PRIMARY_TUMOUR_STUDY_SOURCE_SUPPORTED",
}
treatment = {
    "pretreatment_not_established": "PRETREATMENT_NOT_ESTABLISHED",
    "surgery_no_pretreatment_claim":
        "AT_SURGERY_NO_PRETREATMENT_CLAIM",
}
patient = {
    "curated_patient_key_unverified":
        "CURATED_STUDY_PATIENT_KEY_UNVERIFIED",
    "biosample_proxy": "BIOSAMPLE_PROXY_NOT_PATIENT_ID",
    "stage41_patient_proposal":
        "STAGE41_PATIENT_PROPOSAL_NOT_VERIFIED",
}

output = []
for study in study_ids:
    sub = by[study]
    bc = Counter(r["unsigned_evidence_bucket"] for r in sub)
    gc = Counter(r["historical_tissue_grade"] for r in sub)
    tc = Counter(r["treatment_grade_UNREVIEWED"] for r in sub)
    pc = Counter(r["patient_grade_UNREVIEWED"] for r in sub)

    record = {
        "study": study,
        "n_admin_screen": len(sub),
        "scope_status": (
            "SCREEN_CONTRIBUTOR_UNSIGNED" if sub
            else "NO_SCREEN_RECORDS_NOT_CERTIFIED_HERE"
        ),
    }
    for key, value in buckets.items():
        record["n_" + key] = bc[value]
    for key, value in grades.items():
        record["n_" + key] = gc[value]
    for key, value in treatment.items():
        record["n_" + key] = tc[value]
    for key, value in patient.items():
        record["n_" + key] = pc[value]

    record["clinical_admissions"] = 0
    output.append(record)

assert sum(x["n_admin_screen"] for x in output) == 2055
assert sum(x["n_non_tumour"] for x in output) == 22
assert sum(x["n_primary_exact"] for x in output) == 5
assert sum(x["n_site_hold"] for x in output) == 52
assert sum(x["n_v5_context"] for x in output) == 280
assert sum(x["n_historical"] for x in output) == 1696

# Prioritize unresolved source evidence, not apparent eligibility.
output.sort(key=lambda r: (
    -r["n_primary_not_attested"],
    -r["n_surgery_not_certified"],
    -r["n_site_hold"],
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

print("STUDY_REVIEW_MATRIX_PASS")
print("Full inventory studies:", len(studies))
print("Screen-contributing studies:", len(by))
print("Other studies:", len(studies) - len(by))
print("Screen experiments:", len(data))
print("Source file SHA256:", digest(source))
print("Review matrix SHA256:", digest(target))
print()
print("TOP STUDIES BY UNATTESTED PRIMARY EVIDENCE")
print("study\tN\tunattested\tsurgery\tstudy_supported\tpatient_key_unverified")
for r in output[:25]:
    print(
        r["study"],
        r["n_admin_screen"],
        r["n_primary_not_attested"],
        r["n_surgery_not_certified"],
        r["n_study_primary_supported"],
        r["n_curated_patient_key_unverified"],
        sep="\t",
    )
print()
print("No source files modified")
print("No patients certified")
print("No model analysis executed")
