#!/usr/bin/env python3
"""Read-only membership and source preflight for SRP057793."""

import csv
import io
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
base = root / "26_recount3_external_primary_v1"
target = "SRP057793"

v3 = base / "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"
v5 = base / "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"

stage41 = root / (
    "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
)

candidate = stage41 / "03_candidate_integrated_evidence_3602.tsv"
fullsource = stage41 / "01_full_source_experiment_identity_5221.tsv"

def local_table(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return reader.fieldnames or [], list(reader)

def archive_table(path, suffix):
    with zipfile.ZipFile(path) as z:
        matches = [
            x for x in z.namelist()
            if x.endswith("/" + suffix)
        ]
        assert len(matches) == 1, (path.name, suffix)
        raw = z.read(matches[0]).decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(raw), delimiter="\t")
    return reader.fieldnames or [], list(reader)

def category(value):
    s = str(value or "").strip()
    if not s:
        return "[EMPTY]"
    if re.fullmatch(r"[A-Z0-9_./-]{1,85}", s):
        return s
    return "[FREE_TEXT_SUPPRESSED]"

def counts(rows, field):
    c = Counter(category(r.get(field)) for r in rows)
    return dict(sorted(c.items()))

def ids(rows, field):
    return {
        str(r.get(field) or "").strip()
        for r in rows
        if str(r.get(field) or "").strip()
    }

print("SRP057793_COHORT_MEMBERSHIP_V1")

fc, candidates = local_table(candidate)
ff, original = local_table(fullsource)

pv, primary = archive_table(
    v5, "results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
)
rv, studies = archive_table(
    v3, "results/02_ALL_176_STUDY_SOURCE_REVIEW_UNSIGNED.tsv"
)

assert len(candidates) == 3602
assert len(original) == 5221
assert len(primary) == 2055
assert len(studies) == 176

c = [r for r in candidates if r.get("study") == target]
f = [r for r in original if r.get("study") == target]
p = [r for r in primary if r.get("study") == target]
r = [r for r in studies if r.get("study") == target]

assert len(r) == 1
assert len(p) == 140

cid = ids(c, "library_id")
fid = ids(f, "experiment_acc")
pid = ids(p, "library_id")

print("\nMEMBERSHIP")
print("ORIGINAL_SOURCE_ROWS", len(f))
print("STAGE41_CANDIDATE_ROWS", len(c))
print("V5_PRIMARY_ROWS", len(p))
print("ORIGINAL_DISTINCT_EXPERIMENTS", len(fid))
print("CANDIDATE_DISTINCT_EXPERIMENTS", len(cid))
print("PRIMARY_DISTINCT_EXPERIMENTS", len(pid))
print("CANDIDATE_NOT_IN_ORIGINAL", len(cid - fid))
print("PRIMARY_NOT_IN_ORIGINAL", len(pid - fid))
print("PRIMARY_NOT_IN_CANDIDATES", len(pid - cid))
print("ORIGINAL_NOT_IN_PRIMARY", len(fid - pid))
print("CANDIDATES_NOT_IN_PRIMARY", len(cid - pid))

print("\nHISTORICAL_STUDY_REVIEW")
review = r[0]

for field in (
    "cancer",
    "publication_doi",
    "libraries",
    "runs",
    "patients",
    "patient_key_basis",
    "qa_original_libraries",
    "qa_stage41_primary_screen",
    "qa_missing_pretreatment",
    "qa_clinical_holds",
    "qa_patient_key_repeated_libraries",
    "overall_study_disposition",
    "study_level_source_verification",
    "public_source_review_status",
    "review_action_from_publication",
    "FINAL_STUDY_CLINICAL_APPROVAL_AS_OF_V3"
):
    if field in rv:
        print("STUDY_FIELD", field, category(review.get(field)))

# Show which source evidence is available,
# without exposing individual patient/source values.
for field in (
    "publication_doi",
    "public_source_url",
    "STUDY_EVIDENCE_URL",
    "public_source_tissue_evidence",
    "public_source_treatment_evidence",
    "public_source_patient_evidence",
    "full_source_curation_comment",
    "review_caveat",
    "required_evidence"
):
    if field in rv:
        value = str(review.get(field) or "")
        print(
            "STUDY_EVIDENCE_FIELD", field,
            "POPULATED", bool(value.strip()),
            "CHARACTERS", len(value)
        )

print("\nSTAGE41_CANDIDATE_GATES")
for field in (
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
    if field in fc:
        print("CANDIDATE_STATUS", field, counts(c, field))
    else:
        print("CANDIDATE_FIELD_ABSENT", field)

print("\nV5_PRIMARY_METADATA")
for field in (
    "patient_key",
    "biosample",
    "patient_key_basis",
    "metastatic",
    "treatment",
    "treatment_timing",
    "collection_time_point",
    "source_supported_tissue_context_NOT_SRX_CERTIFIED",
    "source_supported_treatment_context_NOT_SRX_CERTIFIED",
    "INDIVIDUAL_SOURCE_ATTESTATION",
    "FINAL_PRIMARY_RELEASE_v3",
    "v5_final_primary_approval"
):
    if field not in pv:
        print("V5_FIELD_ABSENT", field)
        continue

    values = [
        str(x.get(field) or "").strip()
        for x in p
    ]
    populated = [x for x in values if x]

    print(
        "V5_FIELD", field,
        "POPULATED", len(populated),
        "DISTINCT", len(set(populated))
    )

for field in (
    "patient_key_basis",
    "FINAL_PRIMARY_RELEASE_v3",
    "v5_final_primary_approval"
):
    if field in pv:
        print("V5_STATUS", field, counts(p, field))

# Patient-key multiplicity is descriptive only.
if "patient_key" in pv:
    patient_keys = [
        str(x.get("patient_key") or "").strip()
        for x in p
    ]
    valid = [x for x in patient_keys if x]
    multiplicities = Counter(valid)
    distribution = Counter(multiplicities.values())
    print("PATIENT_KEYS_NONEMPTY", len(valid))
    print("PATIENT_KEYS_DISTINCT", len(multiplicities))
    print(
        "PATIENT_KEY_MULTIPLICITY_DISTRIBUTION",
        dict(sorted(distribution.items()))
    )

print("\nAUDIT_LIMITATIONS")
print("SOURCE_PATIENT_KEYS_NOT_CERTIFIED")
print("NO_PATIENT_OR_EXPERIMENT_VALUES_PRINTED")
print("NO_ORIGINAL_FILES_MODIFIED")
print("NO_CLINICAL_APPROVAL_INFERRED")
print("CLINICAL_ADMISSIONS 0")
print("SRP057793_COHORT_MEMBERSHIP_V1_COMPLETE")
