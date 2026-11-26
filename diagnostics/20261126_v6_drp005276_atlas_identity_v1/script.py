#!/usr/bin/env python3
"""Targeted DRP005276-to-ARCHS4 identifier overlap audit.

Read-only; aggregate stdout only; Python 3.6 compatible.
Negative identifier matches never certify patient independence.
"""

import csv
import gzip
import hashlib
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

raw_path = rec / (
    "raw_studies/DRP005276/metadata/"
    "sra.sra.DRP005276.MD.gz"
)

source_path = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)

v5_path = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_"
    "CROSSWALK_UNSIGNED.tsv"
)

atlas_path = root / (
    "RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/"
    "PROJECT_TREE/"
    "14_projection_use_analysis_ready_atlas_v1_0_3/"
    "03_projection_analysis_ready_core.tsv.gz"
)

EXPECTED = {
    raw_path:
        "cba60a28b225d7852e17e180ddaad2680504e7f1cd5cf94de5d73401db12d36c",
    source_path:
        "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb",
    v5_path:
        "ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a"
}

GSM = re.compile(r"\bGSM[0-9]+\b", re.I)

def norm(x):
    return " ".join(str(x or "").strip().casefold().split())

def file_sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def read_tsv(path, compressed=False):
    opener = gzip.open if compressed else open
    with opener(
        str(path), "rt",
        encoding="utf-8-sig", newline=""
    ) as f:
        return list(csv.DictReader(f, delimiter="\t"))

def gsm_tokens(row):
    """Recover literal GSM aliases from identity-related fields.

    Exclude study abstracts and general descriptive fields;
    incidental citations are not accession crosswalks.
    """
    values = set()

    for field, value in row.items():
        key = field.casefold()

        if not any(x in key for x in (
            "name", "title", "alias", "accession",
            "attribute", "geo", "gsm"
        )):
            continue

        for token in GSM.findall(value or ""):
            values.add(token.upper())

    return values

print("V6_DRP005276_ATLAS_IDENTITY_V1")

# ------------------------------------------------------------
# 1. Verify frozen inputs by their previous SHA256 values.
# ------------------------------------------------------------

for path, checksum in EXPECTED.items():
    assert path.is_file(), "FROZEN_INPUT_MISSING"
    assert file_sha(path) == checksum, "FROZEN_INPUT_CHANGED"

assert atlas_path.is_file(), "FROZEN_ATLAS_MISSING"

print("FROZEN_LOCAL_INPUT_SHA256_PASS", len(EXPECTED))
print("FROZEN_ATLAS_SHA256", file_sha(atlas_path))

raw_rows = read_tsv(raw_path, compressed=True)
source_rows = read_tsv(source_path)
v5_rows = read_tsv(v5_path)
atlas_rows = read_tsv(atlas_path, compressed=True)

assert len(raw_rows) == 103
assert len(source_rows) == 599
assert len(v5_rows) == 5241
assert len(atlas_rows) == 9465

original = {r["experiment_acc"]: r for r in raw_rows}
selected_source = {
    r["experiment_acc"]: r
    for r in source_rows
    if r["study"] == "DRP005276"
}

assert len(original) == 103
assert len(selected_source) == 83
assert set(selected_source) <= set(original)

selected = set(selected_source)
held = set(original) - selected
assert len(held) == 20

v5_drp = {
    r["library_id"]: r
    for r in v5_rows if r["study"] == "DRP005276"
}

assert len(v5_drp) == 103
assert set(v5_drp) == set(original)

# The atlas frozen sample ID axis is GEO GSM accessions.
atlas_ids = [r["sample_id"] for r in atlas_rows]
assert len(set(atlas_ids)) == 9465
assert all(
    re.fullmatch(r"GSM[0-9]+", x)
    for x in atlas_ids
)

atlas_gsm = set(atlas_ids)

print("\nDATA_UNIVERSES")
print("ORIGINAL_103", len(original))
print("SELECTED_83", len(selected))
print("HISTORICAL_QC_HELD_20", len(held))
print("FROZEN_ATLAS_9465", len(atlas_gsm))
print("V5_INVENTORY_5241", len(v5_rows))

# ------------------------------------------------------------
# 2. Literal GSM links: raw SRA versus previously extracted
#    source evidence for the 83 selected experiments.
# ------------------------------------------------------------

alias_by_exp = {}

for acc, raw in original.items():
    tokens = gsm_tokens(raw)

    if acc in selected_source:
        tokens |= gsm_tokens(selected_source[acc])

    alias_by_exp[acc] = tokens

print("\nGSM_ALIAS_EVIDENCE")

for label, ids in (
    ("ALL_103", set(original)),
    ("SELECTED_83", selected),
    ("QC_HELD_20", held)
):
    with_alias = sum(
        bool(alias_by_exp[x]) for x in ids
    )

    direct_atlas = sum(
        bool(alias_by_exp[x] & atlas_gsm)
        for x in ids
    )

    print(
        label,
        "EXPERIMENTS_WITH_RECOVERED_GSM", with_alias,
        "EXACT_FROZEN_ATLAS_GSM_MATCH", direct_atlas
    )

print(
    "DISTINCT_GSM_ALIASES",
    len(set().union(*alias_by_exp.values()))
)

# ------------------------------------------------------------
# 3. Exact experiment reuse across studies in V5.
#
# Experiment accessions are strong technical identifiers.
# Duplicate sample names alone are weaker evidence.
# ------------------------------------------------------------

other = [
    r for r in v5_rows
    if r["study"] != "DRP005276"
]

other_experiments = {
    r["library_id"] for r in other
}

print("\nCROSS_STUDY_EXPERIMENT_ID")

print(
    "DRP_EXPERIMENT_ACCESSIONS_IN_OTHER_STUDIES",
    len(set(original) & other_experiments)
)

# ------------------------------------------------------------
# 4. Cross-study source-sample-name collisions.
#
# This is a screening comparison, NOT a patient identity join.
# Compare only directly named sample-name fields.
# ------------------------------------------------------------

name_columns = [
    field for field in v5_rows[0]
    if field in (
        "v5_original_sra_sample_name",
        "source_sample_name",
        "sample_name"
    )
]

print("\nCROSS_STUDY_SAMPLE_NAME_SCREEN")
print(
    "V5_NAMED_SAMPLE_FIELDS_AVAILABLE",
    len(name_columns)
)

drp_names = {
    acc: norm(row.get("sample_name", ""))
    for acc, row in original.items()
}

other_names = defaultdict(set)

for row in other:
    for field in name_columns:
        value = norm(row.get(field, ""))
        if value:
            other_names[value].add(row["study"])

for label, ids in (
    ("SELECTED_83", selected),
    ("QC_HELD_20", held)
):
    found = [
        acc for acc in ids
        if drp_names[acc]
        and drp_names[acc] in other_names
    ]

    print(
        label,
        "SAMPLE_NAME_COLLISION_EXPERIMENTS",
        len(found)
    )

# ------------------------------------------------------------
# 5. Exact SRA sample-accession screening where V5 provides
#    an explicitly named sample-accession field.
#
# Do not confuse an ordinary sample label with BioSample identity.
# ------------------------------------------------------------

accession_columns = [
    field for field in v5_rows[0]
    if (
        "sample" in field.casefold()
        and any(term in field.casefold()
                for term in ("acc", "accession", "biosample"))
    )
]

print("\nCROSS_STUDY_SAMPLE_ACCESSION_SCREEN")
print(
    "V5_SAMPLE_ACCESSION_FIELD_COUNT",
    len(accession_columns)
)

# Field names are non-patient schema metadata and can be printed.
print(
    "V5_SAMPLE_ACCESSION_FIELDS",
    ";".join(sorted(accession_columns))
    if accession_columns else "NONE"
)

other_accessions = set()

for row in other:
    for field in accession_columns:
        value = norm(row.get(field, ""))
        if value:
            other_accessions.add(value)

for label, ids in (
    ("SELECTED_83", selected),
    ("QC_HELD_20", held)
):
    matched = sum(
        norm(original[acc].get("sample_acc", ""))
        in other_accessions
        for acc in ids
        if norm(original[acc].get("sample_acc", ""))
    )

    print(
        label,
        "EXACT_SAMPLE_ACCESSION_COLLISIONS",
        matched
    )

# ------------------------------------------------------------
# 6. Existing V5 overlap flag availability.
#    Do not interpret unknown enum values as approvals.
# ------------------------------------------------------------

print("\nEXISTING_V5_OVERLAP_FIELDS")

flag_fields = [
    field for field in v5_rows[0]
    if "atlas" in field.casefold()
    or "overlap" in field.casefold()
]

print("V5_ATLAS_OR_OVERLAP_FIELD_COUNT", len(flag_fields))

for field in sorted(flag_fields):
    values = [
        norm(v5_drp[x].get(field, ""))
        for x in selected
    ]

    print(
        "FIELD", field,
        "NONEMPTY_SELECTED",
        sum(bool(v) for v in values),
        "EXPLICIT_YES",
        sum(v == "yes" for v in values),
        "EXPLICIT_NO",
        sum(v == "no" for v in values),
        "OTHER_NONEMPTY",
        sum(v not in ("", "yes", "no") for v in values)
    )

# ------------------------------------------------------------
# 7. Scientific interpretation guardrails.
# ------------------------------------------------------------

print("\nSCIENTIFIC_GATES")
print(
    "EXACT_GSM_SCREEN",
    "ASSESS_FROM_COUNTS_NOT_GLOBAL_INDEPENDENCE"
)
print(
    "SAMPLE_NAME_COLLISIONS",
    "SCREEN_ONLY_NOT_PATIENT_IDENTITY"
)
print(
    "SAMPLE_ACCESSION_COVERAGE",
    "REQUIRES_FIELD_AND_SOURCE_VALIDATION"
)
print(
    "ARCHS4_GLOBAL_PATIENT_INDEPENDENCE",
    "NOT_CERTIFIED"
)
print("PRETREATMENT_STATUS", "NOT_CERTIFIED")
print("CLINICAL_ADMISSIONS", 0)
print("NO_INDIVIDUAL_IDENTIFIERS_PRINTED")
print("NO_FROZEN_SOURCES_MODIFIED")

print("V6_DRP005276_ATLAS_IDENTITY_V1_COMPLETE")
