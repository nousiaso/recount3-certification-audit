#!/usr/bin/env python3
"""DRP005276 full 103-experiment source census. Read-only."""
import csv
import gzip
import hashlib
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"
v5 = rec / "clinical_specimen_patient_v5_unsigned_20261122"

paths = {
    "raw": rec / (
        "raw_studies/DRP005276/metadata/"
        "sra.sra.DRP005276.MD.gz"
    ),
    "source": repo / (
        "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
        "source_evidence_599_LOCAL_ONLY.tsv"
    ),
    "inventory": v5 / (
        "01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_"
        "CROSSWALK_UNSIGNED.tsv"
    ),
    "screen": v5 / (
        "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
    ),
}

expected = {
    "raw":
        "cba60a28b225d7852e17e180ddaad2680504e7f1cd5cf94de5d73401db12d36c",
    "source":
        "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb",
    "inventory":
        "ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a",
    "screen":
        "06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49",
}

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def norm(value):
    return " ".join(str(value or "").strip().casefold().split())

def attributes(value):
    out = defaultdict(set)
    for part in (value or "").split("|"):
        if ";;" in part:
            key, val = part.split(";;", 1)
            if norm(key) and norm(val):
                out[norm(key)].add(norm(val))
    return out

def structure(values):
    values = [x for x in values if x]
    counts = Counter(values)
    return (
        len(values),
        len(counts),
        sum(n > 1 for n in counts.values()),
        sum(n for n in counts.values() if n > 1)
    )

print("V6_DRP005276_FULL103_SOURCE_CENSUS")

for key, path in paths.items():
    assert path.is_file(), "Missing documented input: " + key
    assert digest(path) == expected[key], (
        "Original source SHA256 changed: " + key
    )
print("FROZEN_SOURCE_SHA256_PASS", len(paths))

source_all = tsv(paths["source"])
inventory_all = tsv(paths["inventory"])
screen_all = tsv(paths["screen"])

assert len(source_all) == 599
assert len(inventory_all) == 5241
assert len(screen_all) == 2055

selected = {
    r["experiment_acc"]: r
    for r in source_all
    if r["study"] == "DRP005276"
}

inventory = {
    r["library_id"]: r
    for r in inventory_all
    if r["study"] == "DRP005276"
}

screen = {
    r["library_id"]: r
    for r in screen_all
    if r["study"] == "DRP005276"
}

assert len(selected) == 83

raw = defaultdict(list)
with gzip.open(
    str(paths["raw"]), "rt",
    encoding="utf-8", errors="replace",
    newline=""
) as f:
    reader = csv.DictReader(f, delimiter="\t")
    fields = set(reader.fieldnames or [])
    assert {
        "experiment_acc", "sample_name",
        "sample_acc", "sample_title",
        "sample_attributes"
    } <= fields

    for row in reader:
        if row["experiment_acc"]:
            raw[row["experiment_acc"]].append(row)

assert len(raw) == 103
assert set(selected) <= set(raw)

# Detect metadata disagreement across multiple runs.
conflicts = 0
for exp, rows in raw.items():
    for field in (
        "sample_acc", "sample_name",
        "sample_title", "sample_attributes"
    ):
        if len({r.get(field, "") for r in rows}) > 1:
            conflicts += 1

assert conflicts == 0, "Original source specimen conflict"

first = {k: v[0] for k, v in raw.items()}
outside_ids = set(raw) - set(selected)

assert len(outside_ids) == 20

# Previously verified selected-source identity is preserved.
for acc, src in selected.items():
    r = first[acc]
    assert r["sample_title"] == src["original_sample_title"]
    assert r["sample_attributes"] == src["original_sample_attributes"]
    assert r["sample_acc"] == src["sra_sample_acc"]

print("\nSOURCE_UNIVERSE")
print("ORIGINAL_DISTINCT_EXPERIMENTS", len(raw))
print("ORIGINAL_METADATA_ROWS", sum(map(len, raw.values())))
print("SELECTED_EXPERIMENTS", len(selected))
print("OUTSIDE_SELECTED", len(outside_ids))
print("SPECIMEN_METADATA_CONFLICTS", conflicts)

print("\nFROZEN_V5_MEMBERSHIP")
print("V5_STUDY_INVENTORY", len(inventory))
print("V5_STUDY_SCREEN", len(screen))
print("SELECTED_IN_V5_INVENTORY",
      len(set(selected) & set(inventory)))
print("SELECTED_IN_V5_SCREEN",
      len(set(selected) & set(screen)))
print("OUTSIDE_IN_V5_INVENTORY",
      len(outside_ids & set(inventory)))
print("OUTSIDE_IN_V5_SCREEN",
      len(outside_ids & set(screen)))
print("OUTSIDE_NOT_IN_V5_INVENTORY",
      len(outside_ids - set(inventory)))

print("OUTSIDE_V5_QC_HOLDS", sum(
    inventory[x].get("stage41_partition_NOT_ADMISSION") ==
    "QC_HOLD_NOT_ADMITTED"
    for x in outside_ids & set(inventory)
))

selected_attrs = {
    x: attributes(first[x]["sample_attributes"])
    for x in selected
}
outside_attrs = {
    x: attributes(first[x]["sample_attributes"])
    for x in outside_ids
}

def group_report(name, ids, attrs):
    ids = sorted(ids)

    print("\nGROUP", name, "N", len(ids))

    for field in ("sample_name", "sample_acc", "sample_title"):
        print("FIELD", field, "STRUCTURE",
              structure(norm(first[x].get(field, "")) for x in ids))

    for key in ("isolate", "tissue", "disease"):
        vals = [
            next(iter(attrs[x][key]))
            if len(attrs[x].get(key, set())) == 1
            else ""
            for x in ids
        ]
        print("ATTRIBUTE", key, "STRUCTURE", structure(vals))
        print("ATTRIBUTE", key, "MULTIVALUE",
              sum(len(attrs[x].get(key, set())) > 1 for x in ids))

    print("TISSUE_EXACT_BLADDER", sum(
        "bladder" in attrs[x].get("tissue", set())
        for x in ids
    ))
    print("DISEASE_EXACT_BLADDER_CANCER", sum(
        "bladder cancer" in attrs[x].get("disease", set())
        for x in ids
    ))

    terms = (
        "normal", "healthy", "primary",
        "tumor", "tumour", "metasta",
        "recur", "adjacent", "control",
        "blood", "urine", "cell line",
        "pretreat", "untreated", "chemotherapy"
    )
    for term in terms:
        title_n = sum(
            term in norm(first[x]["sample_title"])
            for x in ids
        )
        attr_n = sum(
            term in norm(first[x]["sample_attributes"])
            for x in ids
        )
        if title_n or attr_n:
            print("TERM_HINT", term,
                  "TITLE", title_n, "ATTRIBUTES", attr_n)

group_report("SELECTED_83", set(selected), selected_attrs)
group_report("OUTSIDE_20", outside_ids, outside_attrs)

print("\nCROSS_GROUP_IDENTIFIER_OVERLAP")

for field in ("sample_name", "sample_acc", "sample_title"):
    left = {
        norm(first[x].get(field))
        for x in selected
        if norm(first[x].get(field))
    }
    right = {
        norm(first[x].get(field))
        for x in outside_ids
        if norm(first[x].get(field))
    }
    print(field, "EXACT_SHARED_VALUES", len(left & right))

for key in ("isolate", "tissue", "disease"):
    left = set().union(
        *(selected_attrs[x].get(key, set()) for x in selected)
    )
    right = set().union(
        *(outside_attrs[x].get(key, set()) for x in outside_ids)
    )
    print("ATTRIBUTE", key, "EXACT_SHARED_VALUES", len(left & right))

print("\nSOURCE_STUDY_DESCRIPTIVE_FIELDS")
for field in (
    "study_title", "study_abstract",
    "experiment_title", "library_strategy",
    "library_source", "library_selection"
):
    if field in fields:
        vals = [norm(first[x].get(field, "")) for x in raw]
        print("FIELD", field,
              "POPULATED", sum(bool(v) for v in vals),
              "DISTINCT", len(set(v for v in vals if v)))
    else:
        print("FIELD", field, "NOT_IN_SOURCE_SCHEMA")

print("\nSCIENTIFIC_DISPOSITION")
print("FULL_SOURCE_MEMBERSHIP", "RECONCILED")
print("OUTSIDE_EXPERIMENTS_SOURCE_CLASS", "PENDING_INTERPRETATION")
print("PATIENT_IDENTITY", "UNSIGNED")
print("PRIMARY_TUMOUR_STATUS", "UNSIGNED")
print("PRETREATMENT_STATUS", "UNSIGNED")
print("CLINICAL_ADMISSIONS", 0)
print("NO_SAMPLE_OR_PATIENT_VALUES_PRINTED")
print("NO_FROZEN_FILES_MODIFIED")
print("V6_DRP005276_FULL103_SOURCE_CENSUS_COMPLETE")
