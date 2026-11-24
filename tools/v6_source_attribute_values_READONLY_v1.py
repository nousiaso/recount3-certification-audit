#!/usr/bin/env python3
"""Study-level value inventory from frozen 599-record source assembly."""

import csv
import hashlib
import sys
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
folder = repo / "diagnostics/20261123_v6_primary_unattested_review_v1_1"
source = folder / "source_evidence_599_LOCAL_ONLY.tsv"
out = repo / "diagnostics/20261123_v6_source_attribute_values"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

assert sha(source) == (
    "8b6f017dd575ec6fc3917b784740d37c8556fe567"
    "474e6a5e69dcad066ef35eb"
)

with source.open(newline="", encoding="utf-8-sig") as f:
    records = list(csv.DictReader(f, delimiter="\t"))

assert len(records) == 599
assert len({r["experiment_acc"] for r in records}) == 599
assert len({r["study"] for r in records}) == 39
assert all(r["clinical_approval"] == "NO_UNSIGNED"
           for r in records)

GROUPS = {
    "specimen": {
        "tissue", "organ", "organism part", "site",
        "sample type", "sample_type", "tissue type",
        "tumor location", "tumour location",
        "tumor in specimen", "tumour in specimen",
        "sample from primary gbm diagnosis",
    },
    "disease": {
        "disease", "diagnosis", "disease status",
        "disease state", "tumor type", "tumour type",
        "histological tumor classification",
        "who histotype", "intrinsic subtype",
        "disease staging", "disease_stage", "stage",
    },
    "treatment_timing": {
        "treatment", "treatment status",
        "first-line therapy", "therapy",
        "time point", "collection timepoint",
        "subject status", "prior treatment",
    },
    "patient_token": {
        "patient", "patient id", "subject id",
        "individual", "sample donor",
    },
}

by_study = defaultdict(list)
value_counts = Counter()

for row in records:
    attributes = row["original_sample_attributes"]
    parsed = set()

    for item in attributes.split("|"):
        if ";;" not in item:
            continue
        key, value = item.split(";;", 1)
        key = " ".join(key.lower().split())
        value = value.strip()
        if key and value:
            parsed.add((key, value))

    by_study[row["study"]].append((row, parsed))

    # Retain actual values only in the local, ignored detail table.
    for key, value in parsed:
        if any(key in keys for keys in GROUPS.values()):
            value_counts[(row["study"], key, value)] += 1

summary = []

for study, items in sorted(by_study.items()):
    entry = {
        "study": study,
        "n_experiments": len(items),
        "n_without_sample_attributes": sum(
            not r["original_sample_attributes"].strip()
            for r, _ in items
        ),
    }

    for group, keys in GROUPS.items():
        entry["n_with_" + group + "_key"] = sum(
            any(k in keys for k, _ in parsed)
            for _, parsed in items
        )
        entry["n_distinct_" + group + "_values"] = len({
            (k, v)
            for _, parsed in items
            for k, v in parsed if k in keys
        })

    summary.append(entry)

assert sum(r["n_experiments"] for r in summary) == 599

out.mkdir(parents=True, exist_ok=True)

safe = out / "study_attribute_coverage_39.tsv"
with safe.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(summary[0]),
        delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(summary)

local = out / "attribute_value_counts_599_LOCAL_ONLY.tsv"
with local.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(["study", "attribute_key", "attribute_value",
                     "n_experiments"])
    for (study, key, value), n in sorted(value_counts.items()):
        writer.writerow([study, key, value, n])

print("VALUE_LEVEL_SOURCE_AUDIT_PASS")
print("Experiments:", len(records))
print("Studies:", len(summary))
print("No automated specimen classifications")
print()
print("study\tN\tspecimen_field\tdisease_field\ttiming_field\tpatient_token_field")

for r in sorted(summary, key=lambda x: -x["n_experiments"]):
    print(
        r["study"],
        r["n_experiments"],
        r["n_with_specimen_key"],
        r["n_with_disease_key"],
        r["n_with_treatment_timing_key"],
        r["n_with_patient_token_key"],
        sep="\t"
    )

print()
print("Aggregate summary SHA256:", sha(safe))
print("Local full-value table SHA256:", sha(local))
print("NO CLINICAL ADMISSIONS")
