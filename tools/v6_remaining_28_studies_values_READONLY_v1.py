#!/usr/bin/env python3
"""Review categorical clinical evidence in the remaining 28 studies."""
import csv, hashlib, re, sys
from pathlib import Path
from collections import defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"

source = repo / (
    "diagnostics/20261123_v6_source_attribute_values/"
    "attribute_value_counts_599_LOCAL_ONLY.tsv"
)
coverage = repo / (
    "diagnostics/20261123_v6_source_attribute_values/"
    "study_attribute_coverage_39.tsv"
)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

assert sha(source) == (
    "400329705cb6fc82c223b7a2384f2c0b24e02650e192b1d5a5988465f35fd452"
)
assert sha(coverage) == (
    "4bbfcf2769fb5cc4108e04c3c3443f8805b2c1f21c29bf506128257ff5eb220e"
)

def read(p):
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

all_studies = read(coverage)
assert len(all_studies) == 39
assert sum(int(r["n_experiments"]) for r in all_studies) == 599

previously_reviewed = {
    "DRP005276", "SRP113755", "SRP221385",
    "SRP113639", "SRP156394", "ERP023222",
    "SRP072494", "ERP010930", "SRP090849",
    "SRP094795", "SRP018008"
}
remaining = {
    r["study"]: int(r["n_experiments"])
    for r in all_studies
    if r["study"] not in previously_reviewed
}
assert len(remaining) == 28
assert sum(remaining.values()) == 217

# Clinical categorical keys only. Patient, age, isolate and
# source_name values are deliberately excluded.
allowed = {
    "tissue", "organ", "organism part", "site",
    "sample type", "sample_type", "tissue type",
    "tumor type", "tumour type",
    "tumor in specimen", "tumour in specimen",
    "tumor location", "tumour location",
    "disease", "diagnosis", "disease status",
    "disease state", "disease staging", "disease_stage",
    "stage", "histological tumor classification",
    "who histotype", "intrinsic subtype",
    "treatment", "treatment status", "prior treatment",
    "first-line therapy", "time point",
    "sample from primary gbm diagnosis",
}
by = defaultdict(list)

for r in read(source):
    if r["study"] in remaining and r["attribute_key"] in allowed:
        by[(r["study"], r["attribute_key"])].append(r)

def safe_display(value):
    value = value.strip()
    if (
        not value or len(value) > 80
        or re.search(r"(?i)(patient|subject|donor|sample)\s*"
                     r"[-_:#]*[A-Za-z]*\d+", value)
        or re.search(r"(?i)\b(?:SRX|SRS|GSM|SAMN)\d+\b", value)
        or "@" in value or "\n" in value
    ):
        return "[VALUE_RETAINED_LOCALLY]"
    return value

print("REMAINING_28_STUDIES_SOURCE_VALUES")
print("Total studies:", len(remaining))
print("Total selected experiments:", sum(remaining.values()))
print("Original SRA sources not reread")
print("No automatic clinical classifications")

for study, n in sorted(remaining.items(), key=lambda x: (-x[1], x[0])):
    print("\nSTUDY", study, "N", n)
    keys = sorted(k for s, k in by if s == study)
    if not keys:
        print("NO_RECOGNIZED_CLINICAL_CATEGORICAL_FIELDS")
        continue

    for key in keys:
        values = sorted(
            by[(study, key)],
            key=lambda r: (
                -int(r["n_experiments"]), r["attribute_value"]
            )
        )
        print("FIELD", key, "N_DISTINCT", len(values))
        for r in values[:12]:
            print(
                "VALUE",
                key,
                repr(safe_display(r["attribute_value"])),
                "n=" + r["n_experiments"],
                sep="\t"
            )
        if len(values) > 12:
            print("OTHER_DISTINCT_VALUES", len(values) - 12)

print("\nFIRST_PASS_REMAINING_28_COMPLETE")
print("No original files modified")
print("No clinical admissions")
