#!/usr/bin/env python3
"""Reconcile source-metastasis evidence and non-V5 tumour exceptions."""
import csv, hashlib, sys
from pathlib import Path
from collections import Counter

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

source = repo / (
    "diagnostics/20261124_v6_osteosarcoma_source_census/"
    "original_source_membership_LOCAL_ONLY.tsv"
)
v5 = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_CROSSWALK_UNSIGNED.tsv"
)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

assert sha(source) == (
    "c9f871261e9a64fd7857d443a9ebc14d86e480e0d538c987c506197229703316"
)
assert sha(v5) == (
    "ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a"
)

def read(p):
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

src = read(source)
historical = [
    r for r in read(v5) if r["study"] == "SRP090849"
]
assert len(src) == 52 and len(historical) == 39

hist_by_id = {r["library_id"]: r for r in historical}
src_by_id = {r["experiment_acc"]: r for r in src}
assert len(hist_by_id) == 39 and len(src_by_id) == 52

# Inspect original historical grades without guessing which
# V5 column contains the METASTATIC_CONTEXT designation.
historic_met = set()
grade_columns = Counter()

for accession, r in hist_by_id.items():
    for field, value in r.items():
        if "METASTATIC_CONTEXT" in (value or ""):
            historic_met.add(accession)
            grade_columns[field] += 1

source_met_in_v5 = {
    r["experiment_acc"] for r in src
    if r["tumor_type_source_value"].strip().lower() == "metastasis"
    and r["in_v5"] == "YES"
}

print("SRP090849_EXCEPTION_RECONCILIATION")
print("Historical metastatic-context field counts:",
      dict(grade_columns))
print("Historical metastatic-context experiments:",
      len(historic_met))
print("SRA metastatic-labelled experiments in V5:",
      len(source_met_in_v5))
print("Exact accession intersection:",
      len(historic_met & source_met_in_v5))
print("Historical only:", len(historic_met - source_met_in_v5))
print("Source only:", len(source_met_in_v5 - historic_met))

if historic_met:
    for acc in sorted(historic_met ^ source_met_in_v5):
        print("METASTASIS_DISCORDANCE", acc)
else:
    print("HISTORICAL_FLAG_FIELD_NOT_LOCATED_IN_V5_VALUES")

print("\nNON-V5 SOURCE RECORDS")
outside = [
    r for r in src if r["in_v5"] == "NO"
]
assert len(outside) == 13

for r in sorted(outside, key=lambda x: x["experiment_acc"]):
    print(
        "NON_V5",
        r["experiment_acc"],
        r["tumor_type_source_value"],
        r["tissue_source_value"],
        sep="\t"
    )

print("\nQC-HELD PRIMARY-LABELLED EXPERIMENTS")
held = [
    r for r in src
    if r["v5_partition"] == "QC_HOLD_NOT_ADMITTED"
    and r["tumor_type_source_value"].lower() == "primary"
]
assert len(held) == 4
for r in sorted(held, key=lambda x: x["experiment_acc"]):
    print("QC_HOLD_PRIMARY_SOURCE_LABEL", r["experiment_acc"])

print("\nNO SOURCE FILES MODIFIED")
print("NO CLINICAL ADMISSIONS")
print("EXCEPTION_REVIEW_COMPLETE")
