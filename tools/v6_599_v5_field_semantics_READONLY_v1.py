#!/usr/bin/env python3
"""Audit semantic evidence content of populated V5 clinical fields."""
import csv, hashlib, re, sys
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
p = repo / (
    "diagnostics/20261124_v6_599_integrated_review/"
    "599_accession_integrated_evidence_LOCAL_ONLY.tsv"
)

def sha(x):
    return hashlib.sha256(x.read_bytes()).hexdigest()

assert sha(p) == (
    "c77599f38ecfcf16be5f8a4648f0386421f65de215a9a3a3c13bcb965f82291a"
)

with p.open(newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 599
assert len({(r["study"], r["experiment_acc"]) for r in rows}) == 599

patient_grade = Counter(r["v5_patient_evidence_grade"] for r in rows)
treatment_grade = Counter(r["v5_treatment_evidence_grade"] for r in rows)
tissue_grade = Counter(r["v5_tissue_evidence_grade"] for r in rows)

patient_groups = defaultdict(list)
treatment_values = Counter()

for r in rows:
    key = r["v5_patient_key_LOCAL_ONLY"].strip()
    timing = r["v5_treatment_timing_LOCAL_ONLY"].strip()
    assert key and timing
    patient_groups[(r["study"], key)].append(r["experiment_acc"])
    treatment_values[timing] += 1

def label(x):
    # Only display simple categorical codes. Do not disclose
    # free-text specimen descriptions or patient-linked values.
    if re.fullmatch(r"[A-Z][A-Z0-9_:/+-]{0,90}", x):
        return x
    return "LOCAL_VALUE_SHA256_" + hashlib.sha256(
        x.encode()).hexdigest()[:12]

print("V5_599_FIELD_SEMANTICS_AUDIT")
print("Experiments:", len(rows))
print("Distinct study-scoped V5 patient keys:", len(patient_groups))
print("Patient-key groups reused by >1 experiment:",
      sum(len(v) > 1 for v in patient_groups.values()))
print("Experiments in reused-key groups:",
      sum(len(v) for v in patient_groups.values() if len(v) > 1))
print("Distinct V5 treatment values:", len(treatment_values))

for heading, counts in [
    ("PATIENT_EVIDENCE_GRADES", patient_grade),
    ("TREATMENT_EVIDENCE_GRADES", treatment_grade),
    ("TISSUE_EVIDENCE_GRADES", tissue_grade)
]:
    print("\n", heading)
    for value, n in counts.most_common():
        print(label(value), n, sep="\t")

print("\nTREATMENT_TIMING_VALUE_FREQUENCIES")
for value, n in treatment_values.most_common(25):
    print(label(value), n, sep="\t")

print("\nREPEATED_V5_PATIENT_KEYS_BY_STUDY")
by_study = Counter()
for (study, key), group in patient_groups.items():
    if len(group) > 1:
        by_study[study] += 1
for study, n in sorted(by_study.items()):
    print(study, n, sep="\t")

print("\nNo patient keys printed")
print("No clinical grades changed")
print("NO CLINICAL ADMISSIONS")
print("V5_FIELD_SEMANTICS_REVIEW_COMPLETE")
