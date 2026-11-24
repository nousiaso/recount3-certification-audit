#!/usr/bin/env python3
"""Integrate frozen V5, V6 and 599 original SRA records.
Evidence only. No clinical classification, approval or model execution.
"""
import csv, hashlib, json, sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

source = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)
v6 = repo / (
    "diagnostics/20261123_v6_evidence_precedence_dryrun/"
    "proposed_evidence_2055_UNSIGNED.tsv"
)
v5 = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
)
out = repo / "diagnostics/20261124_v6_599_integrated_review"

expected = {
    source: "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb",
    v6: "08efe2c8f469edcbdfb98231dc87f705d801434e78dcbea85f57c4b94da88eb0",
    v5: "06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49",
}

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def read(p):
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

for p, checksum in expected.items():
    assert p.is_file() and sha(p) == checksum, f"INPUT_CHANGED: {p}"

original = read(source)
frozen_v6 = read(v6)
frozen_v5 = read(v5)

assert len(original) == 599
assert len(frozen_v6) == len(frozen_v5) == 2055

v6_by_id = {r["library_id"]: r for r in frozen_v6}
v5_by_id = {r["library_id"]: r for r in frozen_v5}
assert len(v6_by_id) == len(v5_by_id) == 2055

selected = {
    r["library_id"]
    for r in frozen_v6
    if r["unsigned_evidence_bucket"] ==
       "HISTORICAL_EVIDENCE_REVIEW_PENDING"
    and r["historical_tissue_grade"] ==
        "PRIMARY_TUMOUR_NOT_ATTESTED"
}

assert len(selected) == 599
assert {r["experiment_acc"] for r in original} == selected
assert len({r["study"] for r in original}) == 39

groups = {
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
relevant_keys = set().union(*groups.values())

def parse_attributes(text):
    result = defaultdict(set)
    for segment in text.split("|"):
        if ";;" not in segment:
            continue
        key, value = segment.split(";;", 1)
        key = " ".join(key.strip().lower().split())
        value = value.strip()
        if key and value:
            result[key].add(value)
    return {k: sorted(v) for k, v in result.items()}

records = []
per_study = defaultdict(lambda: Counter())
srs_conflicts = []

for src in sorted(original, key=lambda r: (
    r["study"], r["experiment_acc"]
)):
    acc = src["experiment_acc"]
    old = v5_by_id[acc]
    current = v6_by_id[acc]
    study = src["study"]

    assert old["study"] == current["study"] == study
    assert current["historical_tissue_grade"] == \
           old["tissue_evidence_grade"]
    assert current["FINAL_APPROVAL"] == "NO_UNSIGNED"
    assert old["FINAL_CLINICAL_APPROVAL"] == "NO_UNSIGNED"
    assert src["clinical_approval"] == "NO_UNSIGNED"
    assert src["raw_source_status"] == "EXACT_EXPERIMENT_MATCH"

    a = parse_attributes(src["original_sample_attributes"])
    clinical = {
        k: a[k] for k in sorted(a) if k in relevant_keys
    }

    sra_sample = src["sra_sample_acc"]
    older_sample = old.get("source_SRA_sample_acc", "")
    if sra_sample and older_sample and sra_sample != older_sample:
        srs_conflicts.append((study, acc))

    flags = {
        category: int(any(k in a for k in keys))
        for category, keys in groups.items()
    }

    counters = per_study[study]
    counters["n"] += 1
    for category, present in flags.items():
        counters["has_" + category] += present
    counters["has_v5_patient_key"] += int(
        bool(old.get("patient_key", "").strip())
    )
    counters["has_v5_treatment_value"] += int(
        bool(old.get("treatment_timing", "").strip())
    )

    records.append({
        "study": study,
        "experiment_acc": acc,
        "sra_sample_acc": sra_sample,
        "v5_sra_sample_acc": older_sample,
        "original_sample_title_LOCAL_ONLY":
            src["original_sample_title"],
        "original_clinical_values_json_LOCAL_ONLY":
            json.dumps(clinical, sort_keys=True),
        "v5_tissue_evidence_grade":
            old["tissue_evidence_grade"],
        "v5_treatment_evidence_grade":
            old["treatment_evidence_grade"],
        "v5_patient_evidence_grade":
            old["patient_evidence_grade"],
        "v5_treatment_timing_LOCAL_ONLY":
            old["treatment_timing"],
        "v5_patient_key_LOCAL_ONLY":
            old["patient_key"],
        "v5_within_study_group_LOCAL_ONLY":
            old.get("v5_patient_group_within_study", ""),
        "v5_stage41_partition":
            old["stage41_partition_NOT_ADMISSION"],
        "v6_bucket":
            current["unsigned_evidence_bucket"],
        "source_file_sha256":
            src["source_file_sha256"],
        "primary_sampled_site_certified": "NO",
        "pretreatment_certified": "NO",
        "patient_independence_certified": "NO",
        "technical_qc_certified_for_release": "NO",
        "atlas_independence_certified": "NO",
        "model_transfer_certified": "NO",
        "clinical_approval": "NO_UNSIGNED",
        "review_status": "EVIDENCE_ASSEMBLED_ADJUDICATION_OPEN",
    })

assert len(records) == 599
assert len(per_study) == 39

totals = Counter()
for counts in per_study.values():
    totals.update(counts)

assert totals["n"] == 599
for category, expected_n in {
    "specimen": 571,
    "disease": 441,
    "treatment_timing": 24,
    "patient_token": 75,
}.items():
    assert totals["has_" + category] == expected_n, (
        category, totals["has_" + category]
    )

out.mkdir(parents=True, exist_ok=True)

local = out / "599_accession_integrated_evidence_LOCAL_ONLY.tsv"
with local.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(records[0]),
        delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(records)

study_rows = []
for study, c in sorted(per_study.items()):
    study_rows.append({
        "study": study,
        "n_experiments": c["n"],
        "n_with_sra_specimen_field": c["has_specimen"],
        "n_with_sra_disease_field": c["has_disease"],
        "n_with_sra_timing_field": c["has_treatment_timing"],
        "n_with_sra_patient_token": c["has_patient_token"],
        "n_with_v5_patient_key": c["has_v5_patient_key"],
        "n_with_v5_treatment_value": c["has_v5_treatment_value"],
        "review_status": "CLINICAL_ADJUDICATION_OPEN",
    })

study_file = out / "39_study_integrated_coverage_UNSIGNED.tsv"
with study_file.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(study_rows[0]),
        delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(study_rows)

summary = {
    "status": (
        "SOURCE_JOIN_CONFLICT_HOLD" if srs_conflicts
        else "INTEGRATED_EVIDENCE_PASS_NOT_CERTIFIED"
    ),
    "n_experiments": 599,
    "n_studies": 39,
    "n_sra_v5_sample_accession_conflicts": len(srs_conflicts),
    "source_hashes": {
        str(p.relative_to(root)): expected[p]
        for p in expected
    },
    "source_attribute_coverage": {
        k: totals["has_" + k] for k in groups
    },
    "local_599_matrix_sha256": sha(local),
    "study_summary_sha256": sha(study_file),
    "clinical_approvals": 0,
}
summary_file = out / "summary.json"
summary_file.write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)

print("V6_INTEGRATED_599_EVIDENCE_REVIEW")
print("Experiments:", len(records))
print("Studies:", len(study_rows))
print("SRA-to-V5 sample-accession conflicts:", len(srs_conflicts))
print("Study coverage totals:", dict(totals))
print("Local 599 matrix SHA256:", sha(local))
print("39-study summary SHA256:", sha(study_file))
print("Summary JSON SHA256:", sha(summary_file))
print("Clinical approvals:", 0)
print(summary["status"])
