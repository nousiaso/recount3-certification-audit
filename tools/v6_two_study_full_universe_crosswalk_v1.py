#!/usr/bin/env python3
"""Reconcile 42 source labels against the existing 75-study-record universe.
No original data modification, approval, run aggregation or model fitting.
"""
import csv, hashlib, json, sys
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

v5 = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_CROSSWALK_UNSIGNED.tsv"
)
folder = repo / "diagnostics/20261124_v6_GBM_osteosarcoma_review"
detail = folder / "42_accession_source_evidence_LOCAL_ONLY.tsv"
prior = folder / "summary.json"
out = repo / "diagnostics/20261124_v6_two_study_full_crosswalk"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

assert sha(v5) == (
    "ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a"
)
assert sha(detail) == (
    "80f22d1bae0e56f985afcee7130b772a5b17816aebb633903bc0736cdb8c0603"
)
assert sha(prior) == (
    "896a1457df47ce656457db417ad5f0c5e427349323cb2954ff913420e65b5160"
)

def read(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

full = [
    r for r in read(v5)
    if r["study"] in {"SRP072494", "SRP090849"}
]
selected = read(detail)

assert Counter(r["study"] for r in full) == {
    "SRP072494": 36, "SRP090849": 39
}
assert Counter(r["study"] for r in selected) == {
    "SRP072494": 20, "SRP090849": 22
}

lookup = {(r["study"], r["experiment_acc"]): r for r in selected}
assert len(lookup) == 42

by_group = defaultdict(list)
joined = []
for r in full:
    k = (r["study"], r["library_id"])
    src = lookup.get(k)
    part = r["stage41_partition_NOT_ADMISSION"]

    if src:
        assert part == "NO_DETECTED_STUDY_ALIAS_OVERLAP_PENDING_GATES"
        assert r["v5_original_sra_sample_name"].startswith("GSM")
        assert r["source_SRA_sample_acc"] == src["sra_sample_acc"]

    by_group[(r["study"], r["patient_key"])].append(r)

    joined.append({
        "study": r["study"],
        "experiment_acc": r["library_id"],
        "stage41_partition": part,
        "GEO_GSM_FROM_ORIGINAL_SRA":
            r["v5_original_sra_sample_name"],
        "original_SRA_sample":
            r["source_SRA_sample_acc"],
        "curated_patient_key_UNVERIFIED":
            r["patient_key"],
        "original_title":
            r["v5_original_sra_sample_title"],
        "source_specimen_evidence_NOT_APPROVAL":
            src["evidence_flag_NOT_APPROVAL"] if src
            else "NOT_IN_42_RECORD_SOURCE_REVIEW",
        "original_treatment_timing_CURATED":
            r["treatment_timing"],
        "clinical_approval": "NO_UNSIGNED"
    })

assert len(joined) == 75
assert len({(r["study"], r["experiment_acc"]) for r in joined}) == 75

gbm = [items for (study, key), items in by_group.items()
       if study == "SRP072494"]
os = [items for (study, key), items in by_group.items()
      if study == "SRP090849"]

assert len(gbm) == 20
assert Counter(map(len, gbm)) == {2: 16, 1: 4}
assert len(os) == 39
assert all(len(x) == 1 for x in os)

assert Counter(
    r["stage41_partition"] for r in joined
    if r["study"] == "SRP072494"
) == {
    "NO_DETECTED_STUDY_ALIAS_OVERLAP_PENDING_GATES": 20,
    "NOT_IN_STAGE41_CANDIDATES": 16,
}
assert Counter(
    r["stage41_partition"] for r in joined
    if r["study"] == "SRP090849"
) == {
    "NO_DETECTED_STUDY_ALIAS_OVERLAP_PENDING_GATES": 22,
    "QC_HOLD_NOT_ADMITTED": 9,
    "NOT_IN_STAGE41_CANDIDATES": 8,
}

out.mkdir(parents=True, exist_ok=True)
target = out / "75_full_study_crosswalk_LOCAL_ONLY.tsv"
with target.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(joined[0]), delimiter="\t",
        lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(sorted(
        joined, key=lambda x: (x["study"], x["experiment_acc"])
    ))

report = {
    "status": "CROSSWALK_PASS_NOT_CLINICAL_CERTIFICATION",
    "full_V5_records": 75,
    "source_reviewed_experiments": 42,
    "gbm_study_records": 36,
    "gbm_baseline_screen": 20,
    "gbm_other_study_records": 16,
    "gbm_curated_two_experiment_groups": 16,
    "gbm_curated_single_experiment_groups": 4,
    "osteosarcoma_study_records": 39,
    "osteosarcoma_screen": 22,
    "osteosarcoma_qc_holds": 9,
    "osteosarcoma_other_records": 8,
    "v5_sha256": sha(v5),
    "prior_42_evidence_sha256": sha(detail),
    "local_crosswalk_sha256": sha(target),
    "clinical_approvals": 0,
}
summary = out / "summary.json"
summary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

print("FULL_STUDY_CROSSWALK_PASS")
for k, v in report.items():
    print(k, v, sep="\t")
print("No original metadata changed")
print("No clinical decisions or model outputs changed")
