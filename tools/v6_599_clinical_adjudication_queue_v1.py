#!/usr/bin/env python3
"""Build reviewer queue from the frozen 599-record evidence matrix.
Does not infer biological eligibility, certify patients or alter sources.
"""
import csv, hashlib, json, sys
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
source = repo / (
    "diagnostics/20261124_v6_599_integrated_review/"
    "599_accession_integrated_evidence_LOCAL_ONLY.tsv"
)
out = repo / "diagnostics/20261124_v6_599_adjudication_queue"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

assert sha(source) == (
    "c77599f38ecfcf16be5f8a4648f0386421f65de215a9a3a3c13bcb965f82291a"
)

with source.open(newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 599
assert len({(r["study"], r["experiment_acc"]) for r in rows}) == 599
assert len({r["study"] for r in rows}) == 39
assert all(r["clinical_approval"] == "NO_UNSIGNED" for r in rows)

P = [
    "CURATED_STUDY_PATIENT_KEY_UNVERIFIED",
    "BIOSAMPLE_PROXY_NOT_PATIENT_ID",
    "STAGE41_PATIENT_PROPOSAL_NOT_VERIFIED",
]
T = [
    "PRETREATMENT_NOT_ESTABLISHED",
    "CURATED_PRETREATMENT_LABEL_REQUIRES_SOURCE_CONFIRMATION",
    "GEO_BASELINE_TEXT_SOURCE_SUPPORTED",
]

pg = Counter(r["v5_patient_evidence_grade"] for r in rows)
tg = Counter(r["v5_treatment_evidence_grade"] for r in rows)

assert pg == dict(zip(P, (467, 90, 42)))
assert tg == dict(zip(T, (552, 24, 23)))

groups = defaultdict(list)
by_study = defaultdict(Counter)
timing_cross = Counter()
queue = []

for r in sorted(rows, key=lambda x: (x["study"], x["experiment_acc"])):
    study = r["study"]
    patient_key = r["v5_patient_key_LOCAL_ONLY"]
    timing = r["v5_treatment_timing_LOCAL_ONLY"]

    assert patient_key.strip() and timing.strip()
    assert r["v6_bucket"] == "HISTORICAL_EVIDENCE_REVIEW_PENDING"
    assert r["v5_tissue_evidence_grade"] == "PRIMARY_TUMOUR_NOT_ATTESTED"

    groups[(study, patient_key)].append(r["experiment_acc"])

    timing_digest = hashlib.sha256(
        timing.encode()
    ).hexdigest()[:12]

    timing_cross[
        (r["v5_treatment_evidence_grade"], timing_digest)
    ] += 1

    c = by_study[study]
    c["n"] += 1
    c["patient:" + r["v5_patient_evidence_grade"]] += 1
    c["treatment:" + r["v5_treatment_evidence_grade"]] += 1

    # Preserve original clinical evidence in LOCAL-ONLY output.
    # These review fields require external verification by a reviewer.
    queue.append({
        "study": study,
        "experiment_acc": r["experiment_acc"],
        "sra_sample_acc": r["sra_sample_acc"],
        "original_sample_title_LOCAL_ONLY":
            r["original_sample_title_LOCAL_ONLY"],
        "source_clinical_values_LOCAL_ONLY":
            r["original_clinical_values_json_LOCAL_ONLY"],
        "historical_tissue_grade":
            r["v5_tissue_evidence_grade"],
        "historical_treatment_grade":
            r["v5_treatment_evidence_grade"],
        "historical_patient_grade":
            r["v5_patient_evidence_grade"],
        "v5_patient_key_LOCAL_ONLY": patient_key,
        "v5_treatment_timing_LOCAL_ONLY": timing,
        "primary_sampled_site_review": "NOT_ADJUDICATED",
        "primary_tumour_origin_review": "NOT_ADJUDICATED",
        "collection_vs_treatment_review": "NOT_ADJUDICATED",
        "patient_identity_review": "NOT_ADJUDICATED",
        "within_study_patient_group_review": "NOT_ADJUDICATED",
        "cross_study_independence_review": "NOT_ADJUDICATED",
        "technical_release_review": "NOT_ADJUDICATED",
        "source_reference_LOCAL_ONLY": "",
        "reviewer_reason_LOCAL_ONLY": "",
        "provisional_disposition": "HOLD_PENDING_CLINICAL_ADJUDICATION",
        "final_clinical_approval": "NO_UNSIGNED",
    })

reused_by_study = Counter()
for (study, key), accessions in groups.items():
    if len(accessions) > 1:
        reused_by_study[study] += 1

assert len(groups) == 591
assert sum(reused_by_study.values()) == 8
assert dict(reused_by_study) == {
    "SRP058107": 4,
    "SRP089805": 3,
    "SRP131744": 1,
}

study_rows = []
for study, c in sorted(by_study.items()):
    record = {
        "study": study,
        "n_experiments": c["n"],
        "n_reused_patient_key_groups": reused_by_study[study],
    }
    for code in P:
        record["n_patient_" + code] = c["patient:" + code]
    for code in T:
        record["n_treatment_" + code] = c["treatment:" + code]
    record["clinical_status"] = "UNSIGNED"
    study_rows.append(record)

assert len(study_rows) == 39
assert sum(x["n_experiments"] for x in study_rows) == 599
assert len({v for _, v in timing_cross}) == 2

out.mkdir(parents=True, exist_ok=True)

def write_tsv(path, records):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(records[0]),
            delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)

local = out / "599_clinical_review_queue_LOCAL_ONLY.tsv"
study_file = out / "39_study_gate_summary_UNSIGNED.tsv"
write_tsv(local, queue)
write_tsv(study_file, study_rows)

summary = {
    "status": "REVIEW_QUEUE_ASSEMBLED_NOT_ADJUDICATED",
    "input_sha256": sha(source),
    "n_experiments": 599,
    "n_studies": 39,
    "patient_evidence_grades": dict(pg),
    "treatment_evidence_grades": dict(tg),
    "n_study_scoped_patient_keys": 591,
    "n_reused_patient_key_groups": 8,
    "treatment_grade_by_timing_value_hash": {
        f"{grade}|{digest}": n
        for (grade, digest), n in sorted(timing_cross.items())
    },
    "local_review_queue_sha256": sha(local),
    "study_summary_sha256": sha(study_file),
    "clinical_approvals": 0,
}

report = out / "summary.json"
report.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

print("V6_599_REVIEW_QUEUE_ASSEMBLED")
print("Experiments:", 599)
print("Studies:", 39)
print("Distinct study-scoped V5 keys:", len(groups))
print("Reused-key groups:", sum(reused_by_study.values()))
print()
print("TREATMENT_GRADE_BY_TIMING_TOKEN")
for (grade, digest), n in sorted(timing_cross.items()):
    print(grade, "SHA256_" + digest, n, sep="\t")

print()
print("39-study summary SHA256:", sha(study_file))
print("599-record local queue SHA256:", sha(local))
print("Aggregate summary SHA256:", sha(report))
print("Clinical approvals: 0")
print("REVIEW_QUEUE_READY_FOR_SOURCE_ADJUDICATION")
