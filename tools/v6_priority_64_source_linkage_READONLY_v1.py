#!/usr/bin/env python3
"""Read-only exact experiment join for four prioritized clinical studies.
No automatic specimen adjudication, approvals or source modifications.
"""
import csv, hashlib, json, re, sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

queue = repo / (
    "diagnostics/20261124_v6_599_adjudication_queue/"
    "599_clinical_review_queue_LOCAL_ONLY.tsv"
)
v5 = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
)
out = repo / "diagnostics/20261124_v6_priority_64_source_linkage"

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

assert sha(queue) == (
    "b5fdf922e4791b4b30c2ffa798c6ffed99b43478a0681c02a1faa2aaa248f3bd"
)
assert sha(v5) == (
    "06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49"
)

def read(p):
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

all_queue = read(queue)
all_v5 = read(v5)

assert len(all_queue) == 599
assert len(all_v5) == 2055

v5_by_id = {r["library_id"]: r for r in all_v5}
assert len(v5_by_id) == 2055

priority = {
    "SRP094795": 23,
    "SRP072494": 20,
    "SRP194626": 4,
    "SRP112751": 17,
}

# Exact GEO sample lists taken from the original series pages.
# Membership alone is not independently verified GSM-to-SRX linkage.
geo_sets = {
    "SRP094795": {
        f"GSM{x}" for x in range(2420393, 2420417)
    },
    "SRP112751": {
        f"GSM{x}" for x in range(2706444, 2706461)
    },
}

selected = [
    r for r in all_queue if r["study"] in priority
]
assert Counter(r["study"] for r in selected) == priority
assert len(selected) == 64
assert len({r["experiment_acc"] for r in selected}) == 64
assert all(
    r["final_clinical_approval"] == "NO_UNSIGNED"
    for r in selected
)

local_rows = []
per_study = defaultdict(Counter)
aliases = defaultdict(list)
keys = defaultdict(list)

for r in sorted(selected, key=lambda x: (
    x["study"], x["experiment_acc"]
)):
    study = r["study"]
    acc = r["experiment_acc"]
    old = v5_by_id[acc]
    assert old["study"] == study
    assert old["FINAL_CLINICAL_APPROVAL"] == "NO_UNSIGNED"

    original_sample = r["sra_sample_acc"]
    v5_sample = old.get("source_SRA_sample_acc", "")

    if original_sample and v5_sample:
        assert original_sample == v5_sample, (
            "SRA sample mismatch", study, acc
        )

    gsm_alias = old.get("v5_original_sra_sample_name", "").strip()
    is_gsm = bool(re.fullmatch(r"GSM[0-9]+", gsm_alias))

    if not gsm_alias:
        alias_grade = "GSM_ALIAS_MISSING"
    elif not is_gsm:
        alias_grade = "ORIGINAL_ALIAS_NOT_GSM"
    elif study in geo_sets:
        alias_grade = (
            "GSM_IN_ORIGINAL_GEO_SAMPLE_LIST"
            if gsm_alias in geo_sets[study]
            else "GSM_NOT_IN_EXPECTED_GEO_SAMPLE_LIST"
        )
    else:
        alias_grade = "GSM_ALIAS_PRESENT_SERIES_LINK_NOT_CHECKED"

    a = json.loads(r["source_clinical_values_LOCAL_ONLY"])
    assert isinstance(a, dict)

    c = per_study[study]
    c["n"] += 1
    c["alias:" + alias_grade] += 1

    for grade_field in (
        "historical_treatment_grade",
        "historical_patient_grade",
    ):
        c[grade_field + ":" + r[grade_field]] += 1

    for field in (
        "sample from primary gbm diagnosis",
        "time point",
        "tissue",
        "histological tumor classification",
    ):
        values = a.get(field, [])
        if len(values) == 1:
            c["source:" + field + ":" + values[0]] += 1

    aliases[study].append(gsm_alias)
    keys[(study, r["v5_patient_key_LOCAL_ONLY"])].append(acc)

    local_rows.append({
        "study": study,
        "experiment_acc": acc,
        "sra_sample_acc": original_sample,
        "v5_GSM_alias_LOCAL_ONLY": gsm_alias,
        "GSM_alias_check_NOT_CERTIFICATION": alias_grade,
        "v5_original_sample_title_LOCAL_ONLY":
            old.get("v5_original_sra_sample_title", ""),
        "submitted_sample_title_LOCAL_ONLY":
            r["original_sample_title_LOCAL_ONLY"],
        "submitted_clinical_values_LOCAL_ONLY":
            r["source_clinical_values_LOCAL_ONLY"],
        "patient_key_LOCAL_ONLY":
            r["v5_patient_key_LOCAL_ONLY"],
        "historical_treatment_grade":
            r["historical_treatment_grade"],
        "historical_patient_grade":
            r["historical_patient_grade"],
        "source_sample_join": (
            "MATCHED" if original_sample and v5_sample
            else "ONE_OR_BOTH_SAMPLE_FIELDS_EMPTY"
        ),
        "specimen_review": "NOT_ADJUDICATED",
        "treatment_review": "NOT_ADJUDICATED",
        "patient_review": "NOT_ADJUDICATED",
        "final_approval": "NO_UNSIGNED",
    })

assert sum(c["n"] for c in per_study.values()) == 64

summary = {
    "status": "SOURCE_ALIAS_REVIEW_COMPLETE_NOT_CERTIFIED",
    "input_queue_sha256": sha(queue),
    "input_v5_sha256": sha(v5),
    "n_priority_experiments": 64,
    "n_priority_studies": 4,
    "studies": {},
    "clinical_approvals": 0,
}

for study, n in priority.items():
    c = per_study[study]
    present = [x for x in aliases[study] if x]
    reused = sum(
        len(v) > 1 for (s, _), v in keys.items()
        if s == study
    )
    summary["studies"][study] = {
        "n_experiments": n,
        "n_nonempty_GSM_aliases": len(present),
        "n_unique_nonempty_GSM_aliases": len(set(present)),
        "n_GSM_in_expected_series_list":
            c["alias:GSM_IN_ORIGINAL_GEO_SAMPLE_LIST"],
        "n_GSM_outside_expected_series_list":
            c["alias:GSM_NOT_IN_EXPECTED_GEO_SAMPLE_LIST"],
        "n_GSM_alias_missing": c["alias:GSM_ALIAS_MISSING"],
        "n_alias_non_GSM": c["alias:ORIGINAL_ALIAS_NOT_GSM"],
        "n_GSM_alias_present_unchecked":
            c["alias:GSM_ALIAS_PRESENT_SERIES_LINK_NOT_CHECKED"],
        "n_reused_study_scoped_V5_keys": reused,
        "source_field_counts": {
            k: v for k, v in sorted(c.items())
            if k.startswith("source:")
        },
    }

out.mkdir(parents=True, exist_ok=True)

local = out / "64_source_linkage_LOCAL_ONLY.tsv"
with local.open("w", newline="") as f:
    w = csv.DictWriter(
        f, fieldnames=list(local_rows[0]),
        delimiter="\t", lineterminator="\n"
    )
    w.writeheader()
    w.writerows(local_rows)

summary["local_detail_sha256"] = sha(local)

aggregate = out / "summary.json"
aggregate.write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)

print("V6_PRIORITY_64_SOURCE_LINKAGE_REVIEW")
for study, d in summary["studies"].items():
    print(
        study, "N", d["n_experiments"],
        "GSM_present", d["n_nonempty_GSM_aliases"],
        "GSM_in_expected_GEO_list",
        d["n_GSM_in_expected_series_list"],
        "GSM_outside_list",
        d["n_GSM_outside_expected_series_list"],
        "reused_V5_key_groups",
        d["n_reused_study_scoped_V5_keys"],
        sep="\t"
    )
    print("SOURCE_FIELDS", study, d["source_field_counts"])

print("LOCAL_DETAIL_SHA256", sha(local))
print("AGGREGATE_SHA256", sha(aggregate))
print("GEO sample-list consistency is NOT independent SRX/GSM certification")
print("NO CLINICAL ADMISSIONS")
print("PRIORITY_64_SOURCE_LINKAGE_COMPLETE")
