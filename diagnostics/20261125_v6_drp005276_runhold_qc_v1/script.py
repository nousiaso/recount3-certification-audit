#!/usr/bin/env python3
"""DRP005276: frozen V5 versus historical run-hold evidence.

Read-only. No accession, specimen or patient identifiers printed.
No clinical admission. Compatible with Python 3.6.
"""
import csv
import hashlib
import math
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
rec = root / "26_recount3_external_primary_v1"
v5 = rec / "clinical_specimen_patient_v5_unsigned_20261122"

paths = {
    "inventory": (
        v5 / (
            "01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_"
            "CROSSWALK_UNSIGNED.tsv"
        ),
        "ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a"
    ),
    "screen": (
        v5 / "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv",
        "06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49"
    ),
    "readiness": (
        rec / "cohort_readiness_v1/preview/01_per_library_readiness_3602.tsv",
        "404c66603c4d6934e42cfbd6a21eaf0233b02ead9edb568193ba26b144b33d0f"
    ),
    "run_holds": (
        root / (
            "29_recount3_adjudication_20261119/"
            "RECOUNT3_CERTIFICATION_ADJUDICATION_20261119/"
            "REVIEW_OUTPUT/05_run_hold_queue_409.tsv"
        ),
        "d275c77591e7b93582708746e45d69ae817c8e2d13b3bb5b7f63e2e56fa56f0c"
    )
}

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def read(name):
    path, checksum = paths[name]
    assert path.is_file(), "MISSING_FROZEN_INPUT_" + name
    assert digest(path) == checksum, "INPUT_HASH_MISMATCH_" + name

    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert reader.fieldnames, "MISSING_HEADER_" + name
        return list(reader)

def study_records(rows, study):
    out = {}
    for r in rows:
        if r.get("study") != study:
            continue
        acc = r.get("library_id", "")
        assert acc and acc not in out, "DUPLICATE_OR_EMPTY_LIBRARY_ID"
        out[acc] = r
    return out

def restricted_codes(value):
    result = []
    for code in (value or "").split(";"):
        code = code.strip()
        if not code:
            continue
        if re.match(r"^[A-Z][A-Z0-9_]{0,100}$", code):
            result.append(code)
        else:
            result.append("UNREDACTED_CODE_WITHHELD")
    return result

def numeric_summary(rows, field):
    nums = []
    for r in rows:
        try:
            x = float(r.get(field, ""))
            if math.isfinite(x):
                nums.append(x)
        except (ValueError, TypeError):
            pass

    nums.sort()
    print("METRIC", field, "ROWS_WITH_FINITE_VALUE", len(nums))

    if not nums:
        return

    print(
        "METRIC", field, "MIN_MEDIAN_MAX",
        "%.6g" % nums[0],
        "%.6g" % nums[len(nums) // 2],
        "%.6g" % nums[-1]
    )

    if field == "top10_coverage_fraction":
        print("METRIC", field, "ROWS_GT_0_50",
              sum(x > 0.50 for x in nums))

    elif field == "average_mapped_length_bp":
        print("METRIC", field, "ROWS_LT_40",
              sum(x < 40 for x in nums))

    elif field == "gene_detection_fraction":
        print("METRIC", field, "ROWS_LT_0_10",
              sum(x < 0.10 for x in nums))

print("V6_DRP005276_ORIGINAL_RUN_HOLD_RECONCILIATION_V1")

D = {name: read(name) for name in paths}

assert [
    len(D[x]) for x in
    ("inventory", "screen", "readiness", "run_holds")
] == [5241, 2055, 3602, 409]

print("FROZEN_INPUT_SHA256_PASS", len(D))

study = "DRP005276"
inventory = study_records(D["inventory"], study)
screen = study_records(D["screen"], study)
readiness = study_records(D["readiness"], study)

assert len(inventory) == 103
assert len(screen) == 83
assert set(screen) <= set(inventory)

selected = set(screen)
held = set(inventory) - selected

assert len(held) == 20
assert all(
    inventory[x].get("stage41_partition_NOT_ADMISSION")
    == "QC_HOLD_NOT_ADMITTED"
    for x in held
)

print("\nFROZEN_MEMBERSHIP")
print("V5_INVENTORY", len(inventory))
print("V5_SCREEN_SELECTED", len(selected))
print("V5_QC_HOLD_OUTSIDE_SCREEN", len(held))
print("READINESS_STUDY_ROWS", len(readiness))
print("READINESS_SELECTED_MATCHES",
      len(selected & set(readiness)))
print("READINESS_HELD_MATCHES",
      len(held & set(readiness)))
print("READINESS_HELD_MISSING",
      len(held - set(readiness)))

# Match original run-hold rows by exact library ID.
# Multiple runs belonging to one experiment remain grouped.
hold_by_id = defaultdict(list)
study_mismatch = 0

for r in D["run_holds"]:
    acc = r.get("library_id", "")

    if acc not in inventory:
        continue

    if r.get("study") and r["study"] != study:
        study_mismatch += 1
        continue

    hold_by_id[acc].append(r)

print("RUN_HOLD_QUEUE_STUDY_MISMATCH_ROWS", study_mismatch)

for label, ids in (
    ("SELECTED_83", selected),
    ("HELD_20", held)
):
    records = [
        r
        for acc in sorted(ids)
        for r in hold_by_id.get(acc, [])
    ]

    with_holds = {
        acc for acc in ids if hold_by_id.get(acc)
    }

    row_reason_counts = Counter()
    experiment_reason_counts = Counter()

    for r in records:
        row_reason_counts.update(
            restricted_codes(r.get("reason_codes", ""))
        )

    for acc in ids:
        reasons = set()
        for r in hold_by_id.get(acc, []):
            reasons.update(
                restricted_codes(r.get("reason_codes", ""))
            )
        experiment_reason_counts.update(reasons)

    print("\nGROUP", label, "EXPERIMENTS", len(ids))

    print("V5_PARTITION_COUNTS", dict(sorted(Counter(
        inventory[acc].get(
            "stage41_partition_NOT_ADMISSION", "<MISSING>"
        )
        for acc in ids
    ).items())))

    print("V5_DISPOSITION_COUNTS", dict(sorted(Counter(
        inventory[acc].get(
            "recommended_provisional_disposition", "<MISSING>"
        )
        for acc in ids
    ).items())))

    print("READINESS_HISTORICAL_RUN_HOLD",
          dict(sorted(Counter(
              readiness[acc].get(
                  "historical_run_hold", "<BLANK>"
              ) if acc in readiness else "<NO_READINESS_ROW>"
              for acc in ids
          ).items())))

    print(
        "READINESS_OFFICIAL_REASON_FIELD_POPULATED",
        sum(
            bool(readiness[acc].get(
                "official_v101_hold_reasons", ""
            ).strip())
            for acc in ids if acc in readiness
        )
    )

    print(
        "ORIGINAL_HOLD_QUEUE_EXPERIMENTS_MATCHED",
        len(with_holds)
    )
    print(
        "ORIGINAL_HOLD_QUEUE_EXPERIMENTS_NOT_MATCHED",
        len(ids - with_holds)
    )
    print("ORIGINAL_HOLD_QUEUE_ROWS", len(records))

    print(
        "ORIGINAL_HOLD_QUEUE_ROWS_PER_EXPERIMENT",
        dict(sorted(Counter(
            len(hold_by_id.get(acc, []))
            for acc in ids
        ).items()))
    )

    print(
        "QC_REASON_EXPERIMENT_COUNTS",
        dict(sorted(experiment_reason_counts.items()))
    )
    print(
        "QC_REASON_ROW_COUNTS",
        dict(sorted(row_reason_counts.items()))
    )

    for field in (
        "gene_detection_fraction",
        "top10_coverage_fraction",
        "average_mapped_length_bp"
    ):
        numeric_summary(records, field)

print("\nINTERPRETATION_GATES")
print("HOLD_MEMBERSHIP FROZEN_V5_CONFIRMED")
print("QC_REASON_COMPLETENESS ASSESS_FROM_MATCH_COUNTS")
print("TECHNICAL_REINSTATEMENTS 0")
print("PATIENT_INDEPENDENCE NOT_CERTIFIED")
print("PRIMARY_TUMOUR_SITE NOT_CERTIFIED")
print("PRETREATMENT_TIMING NOT_CERTIFIED")
print("CLINICAL_ADMISSIONS 0")
print("NO_PATIENT_OR_SAMPLE_VALUES_PRINTED")
print("NO_FROZEN_FILES_MODIFIED")
print("V6_DRP005276_ORIGINAL_RUN_HOLD_RECONCILIATION_V1_COMPLETE")
