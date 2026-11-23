#!/usr/bin/env python3
"""Read-only audit of repeated experiment rows across the 599-case queue."""
import csv, gzip, hashlib, json, sys
from pathlib import Path
from collections import defaultdict, Counter

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"
v6 = repo / (
    "diagnostics/20261123_v6_evidence_precedence_dryrun/"
    "proposed_evidence_2055_UNSIGNED.tsv"
)
expected = (
    "08efe2c8f469edcbdfb98231dc87f705d801434e78"
    "dcbea85f57c4b94da88eb0"
)
assert hashlib.sha256(v6.read_bytes()).hexdigest() == expected

with v6.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

selected = [
    r for r in rows
    if r["unsigned_evidence_bucket"]
       == "HISTORICAL_EVIDENCE_REVIEW_PENDING"
    and r["historical_tissue_grade"]
       == "PRIMARY_TUMOUR_NOT_ATTESTED"
]
assert len(selected) == 599

by_study = defaultdict(set)
for r in selected:
    by_study[r["study"]].add(r["library_id"])

def fingerprint(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True,
                   ensure_ascii=False).encode("utf-8")
    ).hexdigest()

duplicates = 0
conflicts = 0
missing = 0
matched = 0
problem_studies = 0

print("SOURCE_DUPLICATE_AUDIT_V1")
print("Selected experiments:", len(selected))
print("Studies:", len(by_study))

for study, wanted in sorted(by_study.items()):
    path = rec / (
        f"raw_studies/{study}/metadata/sra.sra.{study}.MD.gz"
    )
    if not path.is_file():
        print("MISSING_SOURCE", study, len(wanted), sep="\t")
        missing += len(wanted)
        problem_studies += 1
        continue

    observations = defaultdict(list)
    with gzip.open(path, "rt", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []
        if "experiment_acc" not in fields:
            print("UNSUPPORTED_SCHEMA", study, sep="\t")
            missing += len(wanted)
            problem_studies += 1
            continue

        runfields = [
            x for x in fields if "run" in x.lower()
        ]

        for row in reader:
            acc = row.get("experiment_acc", "")
            if acc not in wanted:
                continue

            specimen = {
                k: row.get(k, "")
                for k in (
                    "sample_acc", "biosample",
                    "sample_title", "sample_attributes"
                )
                if k in fields
            }
            observations[acc].append({
                "full_sha": fingerprint(row),
                "specimen_sha": fingerprint(specimen),
                "sample_acc": row.get("sample_acc", ""),
                "biosample": row.get("biosample", ""),
                "run_columns": {
                    k: row.get(k, "") for k in runfields
                },
            })

    matched += len(observations)
    missing_here = wanted - set(observations)
    missing += len(missing_here)
    repeated = {
        acc: items for acc, items in observations.items()
        if len(items) > 1
    }

    if repeated or missing_here:
        problem_studies += 1
        print(
            "STUDY", study,
            "selected", len(wanted),
            "matched", len(observations),
            "repeated_experiments", len(repeated),
            "missing", len(missing_here),
            sep="\t"
        )

    for acc, items in sorted(repeated.items()):
        duplicates += 1
        n_full = len({x["full_sha"] for x in items})
        n_spec = len({x["specimen_sha"] for x in items})
        if n_spec > 1:
            conflicts += 1

        print(
            "DUPLICATE", study, acc,
            "rows", len(items),
            "unique_full_rows", n_full,
            "unique_specimen_metadata", n_spec,
            sep="\t"
        )
        for i, x in enumerate(items, 1):
            print(
                "ROW_SIGNATURE",
                acc,
                i,
                x["full_sha"][:16],
                x["specimen_sha"][:16],
                "SRS=" + x["sample_acc"],
                "BioSample=" + x["biosample"],
                "run_fields=" + json.dumps(
                    x["run_columns"], sort_keys=True
                ),
                sep="\t"
            )

print("TOTAL_EXACT_EXPERIMENT_MATCHES", matched)
print("TOTAL_UNMATCHED_EXPERIMENTS", missing)
print("STUDIES_WITH_DUPLICATE_OR_MISSING_RECORDS", problem_studies)
print("EXPERIMENTS_WITH_MULTIPLE_SOURCE_ROWS", duplicates)
print("DUPLICATES_WITH_DIFFERENT_SPECIMEN_METADATA", conflicts)
print("NO_DUPLICATES_COLLAPSED")
print("READ_ONLY_COMPLETE")
