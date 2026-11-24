#!/usr/bin/env python3
"""Cross-tab existing 52-record SRA census; no source re-extraction."""
import csv, hashlib, json, sys
from pathlib import Path
from collections import Counter

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
folder = repo / "diagnostics/20261124_v6_osteosarcoma_source_census"
source = folder / "original_source_membership_LOCAL_ONLY.tsv"
out = repo / "diagnostics/20261124_v6_osteosarcoma_membership_crosstab"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

assert sha(source) == (
    "c9f871261e9a64fd7857d443a9ebc14d86e480e0d538c987c506197229703316"
)

with source.open(newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == len({r["experiment_acc"] for r in rows}) == 52
assert Counter(r["in_v5"] for r in rows) == {"YES":39, "NO":13}
assert all(r["specimen_field_conflict"] == "NO" for r in rows)

partitions = [
    "NO_DETECTED_STUDY_ALIAS_OVERLAP_PENDING_GATES",
    "QC_HOLD_NOT_ADMITTED",
    "NOT_IN_STAGE41_CANDIDATES",
    "NOT_IN_FROZEN_V5_UNIVERSE",
]
expected_partitions = [22, 9, 8, 13]
types = ["primary", "metastasis", "unknown", "cell line", "normal bone"]

by_partition = Counter(r["v5_partition"] for r in rows)
assert [by_partition[p] for p in partitions] == expected_partitions

by_type = Counter(
    r["tumor_type_source_value"].strip().lower() for r in rows
)
assert by_type == {
    "primary": 24, "metastasis": 9, "unknown": 11,
    "cell line": 5, "normal bone": 3
}

cross = Counter(
    (r["v5_partition"], r["tumor_type_source_value"].strip().lower())
    for r in rows
)

assert [cross[(partitions[0], t)] for t in types] == [17,0,5,0,0]

report = {
    "status": "EXACT_SOURCE_MEMBERSHIP_CROSSTAB_NOT_CLINICAL_APPROVAL",
    "input_sha256": sha(source),
    "n_original_experiments": 52,
    "n_in_v5": 39,
    "n_outside_v5": 13,
    "specimen_field_conflicts": 0,
    "categories": types,
    "partitions": {
        p: {
            "n": by_partition[p],
            "source_tumor_type_counts": {
                t: cross[(p,t)] for t in types
            }
        }
        for p in partitions
    },
    "clinical_approvals": 0,
}

out.mkdir(parents=True, exist_ok=True)
summary = out / "summary.json"
summary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

print("OSTEOSARCOMA_MEMBERSHIP_CROSSTAB_PASS")
print("SOURCE_PARTITION\t" + "\t".join(types) + "\tTOTAL")
for p in partitions:
    print(p, *[cross[(p,t)] for t in types],
          by_partition[p], sep="\t")

print("SOURCE_TOTAL", *[by_type[t] for t in types], 52, sep="\t")
print("SUMMARY_SHA256", sha(summary))
print("NO CLINICAL ADMISSIONS; NO SOURCE RE-EXTRACTION")
