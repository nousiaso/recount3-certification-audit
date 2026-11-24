#!/usr/bin/env python3
"""Original SRA versus V5 study membership. No clinical decisions."""
import csv, gzip, hashlib, json, sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
rec = root / "26_recount3_external_primary_v1"
repo = root / "RECOUNT3_CERTIFICATION_GIT"

v5 = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_CROSSWALK_UNSIGNED.tsv"
)
sra = rec / (
    "raw_studies/SRP090849/metadata/"
    "sra.sra.SRP090849.MD.gz"
)
out = repo / "diagnostics/20261124_v6_osteosarcoma_source_census"

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

assert sha(v5) == (
    "ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a"
)

with v5.open(newline="") as f:
    v5_rows = [
        r for r in csv.DictReader(f, delimiter="\t")
        if r["study"] == "SRP090849"
    ]
assert len(v5_rows) == 39
v5_by_id = {r["library_id"]: r for r in v5_rows}
assert len(v5_by_id) == 39

source_rows = defaultdict(list)
with gzip.open(sra, "rt", errors="replace", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    assert {
        "experiment_acc", "sample_acc",
        "sample_title", "sample_attributes"
    } <= set(reader.fieldnames or [])

    for row in reader:
        source_rows[row["experiment_acc"]].append(row)

assert source_rows
summary = Counter()
details = []

for experiment, raw_rows in sorted(source_rows.items()):
    # Inspect specimen concordance without collapsing run evidence.
    specimen_fields = (
        "sample_acc", "sample_title", "sample_attributes"
    )
    signatures = {
        tuple(r.get(k, "") for k in specimen_fields)
        for r in raw_rows
    }
    conflict = len(signatures) > 1
    if conflict:
        summary["specimen_field_conflicts"] += 1

    first = raw_rows[0]
    attrs = defaultdict(set)
    for part in first["sample_attributes"].split("|"):
        if ";;" in part:
            key, value = part.split(";;", 1)
            attrs[key.strip().lower()].add(value.strip())

    def value(key):
        vals = attrs.get(key, set())
        return "|".join(sorted(vals))

    in_v5 = experiment in v5_by_id
    status = (
        v5_by_id[experiment]["stage41_partition_NOT_ADMISSION"]
        if in_v5 else "NOT_IN_FROZEN_V5_UNIVERSE"
    )

    summary["raw_unique_experiments"] += 1
    summary["in_v5" if in_v5 else "not_in_v5"] += 1
    summary["partition:" + status] += 1
    summary["tumor_type:" + (value("tumor type") or "MISSING")] += 1
    summary["tissue:" + (value("tissue") or "MISSING")] += 1

    details.append({
        "experiment_acc": experiment,
        "sra_sample": first["sample_acc"],
        "sample_title": first["sample_title"],
        "tumor_type_source_value": value("tumor type"),
        "tissue_source_value": value("tissue"),
        "in_v5": "YES" if in_v5 else "NO",
        "v5_partition": status,
        "n_source_metadata_rows": len(raw_rows),
        "specimen_field_conflict": "YES" if conflict else "NO",
    })

assert summary["in_v5"] == 39, (
    "V5 experiments not all located in source", dict(summary)
)

out.mkdir(parents=True, exist_ok=True)
local = out / "original_source_membership_LOCAL_ONLY.tsv"

with local.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(details[0]),
        delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(details)

report = {
    "status": "READ_ONLY_SOURCE_CENSUS_NOT_CERTIFICATION",
    "study": "SRP090849",
    "source_sra_sha256": sha(sra),
    "v5_sha256": sha(v5),
    "census": dict(sorted(summary.items())),
    "local_detail_sha256": sha(local),
    "original_geo_series": "GSE87624",
    "geo_displayed_samples": 52,
    "clinical_admissions": 0,
}
summary_file = out / "summary.json"
summary_file.write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n"
)

print("OSTEOSARCOMA_SOURCE_CENSUS_COMPLETE")
for k, n in sorted(summary.items()):
    print(k, n, sep="\t")
print("LOCAL_DETAIL_SHA256", sha(local))
print("SUMMARY_SHA256", sha(summary_file))
print("No source metadata or clinical decisions changed")
