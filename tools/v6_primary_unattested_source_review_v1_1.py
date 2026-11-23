#!/usr/bin/env python3
"""Assemble original SRA evidence for 599 historical-review records.
No classification, admission, source modification, or model execution.
"""
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"
out = repo / "diagnostics/20261123_v6_primary_unattested_review_v1_1"

v6 = repo / (
    "diagnostics/20261123_v6_evidence_precedence_dryrun/"
    "proposed_evidence_2055_UNSIGNED.tsv"
)
v5 = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
)

def digest(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def readtsv(p):
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

assert digest(v6) == (
    "08efe2c8f469edcbdfb98231dc87f705d801434e78"
    "dcbea85f57c4b94da88eb0"
)
assert digest(v5) == (
    "06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5"
    "fc09ff9c167d0e2e9cf49e49"
)

rows = readtsv(v6)
old = readtsv(v5)
assert len(rows) == len(old) == 2055

original = {r["library_id"]: r for r in old}
assert len(original) == 2055
assert set(original) == {r["library_id"] for r in rows}

selected = [
    r for r in rows
    if (
        r["unsigned_evidence_bucket"]
        == "HISTORICAL_EVIDENCE_REVIEW_PENDING"
        and r["historical_tissue_grade"]
        == "PRIMARY_TUMOUR_NOT_ATTESTED"
    )
]
assert len(selected) == 599

by_study = defaultdict(list)
for row in selected:
    by_study[row["study"]].append(row)

results = []
summaries = []

for study in sorted(by_study):
    target = {r["library_id"] for r in by_study[study]}
    path = rec / (
        f"raw_studies/{study}/metadata/sra.sra.{study}.MD.gz"
    )

    found = {}
    row_counts = Counter()
    source_status = "READABLE"
    source_sha = ""

    if not path.is_file():
        source_status = "RAW_SOURCE_FILE_MISSING"
    else:
        source_sha = digest(path)

        with gzip.open(
            path, "rt", encoding="utf-8",
            errors="replace", newline=""
        ) as f:
            reader = csv.DictReader(f, delimiter="\t")
            fields = set(reader.fieldnames or [])

            required = {
                "experiment_acc", "sample_title",
                "sample_attributes"
            }

            if not required <= fields:
                source_status = (
                    "UNSUPPORTED_SCHEMA:" +
                    ",".join(sorted(required - fields))
                )
            else:
                for raw in reader:
                    accession = raw["experiment_acc"]
                    if accession not in target:
                        continue
                    if accession in found:
                        previous = found[accession]
                        specimen_fields = (
                            "sample_acc", "biosample",
                            "sample_title", "sample_attributes"
                        )
                        if any(
                            previous.get(k, "") != raw.get(k, "")
                            for k in specimen_fields
                        ):
                            raise RuntimeError(
                                f"SPECIMEN_METADATA_CONFLICT: "
                                f"{study} {accession}"
                            )
                    else:
                        found[accession] = raw
                    row_counts[accession] += 1

    keys = Counter()
    example_titles = []

    for row in by_study[study]:
        accession = row["library_id"]
        raw = found.get(accession, {})
        attributes = raw.get("sample_attributes", "")

        for segment in attributes.split("|"):
            if ";;" in segment:
                key, _ = segment.split(";;", 1)
                keys[key.strip().lower()] += 1

        if raw and len(example_titles) < 3:
            example_titles.append((
                accession,
                raw.get("sample_title", ""),
                attributes[:350],
            ))

        results.append({
            "study": study,
            "experiment_acc": accession,
            "raw_source_status": (
                "EXACT_EXPERIMENT_MATCH" if raw
                else source_status if source_status != "READABLE"
                else "EXPERIMENT_NOT_FOUND_IN_RAW_SOURCE"
            ),
            "sra_sample_acc": raw.get("sample_acc", ""),
            "original_sample_title":
                raw.get("sample_title", ""),
            "original_sample_attributes": attributes,
            "stage41_tissue_source_evidence":
                original[accession].get(
                    "stage41_tissue_source_evidence", ""
                ),
            "historical_tissue_grade":
                row["historical_tissue_grade"],
            "treatment_grade":
                row["treatment_grade_UNREVIEWED"],
            "patient_grade":
                row["patient_grade_UNREVIEWED"],
            "source_file": str(path),
            "source_file_sha256": source_sha,
            "n_original_metadata_rows": row_counts.get(accession, 0),
            "clinical_approval": "NO_UNSIGNED",
        })

    summaries.append({
        "study": study,
        "n_historical_unattested": len(target),
        "n_exact_sra_matches": len(found),
        "n_unmatched": len(target) - len(found),
        "n_multirow_experiments": sum(
            x > 1 for x in row_counts.values()
        ),
        "n_additional_metadata_rows": sum(
            x - 1 for x in row_counts.values()
        ),
        "source_status": source_status,
        "source_sha256": source_sha,
        "observed_attribute_keys":
            ";".join(
                f"{key}:{n}"
                for key, n in keys.most_common(15)
            ),
    })

    print("\nSTUDY", study)
    print("SELECTED", len(target))
    print("EXACT_RAW_MATCHES", len(found))
    print("SOURCE_STATUS", source_status)
    print("ATTRIBUTE_KEYS", keys.most_common(12))

    # Patient/sample values retained only in the local detail file.

# Independent consistency expectations from the completed duplicate audit.
assert len(results) == 599
assert len(summaries) == 39
assert sum(x["n_exact_sra_matches"] for x in summaries) == 599
assert sum(x["n_unmatched"] for x in summaries) == 0
assert sum(x["n_multirow_experiments"] for x in summaries) == 59
assert sum(x["n_additional_metadata_rows"] for x in summaries) == 147

out.mkdir(parents=True, exist_ok=True)

def write(name, records):
    p = out / name
    with p.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=list(records[0]),
            delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)
    return p

detail = write(
    "source_evidence_599_LOCAL_ONLY.tsv", results
)
summary = write(
    "study_source_review_summary.tsv", summaries
)

print("\n=== FINAL AUDIT SUMMARY ===")
print("Historical unattested experiments:", len(results))
print("Studies represented:", len(summaries))
print("Exact raw SRA matches:", sum(
    x["n_exact_sra_matches"] for x in summaries
))
print("Unmatched experiments:", sum(
    x["n_unmatched"] for x in summaries
))
print("Detail SHA256:", digest(detail))
print("Summary SHA256:", digest(summary))
print("No specimen classifications changed")
print("No clinical approvals")
print("Multirow experiments:", sum(
    x["n_multirow_experiments"] for x in summaries
))
print("Additional metadata rows:", sum(
    x["n_additional_metadata_rows"] for x in summaries
))
print("Specimen fields concordant; official run accessions not certified")
print("SOURCE_REVIEW_ASSEMBLY_COMPLETE_V1_1")
