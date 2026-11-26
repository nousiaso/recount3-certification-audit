#!/usr/bin/env python3
"""Inspect original identity-join input availability and schemas.

Only reads file metadata and Stage41 TSV aggregates.
No original historical join or biological inference.
Python 3.6 compatible.
"""

import csv
import gzip
import hashlib
import sys
from pathlib import Path
from collections import Counter

root = Path(sys.argv[1])

bases = [
    ("PROJECT_ROOT", root),
    (
        "FROZEN_PROJECT_TREE",
        root / "RELEASES" /
        "ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0" /
        "PROJECT_TREE"
    ),
    (
        "RECOUNT3_ROOT",
        root / "26_recount3_external_primary_v1"
    ),
]

inputs = [
    (
        "STAGE04",
        "04_candidate_evidence_v1_0_4/"
        "04_candidate_sample_full_metadata.parquet"
    ),
    (
        "STAGE05",
        "05_study_sample_adjudication_v1_0_2_1/"
        "18_complete_adjudication_manifest.parquet"
    ),
    (
        "ATLAS",
        "14_projection_use_analysis_ready_atlas_v1_0_3/"
        "03_projection_analysis_ready_core.tsv.gz"
    ),
    (
        "CANDIDATES",
        "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/"
        "03_candidate_integrated_evidence_3602.tsv"
    ),
    (
        "FULL_SOURCE",
        "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/"
        "01_full_source_experiment_identity_5221.tsv"
    ),
]

try:
    import pyarrow.parquet as pq
    parquet_available = True
except ImportError:
    parquet_available = False

print("ARCHS4_IDENTITY_INPUT_SCHEMA_V1")
print("PYARROW_METADATA_AVAILABLE", parquet_available)
print("ORIGINAL_ATLAS_WIDE_JOIN_EXECUTED NO")

found = Counter()

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def tabular_schema(path, label):
    opener = gzip.open if path.name.endswith(".gz") else open

    with opener(str(path), "rt", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fields = reader.fieldnames or []

        print("TSV_SCHEMA", label, "FIELD_COUNT", len(fields))
        print("TSV_FIELDS", label, ";".join(fields[:70]))

        count = 0
        study_counts = Counter()
        identity_counts = Counter()
        identity_fields = (
            "experiment_acc", "biosample_IDs",
            "biosamples_all_sources", "GEO_GSM_aliases",
            "atlas_exact_matches", "in_candidate_3602",
        )

        for row in reader:
            count += 1

            if row.get("study", "").strip() == "DRP005276":
                study_counts["DRP005276"] += 1

                for field in identity_fields:
                    if field in row and str(row[field]).strip():
                        identity_counts[field] += 1

        print("TSV_ROWS", label, count)
        print(
            "DRP005276_ROWS", label,
            study_counts["DRP005276"]
        )

        for field in identity_fields:
            if field in fields:
                print(
                    "DRP005276_POPULATED",
                    label, field, identity_counts[field]
                )

for base_name, base in bases:
    print("\nBASE", base_name, "EXISTS", base.is_dir())

    if not base.is_dir():
        continue

    for label, relative in inputs:
        path = base / relative

        if not path.is_file():
            print("INPUT_MISSING", base_name, label)
            continue

        found[label] += 1

        print(
            "INPUT_FOUND", base_name, label,
            "BYTES", path.stat().st_size
        )
        print("INPUT_PATH", str(path.relative_to(root)))

        if path.name.endswith(".parquet"):
            if parquet_available:
                metadata = pq.ParquetFile(str(path))
                names = metadata.schema_arrow.names

                print(
                    "PARQUET_SCHEMA", label,
                    "ROWS", metadata.metadata.num_rows,
                    "ROW_GROUPS", metadata.metadata.num_row_groups,
                    "FIELD_COUNT", len(names)
                )
                print(
                    "PARQUET_RELEVANT_FIELDS", label,
                    ";".join(
                        n for n in names
                        if n in (
                            "geo_accession", "sample_id",
                            "relation", "biosample",
                            "experiment_acc"
                        )
                    )
                )
            else:
                print(
                    "PARQUET_SCHEMA_UNAVAILABLE",
                    label, "PYARROW_NOT_INSTALLED"
                )
        else:
            print("INPUT_SHA256", label, digest(path))
            tabular_schema(path, label)

print("\nINPUT_AVAILABILITY_SUMMARY")

for label, _ in inputs:
    print(
        "INPUT", label,
        "LOCATIONS_FOUND", found[label]
    )

print("SCHEMA_DISCOVERY_ONLY YES")
print("ORIGINAL_JOIN_RERUN NO")
print("PATIENT_INDEPENDENCE_CERTIFIED NO")
print("PRETREATMENT_CERTIFIED NO")
print("CLINICAL_ADMISSIONS 0")
print("NO_PATIENT_OR_SAMPLE_VALUES_PRINTED")
print("NO_FROZEN_FILES_MODIFIED")
print("ARCHS4_IDENTITY_INPUT_SCHEMA_V1_COMPLETE")
