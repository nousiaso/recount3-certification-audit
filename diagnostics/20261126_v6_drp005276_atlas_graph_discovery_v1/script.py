#!/usr/bin/env python3
"""DRP005276: V5 atlas flags and archived identity-file discovery.

Read-only. Python 3.6 compatible.
No individual patient, sample or experiment values printed.
"""

import csv
import hashlib
import os
import re
import sys

from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
rec = root / "26_recount3_external_primary_v1"
repo = root / "RECOUNT3_CERTIFICATION_GIT"

v5 = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_"
    "CROSSWALK_UNSIGNED.tsv"
)

source = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)

expected = {
    v5:
        "ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a",
    source:
        "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb"
}

def digest(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(1048576), b""
        ):
            h.update(block)

    return h.hexdigest()

def read_table(path):
    with path.open(
        encoding="utf-8-sig", newline=""
    ) as f:
        reader = csv.DictReader(
            f, delimiter="\t"
        )
        return reader.fieldnames, list(reader)

def safe_state(value):
    """Show short status enums, hash other values.

    The hash is a value-category marker, not a patient ID.
    """
    s = " ".join(
        str(value or "").strip().split()
    )

    if not s:
        return "[EMPTY]"

    # Allow only short text-like categorical values.
    # No accession-shaped strings or free text with digits.
    if (
        len(s) <= 72
        and re.fullmatch(
            r"[A-Za-z][A-Za-z_ -]*", s
        )
    ):
        return s.upper().replace(" ", "_")

    return (
        "[VALUE_SHA256_" +
        hashlib.sha256(
            s.encode("utf-8")
        ).hexdigest()[:12] +
        "]"
    )

print("V6_DRP005276_ATLAS_GRAPH_DISCOVERY_V1")

# ------------------------------------------------------------
# 1. Validate immutable original files.
# ------------------------------------------------------------

for path, sha in expected.items():
    assert path.is_file(), "SOURCE_FILE_MISSING"
    assert digest(path) == sha, \
        "FROZEN_SOURCE_CHECKSUM_CHANGED"

print("FROZEN_INPUT_SHA256_PASS", len(expected))

fields, rows = read_table(v5)
_, queue = read_table(source)

assert len(rows) == 5241
assert len(queue) == 599
assert fields is not None

selected = {
    r["experiment_acc"]
    for r in queue
    if r["study"] == "DRP005276"
}

study = {
    r["library_id"]: r
    for r in rows
    if r["study"] == "DRP005276"
}

assert len(selected) == 83
assert len(study) == 103
assert selected <= set(study)

held = set(study) - selected
assert len(held) == 20

print(
    "V5_STUDY_ROWS", len(study),
    "SELECTED", len(selected),
    "QC_HELD", len(held)
)

# ------------------------------------------------------------
# 2. Inspect actual atlas/overlap flag values.
#
# Previous diagnostic counted YES and NO but did not
# decode nonboolean categorical status values.
# ------------------------------------------------------------

print("\nV5_ATLAS_FLAG_VALUE_SEMANTICS")

flags = sorted(
    f for f in fields
    if "atlas" in f.lower()
    or "overlap" in f.lower()
)

print("FLAG_FIELDS", len(flags))

for field in flags:

    for group, ids in (
        ("SELECTED", selected),
        ("QC_HELD", held)
    ):
        values = Counter(
            safe_state(
                study[x].get(field, "")
            )
            for x in ids
        )

        print(
            "FLAG", field,
            "GROUP", group,
            "STATES",
            dict(sorted(values.items()))
        )

# ------------------------------------------------------------
# 3. Inspect availability of source-identifier fields.
#
# Only column names, populated counts and distinct
# value counts are printed.
# ------------------------------------------------------------

print("\nV5_IDENTITY_FIELD_COVERAGE")

keywords = (
    "biosample",
    "sample_acc",
    "experiment_acc",
    "source_sra",
    "gsm",
    "geo",
    "run_acc"
)

identity_fields = sorted(
    f for f in fields
    if any(k in f.lower() for k in keywords)
)

print(
    "CANDIDATE_IDENTITY_FIELDS",
    len(identity_fields)
)

for field in identity_fields:

    values = [
        str(study[x].get(field, "")).strip()
        for x in selected
    ]

    nonempty = [
        x for x in values if x
    ]

    print(
        "ID_FIELD", field,
        "POPULATED_SELECTED", len(nonempty),
        "DISTINCT_VALUES", len(set(nonempty))
    )

# ------------------------------------------------------------
# 4. Locate historical identity-crosswalk files.
#
# IMPORTANT:
# This is only a filename inventory.
# No candidate file is read or assumed authoritative.
#
# We search known historical roots with bounded depth.
# Large raw-study/matrix directories are skipped.
# ------------------------------------------------------------

print("\nAUTHORITATIVE_CROSSWALK_FILE_DISCOVERY")

bases = [
    (rec, 5),
    (
        root / "27_recount3_final_cert_v1_0_0",
        4
    ),
    (
        root / "28_recount3_final_cert_v1_0_1",
        4
    ),
    (
        root / "29_recount3_adjudication_20261119",
        4
    ),
    (
        root / "RELEASES",
        7
    )
]

pattern = re.compile(
    r"(?:biosample|crosswalk|overlap|"
    r"(?:atlas|archs4).*(?:ident|sample|sra|gsm)|"
    r"(?:ident|gsm|experiment).*(?:atlas|archs4))",
    re.I
)

prune = {
    "raw_studies",
    "__pycache__",
    ".git",
    "matrices",
    "fastq",
    "counts",
    "00_extracted",
    "gene_expression",
    "hdf5"
}

extensions = (
    ".tsv",
    ".tsv.gz",
    ".csv",
    ".csv.gz",
    ".json",
    ".zip",
    ".tar.gz"
)

for base, maxdepth in bases:

    if not base.is_dir():
        print(
            "SEARCH_ROOT_MISSING",
            str(base.relative_to(root))
        )
        continue

    visited = 0
    matches = []
    truncated = False

    for parent, dirs, files in os.walk(
        str(base)
    ):
        depth = len(
            Path(parent).relative_to(base).parts
        )

        dirs[:] = sorted(
            d for d in dirs
            if d not in prune
            and depth < maxdepth
        )

        for name in files:
            visited += 1

            if (
                pattern.search(name)
                and name.lower().endswith(
                    extensions
                )
            ):
                matches.append(
                    str(
                        (
                            Path(parent) / name
                        ).relative_to(root)
                    )
                )

            # Bound the scan to avoid traversing
            # enormous unrelated data trees.
            if visited >= 35000:
                truncated = True
                break

        if truncated:
            break

    print(
        "SEARCH_ROOT",
        str(base.relative_to(root)),
        "FILES_SCANNED", visited,
        "CANDIDATE_FILES", len(matches),
        "SCAN_TRUNCATED", truncated
    )

    for path in sorted(matches)[:60]:
        print("CANDIDATE_PATH", path)

    if len(matches) > 60:
        print(
            "CANDIDATE_PATHS_NOT_PRINTED",
            len(matches) - 60
        )

# ------------------------------------------------------------
# 5. Conservative scientific status.
# ------------------------------------------------------------

print("\nSCIENTIFIC_GATES")

print(
    "IDENTIFIER_GRAPH",
    "SOURCE_FILE_SELECTION_PENDING_SCHEMA_REVIEW"
)

print(
    "ARCHS4_PATIENT_INDEPENDENCE",
    "NOT_CERTIFIED"
)

print(
    "PRETREATMENT_STATUS",
    "NOT_CERTIFIED"
)

print("CLINICAL_ADMISSIONS 0")

print(
    "NO_INDIVIDUAL_PATIENT_SAMPLE_OR_RUN_VALUES_PRINTED"
)

print("NO_FROZEN_DATA_MODIFIED")

print(
    "V6_DRP005276_ATLAS_GRAPH_DISCOVERY_V1_COMPLETE"
)
