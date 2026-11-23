#!/usr/bin/env python3
"""Read-only inventory of historical audit archive contents."""

from pathlib import Path
from collections import Counter
import zipfile
import tarfile
import sys

rec = Path(sys.argv[1])

archives = [
    "ARCHS4_EXECUTION_CONTRACT_AUDIT_V1_1_20261119.zip",
    "RECOUNT3_W34_ATLAS_CANDIDATES.tar.gz",
    "RECOUNT3_SOURCE_EVIDENCE_3602.tar.gz",
    "RECOUNT3_CERTIFICATION_ADJUDICATION_20261119.zip",
    "RECOUNT3_CERTIFICATION_V1_0_1_20261118.zip",
    "RECOUNT3_STAGE41_V101_REVIEW_RESULTS.tar.gz",
    "FUSCC_COHORT_RECONCILIATION_UNSIGNED_20261120.zip",
    "FUSCC_MODEL_COMPATIBILITY_EVIDENCE.tar.gz",
    "FUSCC_FROZEN_W_PROJECTION_METHODS.tar.gz",
    "FUSCC_ORIGINAL_PROJECTION_SOURCE.tar.gz",
    "FUSCC_FIXED_W_SOURCE_REVIEW_AND_PREFLIGHT_20261120.zip",
    "ARCHS4_FROZEN_GENE_PROVENANCE_SOURCE.tar.gz",
]

code_suffixes = (".py", ".r", ".sh", ".rmd")
evidence_suffixes = (
    ".json", ".md", ".txt", ".log",
    ".tsv", ".csv", ".yaml", ".yml"
)

def members(path):
    if path.name.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            return [
                (i.filename, i.file_size)
                for i in z.infolist()
                if not i.is_dir()
            ]

    with tarfile.open(path, "r:*") as t:
        return [
            (i.name, i.size)
            for i in t.getmembers()
            if i.isfile()
        ]

print("=== HISTORICAL ARCHIVE MEMBER INVENTORY ===")
print("No extraction or scientific execution performed")

for name in archives:
    p = rec / name

    print("\n" + "=" * 70)
    print("ARCHIVE:", name)

    if not p.is_file():
        print("STATUS: NOT FOUND")
        continue

    try:
        entries = members(p)
    except Exception as e:
        print("READ ERROR:", repr(e))
        continue

    print("STATUS: READABLE")
    print("FILE COUNT:", len(entries))
    print("UNCOMPRESSED BYTES:", sum(s for _, s in entries))

    suffix_counts = Counter(
        Path(n).suffix.lower() or "NO_EXTENSION"
        for n, _ in entries
    )
    print("FILE TYPES:", dict(suffix_counts))

    code = [
        (n, s) for n, s in entries
        if n.lower().endswith(code_suffixes)
    ]
    evidence = [
        (n, s) for n, s in entries
        if n.lower().endswith(evidence_suffixes)
        and not n.lower().endswith(code_suffixes)
    ]

    print("\nCODE MEMBERS:", len(code))
    for n, s in code[:25]:
        print("  CODE", s, n)

    print("\nEVIDENCE / REPORT MEMBERS:", len(evidence))
    for n, s in evidence[:35]:
        print("  EVIDENCE", s, n)

    other = [
        (n, s) for n, s in entries
        if (n, s) not in code and (n, s) not in evidence
    ]

    print("\nOTHER MEMBERS:", len(other))
    for n, s in other[:12]:
        print("  OTHER", s, n)

    if len(code) > 25 or len(evidence) > 35 or len(other) > 12:
        print("NOTE: Some member names omitted from printout")

print("\nREAD_ONLY_ARCHIVE_INVENTORY_COMPLETE")
