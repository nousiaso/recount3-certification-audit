#!/usr/bin/env python3
"""Recover selected historical audit evidence without scientific execution."""

import csv
import hashlib
import sys
import tarfile
import zipfile
from pathlib import Path

rec = Path(sys.argv[1])
repo = Path(sys.argv[2])

# Archive hashes established by the actual cluster inventory.
SOURCES = {
    "ARCHS4_EXECUTION_CONTRACT_AUDIT_V1_1_20261119.zip": (
        "46e7e06fc553d69f5a2651aaf46cd75c4cef1387585188066f32486555b980ce",
        [
            "EXECUTION_CONTRACT_SOURCE_AUDIT_V1_1.md",
            "STAGE11_SOURCE_ADJUDICATION.json",
            "stage11_reviewed_vs_authoritative_frozen.diff",
        ],
    ),
    "RECOUNT3_W34_ATLAS_CANDIDATES.tar.gz": (
        "737d81fbfe4ca75445d6218f7c43eaca63281dee426ad22d0c347817d74e5873",
        [
            "03_readonly_probe_summary.json",
        ],
    ),
    "RECOUNT3_SOURCE_EVIDENCE_3602.tar.gz": (
        "1f5737e7b10bdfa6ac9c884d246e7d84e949e69fe3a87bfad22a60dd63081788",
        [
            "source_evidence_summary.json",
        ],
    ),
    "RECOUNT3_CERTIFICATION_ADJUDICATION_20261119.zip": (
        "3de638a96cf1f6b851fe2529832898414ea4a9fb17aeb7def9ac7800bda00ea5",
        [
            "probe_W34_and_ARCHS4_atlas_READONLY.py",
            "build_adjudication.py",
            "test_reconciled_review.py",
            "README_CLUSTER_AND_SCIENCE.txt",
            "09_adjudication_summary.json",
        ],
    ),
    "RECOUNT3_CERTIFICATION_V1_0_1_20261118.zip": (
        "0fedf1f32c33f73d2a2bbdf459c98e38779319b9a277c3f9cff6ef7a291e7f61",
        [
            "_v100_reference_readonly.py",
            "certify.py",
            "test_v101.py",
            "README_START_HERE.txt",
            "SCIENCE_AND_RELEASE_POLICY_v1_0_1.txt",
        ],
    ),
    "RECOUNT3_STAGE41_V101_REVIEW_RESULTS.tar.gz": (
        "258c5f0ae5d337f3d78321cd5de9f9aa18bb04722cd5fac8e7f4b7b49751c2f1",
        [
            "README.md",
            "10_summary.json",
        ],
    ),
    "FUSCC_COHORT_RECONCILIATION_UNSIGNED_20261120.zip": (
        "642d372af2fb6adc3d2f7a7033863b4479a67ca898faf990ff9d459a13c3ee4b",
        [
            "06_FUSCC_original_hold_QC_readonly.py",
            "00_FUSCC_SCIENTIFIC_ASSESSMENT.md",
            "05_summary.json",
        ],
    ),
    "FUSCC_MODEL_COMPATIBILITY_EVIDENCE.tar.gz": (
        "442768211e425f698882e6c7eaac90cbcac141c2afb010dd9520301abdb794e1",
        [
            "00_preflight_summary.json",
            "04_final_gate_summary.json",
            "bridge_report.json",
            "gene_axis_diagnostic_report.json",
        ],
    ),
    "FUSCC_FIXED_W_SOURCE_REVIEW_AND_PREFLIGHT_20261120.zip": (
        "9fed98cc9f4b86f2184629d2b1a6006ad9a896c11e7f3e45aa672125ab6dccbf",
        [
            "03_exact_cluster_transfer_contract_READONLY.py",
            "00_source_anchored_scientific_review.md",
            "01_nine_annotation_discrepancy_W_loadings.tsv",
            "02_nine_gene_bridge_observations.tsv",
        ],
    ),
}

def sha(data):
    return hashlib.sha256(data).hexdigest()

def label(name):
    return name.removesuffix(".tar.gz").removesuffix(".zip")

def basename(name):
    return Path(name).name.replace("\\.", ".")

def read_members(path, wanted):
    result = {}
    if path.name.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            for info in z.infolist():
                key = basename(info.filename)
                if info.is_dir() or key not in wanted:
                    continue
                if key in result:
                    raise RuntimeError(f"Duplicate member: {key}")
                if info.file_size > 250000:
                    raise RuntimeError(f"Oversized member: {key}")
                result[key] = (info.filename, z.read(info))
    else:
        with tarfile.open(path, "r:*") as t:
            for info in t.getmembers():
                key = basename(info.name)
                if not info.isfile() or key not in wanted:
                    continue
                if key in result:
                    raise RuntimeError(f"Duplicate member: {key}")
                if info.size > 250000:
                    raise RuntimeError(f"Oversized member: {key}")
                handle = t.extractfile(info)
                if handle is None:
                    raise RuntimeError(f"Unreadable: {key}")
                with handle:
                    result[key] = (info.name, handle.read())
    return result

planned = []

# Preflight all archives and target files before writing.
for archive, (expected_hash, wanted) in SOURCES.items():
    path = rec / archive
    if not path.is_file():
        raise SystemExit(f"MISSING ARCHIVE: {path}")

    archive_hash = sha(path.read_bytes())
    if archive_hash != expected_hash:
        raise SystemExit(f"ARCHIVE HASH MISMATCH: {archive}")

    found = read_members(path, set(wanted))

    for filename in wanted:
        if filename not in found:
            raise SystemExit(
                f"MISSING MEMBER: {archive} / {filename}"
            )

        member_path, data = found[filename]

        is_code = filename.lower().endswith(
            (".py", ".r", ".sh", ".rmd")
        )

        area = (
            "scripts/originals/archived_packages"
            if is_code else
            "evidence/history/archived_packages"
        )

        destination = repo / area / label(archive) / filename

        if destination.exists():
            if sha(destination.read_bytes()) != sha(data):
                raise SystemExit(
                    f"CONFLICTING EXISTING FILE: {destination}"
                )

        planned.append((
            archive, archive_hash, member_path,
            destination, data
        ))

print("PREFLIGHT PASS")
print("Archives verified:", len(SOURCES))
print("Selected files verified:", len(planned))

# Copy only into the audit repository.
rows = []

for archive, archive_hash, member, dest, data in planned:
    dest.parent.mkdir(parents=True, exist_ok=True)

    if not dest.exists():
        dest.write_bytes(data)

    if sha(dest.read_bytes()) != sha(data):
        raise SystemExit(f"COPY FAILURE: {dest}")

    relative = dest.relative_to(repo)

    rows.append([
        archive,
        archive_hash,
        member,
        str(relative),
        len(data),
        sha(data),
        "RECOVERED_NOT_REEXECUTED",
    ])

    print("RECOVERED", relative)

manifest = (
    repo / "docs" /
    "HISTORICAL_ARCHIVE_RECOVERY_MANIFEST.tsv"
)

with manifest.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "archive", "archive_sha256", "original_member",
        "repository_path", "bytes", "sha256",
        "scientific_execution_status",
    ])
    writer.writerows(rows)

print("\nRECOVERY COMPLETE")
print("Manifest:", manifest)
print("Original scientific sources unchanged")
print("No biological or model computations executed")
