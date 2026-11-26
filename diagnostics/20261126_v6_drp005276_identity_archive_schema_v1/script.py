#!/usr/bin/env python3
"""Inspect archived ARCHS4/recount3 identity evidence.

No extraction, no expression reads, no record-level printing.
Python 3.6 compatible.
"""

import csv
import gzip
import hashlib
import io
import json
import re
import sys
import zipfile

from pathlib import Path, PurePosixPath

root = Path(sys.argv[1])
rec = root / "26_recount3_external_primary_v1"

archive = rec / (
    "ARCHS4_RECOUNT3_EXACT_IDENTITY_JOIN_READONLY_20261120.zip"
)

# Existing historical overlap contracts.
# These are inspected, never modified.
previous = [
    root / (
        "27_recount3_final_cert_v1_0_0/"
        "draft/ARCHS4_overlap_status.tsv"
    ),
    root / (
        "28_recount3_final_cert_v1_0_1/"
        "draft/ARCHS4_overlap_status.tsv"
    ),
    root / (
        "27_recount3_final_cert_v1_0_0/"
        "draft/ARCHS4_overlap_contract.json"
    ),
    root / (
        "28_recount3_final_cert_v1_0_1/"
        "draft/ARCHS4_overlap_contract.json"
    ),
]

# Bound decompression and structured-data inspection.
MAX_READ = 48 * 1024 * 1024
MAX_JSON = 4 * 1024 * 1024

def digest_file(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(1048576), b""
        ):
            h.update(block)

    return h.hexdigest()

def safe_label(value):
    """Redact common identifier patterns from metadata labels."""
    s = str(value)

    s = re.sub(
        r"\b(?:DRX|DRR|DRS|SRX|SRR|SRS|"
        r"ERX|ERR|ERS|GSM|SAMN|SAMD|SAMEA)\d+\b",
        "[ACCESSION]", s, flags=re.I
    )

    s = re.sub(
        r"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}",
        "[EMAIL]", s, flags=re.I
    )

    s = re.sub(
        r"\b(?:patient|subject|isolate)[_-]?\d+\b",
        "[SUBJECT_TOKEN]", s, flags=re.I
    )

    if len(s) > 180:
        s = s[:177] + "..."

    return s

def print_schema(fields, prefix):
    if fields is None:
        print(prefix, "NO_TABLE_HEADER")
        return

    print(prefix, "FIELD_COUNT", len(fields))

    for i, field in enumerate(fields[:65], 1):
        print(
            prefix, "FIELD", i,
            safe_label(field)
        )

    if len(fields) > 65:
        print(
            prefix, "FIELDS_NOT_PRINTED",
            len(fields) - 65
        )

def inspect_table(handle, filename, prefix):
    """Read headers and count rows; never print row values."""

    if filename.endswith(".gz"):
        stream = gzip.GzipFile(
            fileobj=handle, mode="rb"
        )
    else:
        stream = handle

    delimiter = (
        "," if filename.endswith(
            (".csv", ".csv.gz")
        ) else "\t"
    )

    with io.TextIOWrapper(
        stream,
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    ) as text:

        reader = csv.DictReader(
            text, delimiter=delimiter
        )

        print_schema(
            reader.fieldnames, prefix
        )

        n = 0
        drp = 0

        for row in reader:
            n += 1

            # Aggregate study count only.
            # No accession or sample field is printed.
            if any(
                (row.get(k) or "").strip() == "DRP005276"
                for k in (
                    "study",
                    "study_acc",
                    "study_accession"
                )
            ):
                drp += 1

        print(prefix, "ROW_COUNT", n)
        print(
            prefix,
            "EXPLICIT_DRP005276_STUDY_ROWS",
            drp
        )

def inspect_json(handle, prefix):
    """Print schema/type information, never JSON values."""

    obj = json.load(
        io.TextIOWrapper(
            handle,
            encoding="utf-8-sig",
            errors="replace"
        )
    )

    if isinstance(obj, dict):
        print(
            prefix,
            "JSON_TOP_LEVEL_KEYS",
            len(obj)
        )

        # Only reveal conventional descriptive field names.
        # Map keys resembling IDs are withheld.
        for k in sorted(obj)[:65]:
            value = obj[k]
            key = str(k)

            if (
                not re.fullmatch(
                    r"[A-Za-z_][A-Za-z0-9_.-]{0,75}",
                    key
                )
                or re.search(
                    r"(?:patient|subject|sample|isolate)"
                    r"[_-]?\d+",
                    key, re.I
                )
            ):
                key = "[KEY_WITHHELD]"

            kind = type(value).__name__

            if isinstance(
                value, (list, dict)
            ):
                kind += "_LEN_" + str(len(value))

            print(
                prefix,
                "JSON_KEY_TYPE",
                safe_label(key), kind
            )

    elif isinstance(obj, list):
        print(
            prefix,
            "JSON_LIST_LENGTH",
            len(obj)
        )

    else:
        print(
            prefix,
            "JSON_ROOT_TYPE",
            type(obj).__name__
        )

# ============================================================
# 1. Original identity archive: file identity.
# ============================================================

print(
    "V6_DRP005276_ARCHS4_IDENTITY_ARCHIVE_SCHEMA_V1"
)

assert archive.is_file(), \
    "HISTORICAL_IDENTITY_ARCHIVE_MISSING"

print(
    "ZIP_SHA256",
    digest_file(archive)
)
print(
    "ZIP_COMPRESSED_BYTES",
    archive.stat().st_size
)

# ============================================================
# 2. Inspect ZIP central directory and bounded member content.
#
# ZIP members are read in memory and never extracted.
# The archive is not declared authoritative merely because
# it exists or its CRC checks pass.
# ============================================================

with zipfile.ZipFile(
    str(archive), "r"
) as zf:

    members = [
        m for m in zf.infolist()
        if not m.is_dir()
    ]

    names = [m.filename for m in members]

    print("ZIP_FILES", len(members))

    print(
        "ZIP_DUPLICATE_MEMBER_NAMES",
        len(names) - len(set(names))
    )

    print(
        "ZIP_TOTAL_UNCOMPRESSED_BYTES",
        sum(m.file_size for m in members)
    )

    assert len(members) <= 150, \
        "TOO_MANY_ZIP_MEMBERS"

    assert len(names) == len(set(names)), \
        "DUPLICATE_ZIP_MEMBER_NAMES"

    for index, info in enumerate(
        members, 1
    ):
        name = info.filename
        p = PurePosixPath(name)

        unsafe = (
            name.startswith("/")
            or ".." in p.parts
            or "\\" in name
        )

        print(
            "ZIP_MEMBER", index,
            "NAME", safe_label(name),
            "UNCOMPRESSED_BYTES",
            info.file_size,
            "CRC32",
            "{:08x}".format(info.CRC),
            "UNSAFE_NAME",
            unsafe
        )

        if unsafe:
            print(
                "MEMBER", index,
                "UNSAFE_PATH_SKIP"
            )
            continue

        # Do not decompress unusually large members.
        if info.file_size > MAX_READ:
            print(
                "MEMBER", index,
                "SKIP_FULL_READ_SIZE_LIMIT"
            )
            continue

        # Read member completely only within the size limit.
        # ZipFile validates its CRC during complete reads.
        h = hashlib.sha256()

        with zf.open(info, "r") as f:
            for block in iter(
                lambda: f.read(1048576),
                b""
            ):
                h.update(block)

        print(
            "MEMBER", index,
            "UNCOMPRESSED_SHA256",
            h.hexdigest(),
            "CRC_READ_PASS"
        )

        filename = name.lower()
        prefix = "MEMBER_" + str(index)

        if filename.endswith((
            ".tsv", ".csv",
            ".tsv.gz", ".csv.gz"
        )):
            print(
                "MEMBER", index,
                "TABULAR_SCHEMA_BEGIN"
            )

            with zf.open(
                info, "r"
            ) as f:
                inspect_table(
                    f, filename, prefix
                )

        elif (
            filename.endswith(".json")
            and info.file_size <= MAX_JSON
        ):
            print(
                "MEMBER", index,
                "JSON_SCHEMA_BEGIN"
            )

            with zf.open(
                info, "r"
            ) as f:
                inspect_json(
                    f, prefix
                )

        elif filename.endswith((
            ".txt", ".md", ".py", ".sh"
        )):
            print(
                "MEMBER", index,
                "TEXT_OR_SCRIPT_PRESENT",
                "CONTENT_NOT_EMITTED"
            )

# ============================================================
# 3. Inspect historical v1.0.0/v1.0.1 overlap contracts.
#
# These files are NOT assumed equivalent to the 2026-11-20 ZIP.
# ============================================================

print("\nOLDER_ATLAS_CONTRACT_FILES")

for path in previous:
    relative = str(
        path.relative_to(root)
    )

    if not path.is_file():
        print(
            "EXTERNAL_FILE_MISSING",
            relative
        )
        continue

    size = path.stat().st_size

    print(
        "EXTERNAL_FILE",
        relative,
        "BYTES", size,
        "SHA256", digest_file(path)
    )

    if size > MAX_READ:
        print(
            "EXTERNAL_FILE_SKIPPED_SIZE_LIMIT",
            relative
        )
        continue

    prefix = (
        "EXTERNAL_" +
        path.parent.parent.name
    )

    if path.suffix == ".tsv":
        with path.open("rb") as f:
            inspect_table(
                f, ".tsv", prefix
            )

    elif (
        path.suffix == ".json"
        and size <= MAX_JSON
    ):
        with path.open("rb") as f:
            inspect_json(
                f, prefix
            )

# ============================================================
# 4. No patient-independence certification here.
# ============================================================

print("\nSCIENTIFIC_GATES")

print(
    "IDENTITY_ZIP",
    "CONTENTS_AND_SCHEMAS_INSPECTED_"
    "NOT_ATTESTED_AUTHORITATIVE"
)

print(
    "DRP005276_ATLAS_SOURCE_ID_JOIN",
    "NOT_YET_EXECUTED"
)

print(
    "ARCHS4_PATIENT_INDEPENDENCE",
    "NOT_CERTIFIED"
)

print(
    "PRETREATMENT_STATUS",
    "NOT_CERTIFIED"
)

print("CLINICAL_ADMISSIONS", 0)
print("NO_SAMPLE_OR_PATIENT_ROW_VALUES_PRINTED")
print("NO_ARCHIVE_MEMBERS_EXTRACTED")
print("NO_FROZEN_FILES_MODIFIED")

print(
    "V6_DRP005276_ARCHS4_IDENTITY_"
    "ARCHIVE_SCHEMA_V1_COMPLETE"
)
