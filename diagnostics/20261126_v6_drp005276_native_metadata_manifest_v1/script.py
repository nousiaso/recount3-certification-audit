#!/usr/bin/env python3
"""Inspect clinical ZIP members and native metadata headers only."""

import gzip
import os
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

root = Path(sys.argv[1])
base = root / "26_recount3_external_primary_v1"

zip_names = {
    "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V2_20261121.zip",
    "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip",
    "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip",
    "RECOUNT3_FULL_CLINICAL_METADATA_ADJUDICATION_V1_20261120.zip",
}

native_names = {
    "sra.recount_pred.DRP005276.MD.gz",
    "sra.recount_project.DRP005276.MD.gz",
    "sra.recount_qc.DRP005276.MD.gz",
    "sra.recount_seq_qc.DRP005276.MD.gz",
    "sra.sra.DRP005276.MD.gz",
}

wanted = zip_names | native_names
found = defaultdict(list)

def redact(value):
    return re.sub(
        r"(?:SRX|ERX|DRX|SRR|ERR|DRR|GSM|SAMN|SAMD|SAMEA)\d+",
        "[ACCESSION]",
        str(value),
        flags=re.I
    )[:180]

print("DRP005276_NATIVE_METADATA_MANIFEST_V1")

for parent, dirs, names in os.walk(str(base)):
    depth = len(Path(parent).relative_to(base).parts)
    dirs[:] = sorted(
        d for d in dirs
        if d not in (".git", "__pycache__")
        and depth < 8
    )

    for name in names:
        if name in wanted:
            found[name].append(Path(parent) / name)

print("EXPECTED_FILENAMES", len(wanted))

for name in sorted(wanted):
    paths = sorted(found[name])

    print(
        "\nTARGET", name,
        "MATCHING_PATHS", len(paths)
    )

    for path in paths[:3]:
        print(
            "LOCATION",
            redact(str(path.relative_to(root))),
            "COMPRESSED_BYTES",
            path.stat().st_size
        )

        if name.endswith(".zip"):
            try:
                with zipfile.ZipFile(str(path)) as z:
                    members = z.infolist()
                    print(
                        "ZIP_MEMBER_COUNT", len(members),
                        "ZIP_UNCOMPRESSED_BYTES",
                        sum(m.file_size for m in members)
                    )

                    relevant = [
                        m for m in members
                        if re.search(
                            r"readme|clinical|patient|specimen|"
                            r"treat|therapy|source|metadata|"
                            r"evidence|sdrf|idf",
                            m.filename, re.I
                        )
                    ]

                    print(
                        "ZIP_RELEVANT_MEMBER_NAMES",
                        len(relevant)
                    )

                    for m in relevant[:18]:
                        print(
                            "ZIP_MEMBER",
                            redact(m.filename),
                            "BYTES", m.file_size
                        )

                    print(
                        "ZIP_MEMBERS_NOT_PRINTED",
                        max(0, len(relevant) - 18)
                    )

            except (OSError, zipfile.BadZipFile) as exc:
                print(
                    "ZIP_INSPECTION_ERROR",
                    type(exc).__name__
                )

        elif name.endswith(".MD.gz"):
            try:
                with gzip.open(str(path), "rt",
                               encoding="utf-8",
                               errors="replace") as f:
                    first = f.readline(32768)

                print(
                    "GZIP_FIRST_LINE_BYTES",
                    len(first.encode("utf-8")),
                    "LINE_COMPLETE",
                    first.endswith("\n")
                )

                if first.endswith("\n"):
                    delimiter = (
                        "\t" if "\t" in first
                        else "," if "," in first
                        else None
                    )

                    print(
                        "HEADER_DELIMITER",
                        "TAB" if delimiter == "\t"
                        else "COMMA" if delimiter == ","
                        else "UNKNOWN"
                    )

                    if delimiter:
                        columns = first.rstrip(
                            "\r\n"
                        ).split(delimiter)

                        hints = [
                            word for word in (
                                "patient", "clinical",
                                "treatment", "therapy",
                                "chemotherapy", "bcg",
                                "pretreatment", "collection",
                                "specimen", "biosample",
                                "sample", "diagnosis",
                                "timepoint", "prior"
                            )
                            if any(
                                word in col.lower()
                                for col in columns
                            )
                        ]

                        print(
                            "POSSIBLE_HEADER_FIELD_COUNT",
                            len(columns)
                        )
                        print(
                            "HEADER_KEYWORD_HINTS",
                            ";".join(hints) or "[NONE]"
                        )

            except (OSError, EOFError, UnicodeError) as exc:
                print(
                    "GZIP_HEADER_ERROR",
                    type(exc).__name__
                )

    if len(paths) > 3:
        print("ADDITIONAL_LOCATIONS", len(paths) - 3)

print("\nSCIENTIFIC_GATES")
print("ZIP_MANIFESTS_ONLY_NO_ARCHIVE_EXTRACTION")
print("NATIVE_METADATA_FIRST_LINES_ONLY")
print("NO_PATIENT_OR_SAMPLE_ROW_VALUES_PRINTED")
print("NO_HISTORICAL_JOIN_EXECUTED")
print("PATIENT_INDEPENDENCE_UNSIGNED")
print("PRETREATMENT_UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("NO_SOURCE_FILES_MODIFIED")
print("DRP005276_NATIVE_METADATA_MANIFEST_V1_COMPLETE")
