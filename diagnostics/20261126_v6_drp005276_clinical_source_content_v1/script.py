#!/usr/bin/env python3
"""Inspect clinical-source structure without publishing patient values."""

import csv
import gzip
import hashlib
import io
import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
base = root / "26_recount3_external_primary_v1"

archives = {
    "V2": (
        "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V2_20261121.zip",
        [
            "results/04_PUBLICATION_REVIEWED_STUDIES.tsv",
            "results/02_ALL_176_STUDY_SOURCE_REVIEW_UNSIGNED.tsv",
        ],
    ),
    "V3": (
        "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip",
        [
            "results/04_PUBLICATION_REVIEWED_STUDIES.tsv",
            "results/02_ALL_176_STUDY_SOURCE_REVIEW_UNSIGNED.tsv",
            "inputs/01_public_source_study_evidence_reviewed.json",
        ],
    ),
    "V5": (
        "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip",
        [
            "results/07_NINE_STUDY_SOURCE_CROSSWALK_SUMMARY.tsv",
        ],
    ),
    "V1": (
        "RECOUNT3_FULL_CLINICAL_METADATA_ADJUDICATION_V1_20261120.zip",
        [
            "results_unsigned/03_STUDIES_176_prioritized_evidence_review_UNSIGNED.tsv",
        ],
    ),
}

native = base / "raw_studies/DRP005276/metadata"

native_names = [
    "sra.recount_pred.DRP005276.MD.gz",
    "sra.recount_project.DRP005276.MD.gz",
    "sra.recount_qc.DRP005276.MD.gz",
    "sra.recount_seq_qc.DRP005276.MD.gz",
    "sra.sra.DRP005276.MD.gz",
]

study_pattern = re.compile(r"(?<!\w)DRP005276(?!\w)", re.I)

clinical_pattern = re.compile(
    r"pretreat|pre[- ]?operat|treat(?:ment|ed)?|"
    r"neoadjuv|chemotherap|radiotherap|"
    r"intravesic|immunotherap|"
    r"\bBCG\b|prior.?therap|previous.?therap",
    re.I,
)

timing_pattern = re.compile(
    r"collect|sampl(?:ing|ed)|surg|resect|"
    r"diagnos|timepoint|time_point",
    re.I,
)

def digest(data):
    return hashlib.sha256(
        str(data).encode("utf-8")
    ).hexdigest()[:16]

def print_study_table(label, member, raw):
    text = raw.decode("utf-8-sig")
    reader = csv.DictReader(
        io.StringIO(text),
        delimiter="\t"
    )
    fields = reader.fieldnames or []
    matches = []

    for row in reader:
        if any(
            study_pattern.search(str(v or ""))
            for v in row.values()
        ):
            matches.append(row)

    print(
        "STUDY_TABLE", label, Path(member).name,
        "FIELD_COUNT", len(fields),
        "DRP005276_ROWS", len(matches)
    )

    for field in fields:
        if clinical_pattern.search(field) or timing_pattern.search(field):
            populated = sum(
                bool((r.get(field) or "").strip())
                for r in matches
            )
            print(
                "RELEVANT_FIELD",
                label,
                field[:100],
                "POPULATED_STUDY_ROWS", populated
            )

    clinical_hits = []
    timing_hits = []

    for row in matches:
        for field, value in row.items():
            value = str(value or "")
            if clinical_pattern.search(value):
                clinical_hits.append(field)
            if timing_pattern.search(value):
                timing_hits.append(field)

    print(
        "STUDY_ROW_CLINICAL_TEXT_FIELDS",
        label,
        sorted(set(clinical_hits))[:25]
    )
    print(
        "STUDY_ROW_TIMING_TEXT_FIELDS",
        label,
        sorted(set(timing_hits))[:25]
    )
    print(
        "STUDY_ROW_CLINICAL_TEXT_FIELD_HITS",
        label, len(clinical_hits)
    )

print("DRP005276_CLINICAL_SOURCE_CONTENT_V1")

# Previous derived clinical reviews: source presence and fields.
for label, (filename, wanted_members) in archives.items():
    path = base / filename
    assert path.is_file()

    with zipfile.ZipFile(str(path)) as archive:
        names = archive.namelist()

        for suffix in wanted_members:
            matches = [
                n for n in names
                if n.endswith("/" + suffix)
            ]
            assert len(matches) == 1, (label, suffix)

            member = matches[0]
            raw = archive.read(member)

            if suffix.endswith(".tsv"):
                print_study_table(label, member, raw)
            else:
                obj = json.loads(raw.decode("utf-8-sig"))
                text = json.dumps(
                    obj, ensure_ascii=False
                )

                print(
                    "STUDY_EVIDENCE_JSON",
                    label,
                    "HAS_DRP005276",
                    bool(study_pattern.search(text)),
                    "HAS_TREATMENT_TERMS",
                    bool(clinical_pattern.search(text)),
                    "HAS_TIMING_TERMS",
                    bool(timing_pattern.search(text)),
                    "SHA256_PREFIX", digest(text)
                )

# Native recount3 metadata: actual column names plus
# counts of treatment/timing terminology, never values.
for filename in native_names:
    path = native / filename
    assert path.is_file(), filename

    with gzip.open(
        str(path), "rt",
        encoding="utf-8-sig",
        newline=""
    ) as handle:
        reader = csv.DictReader(
            handle, delimiter="\t"
        )
        fields = reader.fieldnames or []

        nrows = 0
        clinical_rows = 0
        timing_rows = 0
        clinical_fields = Counter()
        timing_fields = Counter()

        for row in reader:
            nrows += 1

            c = set()
            t = set()

            for field, value in row.items():
                value = str(value or "")

                if clinical_pattern.search(value):
                    c.add(field)

                if timing_pattern.search(value):
                    t.add(field)

            if c:
                clinical_rows += 1
                clinical_fields.update(c)

            if t:
                timing_rows += 1
                timing_fields.update(t)

    print(
        "NATIVE_METADATA", filename,
        "FIELDS", len(fields),
        "ROWS", nrows,
        "CLINICAL_TERM_ROWS", clinical_rows,
        "TIMING_TERM_ROWS", timing_rows
    )

    print(
        "NATIVE_CLINICAL_HEADER_FIELDS",
        filename,
        [
            x[:100] for x in fields
            if clinical_pattern.search(x)
        ][:30]
    )

    print(
        "NATIVE_CLINICAL_VALUE_FIELDS",
        filename,
        sorted(clinical_fields)[:30]
    )

    print(
        "NATIVE_TIMING_VALUE_FIELDS",
        filename,
        sorted(timing_fields)[:30]
    )

print("SOURCE_CONTENT_TRIAGE_ONLY")
print("DERIVED_ARCHIVES_NOT_INDEPENDENT_CLINICAL_ATTESTATION")
print("NO_INDIVIDUAL_SAMPLE_VALUES_PRINTED")
print("NO_PATIENT_TREATMENT_CERTIFICATION")
print("PATIENT_INDEPENDENCE_UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("NO_SOURCE_FILES_MODIFIED")
print("DRP005276_CLINICAL_SOURCE_CONTENT_V1_COMPLETE")
