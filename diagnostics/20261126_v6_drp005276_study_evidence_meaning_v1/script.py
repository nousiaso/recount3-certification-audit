#!/usr/bin/env python3
"""Review the DRP005276 study-level clinical text, not patient rows."""

import csv
import hashlib
import io
import json
import re
import sys
import zipfile
from pathlib import Path

base = Path(sys.argv[1]) / "26_recount3_external_primary_v1"

archives = (
    ("V2", "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V2_20261121.zip"),
    ("V3", "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"),
)

fields_of_interest = (
    "on_or_after_treatment_libraries",
    "collection_procedures",
    "qa_missing_pretreatment",
    "public_source_treatment_evidence",
    "full_source_curation_comment",
    "review_action_from_publication",
    "review_caveat",
    "required_evidence",
)

study = re.compile(r"(?<!\w)DRP005276(?!\w)", re.I)
accession = re.compile(
    r"\b(?:SRX|DRX|ERX|SRR|DRR|ERR|"
    r"GSM|SAMN|SAMD|SAMEA)\d+\b", re.I
)
email = re.compile(
    r"\b[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}\b"
)
patient_label = re.compile(
    r"\b(?:patient|subject|participant)[-_ ]*"
    r"(?:id|number|no\.?)?[: #_-]*[A-Za-z]*\d{2,}\b",
    re.I
)

def clean(value):
    s = re.sub(r"\s+", " ", str(value or "")).strip()
    s = accession.sub("[ACCESSION]", s)
    s = email.sub("[EMAIL]", s)
    s = patient_label.sub("[PATIENT_ID]", s)
    return s[:500] + (" [TRUNCATED]" if len(s) > 500 else "")

print("DRP005276_STUDY_EVIDENCE_MEANING_V1")

for label, filename in archives:
    path = base / filename
    assert path.is_file()

    with zipfile.ZipFile(str(path)) as z:
        members = [
            n for n in z.namelist()
            if n.endswith(
                "/results/04_PUBLICATION_REVIEWED_STUDIES.tsv"
            )
        ]
        assert len(members) == 1

        data = z.read(members[0])
        reader = csv.DictReader(
            io.StringIO(data.decode("utf-8-sig")),
            delimiter="\t"
        )

        columns = reader.fieldnames or []
        matched = [
            row for row in reader
            if any(study.search(str(v or ""))
                   for v in row.values())
        ]
        assert len(matched) == 1

        row = matched[0]

        print("\nSOURCE", label)
        print("STUDY_ROWS", len(matched))
        print("ALL_FIELD_COUNT", len(columns))

        source_fields = [
            f for f in columns
            if re.search(
                r"source|doi|pmid|publication|reference|"
                r"caveat|evidence|treatment|collection",
                f, re.I
            )
        ]

        print(
            "SOURCE_AND_CLINICAL_FIELD_NAMES",
            ";".join(source_fields)
        )

        for field in fields_of_interest:
            if field not in row:
                print("FIELD_ABSENT", field)
                continue

            val = str(row[field] or "")
            print(
                "FIELD", field,
                "VALUE", clean(val),
                "SHA256",
                hashlib.sha256(
                    val.encode("utf-8")
                ).hexdigest()[:16]
            )

        if label == "V3":
            json_members = [
                n for n in z.namelist()
                if n.endswith(
                    "/inputs/01_public_source_study_evidence_reviewed.json"
                )
            ]
            assert len(json_members) == 1
            obj = json.loads(z.read(json_members[0]))

            matches = []

            def visit(node, path="$"):
                if isinstance(node, dict):
                    for k, v in node.items():
                        visit(v, path + "." + str(k))
                elif isinstance(node, list):
                    for i, v in enumerate(node):
                        visit(v, path + "[%d]" % i)
                elif isinstance(node, str):
                    if study.search(node):
                        matches.append(path)

            visit(obj)

            print("V3_JSON_STUDY_MATCH_PATHS", len(matches))
            for path_text in matches[:12]:
                print("JSON_STUDY_PATH", clean(path_text))
            print(
                "JSON_STUDY_PATHS_NOT_PRINTED",
                max(0, len(matches) - 12)
            )

print("STUDY_LEVEL_REVIEW_TEXT_ONLY")
print("NO_INDIVIDUAL_PATIENT_ROWS_PRINTED")
print("DERIVED_REVIEW_NOT_INDEPENDENT_SOURCE_PROOF")
print("PATIENT_INDEPENDENCE_UNSIGNED")
print("PRETREATMENT_UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("NO_SOURCE_FILES_MODIFIED")
print("DRP005276_STUDY_EVIDENCE_MEANING_V1_COMPLETE")
