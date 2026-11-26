#!/usr/bin/env python3
"""DRP005276 S1 patient-ID header and treatment evidence audit.

Read-only source inspection. Python 3.6 compatible.
No individual patient, sample or experiment values printed.
"""
import csv
import gzip
import hashlib
import io
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"

raw = root / (
    "26_recount3_external_primary_v1/raw_studies/"
    "DRP005276/metadata/sra.sra.DRP005276.MD.gz"
)
source = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)

checks = {
    raw:
        "cba60a28b225d7852e17e180ddaad2680504e7f1cd5cf94de5d73401db12d36c",
    source:
        "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb",
}

s1_hash = (
    "d43789f6dab336a047dc10690600f0"
    "b57c48de61156c9d72b59ae4e9a1508270"
)

s1_url = (
    "https://journals.plos.org/plosone/article/file"
    "?id=10.1371/journal.pone.0220173.s001&type=supplementary"
)

def norm(s):
    return " ".join(str(s or "").strip().casefold().split())

def hashfile(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def isolate_value(attrs):
    values = set()
    for item in (attrs or "").split("|"):
        if ";;" not in item:
            continue
        k, v = item.split(";;", 1)
        if norm(k) == "isolate" and norm(v):
            values.add(norm(v))
    assert len(values) == 1, "ISOLATE_UNRESOLVED"
    return next(iter(values))

def localtag(tag):
    return tag.rsplit("}", 1)[-1]

def coord(ref):
    m = re.fullmatch(r"([A-Z]+)([1-9][0-9]*)", ref or "")
    assert m, "INVALID_CELL_REFERENCE"
    col = 0
    for letter in m.group(1):
        col = col * 26 + ord(letter) - 64
    return int(m.group(2)), col

# Classify headers, not patient-specific values.
# Unknown labels are withheld rather than printed.
def header_category(text):
    s = norm(text)
    if not s or len(s) > 100:
        return None

    if re.search(
        r"(?:patient|subject|case)\s*(?:id|no\.?|number)\b", s
    ):
        return "PATIENT_IDENTIFIER"

    if re.search(r"\b(?:stage|t\s*-?\s*stage|pt)\b", s):
        return "STAGE"

    if re.search(r"\b(?:grade|grading)\b", s):
        return "GRADE"

    if re.search(r"\b(?:sex|gender)\b", s):
        return "SEX"

    if re.search(r"\bage\b", s):
        return "AGE"

    if "smok" in s:
        return "SMOKING"

    if re.search(
        r"\b(?:pre\s*-?treat|pretreat|neoadjuvant|"
        r"prior\s+therapy|prior\s+treatment)\b", s
    ):
        return "PRE_COLLECTION_TREATMENT_CANDIDATE"

    if re.search(
        r"\b(?:bcg|intravesical|chemotherapy|radiotherapy|"
        r"immunotherapy|treatment|therapy|treated|medication)\b",
        s
    ):
        return "TREATMENT_UNTIMED_CANDIDATE"

    if re.search(
        r"\b(?:collection|biopsy|resection|cystectomy|"
        r"turbt|surgery)\b", s
    ):
        return "COLLECTION_PROCEDURE_CANDIDATE"

    if re.search(
        r"\b(?:recurrence|relapse|follow.?up|survival|outcome)\b",
        s
    ):
        return "OUTCOME_OR_FOLLOWUP"

    return None

print("V6_DRP005276_S1_PATIENT_HEADER_TREATMENT_V1")

# ------------------------------------------------------------
# 1. Verify original files without changing them.
# ------------------------------------------------------------

for path, expected in checks.items():
    assert path.is_file()
    assert hashfile(path) == expected, "FROZEN_SOURCE_CHANGED"

print("FROZEN_INPUT_SHA256_PASS", len(checks))

with gzip.open(
    str(raw), "rt", encoding="utf-8",
    errors="replace", newline=""
) as f:
    localrows = list(csv.DictReader(f, delimiter="\t"))

local = {
    r["experiment_acc"]: isolate_value(r["sample_attributes"])
    for r in localrows
}

assert len(localrows) == 103
assert len(local) == 103
assert len(set(local.values())) == 103

with source.open(encoding="utf-8-sig", newline="") as f:
    sr = list(csv.DictReader(f, delimiter="\t"))

selected = {
    r["experiment_acc"]
    for r in sr if r["study"] == "DRP005276"
}

assert len(selected) == 83
assert selected <= set(local)

print("LOCAL_EXPERIMENTS", len(local))
print("SELECTED_EXPERIMENTS", len(selected))

# ------------------------------------------------------------
# 2. Fetch original PLOS S1 into memory only.
# ------------------------------------------------------------

request = Request(
    s1_url,
    headers={"User-Agent": "clinical-provenance-readonly/1.0"}
)

with urlopen(request, timeout=60) as response:
    data = response.read(20000001)

assert 0 < len(data) <= 20000000
assert hashlib.sha256(data).hexdigest() == s1_hash, \
    "PLOS_S1_SOURCE_CHANGED"

print("PLOS_S1_SHA256_PASS", s1_hash)
assert zipfile.is_zipfile(io.BytesIO(data))

# ------------------------------------------------------------
# 3. Parse worksheet cell coordinates and values.
# ------------------------------------------------------------

with zipfile.ZipFile(io.BytesIO(data)) as z:
    names = set(z.namelist())

    sheets = [
        p for p in names
        if re.fullmatch(
            r"xl/worksheets/sheet\d+\.xml", p
        )
    ]

    assert len(sheets) == 1, "S1_SHEET_COUNT_CHANGED"

    shared = []

    if "xl/sharedStrings.xml" in names:
        ss = ET.fromstring(z.read("xl/sharedStrings.xml"))

        for item in ss.iter():
            if localtag(item.tag) == "si":
                shared.append(
                    "".join(
                        t.text or ""
                        for t in item.iter()
                        if localtag(t.tag) == "t"
                    )
                )

    xml = ET.fromstring(z.read(sheets[0]))
    cells = {}

    for el in xml.iter():
        if localtag(el.tag) != "c":
            continue

        position = coord(el.get("r", ""))
        vn = next(
            (n for n in el if localtag(n.tag) == "v"),
            None
        )

        if el.get("t") == "s" and vn is not None \
                and vn.text is not None:
            value = shared[int(vn.text)]

        elif el.get("t") == "inlineStr":
            value = "".join(
                n.text or ""
                for n in el.iter()
                if localtag(n.tag) == "t"
            )

        else:
            value = (vn.text or "") if vn is not None else ""

        if norm(value):
            assert position not in cells
            cells[position] = value

# ------------------------------------------------------------
# 4. Resolve the physical isolate column.
# ------------------------------------------------------------

reverse = {v: k for k, v in local.items()}
hits = defaultdict(list)

for position, value in cells.items():
    if norm(value) in reverse:
        hits[norm(value)].append(position)

assert len(hits) == 103
assert all(len(hits[v]) == 1 for v in reverse), \
    "ISOLATE_CELL_NONUNIQUE"

positions = {
    reverse[v]: coordinates[0]
    for v, coordinates in hits.items()
}

patient_rows = {p[0] for p in positions.values()}
patient_cols = {p[1] for p in positions.values()}

assert len(patient_rows) == 103
assert len(patient_cols) == 1, "PATIENT_LAYOUT_CHANGED"

identifier_col = next(iter(patient_cols))
first_patient_row = min(patient_rows)

assert first_patient_row > 1, "NO_HEADER_ROWS"

print("\nPATIENT_IDENTITY_COLUMN")
print("ISOLATE_POSITION_UNIQUE", 103)
print("ISOLATE_COLUMN_NUMBER", identifier_col)
print("PATIENT_DISTINCT_ROWS", len(patient_rows))
print("S1_FIRST_PATIENT_ROW", first_patient_row)

# ------------------------------------------------------------
# 5. Inspect headers above patient rows.
# ------------------------------------------------------------

header_rows = list(range(1, first_patient_row))
all_columns = sorted({c for r, c in cells})

by_column = {}

for col in all_columns:
    labels = [
        header_category(cells.get((r, col), ""))
        for r in header_rows
    ]
    by_column[col] = {v for v in labels if v}

id_columns = [
    col for col in all_columns
    if "PATIENT_IDENTIFIER" in by_column[col]
]

patient_header_pass = (
    len(id_columns) == 1
    and id_columns[0] == identifier_col
)

print("PATIENT_ID_HEADER_COLUMNS", id_columns)
print(
    "ISOLATE_COLUMN_EQUALS_UNIQUE_PATIENT_ID_HEADER",
    patient_header_pass
)
print(
    "PATIENT_ID_HEADER_IN_ISOLATE_COLUMN",
    "PATIENT_IDENTIFIER" in by_column[identifier_col]
)

# ------------------------------------------------------------
# 6. Inventory potential clinical/treatment fields.
#
# Print only categories, column indices and aggregate coverage.
# Do not interpret treatment after diagnosis as pretreatment.
# ------------------------------------------------------------

print("\nCLINICAL_HEADER_COVERAGE")

categories = (
    "STAGE",
    "GRADE",
    "SEX",
    "AGE",
    "SMOKING",
    "PRE_COLLECTION_TREATMENT_CANDIDATE",
    "TREATMENT_UNTIMED_CANDIDATE",
    "COLLECTION_PROCEDURE_CANDIDATE",
    "OUTCOME_OR_FOLLOWUP",
)

for category in categories:
    columns = sorted(
        col for col in all_columns
        if category in by_column[col]
    )

    print(
        "HEADER_CATEGORY", category,
        "COLUMNS", len(columns)
    )

    for col in columns:
        populated = sum(
            bool(norm(cells.get((positions[acc][0], col), "")))
            for acc in local
        )

        print(
            "HEADER_CATEGORY_COVERAGE",
            category,
            "COLUMN", col,
            "POPULATED", populated
        )

print(
    "WORKBOOK_HEADER_CLASSIFIED_COLUMNS",
    sum(bool(by_column[c]) for c in all_columns)
)

print(
    "WORKBOOK_HEADER_UNCLASSIFIED_COLUMNS",
    sum(not by_column[c] for c in all_columns)
)

# ------------------------------------------------------------
# 7. Final interpretation gates. No clinical admission.
# ------------------------------------------------------------

print("\nSCIENTIFIC_GATES")

print(
    "WITHIN_STUDY_PATIENT_HEADER_LINK",
    "SUPPORTED" if patient_header_pass else "UNRESOLVED"
)

print("PRETREATMENT_STATUS NOT_CERTIFIED")
print("CROSS_STUDY_INDEPENDENCE NOT_CERTIFIED")
print("CLINICAL_ADMISSIONS 0")

print("NO_INDIVIDUAL_PATIENT_OR_SAMPLE_VALUES_PRINTED")
print("NO_FROZEN_FILES_MODIFIED")

print(
    "V6_DRP005276_S1_PATIENT_HEADER_TREATMENT_V1_COMPLETE"
)
