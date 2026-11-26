#!/usr/bin/env python3
"""Inspect DRP005276 isolate positions and clinical fields in PLOS S1.

Read-only; Python 3.6 compatible.
No individual patient, experiment or sample identifiers printed.
"""

import csv
import gzip
import hashlib
import io
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

raw_path = rec / (
    "raw_studies/DRP005276/metadata/"
    "sra.sra.DRP005276.MD.gz"
)
source_path = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)

SHA = {
    raw_path:
        "cba60a28b225d7852e17e180ddaad2680504e7f1cd5cf94de5d73401db12d36c",
    source_path:
        "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb"
}

S1_SHA = (
    "d43789f6dab336a047dc10690600f0"
    "b57c48de61156c9d72b59ae4e9a1508270"
)
S1_URL = (
    "https://journals.plos.org/plosone/article/file"
    "?id=10.1371/journal.pone.0220173.s001&type=supplementary"
)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def sha_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def norm(value):
    return " ".join(str(value or "").casefold().strip().split())

def lname(tag):
    return tag.rsplit("}", 1)[-1]

def attr_isolate(text):
    values = set()

    for part in (text or "").split("|"):
        if ";;" not in part:
            continue
        key, value = part.split(";;", 1)
        if norm(key) == "isolate" and norm(value):
            values.add(norm(value))

    assert len(values) == 1, "ISOLATE_MISSING_OR_MULTIVALUED"
    return next(iter(values))

def colnum(ref):
    m = re.fullmatch(r"([A-Z]+)([1-9][0-9]*)", ref or "")
    assert m, "INVALID_SHEET_COORDINATE"

    col = 0
    for char in m.group(1):
        col = 26 * col + ord(char) - ord("A") + 1

    return int(m.group(2)), col

# ------------------------------------------------------------
# Controlled metadata-label classification.
# Do not print raw workbook labels or clinical values.
# ------------------------------------------------------------

def classify_label(value):
    s = norm(value)

    if not s or len(s) > 100 or "@" in s:
        return None

    if re.search(
        r"(?:patient|subject|case)\s*(?:id|no|number)\b", s
    ):
        return "PATIENT_IDENTIFIER_LABEL"

    if re.search(
        r"\b(?:tumou?r|patholog(?:ical|ic)|clinical)?\s*"
        r"(?:t\s*-?\s*stage|stage)\b", s
    ) or s in ("pt", "t stage", "stage", "pathological t"):
        return "STAGE_LABEL"

    if re.search(
        r"\b(?:nmibc|mibc|muscle\s*-?\s*invas(?:ive|ion))\b",
        s
    ):
        return "MUSCLE_INVASION_LABEL"

    if "grade" in s:
        return "GRADE_LABEL"

    if re.search(r"\b(?:sex|gender)\b", s):
        return "SEX_LABEL"

    if re.search(r"\bage\b", s):
        return "AGE_LABEL"

    if "smok" in s:
        return "SMOKING_LABEL"

    if re.search(
        r"\b(?:treat(?:ment|ed)?|therapy|bcg|neoadjuvant|"
        r"chemo(?:therapy)?|radiotherapy|cystectomy|"
        r"tur(?:bt)?|surgery)\b", s
    ):
        return "TREATMENT_OR_PROCEDURE_LABEL"

    if "recurr" in s:
        return "RECURRENCE_LABEL"

    if "primary" in s or "metasta" in s:
        return "PRIMARY_OR_METASTASIS_LABEL"

    return None

# Conservative classification using explicit stage notation.
# Unknown/ambiguous codes remain UNKNOWN, never inferred.
def classify_stage(value):
    s = norm(value)
    compact = re.sub(r"[\s_\-/]+", "", s)

    if compact in (
        "nmibc",
        "nonmuscleinvasive",
        "nonmuscleinvasivebladdercancer"
    ):
        return "NMIBC"

    if compact in (
        "mibc",
        "muscleinvasive",
        "muscleinvasivebladdercancer"
    ):
        return "MIBC"

    m = re.fullmatch(
        r"(?:p|c|yp)?t(is|a|[0-4])(?:[a-c])?",
        compact
    )

    if m:
        return (
            "NMIBC" if m.group(1) in ("is", "a", "0", "1")
            else "MIBC"
        )

    return "UNKNOWN"

print("V6_DRP005276_S1_PATIENT_AXIS_STAGE_DIAGNOSTIC_V1")

# ------------------------------------------------------------
# 1. Verify immutable local sources.
# ------------------------------------------------------------

for path, expected in SHA.items():
    assert path.is_file()
    assert sha_file(path) == expected, "LOCAL_SOURCE_SHA_MISMATCH"

print("FROZEN_INPUT_SHA256_PASS", len(SHA))

with gzip.open(
    str(raw_path), "rt",
    encoding="utf-8",
    errors="replace",
    newline=""
) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

local = {r["experiment_acc"]: r for r in rows}
assert len(rows) == len(local) == 103

with source_path.open(
    encoding="utf-8-sig", newline=""
) as f:
    src = list(csv.DictReader(f, delimiter="\t"))

selected = {
    r["experiment_acc"]
    for r in src if r["study"] == "DRP005276"
}

assert len(selected) == 83
assert selected <= set(local)

# Exactly one original source isolate per experiment.
isolate = {
    acc: attr_isolate(r["sample_attributes"])
    for acc, r in local.items()
}

assert len(set(isolate.values())) == 103
reverse = {token: acc for acc, token in isolate.items()}

print("ORIGINAL_EXPERIMENTS", len(local))
print("SELECTED_EXPERIMENTS", len(selected))
print("QC_HELD_EXPERIMENTS", len(local) - len(selected))
print("DISTINCT_ORIGINAL_ISOLATES", len(reverse))

# ------------------------------------------------------------
# 2. Fetch the exact published S1 workbook in memory.
#    Verify its original SHA256 before reading any cell.
# ------------------------------------------------------------

req = Request(
    S1_URL,
    headers={
        "User-Agent": "recount3-v6-S1-layout-review/1.0"
    }
)

with urlopen(req, timeout=60) as response:
    payload = response.read(20000001)

assert 0 < len(payload) <= 20000000
assert sha(payload) == S1_SHA, "PLOS_S1_CONTENT_HASH_CHANGED"
assert zipfile.is_zipfile(io.BytesIO(payload))

print("PLOS_S1_SHA256_PASS", sha(payload))

# ------------------------------------------------------------
# 3. Parse XLSX directly from its OOXML structure.
#    No workbook contents or identifiers are written to disk.
# ------------------------------------------------------------

with zipfile.ZipFile(io.BytesIO(payload)) as z:
    sheets = sorted(
        name for name in z.namelist()
        if re.fullmatch(
            r"xl/worksheets/sheet\d+\.xml", name
        )
    )

    assert len(sheets) == 1, "S1_SHEET_COUNT_CHANGED"

    strings = []

    if "xl/sharedStrings.xml" in z.namelist():
        root_xml = ET.fromstring(z.read("xl/sharedStrings.xml"))

        for item in root_xml.iter():
            if lname(item.tag) == "si":
                strings.append(
                    "".join(
                        t.text or ""
                        for t in item.iter()
                        if lname(t.tag) == "t"
                    )
                )

    sheet = ET.fromstring(z.read(sheets[0]))
    grid = {}

    for cell in sheet.iter():
        if lname(cell.tag) != "c":
            continue

        coord = colnum(cell.get("r", ""))
        vn = next(
            (v for v in cell if lname(v.tag) == "v"),
            None
        )
        typ = cell.get("t", "")

        if typ == "s" and vn is not None and vn.text is not None:
            value = strings[int(vn.text)]
        elif typ == "inlineStr":
            value = "".join(
                t.text or ""
                for t in cell.iter()
                if lname(t.tag) == "t"
            )
        elif vn is not None:
            value = vn.text or ""
        else:
            value = ""

        if norm(value):
            assert coord not in grid, "DUPLICATE_CELL_COORDINATE"
            grid[coord] = value

print("\nS1_WORKBOOK_STRUCTURE")
print("S1_POPULATED_CELLS", len(grid))
print(
    "S1_DISTINCT_POPULATED_ROWS",
    len({r for r, c in grid})
)
print(
    "S1_DISTINCT_POPULATED_COLUMNS",
    len({c for r, c in grid})
)

# ------------------------------------------------------------
# 4. Locate each isolate in the published worksheet.
#
# A 103/103 value match is insufficient by itself.
# We require one exact physical cell position per isolate.
# ------------------------------------------------------------

hits = defaultdict(list)

for (r, c), value in grid.items():
    token = norm(value)

    if token in reverse:
        hits[token].append((r, c))

print("\nISOLATE_CELL_MAPPING")
print("S1_ISOLATES_MATCHED", len(hits))

occurrences = Counter(
    len(hits.get(token, []))
    for token in reverse
)

print(
    "S1_ISOLATE_OCCURRENCE_HISTOGRAM",
    dict(sorted(occurrences.items()))
)

unique = (
    len(hits) == 103
    and all(len(hits[token]) == 1 for token in reverse)
)

print("S1_EXACTLY_ONE_CELL_PER_ISOLATE", unique)

if not unique:
    print("PATIENT_AXIS_GATE STOP_NONUNIQUE")
    print("CLINICAL_ADMISSIONS 0")
    print("NO_INDIVIDUAL_IDENTIFIERS_PRINTED")
    print(
        "V6_DRP005276_S1_PATIENT_AXIS_STAGE_"
        "DIAGNOSTIC_V1_COMPLETE"
    )
    sys.exit(0)

pos = {
    reverse[token]: positions[0]
    for token, positions in hits.items()
}

match_rows = {r for r, c in pos.values()}
match_cols = {c for r, c in pos.values()}

print("S1_ISOLATE_DISTINCT_ROWS", len(match_rows))
print("S1_ISOLATE_DISTINCT_COLUMNS", len(match_cols))

# ------------------------------------------------------------
# 5. Determine whether isolates occupy columns or rows.
#    Inspect labels on the opposite metadata axis.
# ------------------------------------------------------------

if len(match_rows) == 1 and len(match_cols) == 103:
    axis = "PATIENT_COLUMNS"
    pivot = min(match_cols)

    items = sorted({r for r, c in grid})

    def candidate_label(index):
        return [
            grid.get((index, c), "")
            for c in range(1, min(pivot, 11))
        ]

    def value_at(index, exp):
        return grid.get((index, pos[exp][1]), "")

elif len(match_cols) == 1 and len(match_rows) == 103:
    axis = "PATIENT_ROWS"
    pivot = min(match_rows)

    items = sorted({c for r, c in grid})

    def candidate_label(index):
        return [
            grid.get((r, index), "")
            for r in range(1, min(pivot, 11))
        ]

    def value_at(index, exp):
        return grid.get((pos[exp][0], index), "")

else:
    axis = "MIXED_UNRESOLVED"

print("S1_PATIENT_AXIS_STRUCTURE", axis)

if axis == "MIXED_UNRESOLVED":
    print("PATIENT_AXIS_GATE STOP_MIXED")
    print("CLINICAL_ADMISSIONS 0")
    print("NO_INDIVIDUAL_IDENTIFIERS_PRINTED")
    print(
        "V6_DRP005276_S1_PATIENT_AXIS_STAGE_"
        "DIAGNOSTIC_V1_COMPLETE"
    )
    sys.exit(0)

# ------------------------------------------------------------
# 6. Identify labelled clinical fields.
#    Print categories and coverage, never original cell values.
# ------------------------------------------------------------

found_labels = Counter()
stage_candidates = []

print("\nPUBLISHED_METADATA_FIELDS")

for index in items:
    classes = {
        classify_label(value)
        for value in candidate_label(index)
    } - {None}

    if not classes:
        continue

    for cls in sorted(classes):
        found_labels[cls] += 1

        populated = sum(
            bool(norm(value_at(index, exp)))
            for exp in local
        )

        print(
            "ANNOTATION_FIELD",
            cls,
            "INDEX", index,
            "POPULATED_OF_103", populated
        )

        if cls in (
            "STAGE_LABEL",
            "MUSCLE_INVASION_LABEL"
        ):
            stage_candidates.append((index, cls))

print(
    "S1_METADATA_LABEL_CLASS_COUNTS",
    dict(sorted(found_labels.items()))
)
print(
    "S1_CANDIDATE_STAGE_FIELDS",
    len(stage_candidates)
)

# ------------------------------------------------------------
# 7. Attempt conservative tumour-stage classification.
#    Do not silently use unrecognized or ambiguous codes.
# ------------------------------------------------------------

print("\nSTAGE_MAPPING_CANDIDATES")

for index, cls in stage_candidates:
    labels = {
        exp: classify_stage(value_at(index, exp))
        for exp in local
    }

    full = Counter(labels.values())
    retained = Counter(labels[x] for x in selected)
    held = Counter(
        labels[x] for x in set(local) - selected
    )

    print(
        "STAGE_FIELD_INDEX", index,
        "LABEL_CLASS", cls
    )
    print(
        "STAGE_FIELD_ALL_103",
        dict(sorted(full.items()))
    )
    print(
        "STAGE_FIELD_SELECTED_83",
        dict(sorted(retained.items()))
    )
    print(
        "STAGE_FIELD_QC_HELD_20",
        dict(sorted(held.items()))
    )

    matches_paper = (
        full["NMIBC"] == 59
        and full["MIBC"] == 44
        and full["UNKNOWN"] == 0
    )

    print(
        "STAGE_FULL_PUBLISHED_59_44_CHECK",
        "PASS" if matches_paper else "NOT_CONFIRMED"
    )

# ------------------------------------------------------------
# 8. No clinical upgrade without inspecting semantics.
# ------------------------------------------------------------

print("\nSCIENTIFIC_GATES")
print(
    "PATIENT_EXPERIMENT_MAPPING",
    "ISOLATE_TO_UNIQUE_S1_POSITION_PASS"
)
print(
    "CLINICAL_PATIENT_INDEPENDENCE",
    "REQUIRES_HEADER_SEMANTIC_REVIEW"
)
print(
    "STAGE_CLASSIFICATION",
    "REVIEW_CANDIDATE_FIELDS"
)
print("PRETREATMENT_HISTORY UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("NO_INDIVIDUAL_IDENTIFIERS_PRINTED")
print("NO_SOURCE_DATA_WRITTEN")
print(
    "V6_DRP005276_S1_PATIENT_AXIS_STAGE_"
    "DIAGNOSTIC_V1_COMPLETE"
)
