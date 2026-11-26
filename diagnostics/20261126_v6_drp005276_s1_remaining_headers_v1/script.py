#!/usr/bin/env python3
"""DRP005276: S1 header-semantics review.

Read-only; Python 3.6 compatible.
Print workbook headers only, never patient values.
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

raw = root / (
    "26_recount3_external_primary_v1/raw_studies/"
    "DRP005276/metadata/sra.sra.DRP005276.MD.gz"
)

source = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)

expected = {
    raw:
        "cba60a28b225d7852e17e180ddaad2680504e7f1cd5cf94de5d73401db12d36c",
    source:
        "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb"
}

s1_sha = (
    "d43789f6dab336a047dc10690600f0"
    "b57c48de61156c9d72b59ae4e9a1508270"
)

s1_url = (
    "https://journals.plos.org/plosone/article/file"
    "?id=10.1371/journal.pone.0220173.s001&type=supplementary"
)

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1048576), b""):
            h.update(part)
    return h.hexdigest()

def norm(s):
    return " ".join(str(s or "").strip().casefold().split())

def tag(node):
    return node.tag.rsplit("}", 1)[-1]

def column(ref):
    m = re.fullmatch(
        r"([A-Z]+)([1-9][0-9]*)", ref or ""
    )
    assert m, "BAD_XLSX_CELL_COORDINATE"

    num = 0
    for c in m.group(1):
        num = 26 * num + ord(c) - 64

    return int(m.group(2)), num

def isolate(attrs):
    vals = set()

    for part in (attrs or "").split("|"):
        if ";;" in part:
            k, v = part.split(";;", 1)
            if norm(k) == "isolate" and norm(v):
                vals.add(norm(v))

    assert len(vals) == 1, "ISOLATE_AMBIGUOUS"
    return next(iter(vals))

# Only header rows are eligible for printing.
# Unexpected identifier-like strings are withheld.
def safe_header(value, isolate_set):
    s = " ".join(str(value or "").strip().split())

    if not s:
        return "[EMPTY]"

    if len(s) > 100 or norm(s) in isolate_set:
        return "[WITHHELD]"

    if re.search(
        r"\b(?:DRX|DRR|DRS|SRX|SRR|SRS|SAMN|"
        r"SAMD|GSM)\d+\b", s, re.I
    ):
        return "[WITHHELD]"

    if re.search(
        r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b",
        s, re.I
    ):
        return "[WITHHELD]"

    if re.search(
        r"\b[A-Z]{2}-[A-Z]{3}-\d{3}\b", s, re.I
    ):
        return "[WITHHELD]"

    if re.search(
        r"[^a-zA-Z0-9 ,./:;()_+%&\-]", s
    ):
        return "[WITHHELD_NONASCII_OR_UNEXPECTED]"

    return s

# Classification is descriptive; ambiguous headers
# remain unresolved rather than being forced into a type.
def meaning(labels):
    s = " ".join(labels).casefold()

    if re.search(
        r"\b(?:pre.?treat|neoadjuvant|prior.?therap|"
        r"bcg|intravesical|chemotherap|radiotherap|"
        r"immunotherap|treatment|therapy|treated|"
        r"surgery|turbt|resection|biopsy|"
        r"collection.?date|recurre)\w*\b", s
    ):
        return "POTENTIAL_THERAPY_TIMING_OR_COLLECTION"

    if re.search(
        r"\b(?:mutation|mutant|somatic|germline|"
        r"polymorphism|snp|genotype|allele|variant|"
        r"codon|wild.?type|rs1042522)\b", s
    ):
        return "GENOTYPE_OR_MUTATION"

    if re.search(
        r"\b(?:fgfr[1234]|tp53|hras|kras|nras|"
        r"pik3ca|braf|akt1|ctnnb1|rb1|atm|"
        r"egfr|erbb2)\b", s
    ):
        return "GENE_RELATED"

    if re.search(
        r"\b(?:patient|subject|case)\s*"
        r"(?:id|no\.?|number)\b", s
    ):
        return "PATIENT_ID"

    if re.search(
        r"\b(?:stage|grade|age|sex|gender|smoking)\b", s
    ):
        return "CLINICAL_NON_TREATMENT"

    return "UNRESOLVED_HEADER_SEMANTICS"

print("V6_DRP005276_S1_REMAINING_HEADERS_V1")

# ------------------------------------------------------------
# 2. Verify frozen local sources.
# ------------------------------------------------------------

for path, checksum in expected.items():
    assert path.is_file()
    assert digest(path) == checksum, "FROZEN_SOURCE_CHANGED"

print("FROZEN_LOCAL_SHA256_PASS", len(expected))

with gzip.open(
    str(raw), "rt",
    encoding="utf-8",
    errors="replace",
    newline=""
) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 103

by_exp = {
    r["experiment_acc"]: isolate(r["sample_attributes"])
    for r in rows
}

assert len(by_exp) == 103
assert len(set(by_exp.values())) == 103

with source.open(
    encoding="utf-8-sig", newline=""
) as f:
    all_source = list(csv.DictReader(f, delimiter="\t"))

selected = {
    r["experiment_acc"]
    for r in all_source
    if r["study"] == "DRP005276"
}

assert len(selected) == 83
assert selected <= set(by_exp)

isolates = set(by_exp.values())

print(
    "SOURCE_EXPERIMENTS", len(by_exp),
    "SELECTED", len(selected),
    "QC_HELD", len(by_exp) - len(selected)
)

# ------------------------------------------------------------
# 3. Retrieve the exact previously verified S1 workbook.
#    No downloadable data are saved into the repository.
# ------------------------------------------------------------

request = Request(
    s1_url,
    headers={
        "User-Agent": "recount3-clinical-audit/1.0"
    }
)

with urlopen(request, timeout=60) as response:
    data = response.read(20000001)

assert 0 < len(data) <= 20000000

assert hashlib.sha256(data).hexdigest() == s1_sha, \
    "S1_SHA_CHANGED"

print("PLOS_S1_SHA256_PASS", s1_sha)

# ------------------------------------------------------------
# 4. Parse XLSX XML with the Python standard library.
# ------------------------------------------------------------

with zipfile.ZipFile(io.BytesIO(data)) as z:
    names = set(z.namelist())

    sheets = [
        p for p in names
        if re.fullmatch(
            r"xl/worksheets/sheet\d+\.xml", p
        )
    ]

    assert len(sheets) == 1, "SHEET_COUNT_CHANGED"

    shared = []

    if "xl/sharedStrings.xml" in names:
        root_ss = ET.fromstring(
            z.read("xl/sharedStrings.xml")
        )

        for node in root_ss.iter():
            if tag(node) == "si":
                shared.append(
                    "".join(
                        x.text or ""
                        for x in node.iter()
                        if tag(x) == "t"
                    )
                )

    xml = ET.fromstring(z.read(sheets[0]))
    grid = {}

    for cell in xml.iter():
        if tag(cell) != "c":
            continue

        coord = column(cell.get("r", ""))

        value_node = next(
            (x for x in cell if tag(x) == "v"),
            None
        )

        if cell.get("t") == "s" \
                and value_node is not None:
            value = shared[int(value_node.text)]

        elif cell.get("t") == "inlineStr":
            value = "".join(
                x.text or ""
                for x in cell.iter()
                if tag(x) == "t"
            )

        else:
            value = (
                value_node.text
                if value_node is not None else ""
            )

        if norm(value):
            assert coord not in grid, \
                "DUPLICATE_COORDINATE"
            grid[coord] = value

# ------------------------------------------------------------
# 5. Confirm the previously locked patient axis.
# ------------------------------------------------------------

hits = defaultdict(list)

for (r, c), value in grid.items():
    if norm(value) in isolates:
        hits[norm(value)].append((r, c))

assert len(hits) == 103
assert all(len(x) == 1 for x in hits.values())

patient_rows = {
    x[0][0] for x in hits.values()
}
patient_cols = {
    x[0][1] for x in hits.values()
}

assert len(patient_rows) == 103
assert patient_cols == {1}
assert min(patient_rows) == 3

patient_heading = (
    str(grid.get((1, 1), "")) + " " +
    str(grid.get((2, 1), ""))
)

assert "patient" in norm(patient_heading)

print(
    "PATIENT_AXIS_RECONFIRMED",
    len(patient_rows),
    "ROWS_COLUMN_1"
)

# ------------------------------------------------------------
# 6. Inspect ONLY header rows 1 and 2.
#
# The 103 patient-level rows are never printed.
# ------------------------------------------------------------

cols = sorted({c for r, c in grid})
assert len(cols) == 18

print("WORKBOOK_COLUMNS", len(cols))

counts = Counter()

for col in cols:
    header1 = safe_header(
        grid.get((1, col), ""),
        isolates
    )

    header2 = safe_header(
        grid.get((2, col), ""),
        isolates
    )

    category = meaning([header1, header2])
    counts[category] += 1

    # This is a count only; patient values remain private.
    filled = sum(
        bool(norm(grid.get((r, col), "")))
        for r in patient_rows
    )

    print(
        "COLUMN", col,
        "HEADER_ROW1", repr(header1),
        "HEADER_ROW2", repr(header2),
        "SEMANTICS", category,
        "POPULATED_PATIENT_ROWS", filled
    )

print(
    "COLUMN_SEMANTIC_COUNTS",
    dict(sorted(counts.items()))
)

print(
    "TREATMENT_TIMING_HEADER_CANDIDATES",
    counts["POTENTIAL_THERAPY_TIMING_OR_COLLECTION"]
)

# ------------------------------------------------------------
# 7. No clinical upgrade based on header screening.
# ------------------------------------------------------------

print("TREATMENT_NAIVE_STATUS NOT_ESTABLISHED")
print("CROSS_STUDY_INDEPENDENCE UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("NO_PATIENT_LEVEL_VALUES_PRINTED")
print("NO_FROZEN_INPUTS_MODIFIED")

print(
    "V6_DRP005276_S1_REMAINING_HEADERS_V1_COMPLETE"
)
