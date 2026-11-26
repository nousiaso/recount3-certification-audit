#!/usr/bin/env python3
"""Read-only OOXML inventory of JCI124108 supplemental tables."""

import hashlib
import posixpath
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

path = Path(sys.argv[1])

N = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "p": "http://schemas.openxmlformats.org/package/2006/relationships",
}

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()

def text_of(node):
    return "".join(
        x.text or "" for x in node.findall(".//m:t", N)
    )

def cell_value(cell, shared):
    typ = cell.get("t", "")
    if typ == "inlineStr":
        return text_of(cell)
    node = cell.find("m:v", N)
    if node is None:
        return ""
    value = node.text or ""
    if typ == "s":
        idx = int(value)
        return shared[idx] if idx < len(shared) else ""
    return value

header_word = re.compile(
    r"patient|sample|specimen|tumou?r|treat|therapy|"
    r"status|location|site|origin|primary|metasta|"
    r"genotyp|mutat|gastric|diagnos|response|"
    r"anatom|recur|resect|histolog|kit|pdgfra",
    re.I
)

safe_header = re.compile(r"^[A-Za-z][A-Za-z _/-]{2,65}$")

print("SRP185703_JCI_SUPPLEMENT_SCHEMA_V1")
print("SOURCE_DOI 10.1172/JCI124108")
print("PUBLIC_XLSX_SHA256", sha256(path))
print("PUBLIC_XLSX_BYTES", path.stat().st_size)

with zipfile.ZipFile(str(path)) as z:
    assert "xl/workbook.xml" in z.namelist()

    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
        shared = [
            text_of(si)
            for si in ss.findall("m:si", N)
        ]

    book = ET.fromstring(z.read("xl/workbook.xml"))
    relationships = ET.fromstring(
        z.read("xl/_rels/workbook.xml.rels")
    )

    relmap = {
        r.get("Id"): r.get("Target")
        for r in relationships.findall("p:Relationship", N)
    }

    sheets = book.findall("m:sheets/m:sheet", N)
    print("WORKSHEET_COUNT", len(sheets))
    print("SHARED_STRING_COUNT", len(shared))

    for sheet in sheets:
        name = sheet.get("name", "")
        rid = sheet.get("{" + N["r"] + "}id")
        target = relmap.get(rid)

        assert target, "MISSING_WORKSHEET_RELATIONSHIP"

        if target.startswith("/"):
            member = target.lstrip("/")
        else:
            member = posixpath.normpath(
                posixpath.join("xl", target)
            )

        assert member.startswith("xl/worksheets/")
        assert member in z.namelist()

        ws = ET.fromstring(z.read(member))
        rows = ws.findall(".//m:sheetData/m:row", N)

        hints = []
        labels = Counter()
        filled = 0

        for row in rows:
            row_number = int(row.get("r", "0"))

            for cell in row.findall("m:c", N):
                value = cell_value(cell, shared).strip()

                if not value:
                    continue

                filled += 1

                if value.lower() in (
                    "untreated", "sensitive",
                    "resistant", "primary",
                    "metastatic", "recurrent"
                ):
                    labels[value.lower()] += 1

                if row_number > 12:
                    continue

                if not safe_header.fullmatch(value):
                    continue

                if header_word.search(value):
                    hints.append(
                        (row_number, cell.get("r"), value)
                    )

        print("\nSHEET", name)
        print("ROW_COUNT", len(rows))
        print("NONEMPTY_CELLS", filled)
        print("EARLY_HEADER_HINT_COUNT", len(hints))

        for rowno, coord, label in hints[:35]:
            print(
                "EARLY_HEADER_HINT",
                "ROW", rowno,
                "CELL", coord,
                "NAME", label
            )

        print(
            "CLASS_LABEL_CELL_COUNTS",
            repr(dict(sorted(labels.items())))
        )

print("\nIMPORTANT_LIMITS")
print("COLUMN_HINTS_NOT_YET_VALIDATED_AS_HEADERS")
print("LABEL_COUNTS_NOT_PATIENT_COUNTS")
print("NO_PATIENT_IDS_OR_CLINICAL_ROW_VALUES_PRINTED")
print("NO_SRA_CLINICAL_CROSSWALK_PERFORMED")
print("NO_FROZEN_DATA_MODIFIED")
print("CLINICAL_ADMISSIONS 0")
print("SRP185703_JCI_SUPPLEMENT_SCHEMA_V1_COMPLETE")
