import csv
import gzip
import hashlib
import io
import posixpath
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

root = Path(sys.argv[1])
xlsx = Path(sys.argv[2])
base = root / "26_recount3_external_primary_v1"

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def norm(x):
    return str(x or "").strip().upper()

print("SRP185703_PUBLISHED_SRA_LINKABILITY_V1")

expected = (
    "0be42e4d7ab9b66e3934c139cd57a1cba13a52ba2404c318cb989073610d75b5"
)
assert sha256(xlsx) == expected
print("JCI_WORKBOOK_HASH_PASS")

v5 = base / (
    "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"
)

with zipfile.ZipFile(v5) as z:
    names = [
        n for n in z.namelist()
        if n.endswith(
            "/results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
        )
    ]
    assert len(names) == 1
    reader = csv.DictReader(
        io.StringIO(z.read(names[0]).decode("utf-8-sig")),
        delimiter="\t"
    )
    rows = list(reader)

assert len(rows) == 2055
study = [r for r in rows if r.get("study") == "SRP185703"]
assert len(study) == 75
assert len({r["library_id"] for r in study}) == 75

print("V5_PRIMARY_MEMBERSHIP_PASS 75")

fields = (
    "library_id",
    "source_SRA_sample_acc",
    "biosample",
    "v5_original_sra_sample_name"
)

ids = {}

for field in fields:
    assert field in study[0]
    values = {
        norm(r.get(field)) for r in study
        if len(norm(r.get(field))) >= 4
    }
    ids[field] = values
    print("SOURCE_FIELD", field, "MATCHABLE_VALUES", len(values))

M = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
P = "http://schemas.openxmlformats.org/package/2006/relationships"

with zipfile.ZipFile(xlsx) as z:
    shared = []

    if "xl/sharedStrings.xml" in z.namelist():
        xml = ET.fromstring(z.read("xl/sharedStrings.xml"))
        shared = [
            "".join(t.text or "" for t in item.iter("{%s}t" % M))
            for item in xml.findall("{%s}si" % M)
        ]

    book = ET.fromstring(z.read("xl/workbook.xml"))
    relationships = ET.fromstring(
        z.read("xl/_rels/workbook.xml.rels")
    )

    targets = {
        r.attrib["Id"]: r.attrib["Target"]
        for r in relationships.findall("{%s}Relationship" % P)
    }

    sheets = book.findall(".//{%s}sheet" % M)
    assert len(sheets) == 7

    for sheet in sheets:
        name = sheet.attrib["name"]
        target = targets[sheet.attrib["{%s}id" % R]]

        member = (
            target.lstrip("/")
            if target.startswith("/")
            else posixpath.normpath(
                posixpath.join("xl", target)
            )
        )
        assert member in z.namelist()

        xml = ET.fromstring(z.read(member))
        values = []

        for cell in xml.findall(
            ".//{%s}sheetData/{%s}row/{%s}c" % (M, M, M)
        ):
            typ = cell.get("t", "")
            node = cell.find("{%s}v" % M)

            if typ == "inlineStr":
                value = "".join(
                    t.text or ""
                    for t in cell.iter("{%s}t" % M)
                )
            elif node is None:
                continue
            elif typ == "s":
                value = shared[int(node.text)]
            else:
                value = node.text or ""

            if str(value).strip():
                values.append(norm(value))

        cellset = set(values)

        print("SHEET", name, "POPULATED_CELLS", len(values))

        for field in fields:
            print(
                "EXACT_MATCH",
                name,
                field,
                len(ids[field] & cellset)
            )

# Inspect the native SRA metadata only for clinical-text
# availability. Do not publish individual sample values.

native = base / "raw_studies/SRP185703/metadata"
print("NATIVE_METADATA_DIR_PRESENT", native.is_dir())

pattern = re.compile(
    r"untreated|treatment.?na[iï]ve|imatinib|"
    r"neoadjuvant|metastatic|prior.?treat|"
    r"pre.?treat|primary.?tumou?r",
    re.I
)

files = sorted(native.glob("*.MD.gz")) if native.is_dir() else []
print("NATIVE_METADATA_FILES", len(files))

for path in files:
    counts = Counter()
    n = 0

    with gzip.open(
        path, "rt",
        encoding="utf-8-sig",
        errors="replace",
        newline=""
    ) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            n += 1
            for key, value in row.items():
                if value is not None and pattern.search(str(value)):
                    counts[key] += 1

    print(
        "NATIVE_FILE",
        path.name,
        "ROWS", n,
        "CLINICAL_HIT_COLUMNS", len(counts)
    )

    for key in sorted(counts):
        if re.fullmatch(r"[A-Za-z0-9_./-]{1,80}", key):
            print(
                "CLINICAL_HIT_FIELD",
                key,
                "ROWS", counts[key]
            )

print("EXACT_MATCHES_ARE_LINK_CANDIDATES_NOT_CLINICAL_CERTIFICATION")
print("NO_INDIVIDUAL_IDENTIFIERS_PRINTED")
print("PRETREATMENT_UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("SRP185703_PUBLISHED_SRA_LINKABILITY_V1_COMPLETE")
