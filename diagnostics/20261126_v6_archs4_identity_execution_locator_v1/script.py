#!/usr/bin/env python3
"""Inspect original ARCHS4 identity script and README as inert text.

Do not import or execute historical program.
Do not read any biological sample records.
Python 3.6 compatible.
"""

import ast
import hashlib
import re
import sys
import zipfile

from pathlib import Path

root = Path(sys.argv[1])

archive = root / (
    "26_recount3_external_primary_v1/"
    "ARCHS4_RECOUNT3_EXACT_IDENTITY_JOIN_READONLY_20261120.zip"
)

assert archive.is_file()

with zipfile.ZipFile(str(archive)) as z:

    names = set(z.namelist())

    pyname = "archs4_recount3_exact_identity_join_READONLY.py"
    rdname = "ARCHS4_RECOUNT3_EXACT_IDENTITY_JOIN_README.md"

    assert {
        pyname, rdname, "SHA256SUMS.txt"
    } == names

    code = z.read(pyname).decode("utf-8-sig")
    readme = z.read(rdname).decode("utf-8-sig")
    manifest = z.read("SHA256SUMS.txt").decode("utf-8-sig")

print("ARCHS4_IDENTITY_EXECUTION_LOCATOR_V1")

zip_sha = hashlib.sha256(archive.read_bytes()).hexdigest()

print("ORIGINAL_ARCHIVE_SHA256", zip_sha)

assert zip_sha == (
    "20c22c10e07ec8d9bc5847b42c65eae3"
    "786614eabdf675afd46ad21b48502c05"
)

# Static syntax validation; historical program is not executed.
ast.parse(code)

print("SOURCE_SCRIPT_AST_PASS")
print(
    "SCRIPT_SHA256",
    hashlib.sha256(code.encode("utf-8")).hexdigest()
)
print(
    "README_SHA256",
    hashlib.sha256(readme.encode("utf-8")).hexdigest()
)

# Compare the source package's embedded SHA256 manifest.
validated = 0

for line in manifest.splitlines():

    match = re.fullmatch(
        r"([a-f0-9]{64})\s+\*?(\S+)",
        line.strip()
    )

    if not match:
        print("MANIFEST_UNPARSED_LINE")
        continue

    digest, name = match.groups()

    if name not in (pyname, rdname):
        print("MANIFEST_ENTRY_OTHER_FILE")
        continue

    data = (
        code.encode("utf-8")
        if name == pyname
        else readme.encode("utf-8")
    )

    assert hashlib.sha256(data).hexdigest() == digest, \
        "PACKAGE_MANIFEST_MISMATCH"

    validated += 1

print("INTERNAL_MANIFEST_VALIDATED_MEMBERS", validated)

# Redact accession-like identifiers and email addresses.
# Only documentation and static program text are displayed.
def safe(line):

    line = re.sub(
        r"\b(?:DRX|DRR|DRS|SRX|SRR|SRS|"
        r"ERX|ERR|ERS|GSM|SAMN|SAMD|SAMEA)\d+\b",
        "[ACCESSION]",
        line,
        flags=re.I
    )

    line = re.sub(
        r"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}",
        "[EMAIL]",
        line,
        flags=re.I
    )

    return line[:250]

# ------------------------------------------------------------
# 1. Show the original README.
# Small document; no patient-level source tables.
# ------------------------------------------------------------

print("\nORIGINAL_README_LINES")

for number, line in enumerate(readme.splitlines(), 1):
    print(
        "README_LINE",
        number,
        safe(line)
    )

# ------------------------------------------------------------
# 2. Inspect code lines related to input/output provenance.
#
# No import, execution or source table processing.
# ------------------------------------------------------------

pattern = re.compile(
    r"add_argument|--out|--root|"
    r"\.tsv|\.csv|\.json|sha256|"
    r"source|output|read_|write_|"
    r"atlas|GSM|biosample",
    re.I
)

selected = [
    (number, line)
    for number, line in enumerate(code.splitlines(), 1)
    if pattern.search(line)
]

print(
    "\nSOURCE_CODE_RELEVANT_LINE_COUNT",
    len(selected)
)

for number, line in selected[:105]:
    print(
        "SOURCE_LINE",
        number,
        safe(line)
    )

if len(selected) > 105:
    print(
        "SOURCE_LINES_OMITTED",
        len(selected) - 105
    )

print("ORIGINAL_CROSSWALK_PROGRAM_EXECUTED NO")
print("ORIGINAL_RESULTS_RECOVERED NOT_YET")
print("NO_BIOLOGICAL_DATA_READ")
print("NO_SOURCE_FILES_MODIFIED")

print("ARCHS4_IDENTITY_EXECUTION_LOCATOR_V1_COMPLETE")
