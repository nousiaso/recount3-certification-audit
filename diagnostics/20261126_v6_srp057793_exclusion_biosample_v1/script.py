#!/usr/bin/env python3
"""Read-only SRP057793 exclusions and BioSample identity audit."""

import csv
import gzip
import io
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
study = "SRP057793"

stage41 = root / (
    "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"
)

fullpath = stage41 / "01_full_source_experiment_identity_5221.tsv"
candpath = stage41 / "03_candidate_integrated_evidence_3602.tsv"

v5 = root / (
    "26_recount3_external_primary_v1/"
    "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"
)

def local(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        return rd.fieldnames or [], list(rd)

ff, original = local(fullpath)
cf, candidates = local(candpath)

with zipfile.ZipFile(v5) as z:
    names = [
        n for n in z.namelist()
        if n.endswith(
            "/results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
        )
    ]
    assert len(names) == 1
    rd = csv.DictReader(
        io.StringIO(z.read(names[0]).decode("utf-8-sig")),
        delimiter="\t"
    )
    primary = list(rd)

f = [r for r in original if r.get("study") == study]
c = [r for r in candidates if r.get("study") == study]
p = [r for r in primary if r.get("study") == study]

assert (len(f), len(c), len(p)) == (149, 148, 140)

fid = {r["experiment_acc"] for r in f}
cid = {r["library_id"] for r in c}
pid = {r["library_id"] for r in p}

held = [r for r in c if r.get("technical_hold") == "YES"]
nonheld = [r for r in c if r.get("technical_hold") == "NO"]

hid = {r["library_id"] for r in held}
nid = {r["library_id"] for r in nonheld}

assert (len(fid), len(cid), len(pid), len(hid), len(nid)) == (
    149, 148, 140, 8, 140
)

print("SRP057793_EXCLUSION_BIOSAMPLE_V1")

print("\nEXACT_MEMBERSHIP")
print("ORIGINAL_ONLY_EXPERIMENTS", len(fid - cid))
print("CANDIDATES_NOT_IN_ORIGINAL", len(cid - fid))
print("HELD_CANDIDATES", len(hid))
print("NONHELD_CANDIDATES", len(nid))
print("NONHELD_EQUALS_V5_PRIMARY", nid == pid)
print("HELD_INTERSECTION_PRIMARY", len(hid & pid))
print("ORIGINAL_MINUS_PRIMARY", len(fid - pid))

# Identify the existing QC evidence columns.
qc_pattern = re.compile(
    r"qc|hold|reason|coverage|mapped|outlier|"
    r"technical|filter|exclude|detect|selection",
    re.I
)

qc_columns = [x for x in cf if qc_pattern.search(x)]

print("\nSTAGE41_QC_SCHEMA")
print("QC_RELEVANT_COLUMNS", ";".join(qc_columns))

for field in qc_columns:
    held_values = [
        str(r.get(field) or "").strip() for r in held
    ]
    populated = [x for x in held_values if x]

    if not populated:
        continue

    print(
        "HELD_QC_FIELD",
        field,
        "POPULATED", len(populated),
        "DISTINCT", len(set(populated))
    )

# Keyword counts only: these are discovery hints, not
# authoritative exclusion-reason assignments.
patterns = {
    "TOP10_COVERAGE": r"top10|top.?10.?coverage",
    "MAPPED_LENGTH": r"mapped.?len|mapped.?length",
    "STUDY_OUTLIER": r"outlier",
    "GENE_DETECTION": r"gene.?detect|detection.?fraction",
    "LOW_MAPPING": r"low.?map|mapping.?rate",
    "SHORT_TAG": r"short.?tag|dge",
}

print("\nQC_TERM_DISCOVERY")

for label, pattern in patterns.items():
    rx = re.compile(pattern, re.I)
    n = sum(
        any(rx.search(str(r.get(field) or ""))
            for field in qc_columns)
        for r in held
    )
    print("HELD_ROWS_WITH_TERM", label, n)

# The source-only experiment's metadata availability;
# never print its accession or sample-level attributes.
source_only = [
    r for r in f
    if r["experiment_acc"] in fid - cid
]
assert len(source_only) == 1

print("\nSOURCE_ONLY_METADATA")
print("SOURCE_ONLY_ROWS", len(source_only))

for field in ff:
    if re.search(
        r"biosample|sample|run|source|qc|reason|"
        r"selection|status|filter|hold",
        field, re.I
    ):
        print(
            "SOURCE_ONLY_FIELD",
            field,
            "POPULATED",
            bool(str(source_only[0].get(field) or "").strip())
        )

# Report BioSample cardinalities without identifiers.
print("\nBIOSAMPLE_CARDINALITIES")

for name, rows, fields in (
    ("ORIGINAL", f, ff),
    ("STAGE41", c, cf),
    ("V5_PRIMARY", p, list(p[0]))
):
    possible = [
        x for x in fields
        if re.search(r"biosample", x, re.I)
    ]

    print("SCOPE", name, "BIOSAMPLE_COLUMNS", ";".join(possible))

    for field in possible:
        values = [
            str(r.get(field) or "").strip()
            for r in rows
        ]
        values = [x for x in values if x]

        cts = Counter(values)
        repeated_groups = sum(v > 1 for v in cts.values())

        print(
            "BIOSAMPLE_FIELD", name, field,
            "POPULATED", len(values),
            "DISTINCT_RAW_VALUES", len(cts),
            "REPEATED_RAW_VALUE_GROUPS", repeated_groups
        )

# Original native SRA schema, used only for determining
# whether a more precise accession-to-BioSample join is possible.
native = root / (
    "26_recount3_external_primary_v1/"
    "raw_studies/SRP057793/metadata/"
    "sra.sra.SRP057793.MD.gz"
)

print("\nNATIVE_SRA_METADATA_PRESENT", native.is_file())

if native.is_file():
    with gzip.open(
        native, "rt", encoding="utf-8-sig",
        newline=""
    ) as handle:
        rd = csv.DictReader(handle, delimiter="\t")
        headers = rd.fieldnames or []
        native_rows = list(rd)

    print("NATIVE_SRA_ROWS", len(native_rows))
    print("NATIVE_SRA_FIELDS", len(headers))

    for field in headers:
        if re.search(
            r"biosample|experiment|sample.?acc",
            field, re.I
        ):
            values = [
                str(r.get(field) or "").strip()
                for r in native_rows
            ]
            values = [x for x in values if x]
            print(
                "NATIVE_ID_FIELD", field,
                "POPULATED", len(values),
                "DISTINCT", len(set(values))
            )

print("\nNO_INDIVIDUAL_IDENTIFIERS_PRINTED")
print("NO_QC_DECISIONS_MODIFIED")
print("NO_FROZEN_FILES_MODIFIED")
print("PATIENT_INDEPENDENCE_UNSIGNED")
print("PRETREATMENT_UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("SRP057793_EXCLUSION_BIOSAMPLE_V1_COMPLETE")
