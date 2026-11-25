#!/usr/bin/env python3
"""DRP005276: original study description and patient-key provenance.

Read-only. Python 3.6 compatible.
No individual patient/sample identifiers printed.
"""
import csv
import gzip
import hashlib
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"

raw_path = rec / (
    "raw_studies/DRP005276/metadata/"
    "sra.sra.DRP005276.MD.gz"
)
queue_path = repo / (
    "diagnostics/20261124_v6_599_adjudication_queue/"
    "599_clinical_review_queue_LOCAL_ONLY.tsv"
)

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

assert digest(raw_path) == (
    "cba60a28b225d7852e17e180ddaad2680504e7f1cd5cf94de5d73401db12d36c"
)
assert digest(queue_path) == (
    "b5fdf922e4791b4b30c2ffa798c6ffed99b43478a0681c02a1faa2aaa248f3bd"
)

def norm(value):
    return " ".join(str(value or "").strip().casefold().split())

def attrs(value):
    result = defaultdict(set)
    for part in (value or "").split("|"):
        if ";;" in part:
            k, v = part.split(";;", 1)
            if norm(k) and norm(v):
                result[norm(k)].add(norm(v))
    return result

def redact(text):
    text = " ".join(str(text or "").split())
    text = re.sub(
        r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b",
        "[EMAIL_REDACTED]", text, flags=re.I
    )
    text = re.sub(
        r"\b(?:DRX|DRR|DRS|SRX|SRR|SRS|ERX|ERR|ERS|"
        r"SAMN|SAMD|SAMEA|GSM)\d+\b",
        "[ACCESSION_REDACTED]", text
    )
    return text

def token_shape(value):
    s = str(value or "")
    if not s:
        return "<EMPTY>"
    shape = []
    for char in s:
        if char.isalpha():
            shape.append("A")
        elif char.isdigit():
            shape.append("9")
        elif char in "-_./":
            shape.append(char)
        else:
            shape.append("?")
    return "".join(shape)[:100]

def summarize_shapes(label, values):
    values = list(values)
    nonempty = [x for x in values if x]
    counts = Counter(token_shape(x) for x in nonempty)
    print(
        label, "POPULATED", len(nonempty),
        "DISTINCT", len(set(nonempty)),
        "SHAPE_CLASS_COUNT", len(counts)
    )
    for pattern, n in counts.most_common(6):
        print(label, "SHAPE", pattern, "COUNT", n)

print("V6_DRP005276_STUDY_CONTEXT_V1")
print("FROZEN_INPUT_SHA256_PASS 2")

with gzip.open(
    str(raw_path), "rt",
    encoding="utf-8", errors="replace",
    newline=""
) as f:
    reader = csv.DictReader(f, delimiter="\t")
    fields = set(reader.fieldnames or [])
    assert {
        "experiment_acc", "study_title", "study_abstract",
        "sample_name", "sample_acc",
        "sample_title", "sample_attributes"
    } <= fields
    raw = list(reader)

assert len(raw) == 103
by_exp = {r["experiment_acc"]: r for r in raw}
assert len(by_exp) == 103

with queue_path.open(
    encoding="utf-8-sig", newline=""
) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

assert len(rows) == 599

selected = {
    r["experiment_acc"]: r
    for r in rows if r["study"] == "DRP005276"
}
assert len(selected) == 83
assert set(selected) <= set(by_exp)

print("\nORIGINAL_STUDY_DESIGN")

for field in ("study_title", "study_abstract"):
    vals = {
        " ".join(r[field].split())
        for r in raw if r.get(field, "").strip()
    }
    print("FIELD", field, "POPULATED",
          sum(bool(r.get(field, "").strip()) for r in raw))
    print("FIELD", field, "DISTINCT", len(vals))
    assert len(vals) == 1

    value = redact(next(iter(vals)))
    print("FIELD", field, "TEXT_START")
    print(value[:3000])
    print("FIELD", field, "TRUNCATED", len(value) > 3000)
    print("FIELD", field, "TEXT_END")

print("\nSOURCE_PUBLICATION_METADATA_FIELDS")
for field in (
    "study_accession", "secondary_study_accession",
    "bioproject", "study_alias",
    "study_description", "publication",
    "pubmed_id", "doi", "center_name"
):
    if field in fields:
        vals = {
            r.get(field, "").strip()
            for r in raw if r.get(field, "").strip()
        }
        print("FIELD", field, "POPULATED",
              sum(bool(r.get(field, "").strip()) for r in raw),
              "DISTINCT", len(vals))
    else:
        print("FIELD", field, "ABSENT")

print("\nPATIENT_KEY_PROVENANCE_SELECTED_83")

def single_isolate(row):
    values = attrs(row["sample_attributes"]).get("isolate", set())
    return next(iter(values)) if len(values) == 1 else ""

key_values = [
    norm(r.get("v5_patient_key_LOCAL_ONLY", ""))
    for r in selected.values()
]
summarize_shapes("CURATED_PATIENT_KEY", key_values)
summarize_shapes(
    "SRA_SAMPLE_NAME",
    [norm(by_exp[a]["sample_name"]) for a in selected]
)
summarize_shapes(
    "SRA_ISOLATE",
    [norm(single_isolate(by_exp[a])) for a in selected]
)

for label, field in (
    ("EXPERIMENT_ACCESSION", "experiment_acc"),
    ("SRA_SAMPLE_ACCESSION", "sample_acc"),
    ("SRA_SAMPLE_NAME", "sample_name"),
    ("SRA_SAMPLE_TITLE", "sample_title")
):
    matching = sum(
        norm(r.get("v5_patient_key_LOCAL_ONLY", ""))
        == norm(a if field == "experiment_acc"
                else by_exp[a].get(field, ""))
        for a, r in selected.items()
        if norm(r.get("v5_patient_key_LOCAL_ONLY", ""))
    )
    print("CURATED_KEY_EQUALS", label, matching)

print("CURATED_KEY_EQUALS_SOURCE_ISOLATE", sum(
    norm(r.get("v5_patient_key_LOCAL_ONLY", ""))
    == norm(single_isolate(by_exp[a]))
    for a, r in selected.items()
    if norm(r.get("v5_patient_key_LOCAL_ONLY", ""))
))

print("\nSCIENTIFIC_STATUS")
print("SOURCE_STUDY_IDENTITY", "DESCRIPTION_RETRIEVED")
print("PUBLICATION_CROSSWALK", "REQUIRES_VERIFICATION")
print("BIOLOGICAL_PATIENT_INDEPENDENCE", "UNSIGNED")
print("PRIMARY_TUMOUR_STATUS", "UNSIGNED")
print("PRETREATMENT_STATUS", "UNSIGNED")
print("CLINICAL_ADMISSIONS 0")
print("NO_INDIVIDUAL_PATIENT_OR_SAMPLE_VALUES_PRINTED")
print("NO_SOURCE_FILES_MODIFIED")
print("V6_DRP005276_STUDY_CONTEXT_V1_COMPLETE")
