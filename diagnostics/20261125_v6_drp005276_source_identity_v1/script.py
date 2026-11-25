#!/usr/bin/env python3
"""V6 DRP005276 original-SRA specimen and patient semantics.

Read-only scientific diagnostic. No clinical admissions.
No patient/sample identifier values printed.
Compatible with Python >=3.6.
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
study = "DRP005276"

source_path = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)
queue_path = repo / (
    "diagnostics/20261124_v6_599_adjudication_queue/"
    "599_clinical_review_queue_LOCAL_ONLY.tsv"
)
raw_path = rec / (
    "raw_studies/DRP005276/metadata/"
    "sra.sra.DRP005276.MD.gz"
)

EXPECTED = {
    source_path:
        "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb",
    queue_path:
        "b5fdf922e4791b4b30c2ffa798c6ffed99b43478a0681c02a1faa2aaa248f3bd",
}

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        assert reader.fieldnames is not None
        return list(reader)

def norm(x):
    return " ".join(str(x or "").strip().casefold().split())

def parse_attrs(text):
    out = defaultdict(set)
    for fragment in (text or "").split("|"):
        if ";;" not in fragment:
            continue
        key, value = fragment.split(";;", 1)
        key, value = norm(key), norm(value)
        if key and value:
            out[key].add(value)
    return out

def structure(values):
    vals = [v for v in values if v]
    counts = Counter(vals)
    return {
        "populated": len(vals),
        "distinct": len(counts),
        "repeated_groups": sum(n > 1 for n in counts.values()),
        "experiments_in_repeated_groups":
            sum(n for n in counts.values() if n > 1),
    }

print("V6_DRP005276_SOURCE_IDENTITY_DIAGNOSTIC")
print("STUDY", study)

# ----------------------------------------------------------
# 1. Verify original local V6 evidence, without rebuilding.
# ----------------------------------------------------------

for path, expected in EXPECTED.items():
    assert path.is_file(), "Required frozen source missing"
    assert sha256(path) == expected, "Frozen source SHA mismatch"

source_all = read_tsv(source_path)
queue_all = read_tsv(queue_path)

assert len(source_all) == 599
assert len(queue_all) == 599

source = {
    r["experiment_acc"]: r
    for r in source_all if r["study"] == study
}
queue = {
    r["experiment_acc"]: r
    for r in queue_all if r["study"] == study
}

assert len(source) == len(queue) == 83
assert set(source) == set(queue)

assert all(
    r["raw_source_status"] == "EXACT_EXPERIMENT_MATCH"
    and r["clinical_approval"] == "NO_UNSIGNED"
    for r in source.values()
)
assert all(
    r["final_clinical_approval"] == "NO_UNSIGNED"
    for r in queue.values()
)

print("\nFROZEN_EVIDENCE")
print("SOURCE_SHA256_CHECKS PASS")
print("SOURCE_SELECTED_EXPERIMENTS", len(source))
print("QUEUE_SELECTED_EXPERIMENTS", len(queue))
print("SOURCE_QUEUE_EXACT_ID_SET_MATCH", True)

# ----------------------------------------------------------
# 2. Original compressed SRA source, including sample_name.
# ----------------------------------------------------------

assert raw_path.is_file(), "Original SRA study metadata missing"

source_hashes = {
    r["source_file_sha256"]
    for r in source.values()
}
assert len(source_hashes) == 1
assert sha256(raw_path) == next(iter(source_hashes)), (
    "Original SRA checksum differs from frozen V6 evidence"
)

raw_by_exp = defaultdict(list)

with gzip.open(
    str(raw_path), "rt",
    encoding="utf-8",
    errors="replace",
    newline=""
) as f:
    reader = csv.DictReader(f, delimiter="\t")
    raw_fields = set(reader.fieldnames or [])
    assert {
        "experiment_acc", "sample_title",
        "sample_attributes"
    } <= raw_fields, "Unexpected original SRA schema"

    for r in reader:
        acc = r["experiment_acc"]
        raw_by_exp[acc].append(r)

assert set(source) <= set(raw_by_exp)

raw = {}
source_multirows = Counter()

for acc, src in source.items():
    original = raw_by_exp[acc]
    source_multirows[len(original)] += 1

    for row in original:
        assert row["sample_title"] == src["original_sample_title"]
        assert row["sample_attributes"] == src["original_sample_attributes"]
        if "sample_acc" in raw_fields:
            assert row["sample_acc"] == src["sra_sample_acc"]

    # Multiple sequencing runs are not multiple specimens.
    for field in (
        "sample_acc", "sample_name", "biosample",
        "sample_title", "sample_attributes"
    ):
        if field in raw_fields:
            assert len({
                r.get(field, "") for r in original
            }) == 1, "Original source specimen field conflict"

    raw[acc] = original[0]

print("\nORIGINAL_SRA_RECONCILIATION")
print("ORIGINAL_SRA_SHA256", sha256(raw_path))
print("RAW_STUDY_DISTINCT_EXPERIMENTS", len(raw_by_exp))
print("SELECTED_EXACT_SOURCE_MATCHES", len(raw))
print("SELECTED_METADATA_ROWS_PER_EXPERIMENT",
      dict(sorted(source_multirows.items())))
print("ORIGINAL_SAMPLE_NAME_FIELD_PRESENT",
      "sample_name" in raw_fields)
print("ORIGINAL_SAMPLE_ACC_FIELD_PRESENT",
      "sample_acc" in raw_fields)
print("ORIGINAL_BIOSAMPLE_FIELD_PRESENT",
      "biosample" in raw_fields)

# ----------------------------------------------------------
# 3. Source sample_name and isolate are separate concepts.
# ----------------------------------------------------------

attrs = {
    acc: parse_attrs(src["original_sample_attributes"])
    for acc, src in source.items()
}

sample_names = {
    acc: norm(row.get("sample_name", ""))
    for acc, row in raw.items()
}
sample_accs = {
    acc: norm(row.get("sample_acc", ""))
    for acc, row in raw.items()
}
biosamples = {
    acc: norm(row.get("biosample", ""))
    for acc, row in raw.items()
}

isolates = {}
for acc, a in attrs.items():
    vals = a.get("isolate", set())
    isolates[acc] = next(iter(vals)) if len(vals) == 1 else ""

print("\nIDENTIFIER_STRUCTURE")
print("SRA_SAMPLE_NAMES", structure(sample_names.values()))
print("SRA_SAMPLE_ACCESSIONS", structure(sample_accs.values()))
print("BIOSAMPLE_VALUES", structure(biosamples.values()))
print("SOURCE_ISOLATES", structure(isolates.values()))
print("ISOLATE_MULTIVALUE",
      sum(len(a.get("isolate", set())) > 1
          for a in attrs.values()))

for label, left, right in (
    ("SAMPLE_NAME_VS_ISOLATE", sample_names, isolates),
    ("SAMPLE_NAME_VS_SRA_SAMPLE", sample_names, sample_accs),
    ("ISOLATE_VS_SRA_SAMPLE", isolates, sample_accs),
):
    comparable = [
        acc for acc in source if left[acc] and right[acc]
    ]
    matches = sum(
        left[acc] == right[acc] for acc in comparable
    )
    print(label,
          "COMPARABLE", len(comparable),
          "EXACT_MATCHES", matches)

# ----------------------------------------------------------
# 4. Specimen, diagnosis, sampled site and timing signals.
# ----------------------------------------------------------

keys_of_interest = (
    "tissue", "disease", "sample type", "sample_type",
    "specimen type", "tumor type", "tumour type",
    "tumor location", "tumour location",
    "primary site", "metastatic site",
    "patient", "patient id", "patient_id",
    "subject id", "subject_id",
    "donor", "individual", "case id",
    "treatment", "treatment status",
    "therapy", "prior treatment",
    "time point", "collection timepoint",
    "isolate"
)

print("\nSOURCE_ATTRIBUTE_COVERAGE")
for key in keys_of_interest:
    populated = sum(
        bool(a.get(key)) for a in attrs.values()
    )
    if populated:
        multi = sum(
            len(a.get(key, set())) > 1
            for a in attrs.values()
        )
        print("KEY", repr(key),
              "POPULATED", populated,
              "MULTIVALUE", multi)

print("\nEXACT_TISSUE_DISEASE_TERMS")
print("TISSUE_EQUALS_BLADDER", sum(
    "bladder" in a.get("tissue", set())
    for a in attrs.values()
))
print("DISEASE_EQUALS_BLADDER_CANCER", sum(
    "bladder cancer" in a.get("disease", set())
    for a in attrs.values()
))

terms = (
    "primary", "normal", "tumor", "tumour",
    "metasta", "recur", "pretreat", "untreated",
    "neoadjuvant", "chemotherapy", "surgery",
    "cystectomy", "biopsy"
)

print("\nSOURCE_TEXT_TERM_COVERAGE_NOT_CERTIFICATION")

for term in terms:
    title_n = sum(
        term in norm(r["original_sample_title"])
        for r in source.values()
    )
    attr_n = sum(
        term in norm(r["original_sample_attributes"])
        for r in source.values()
    )
    print("TERM", term,
          "TITLE_EXPERIMENTS", title_n,
          "ATTRIBUTE_EXPERIMENTS", attr_n)

# ----------------------------------------------------------
# 5. Existing curated patient and treatment grades.
# ----------------------------------------------------------

print("\nFROZEN_CLINICAL_EVIDENCE_GRADES")

for field in (
    "historical_tissue_grade",
    "historical_patient_grade",
    "historical_treatment_grade",
    "provisional_disposition",
    "final_clinical_approval"
):
    counts = Counter(r[field] for r in queue.values())
    print(field, dict(sorted(counts.items())))

patient_keys = {
    acc: norm(r.get("v5_patient_key_LOCAL_ONLY", ""))
    for acc, r in queue.items()
}

print("\nCURATED_PATIENT_KEY_STRUCTURE_NOT_CERTIFIED")
print("PATIENT_KEYS", structure(patient_keys.values()))

print("PATIENT_KEY_EXACTLY_EQUALS_ISOLATE", sum(
    patient_keys[acc] != ""
    and patient_keys[acc] == isolates[acc]
    for acc in source
))
print("PATIENT_KEY_EXACTLY_EQUALS_SAMPLE_NAME", sum(
    patient_keys[acc] != ""
    and patient_keys[acc] == sample_names[acc]
    for acc in source
))

timing = [
    norm(r.get("v5_treatment_timing_LOCAL_ONLY", ""))
    for r in queue.values()
]

print("CURATED_TREATMENT_TIMING_POPULATED",
      sum(bool(x) for x in timing))
print("CURATED_TREATMENT_TIMING_DISTINCT_TOKENS",
      len(set(x for x in timing if x)))

# ----------------------------------------------------------
# 6. Study-scope classification; no auto admissions.
# ----------------------------------------------------------

print("\nSCIENTIFIC_GATES")
print("EXACT_SRA_JOIN", "PASS")
print("INDEPENDENT_BIOLOGICAL_SPECIMENS", "NOT_CERTIFIED")
print("INDEPENDENT_PATIENTS", "NOT_CERTIFIED")
print("PRIMARY_TUMOUR_SAMPLED_SITE", "NOT_CERTIFIED")
print("PRETREATMENT_TIMING", "NOT_CERTIFIED")
print("CROSS_STUDY_ATLAS_INDEPENDENCE", "NOT_CERTIFIED")
print("MODEL_TRANSFER", "NOT_PERFORMED")
print("CLINICAL_ADMISSIONS", 0)

print("\nNO_PATIENT_OR_SAMPLE_VALUES_PRINTED")
print("NO_SOURCE_FILES_WRITTEN")
print("NO_FROZEN_DATA_MODIFIED")
print("V6_DRP005276_SOURCE_IDENTITY_DIAGNOSTIC_COMPLETE")
