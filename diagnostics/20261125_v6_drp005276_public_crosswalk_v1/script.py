#!/usr/bin/env python3
"""DRA008030 versus frozen DRP005276, plus PLOS S1 supplementary metadata.

Read-only diagnostic. External public sources fetched to memory only.
No individual patient, sample, run or experiment identifiers printed.
Python 3.6 compatible. No clinical approvals.
"""
import csv
import gzip
import hashlib
import io
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

root = Path(sys.argv[1])
rec = root / "26_recount3_external_primary_v1"
repo = root / "RECOUNT3_CERTIFICATION_GIT"
raw_path = rec / "raw_studies/DRP005276/metadata/sra.sra.DRP005276.MD.gz"
source_path = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)
expected = {
    raw_path: "cba60a28b225d7852e17e180ddaad2680504e7f1cd5cf94de5d73401db12d36c",
    source_path: "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb",
}

def digest(data):
    return hashlib.sha256(data).hexdigest()

def file_digest(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def normal(v):
    return " ".join(str(v or "").strip().casefold().split())

def local_name(tag):
    return tag.rsplit("}", 1)[-1]

def nodes(element, name):
    return [x for x in element.iter() if local_name(x.tag) == name]

def attr_values(text, key):
    vals = set()
    for part in (text or "").split("|"):
        if ";;" in part:
            k, v = part.split(";;", 1)
            if normal(k) == key and normal(v):
                vals.add(normal(v))
    return vals

def fetch(url, maximum):
    request = Request(url, headers={
        "User-Agent": "recount3-public-accession-provenance-audit/1.0",
        "Accept": "application/xml,application/octet-stream,*/*",
    })
    with urlopen(request, timeout=55) as response:
        payload = response.read(maximum + 1)
    assert 0 < len(payload) <= maximum, "REMOTE_SIZE_LIMIT_OR_EMPTY"
    return payload

print("V6_DRP005276_PUBLIC_DRA_PUBLICATION_CROSSWALK_V1")
for path, sha in expected.items():
    assert path.is_file() and file_digest(path) == sha, "FROZEN_LOCAL_SOURCE_CHANGED"
print("FROZEN_LOCAL_SHA256_PASS", len(expected))

with gzip.open(str(raw_path), "rt", encoding="utf-8", errors="replace", newline="") as fh:
    local_rows = list(csv.DictReader(fh, delimiter="\t"))
local = {r["experiment_acc"]: r for r in local_rows}
assert len(local) == len(local_rows) == 103
with source_path.open(encoding="utf-8-sig", newline="") as fh:
    source_rows = list(csv.DictReader(fh, delimiter="\t"))
selected = {r["experiment_acc"] for r in source_rows if r["study"] == "DRP005276"}
assert len(selected) == 83 and selected <= set(local)
print("LOCAL_ORIGINAL_EXPERIMENTS", len(local))
print("LOCAL_SELECTED_EXPERIMENTS", len(selected))
print("LOCAL_HISTORICAL_QC_HOLDS", len(set(local) - selected))

base = "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA008/DRA008030/"
files = (
    ("experiment", "DRA008030.experiment.xml", 3000000),
    ("sample", "DRA008030.sample.xml", 3000000),
    ("study", "DRA008030.study.xml", 300000),
    ("submission", "DRA008030.submission.xml", 300000),
)
roots = {}
print("\nPUBLIC_DRA_XML")
for kind, filename, maximum in files:
    try:
        data = fetch(base + filename, maximum)
        roots[kind] = ET.fromstring(data)
        print("XML", kind, "FETCH_PASS", "BYTES", len(data), "SHA256", digest(data))
    except Exception as exc:
        # Print exception class only, never a remote response body.
        print("XML", kind, "FETCH_OR_PARSE_FAILED", type(exc).__name__)
if set(roots) != {x[0] for x in files}:
    print("PUBLIC_DRA_XML_INCOMPLETE_NO_CROSSWALK_ASSERTIONS")
else:
    exps = nodes(roots["experiment"], "EXPERIMENT")
    samples = nodes(roots["sample"], "SAMPLE")
    studies = nodes(roots["study"], "STUDY")
    submissions = nodes(roots["submission"], "SUBMISSION")
    exp_by_acc = {e.get("accession", ""): e for e in exps}
    sample_by_acc = {s.get("accession", ""): s for s in samples}
    assert len(exp_by_acc) == len(exps) and "" not in exp_by_acc
    assert len(sample_by_acc) == len(samples) and "" not in sample_by_acc
    print("PUBLIC_DRA_EXPERIMENTS", len(exps))
    print("PUBLIC_DRA_SAMPLES", len(samples))
    print("PUBLIC_DRA_STUDY_ELEMENTS", len(studies))
    print("PUBLIC_DRA_SUBMISSION_ELEMENTS", len(submissions))
    print("PUBLIC_STUDY_ACCESSION_DRP005276", sum(
        s.get("accession") == "DRP005276" for s in studies))
    print("PUBLIC_SUBMISSION_ACCESSION_DRA008030", sum(
        s.get("accession") == "DRA008030" for s in submissions))
    print("PUBLIC_EXPERIMENT_EXACT_MATCH_LOCAL", len(set(exp_by_acc) & set(local)))
    print("PUBLIC_EXPERIMENT_ONLY", len(set(exp_by_acc) - set(local)))
    print("LOCAL_EXPERIMENT_ONLY", len(set(local) - set(exp_by_acc)))
    print("PUBLIC_EXPERIMENT_EXACT_MATCH_SELECTED", len(set(exp_by_acc) & selected))
    ref_counts = Counter()
    exp_to_sample = {}
    for exp, e in exp_by_acc.items():
        study_refs = nodes(e, "STUDY_REF")
        sample_refs = nodes(e, "SAMPLE_DESCRIPTOR")
        if len(study_refs) == 1:
            ref_counts[study_refs[0].get("accession", "<NO_ACCESSION>")] += 1
        if len(sample_refs) == 1:
            exp_to_sample[exp] = sample_refs[0].get("accession", "")
    print("PUBLIC_EXPERIMENT_DRP005276_STUDY_REFS", ref_counts.get("DRP005276", 0))
    print("PUBLIC_EXPERIMENT_SAMPLE_DESCRIPTOR_PRESENT", sum(bool(v) for v in exp_to_sample.values()))
    matched_exp = sorted(set(local) & set(exp_by_acc))
    print("PUBLIC_EXP_LINKED_SAMPLE_EQUALS_LOCAL_SRA_SAMPLE", sum(
        bool(exp_to_sample.get(x)) and exp_to_sample[x] == local[x].get("sample_acc", "")
        for x in matched_exp))
    print("PUBLIC_EXP_LINKED_SAMPLE_IN_PUBLIC_SAMPLE_SET", sum(
        exp_to_sample.get(x) in sample_by_acc for x in matched_exp))
    print("LOCAL_SAMPLE_ACCESSION_IN_PUBLIC_SAMPLE_SET", sum(
        bool(local[x].get("sample_acc")) and local[x]["sample_acc"] in sample_by_acc
        for x in matched_exp))
    print("PUBLIC_SAMPLE_ALIAS_EQUALS_LOCAL_SAMPLE_NAME", sum(
        bool(exp_to_sample.get(x)) and exp_to_sample[x] in sample_by_acc
        and normal(sample_by_acc[exp_to_sample[x]].get("alias", ""))
        == normal(local[x].get("sample_name", ""))
        for x in matched_exp))
    public_keys = Counter()
    for sample in samples:
        for entry in nodes(sample, "SAMPLE_ATTRIBUTE"):
            keys = nodes(entry, "TAG")
            if len(keys) == 1 and keys[0].text:
                public_keys[normal(keys[0].text)] += 1
    for key in ("tissue", "disease", "isolate", "patient", "patient id", "subject id", "treatment"):
        print("PUBLIC_SAMPLE_ATTRIBUTE", key.replace(" ", "_"), "COUNT", public_keys.get(key, 0))
    if len(exp_by_acc) == 103 and set(exp_by_acc) == set(local):
        print("PUBLIC_EXPERIMENT_ACCESSION_SET", "EXACT_103_PASS")
    else:
        print("PUBLIC_EXPERIMENT_ACCESSION_SET", "NOT_EXACT_REQUIRES_REVIEW")

print("\nPLOS_S1_SUPPLEMENTARY_CHECK")
s1_url = ("https://journals.plos.org/plosone/article/file"
          "?id=10.1371/journal.pone.0220173.s001&type=supplementary")
try:
    payload = fetch(s1_url, 20000000)
    assert zipfile.is_zipfile(io.BytesIO(payload)), "NOT_XLSX_RESPONSE"
    with zipfile.ZipFile(io.BytesIO(payload)) as zf:
        names = set(zf.namelist())
        assert "xl/workbook.xml" in names, "NO_WORKBOOK_XML"
        sheet_names = sorted(p for p in names if re.match(r"xl/worksheets/sheet\d+\.xml$", p))
        assert sheet_names, "NO_SHEET_XML"
        shared = []
        if "xl/sharedStrings.xml" in names:
            ss = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in nodes(ss, "si"):
                shared.append("".join(t.text or "" for t in nodes(si, "t")))
        seen_values = set()
        row_counts = []
        for name in sheet_names:
            sheet = ET.fromstring(zf.read(name))
            rows = nodes(sheet, "row")
            row_counts.append(len(rows))
            for cell in nodes(sheet, "c"):
                vnodes = nodes(cell, "v")
                inline = nodes(cell, "is")
                if cell.get("t") == "s" and len(vnodes) == 1 and vnodes[0].text:
                    idx = int(vnodes[0].text)
                    value = shared[idx] if 0 <= idx < len(shared) else ""
                elif inline:
                    value = "".join(t.text or "" for t in nodes(inline[0], "t"))
                elif len(vnodes) == 1:
                    value = vnodes[0].text or ""
                else:
                    value = ""
                if value:
                    seen_values.add(normal(value))
        print("PLOS_S1_FETCH_AND_XLSX_PARSE PASS")
        print("PLOS_S1_BYTES", len(payload))
        print("PLOS_S1_SHA256", digest(payload))
        print("PLOS_S1_SHEET_COUNT", len(sheet_names))
        print("PLOS_S1_ROWS_PER_SHEET", row_counts)
        print("PLOS_S1_DISTINCT_NONEMPTY_CELL_VALUES", len(seen_values))
        for key, values in (
            ("LOCAL_SAMPLE_NAME", [r.get("sample_name", "") for r in local.values()]),
            ("LOCAL_SAMPLE_ACCESSION", [r.get("sample_acc", "") for r in local.values()]),
            ("LOCAL_SOURCE_ISOLATE", [next(iter(a)) if len(a) == 1 else ""
                for r in local.values() for a in [attr_values(r.get("sample_attributes"), "isolate")]]),
        ):
            tokens = {normal(v) for v in values if normal(v)}
            print("PLOS_S1_EXACT_VALUES_MATCH", key, len(tokens & seen_values))
except Exception as exc:
    print("PLOS_S1_UNAVAILABLE_OR_UNPARSEABLE", type(exc).__name__)
    print("PLOS_S1_JOIN", "NOT_ASSESSED")

print("\nSCIENTIFIC_GATES")
print("PUBLIC_DRA_LINK", "REVIEW_PRINTED_COUNTS")
print("PUBLISHED_PATIENT_TO_EXPERIMENT_IDENTITY", "UNSIGNED")
print("SELECTED_83_NMIBC_MIBC_MAPPING", "UNSIGNED")
print("PRIMARY_TUMOUR_COHORT_SOURCE", "PUBLICATION_LEVEL_SUPPORTED")
print("PRETREATMENT_HISTORY", "UNSIGNED")
print("CROSS_STUDY_PATIENT_INDEPENDENCE", "UNSIGNED")
print("CLINICAL_ADMISSIONS", 0)
print("NO_PATIENT_SAMPLE_OR_RUN_VALUES_PRINTED")
print("NO_FROZEN_DATA_MODIFIED")
print("V6_DRP005276_PUBLIC_DRA_PUBLICATION_CROSSWALK_V1_COMPLETE")
