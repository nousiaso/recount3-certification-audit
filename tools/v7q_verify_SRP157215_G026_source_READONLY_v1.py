#!/usr/bin/env python3
"""Validate SRP157215 G026 source before accepting it; no expression extraction."""
import csv
import gzip
import hashlib
import re
import sys
from pathlib import Path

root, raw = map(Path, sys.argv[1:])
s5 = root / "26_recount3_external_primary_v1/stage5_output"
def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))
bridge = s5 / "03_bridge_165_availability.tsv"
mapping = s5 / "00_gene_row_mapping.tsv"
assert sha(bridge) == "c7d6d7c754e3ca77d62424f26627a40ceaec10fac98a0741c945773f3d3b1b67"
assert sha(mapping) == "68be266dd6bd30bfa79a2bd8ac79f0e49c6f5d0f19d2d083e2ca57850ef17391"
b = read(bridge)
m = read(mapping)
assert len(b) == 165 and len(m) == 63856
assert [int(r["row_zero_based"]) for r in m] == list(range(len(m)))
pairs = [r for r in b if r["study"] == "SRP157215"]
assert len(pairs) == 22
assert {r["bridge_class"] for r in pairs} == {"CALIBRATION_ONLY_NOT_IN_PRIMARY"}
for k in ("library_id", "ARCHS4_GSM", "biosample", "ARCHS4_expression_column_index"):
    assert len({r[k] for r in pairs if r[k]}) == 22
run = re.compile(r"(?:SRR|ERR|DRR)\d+", re.I)
wanted = set()
for r in pairs:
    ids = run.findall(r["recount3_run_ids"])
    assert len(ids) == int(r["n_runs"]) == 1
    wanted.add(ids[0].upper())
assert len(wanted) == 22
with gzip.open(raw, "rt", encoding="utf-8-sig", newline="") as f:
    header = None
    for i in range(20):
        line = f.readline()
        if not line:
            break
        cells = line.rstrip("\r\n").split("\t")
        ids = [x.strip().strip('"').upper() for x in cells[1:]]
        if ids and all(run.fullmatch(x) for x in ids):
            header = ids
            print("HEADER_LINE_ZERO_BASED", i)
            break
    assert header is not None, "No G026 run header"
    assert len(header) == len(set(header)), "Duplicate run axis"
    assert wanted <= set(header), "Missing matched runs"
    print("RAW_RUN_COLUMNS", len(header))
    print("MATCHED_RUNS_PRESENT", len(wanted & set(header)))
    n = 0
    for line in f:
        assert n < len(m), "Extra gene row"
        cells = line.rstrip("\r\n").split("\t")
        assert len(cells) == len(header) + 1, f"Width mismatch at row {n}"
        gene = cells[0].strip().strip('"')
        assert gene == m[n]["source_gene_id"] or gene.split(".", 1)[0] == m[n]["gene_base"], f"Gene-axis mismatch at row {n}"
        n += 1
    assert n == len(m), "Missing gene rows"
print("FULL_GENE_AXIS_MATCHED_ROWS", n)
print("RAW_BYTES", raw.stat().st_size)
print("RAW_SHA256", sha(raw))
print("V7Q_SRP157215_RAW_SOURCE_AXIS_PASS")
print("NO_EXPRESSION_MATRIX_EXTRACTED_NO_MODEL_EXECUTION_CLINICAL_ADMISSIONS_0")
