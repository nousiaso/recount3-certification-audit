#!/usr/bin/env python3
"""Reconstructed from the executed inline diagnostic, 2026-11-27.
Original executed source bytes were not captured; aggregate-only review.
No cohort decisions or file writes.
"""
import csv
import sys
from pathlib import Path
from collections import Counter

root = Path(sys.argv[1])
rec = root / "26_recount3_external_primary_v1"
s41 = root / "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120"

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

study = "SRP057793"
source = [r for r in read(s41 / "01_full_source_experiment_identity_5221.tsv") if r["study"] == study]
candidates = [r for r in read(s41 / "03_candidate_integrated_evidence_3602.tsv") if r["study"] == study]
holds = {r["library_id"] for r in candidates if r["technical_hold"] == "YES"}
retained = {r["library_id"] for r in candidates if r["technical_hold"] == "NO"}
assert len(source) == 149 and len(candidates) == 148
assert len(holds) == 8 and len(retained) == 140 and not holds.intersection(retained)
queue = read(rec / "cohort_readiness_v1/source_snapshot/05_run_hold_queue_409.tsv")
matched = [r for r in queue if r["study"] == study and r["library_id"] in holds]
assert len(matched) == 8 and {r["library_id"] for r in matched} == holds
print("=== EXACT HISTORICAL HOLD REASONS ===")
for field in ("reason_codes","automatic_proposal","decision","triage_priority","proposed_action_not_applied"):
    print(field, dict(Counter(r[field] for r in matched)))
# Unlike the original inline display, retain only aggregate ranges, not eight individual metric vectors.
print("=== HELD EXPERIMENT QC RANGES ===")
for field in ("gene_detection_fraction","top10_coverage_fraction","average_mapped_length_bp","gene_sum_coverage_bp"):
    values = [float(r[field]) for r in matched]
    print(field, "MIN", min(values), "MAX", max(values))
print("=== RETAINED CLINICAL EVIDENCE SCHEMA ===")
rows = [r for r in candidates if r["library_id"] in retained]
for field in rows[0]:
    if any(t in field.lower() for t in ("patient","treat","material","specimen","primary","assay","independ","protocol","overlap","clinical")):
        values = [r[field] for r in rows]
        print(field,"POPULATED",sum(v not in ("","NA","None") for v in values),"DISTINCT",len(set(values)))
print("SRP057793_QC_REASON_AND_CLINICAL_SCHEMA_COMPLETE")
print("NO_FILES_MODIFIED_CLINICAL_ADMISSIONS_0")
