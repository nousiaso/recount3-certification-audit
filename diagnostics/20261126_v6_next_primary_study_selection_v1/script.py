import csv
import io
import re
import sys
import zipfile
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path

root = Path(sys.argv[1])
base = root / "26_recount3_external_primary_v1"

v3 = base / "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"
v5 = base / "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"

def readzip(path, suffix):
    with zipfile.ZipFile(path) as z:
        matches = [n for n in z.namelist()
                   if n.endswith("/" + suffix)]
        assert len(matches) == 1
        raw = z.read(matches[0]).decode("utf-8-sig")
    return list(csv.DictReader(io.StringIO(raw), delimiter="\t"))

def number(value):
    try:
        d = Decimal(str(value or "").strip())
        return d if d.is_finite() else None
    except InvalidOperation:
        return None

def safe(value):
    s = str(value or "").strip()
    return s if re.fullmatch(r"[A-Z0-9_/-]{1,90}", s) \
        else "[UNPRINTED]"

queue = readzip(
    v3, "results/06_REMAINING_STUDIES_PRIORITY_QUEUE.tsv"
)
primary = readzip(
    v5, "results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
)

assert len(queue) == 147
assert len(primary) == 2055

counts = Counter(r["study"] for r in primary)

print("RECOUNT3_NEXT_PRIMARY_STUDY_SELECTION_V1")
print("V3_QUEUE_STUDIES", len(queue))
print("V5_PRIMARY_ROWS", len(primary))
print("V5_PRIMARY_STUDIES", len(counts))

candidates = [
    r for r in queue
    if counts.get(r["study"], 0) > 0
]

print("QUEUED_STUDIES_CONTRIBUTING_TO_PRIMARY", len(candidates))
print("PRIMARY_ROWS_IN_QUEUED_STUDIES",
      sum(counts[r["study"]] for r in candidates))

mismatch = 0

for r in candidates:
    reported = number(r.get("qa_stage41_primary_screen"))
    if reported is None or reported != counts[r["study"]]:
        mismatch += 1

print("V3_V5_PRIMARY_COUNT_MISMATCH_STUDIES", mismatch)

def show(label, records):
    print("\n" + label)
    for r in records[:15]:
        study = r["study"]
        print(
            "STUDY", study,
            "PRIMARY2055", counts[study],
            "V3_REVIEW_SCORE",
            str(number(r.get("review_priority_score"))),
            "V3_QA_SCORE",
            str(number(r.get("qa_priority_score"))),
            "MISSING_PRETREATMENT",
            str(number(r.get("qa_missing_pretreatment"))),
            "STRICT_ATLAS_OVERLAP",
            str(number(r.get("qa_strict_atlas_overlap"))),
            "REVIEW_ACTION",
            safe(r.get("qa_action"))
        )

show(
    "TOP_BY_PRIMARY_SCREEN_CONTRIBUTION",
    sorted(
        candidates,
        key=lambda r: (-counts[r["study"]], r["study"])
    )
)

show(
    "TOP_BY_RECORDED_REVIEW_PRIORITY_SCORE",
    sorted(
        [r for r in candidates
         if number(r.get("review_priority_score")) is not None],
        key=lambda r: (
            -number(r.get("review_priority_score")),
            -counts[r["study"]],
            r["study"]
        )
    )
)

print("\nSCORE_ORDER_IS_DESCRIPTIVE_NOT_CLINICAL_APPROVAL")
print("NO_PATIENT_OR_SAMPLE_IDENTIFIERS_PRINTED")
print("NO_COHORT_FILES_MODIFIED")
print("CLINICAL_ADMISSIONS 0")
print("RECOUNT3_NEXT_PRIMARY_STUDY_SELECTION_V1_COMPLETE")
