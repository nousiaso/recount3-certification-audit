import csv
import io
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

base = Path(sys.argv[1]) / "26_recount3_external_primary_v1"

v3 = base / "RECOUNT3_CLINICAL_REVIEW_EVIDENCE_V3_20261121.zip"
v5 = base / "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"

def table(archive, suffix):
    with zipfile.ZipFile(archive) as z:
        paths = [p for p in z.namelist()
                 if p.endswith("/" + suffix)]
        assert len(paths) == 1
        raw = z.read(paths[0]).decode("utf-8-sig")
    return list(csv.DictReader(io.StringIO(raw), delimiter="\t"))

review = table(v3, "results/02_ALL_176_STUDY_SOURCE_REVIEW_UNSIGNED.tsv")
reviewed = table(v3, "results/04_PUBLICATION_REVIEWED_STUDIES.tsv")
queue = table(v3, "results/06_REMAINING_STUDIES_PRIORITY_QUEUE.tsv")
primary = table(v5, "results/03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv")

assert (len(review), len(reviewed), len(queue), len(primary)) == (
    176, 29, 147, 2055
)

rmap = {r["study"]: r for r in review}
qset = {r["study"] for r in queue}
pset = {r["study"] for r in reviewed}
counts = Counter(r["study"] for r in primary)

assert len(rmap) == 176
assert len(qset) == 147
assert len(pset) == 29
assert not (qset & pset)
assert qset | pset == set(rmap)
assert len(counts) == 77
assert set(counts).issubset(rmap)

queued = set(counts) & qset
published = set(counts) & pset

def safe(x):
    x = str(x or "").strip()
    if not x:
        return "EMPTY"
    return x if re.fullmatch(r"[A-Z0-9_./-]{1,100}", x) else "TEXT_SUPPRESSED"

def show(label, studies, n):
    print("\n" + label, "COUNT", len(studies))
    ordered = sorted(studies, key=lambda s: (-counts[s], s))

    for s in ordered[:n]:
        r = rmap[s]
        print(
            "STUDY", s,
            "PRIMARY", counts[s],
            "CANCER", safe(r.get("cancer")),
            "REVIEW_SCORE", safe(r.get("review_priority_score")),
            "PRETREATMENT_MISSING_STUDY",
            safe(r.get("qa_missing_pretreatment")),
            "PUBLIC_REVIEW", safe(r.get("public_source_review_status")),
            "DISPOSITION", safe(r.get("overall_study_disposition")),
            "SOURCE_ACTION", safe(r.get("review_action_from_publication")),
            "FINAL_APPROVAL", safe(r.get("FINAL_STUDY_CLINICAL_APPROVAL_AS_OF_V3"))
        )

print("RECOUNT3_V6_77_STUDY_CLINICAL_TRIAGE_V1")
print("ALL_STUDIES", len(rmap))
print("PUBLICATION_REVIEWED", len(pset))
print("REMAINING_QUEUE", len(qset))
print("PRIMARY_STUDIES", len(counts))
print("PRIMARY_TOTAL", sum(counts.values()))
print("PRIMARY_PUBLICATION_REVIEWED_STUDIES", len(published))
print("PRIMARY_PUBLICATION_REVIEWED_LIBRARIES",
      sum(counts[s] for s in published))
print("PRIMARY_QUEUED_STUDIES", len(queued))
print("PRIMARY_QUEUED_LIBRARIES",
      sum(counts[s] for s in queued))

show("TOP_ALREADY_PUBLICATION_REVIEWED", published, 20)
show("TOP_IN_REMAINING_REVIEW_QUEUE", queued, 15)

for label, studies in (
    ("REVIEWED", published),
    ("QUEUED", queued)
):
    for field in (
        "overall_study_disposition",
        "FINAL_STUDY_CLINICAL_APPROVAL_AS_OF_V3"
    ):
        c = Counter(safe(rmap[s].get(field)) for s in studies)
        print(
            "STATUS_COUNTS", label, field,
            repr(dict(sorted(c.items())))
        )

print("\nKEY_CANDIDATE_CHECK")
for s in ("SRP124458", "SRP002628", "SRP109781"):
    r = rmap[s]
    print(
        "CANDIDATE", s,
        "PRIMARY", counts[s],
        "REVIEW_GROUP", "QUEUE" if s in qset else "REVIEWED",
        "CANCER", safe(r.get("cancer")),
        "DISPOSITION", safe(r.get("overall_study_disposition")),
        "REVIEW_ACTION", safe(r.get("review_action_from_publication"))
    )

print("REVIEW_STATUS_NOT_CLINICAL_CERTIFICATION")
print("NO_PATIENT_OR_SAMPLE_IDENTIFIERS_PRINTED")
print("NO_COHORT_FILES_MODIFIED")
print("CLINICAL_ADMISSIONS 0")
print("RECOUNT3_V6_77_STUDY_CLINICAL_TRIAGE_V1_COMPLETE")
