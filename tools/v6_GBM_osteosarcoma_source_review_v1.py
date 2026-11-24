#!/usr/bin/env python3
"""Read-only evidence interpretation for 20 GBM and 22 osteosarcoma experiments.
No clinical admission, specimen exclusion or official run verification.
"""
import csv, hashlib, json, re, sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
source = repo / (
    "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
    "source_evidence_599_LOCAL_ONLY.tsv"
)
out = repo / "diagnostics/20261124_v6_GBM_osteosarcoma_review"
expected = (
    "8b6f017dd575ec6fc3917b784740d37c8556fe567"
    "474e6a5e69dcad066ef35eb"
)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

assert sha(source) == expected

with source.open(newline="", encoding="utf-8-sig") as f:
    all_rows = list(csv.DictReader(f, delimiter="\t"))

assert len(all_rows) == 599

targets = [
    r for r in all_rows
    if r["study"] in {"SRP072494", "SRP090849"}
]
assert Counter(r["study"] for r in targets) == {
    "SRP072494": 20,
    "SRP090849": 22,
}
assert len({r["experiment_acc"] for r in targets}) == 42

def parse_attributes(raw):
    parts = defaultdict(set)
    for item in raw.split("|"):
        if ";;" in item:
            key, value = item.split(";;", 1)
            key = " ".join(key.lower().split())
            parts[key].add(value.strip())
    return parts

def unique_value(attributes, key):
    values = {x for x in attributes.get(key, set()) if x}
    assert len(values) <= 1, (
        "Conflicting attribute values", key, values
    )
    return next(iter(values)) if values else ""

records = []
subject_groups = defaultdict(list)

for r in targets:
    study = r["study"]
    accession = r["experiment_acc"]
    a = parse_attributes(r["original_sample_attributes"])

    primary_diagnosis = ""
    tumour_type = ""
    timepoint = ""
    therapy = ""
    subject = ""

    if study == "SRP072494":
        primary_diagnosis = unique_value(
            a, "sample from primary gbm diagnosis"
        ).lower()
        timepoint = unique_value(a, "time point")
        therapy = unique_value(a, "first-line therapy")
        subject = unique_value(a, "subject id")

        assert primary_diagnosis in {"yes", "no"}
        assert timepoint == (
            "before bevacizumab combination treatment"
        )
        assert subject
        subject_groups[subject].append(accession)

        evidence_flag = (
            "PRIMARY_DIAGNOSIS_SOURCE_YES"
            if primary_diagnosis == "yes"
            else "NOT_PRIMARY_DIAGNOSIS_SOURCE_NO"
        )
    else:
        tumour_type = unique_value(
            a, "tumor type"
        ).lower()
        assert tumour_type in {"primary", "unknown"}

        evidence_flag = (
            "PRIMARY_TUMOUR_TYPE_SOURCE_LABEL"
            if tumour_type == "primary"
            else "TUMOUR_TYPE_UNKNOWN"
        )

    title = r["original_sample_title"]
    source_name = unique_value(a, "source_name")
    combined = " ".join([
        title, source_name, unique_value(a, "tissue")
    ])

    potential_culture_keyword = bool(re.search(
        r"(?i)(cell[\s_-]?line|cultured|culture|atcc|immortalized)",
        combined
    ))

    records.append({
        "study": study,
        "experiment_acc": accession,
        "sra_sample_acc": r["sra_sample_acc"],
        "evidence_flag_NOT_APPROVAL": evidence_flag,
        "primary_diagnosis_value": primary_diagnosis,
        "tumour_type_value": tumour_type,
        "time_point_value": timepoint,
        "first_line_therapy_value": therapy,
        "subject_token_LOCAL_ONLY": subject,
        "sample_title_LOCAL_ONLY": title,
        "source_name_LOCAL_ONLY": source_name,
        "possible_culture_keyword_NOT_PROOF":
            str(potential_culture_keyword),
        "clinical_approval": "NO_UNSIGNED",
    })

counts = Counter(
    (r["study"], r["evidence_flag_NOT_APPROVAL"])
    for r in records
)
assert counts == {
    ("SRP072494", "PRIMARY_DIAGNOSIS_SOURCE_YES"): 10,
    ("SRP072494", "NOT_PRIMARY_DIAGNOSIS_SOURCE_NO"): 10,
    ("SRP090849", "PRIMARY_TUMOUR_TYPE_SOURCE_LABEL"): 17,
    ("SRP090849", "TUMOUR_TYPE_UNKNOWN"): 5,
}

subject_multiplicity = Counter(
    len(v) for v in subject_groups.values()
)

out.mkdir(parents=True, exist_ok=True)

local = out / "42_accession_source_evidence_LOCAL_ONLY.tsv"
with local.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(records[0]),
        delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(
        sorted(records, key=lambda r: (
            r["study"], r["experiment_acc"]
        ))
    )

summary = {
    "status": "SOURCE_EVIDENCE_ONLY_NO_APPROVAL",
    "source_599_sha256": expected,
    "n_experiments": 42,
    "source_label_counts": {
        f"{study}:{flag}": n
        for (study, flag), n in sorted(counts.items())
    },
    "gbm_distinct_submitted_subject_tokens_NOT_VERIFIED":
        len(subject_groups),
    "gbm_subject_token_multiplicity":
        dict(sorted(subject_multiplicity.items())),
    "osteosarcoma_possible_culture_keyword_hits_NOT_PROOF":
        sum(
            r["possible_culture_keyword_NOT_PROOF"] == "True"
            for r in records if r["study"] == "SRP090849"
        ),
    "local_table_sha256": sha(local),
    "limitations": [
        "Primary-diagnosis label is not treatment-naive certification",
        "Primary tumour-type label is not cell-line exclusion proof",
        "Submitted subject IDs are not verified patient identities",
        "GSE79671 recurrent GBM treatment context requires review",
        "GSE87624 includes human tumours and cell lines",
        "No official run QC or clinical admission was performed",
    ],
}

report = out / "summary.json"
report.write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)

print("TWO_STUDY_SOURCE_REVIEW_PASS")
print("42/42 existing source records examined")
for (study, flag), n in sorted(counts.items()):
    print(study, flag, n, sep="\t")

print(
    "GBM distinct submitted subject tokens "
    "(NOT certified patients):",
    len(subject_groups)
)
print(
    "GBM samples per submitted subject token:",
    dict(sorted(subject_multiplicity.items()))
)
print(
    "Osteosarcoma culture-keyword hits "
    "(NOT a cell-line exclusion test):",
    summary[
        "osteosarcoma_possible_culture_keyword_hits_NOT_PROOF"
    ]
)
print("LOCAL_DETAIL_SHA256:", sha(local))
print("SUMMARY_SHA256:", sha(report))
print("NO CLINICAL ADMISSIONS")
