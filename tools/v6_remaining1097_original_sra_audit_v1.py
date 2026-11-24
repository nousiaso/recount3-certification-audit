#!/usr/bin/env python3
"""One-pass original-SRA evidence audit: 1097 experiments, 28 studies.
Frozen inputs read-only. New diagnostics only. No clinical admission.
"""
import csv, gzip, hashlib, json, sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
rec = root / "26_recount3_external_primary_v1"
out = repo / "diagnostics/20261124_v6_remaining1097_source_v1"

v6 = repo / (
    "diagnostics/20261123_v6_evidence_precedence_dryrun/"
    "proposed_evidence_2055_UNSIGNED.tsv"
)
reviewed = repo / (
    "diagnostics/20261124_v6_599_adjudication_queue/"
    "599_clinical_review_queue_LOCAL_ONLY.tsv"
)
v5 = rec / (
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
)

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

assert sha(v6) == (
    "08efe2c8f469edcbdfb98231dc87f705d801434e78dcbea85f57c4b94da88eb0"
)
assert sha(reviewed) == (
    "b5fdf922e4791b4b30c2ffa798c6ffed99b43478a0681c02a1faa2aaa248f3bd"
)
assert sha(v5) == (
    "06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49"
)

def read(p):
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

all_v6 = read(v6)
all_v5 = read(v5)
prior = read(reviewed)

assert len(all_v6) == len(all_v5) == 2055
assert len(prior) == 599

prior_ids = {r["experiment_acc"] for r in prior}
assert len(prior_ids) == 599

historical = [
    r for r in all_v6
    if r["unsigned_evidence_bucket"] ==
       "HISTORICAL_EVIDENCE_REVIEW_PENDING"
]
assert len(historical) == 1696
assert prior_ids <= {r["library_id"] for r in historical}

remaining = [
    r for r in historical
    if r["library_id"] not in prior_ids
]
assert len(remaining) == 1097
assert len({r["library_id"] for r in remaining}) == 1097
assert len({r["study"] for r in remaining}) == 28

assert Counter(r["historical_tissue_grade"] for r in remaining) == {
    "PRIMARY_SURGERY_CONTEXT_NOT_PATHOLOGY_CERTIFIED": 603,
    "PRIMARY_TUMOUR_STUDY_SOURCE_SUPPORTED": 494,
}
assert Counter(r["treatment_grade_UNREVIEWED"] for r in remaining) == {
    "AT_SURGERY_NO_PRETREATMENT_CLAIM": 674,
    "COHORT_PRETREATMENT_SOURCE_SUPPORTED": 334,
    "PRETREATMENT_NOT_ESTABLISHED": 75,
    "GEO_BASELINE_TEXT_SOURCE_SUPPORTED": 14,
}
assert Counter(r["patient_grade_UNREVIEWED"] for r in remaining) == {
    "CURATED_STUDY_PATIENT_KEY_UNVERIFIED": 457,
    "SOURCE_CROSS_ASSAY_PATIENT_TOKEN_PROPOSAL": 364,
    "BIOSAMPLE_PROXY_NOT_PATIENT_ID": 221,
    "STUDY_UNIQUE_PATIENT_COHORT_WITH_RNA_ISOLATE_ONLY": 55,
}

v5_by_id = {r["library_id"]: r for r in all_v5}
assert len(v5_by_id) == 2055

by_study = defaultdict(list)
for r in remaining:
    assert r["FINAL_APPROVAL"] == "NO_UNSIGNED"
    assert v5_by_id[r["library_id"]]["study"] == r["study"]
    by_study[r["study"]].append(r)

GROUPS = {
    "specimen": {
        "tissue", "organ", "organism part", "site",
        "sample type", "sample_type", "tissue type",
        "tumor location", "tumour location",
        "tumor in specimen", "tumour in specimen",
    },
    "disease": {
        "disease", "diagnosis", "disease state", "disease status",
        "tumor type", "tumour type", "stage", "disease staging",
        "histological tumor classification", "who histotype",
    },
    "treatment": {
        "treatment", "treatment status", "prior treatment",
        "first-line therapy", "therapy", "time point",
        "collection timepoint", "subject status",
    },
    "patient": {
        "patient", "patient id", "subject id",
        "individual", "sample donor",
    },
}

def attributes(text):
    found = set()
    for item in text.split("|"):
        if ";;" not in item:
            continue
        k, v = item.split(";;", 1)
        k = " ".join(k.lower().split())
        v = v.strip()
        if k and v:
            found.add((k, v))
    return found

records = []
studies = []
values = Counter()
grade_cross = Counter()

for study in sorted(by_study):
    target = {r["library_id"] for r in by_study[study]}
    path = rec / (
        f"raw_studies/{study}/metadata/sra.sra.{study}.MD.gz"
    )

    original = {}
    row_count = Counter()
    conflicts = set()
    source_status = "READABLE"
    source_hash = ""

    if not path.is_file():
        source_status = "SOURCE_FILE_MISSING"
    else:
        source_hash = sha(path)
        try:
            with gzip.open(
                path, "rt", encoding="utf-8",
                errors="replace", newline=""
            ) as f:
                reader = csv.DictReader(f, delimiter="\t")
                required = {
                    "experiment_acc", "sample_title",
                    "sample_attributes"
                }
                missing = required - set(reader.fieldnames or [])
                if missing:
                    source_status = "SOURCE_SCHEMA_MISSING_FIELDS"
                else:
                    for raw in reader:
                        acc = raw["experiment_acc"]
                        if acc not in target:
                            continue

                        row_count[acc] += 1
                        specimen_fields = (
                            "sample_acc", "biosample",
                            "sample_title", "sample_attributes"
                        )

                        if acc not in original:
                            original[acc] = {
                                k: raw.get(k, "")
                                for k in specimen_fields
                            }
                        elif any(
                            original[acc][k] != raw.get(k, "")
                            for k in specimen_fields
                        ):
                            conflicts.add(acc)
        except (OSError, EOFError, UnicodeError, csv.Error):
            source_status = "SOURCE_READ_ERROR"
            original.clear()
            row_count.clear()
            conflicts.clear()

    coverage = Counter()
    grade_counts = Counter()
    status_counts = Counter()

    for r in sorted(
        by_study[study], key=lambda r: r["library_id"]
    ):
        acc = r["library_id"]
        raw = original.get(acc, {})
        old = v5_by_id[acc]

        if source_status != "READABLE":
            status = source_status
        elif acc not in original:
            status = "EXPERIMENT_ABSENT_FROM_SOURCE"
        elif acc in conflicts:
            status = "CONFLICTING_ORIGINAL_SPECIMEN_ROWS"
        elif (
            raw.get("sample_acc") and
            old.get("source_SRA_sample_acc") and
            raw["sample_acc"] != old["source_SRA_sample_acc"]
        ):
            status = "SRA_V5_SAMPLE_ACCESSION_CONFLICT"
        else:
            status = "EXACT_EXPERIMENT_SPECIMEN_MATCH"

        status_counts[status] += 1

        parsed = attributes(raw.get("sample_attributes", ""))
        keys = {k for k, _ in parsed}

        if not parsed:
            coverage["n_without_parsed_attributes"] += 1

        for group, groupkeys in GROUPS.items():
            coverage["n_with_" + group + "_key"] += int(
                bool(keys & groupkeys)
            )

        for key, value in parsed:
            values[(study, key, value)] += 1

        grade_cross[(
            r["historical_tissue_grade"],
            r["treatment_grade_UNREVIEWED"],
            r["patient_grade_UNREVIEWED"]
        )] += 1

        grade_counts[r["historical_tissue_grade"]] += 1

        records.append({
            "study": study,
            "experiment_acc": acc,
            "source_status": status,
            "source_sample_acc_LOCAL_ONLY":
                raw.get("sample_acc", ""),
            "v5_sample_acc_LOCAL_ONLY":
                old.get("source_SRA_sample_acc", ""),
            "source_biosample_LOCAL_ONLY":
                raw.get("biosample", ""),
            "source_sample_title_LOCAL_ONLY":
                raw.get("sample_title", ""),
            "source_sample_attributes_LOCAL_ONLY":
                raw.get("sample_attributes", ""),
            "n_original_metadata_rows": row_count[acc],
            "historical_tissue_grade":
                r["historical_tissue_grade"],
            "historical_treatment_grade":
                r["treatment_grade_UNREVIEWED"],
            "historical_patient_grade":
                r["patient_grade_UNREVIEWED"],
            "stage41_tissue_source_evidence_LOCAL_ONLY":
                old.get("stage41_tissue_source_evidence", ""),
            "source_file_sha256": source_hash,
            "clinical_approval": "NO_UNSIGNED",
        })

    studies.append({
        "study": study,
        "n_experiments": len(target),
        "n_exact_source_matches":
            status_counts["EXACT_EXPERIMENT_SPECIMEN_MATCH"],
        "n_source_holds":
            len(target) -
            status_counts["EXACT_EXPERIMENT_SPECIMEN_MATCH"],
        "n_multirow_experiments": sum(
            row_count[acc] > 1 for acc in target
        ),
        "n_extra_metadata_rows": sum(
            max(0, row_count[acc] - 1) for acc in target
        ),
        "n_without_parsed_attributes":
            coverage["n_without_parsed_attributes"],
        **{
            "n_with_" + group + "_key":
                coverage["n_with_" + group + "_key"]
            for group in GROUPS
        },
        "n_surgery_context": grade_counts[
            "PRIMARY_SURGERY_CONTEXT_NOT_PATHOLOGY_CERTIFIED"
        ],
        "n_study_source_primary": grade_counts[
            "PRIMARY_TUMOUR_STUDY_SOURCE_SUPPORTED"
        ],
        "source_file_status": source_status,
        "source_file_sha256": source_hash,
        "clinical_status": "NO_UNSIGNED",
    })

assert len(records) == 1097
assert len(studies) == 28

n_exact = sum(x["n_exact_source_matches"] for x in studies)
n_holds = 1097 - n_exact

status = (
    "SOURCE_JOIN_1097_PASS_NOT_CERTIFIED"
    if n_holds == 0 else
    "SOURCE_JOIN_HOLDS_REQUIRE_INSPECTION"
)

def write_tsv(path, rows, columns=None):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=columns or list(rows[0]),
            delimiter="\t", lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)

out.mkdir(parents=True, exist_ok=True)
local = out / "1097_original_source_evidence_LOCAL_ONLY.tsv"
values_file = out / "1097_attribute_values_LOCAL_ONLY.tsv"
study_file = out / "28_study_source_summary_UNSIGNED.tsv"
cross_file = out / "historical_grade_crosstab_UNSIGNED.tsv"
summary_file = out / "summary.json"

assert not any(
    p.exists() for p in (
        local, values_file, study_file,
        cross_file, summary_file
    )
), "Outputs already exist; inspect before rerunning"

write_tsv(local, records)
write_tsv(study_file, studies)

write_tsv(values_file, [
    {
        "study": s, "attribute_key": k,
        "attribute_value_LOCAL_ONLY": v,
        "n_experiments": n
    }
    for (s, k, v), n in sorted(values.items())
], ["study", "attribute_key",
    "attribute_value_LOCAL_ONLY", "n_experiments"])

write_tsv(cross_file, [
    {
        "tissue_grade": t,
        "treatment_grade": tr,
        "patient_grade": p,
        "n_experiments": n
    }
    for (t, tr, p), n in sorted(grade_cross.items())
])

summary = {
    "status": status,
    "n_experiments": 1097,
    "n_studies": 28,
    "n_exact_source_matches": n_exact,
    "n_source_holds": n_holds,
    "study_counts": {
        r["study"]: r["n_experiments"] for r in studies
    },
    "attribute_coverage_totals": {
        key: sum(r[key] for r in studies)
        for key in (
            "n_with_specimen_key", "n_with_disease_key",
            "n_with_treatment_key", "n_with_patient_key",
            "n_without_parsed_attributes",
        )
    },
    "input_hashes": {
        "V6_2055": sha(v6),
        "prior_599": sha(reviewed),
        "V5_2055": sha(v5),
    },
    "outputs_sha256": {
        "local_source": sha(local),
        "local_values": sha(values_file),
        "study_summary": sha(study_file),
        "grade_cross": sha(cross_file),
    },
    "clinical_approvals": 0,
}
summary_file.write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)

print("V6_REMAINING_1097_ORIGINAL_SRA_AUDIT")
print("Experiments:", len(records))
print("Studies:", len(studies))
print("Exact original-SRA matches:", n_exact)
print("Source holds:", n_holds)
print("Attribute coverage:", summary["attribute_coverage_totals"])
print()

print("STUDY_SOURCE_QC")
print("study\tN\texact\tholds\tmultirow\tsource_status")
for r in sorted(
    studies, key=lambda r: -r["n_experiments"]
):
    print(
        r["study"], r["n_experiments"],
        r["n_exact_source_matches"], r["n_source_holds"],
        r["n_multirow_experiments"],
        r["source_file_status"], sep="\t"
    )

print("\nLocal source SHA256:", sha(local))
print("Local values SHA256:", sha(values_file))
print("28-study summary SHA256:", sha(study_file))
print("Grade cross-tab SHA256:", sha(cross_file))
print("Aggregate summary SHA256:", sha(summary_file))
print("STATUS:", status)
print("CLINICAL APPROVALS: 0")
