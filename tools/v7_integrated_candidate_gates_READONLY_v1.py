#!/usr/bin/env python3
"""V7 first full-population gate preflight (read-only inputs, local-only outputs).

Historical grades are evidence labels, NEVER automatic clinical certification.
Needs the already frozen V5 TSV, V6 evidence-precedence TSV, and Stage41 TSV.
No SpatialEcoTyper projection; no modifications to original datasets.
"""
import argparse
import csv
import hashlib
import io
import json
import os
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

EXPECTED = {
    "v5": "06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49",
    "v6": "08efe2c8f469edcbdfb98231dc87f705d801434e78dcbea85f57c4b94da88eb0",
    "stage41": "aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81",
}
BUCKETS = {
    "EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE": 22,
    "PRIMARY_TUMOUR_SPECIMEN_ONLY_OTHER_GATES_OPEN": 5,
    "HOLD_SAMPLED_SITE_METASTASIS_DISEASE_FLAG": 52,
    "V5_TUMOUR_OR_CONTEXT_SPECIMEN_GATES_OPEN": 280,
    "HISTORICAL_EVIDENCE_REVIEW_PENDING": 1696,
}

def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()

def read_bytes(path):
    with open(path, "rb") as f:
        return f.read()

def read_v5(rec):
    name = "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"
    plain = rec / "clinical_specimen_patient_v5_unsigned_20261122" / name
    archive = rec / "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip"
    if plain.is_file():
        return read_bytes(plain), str(plain)
    if archive.is_file():
        with zipfile.ZipFile(archive) as z:
            names = [n for n in z.namelist() if n.endswith("/results/" + name)]
            if len(names) != 1:
                raise ValueError("V5 archive member ambiguous/missing: " + repr(names))
            return z.read(names[0]), str(archive) + "::" + names[0]
    raise FileNotFoundError("Frozen V5 primary-screen TSV or archive missing")

def table(raw, label):
    rd = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    columns = rd.fieldnames or []
    if not columns:
        raise ValueError(label + ": no columns")
    rows = list(rd)
    if any(None in x for x in rows):
        raise ValueError(label + ": malformed TSV row")
    return columns, rows

def require(columns, names, label):
    missing = set(names) - set(columns)
    if missing:
        raise ValueError(label + " missing " + ",".join(sorted(missing)))

def exact_index(rows, key, label):
    ids = [r[key] for r in rows]
    if not all(ids) or len(ids) != len(set(ids)):
        raise ValueError(label + ": missing or duplicate " + key)
    return dict(zip(ids, rows))

def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path, help="ARCHS4_V2_6 root")
    ap.add_argument("--output", type=Path, required=True, help="new LOCAL-ONLY output folder, ideally repo/private_data/")
    args = ap.parse_args()
    root = args.root.resolve()
    rec = root / "26_recount3_external_primary_v1"
    repo = root / "RECOUNT3_CERTIFICATION_GIT"
    src_v6 = repo / "diagnostics/20261123_v6_evidence_precedence_dryrun/proposed_evidence_2055_UNSIGNED.tsv"
    src_s41 = root / "41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/03_candidate_integrated_evidence_3602.tsv"
    raw = {"v5": read_v5(rec)[0], "v6": read_bytes(src_v6), "stage41": read_bytes(src_s41)}
    source_names = {"v5": read_v5(rec)[1], "v6": str(src_v6), "stage41": str(src_s41)}
    for label, data in raw.items():
        observed = sha_bytes(data)
        if observed != EXPECTED[label]:
            raise ValueError("FROZEN_INPUT_SHA_MISMATCH " + label + " " + observed)

    cv5, v5 = table(raw["v5"], "V5")
    cv6, v6 = table(raw["v6"], "V6")
    cs41, stage = table(raw["stage41"], "Stage41")
    require(cv5, ["library_id", "study", "FINAL_CLINICAL_APPROVAL", "tissue_evidence_grade", "treatment_evidence_grade", "patient_evidence_grade", "patient_key", "biosample"], "V5")
    require(cv6, ["library_id", "study", "unsigned_evidence_bucket", "FINAL_APPROVAL"], "V6")
    require(cs41, ["library_id", "study", "technical_hold", "cohort_approved", "primary_tumour_certified", "pretreatment_certified", "patient_independence_certified", "atlas_independence_certified"], "Stage41")
    if (len(v5), len(v6), len(stage)) != (2055, 2055, 3602):
        raise ValueError("FROZEN_POPULATION_SIZE_MISMATCH")
    by5 = exact_index(v5, "library_id", "V5")
    by6 = exact_index(v6, "library_id", "V6")
    by41 = exact_index(stage, "library_id", "Stage41")
    if set(by5) != set(by6) or not set(by5) <= set(by41):
        raise ValueError("FROZEN_EXPERIMENT_MEMBERSHIP_MISMATCH")
    if len({r["study"] for r in v5}) != 77:
        raise ValueError("EXPECTED_77_STUDIES")
    actual_buckets = Counter(r["unsigned_evidence_bucket"] for r in v6)
    if actual_buckets != BUCKETS:
        raise ValueError("V6_BUCKETS_CHANGED " + repr(actual_buckets))

    output = []
    bystudy = defaultdict(list)
    patient_keys = defaultdict(list)
    for acc, old in sorted(by5.items(), key=lambda kv: (kv[1]["study"], kv[0])):
        evidence, st = by6[acc], by41[acc]
        if old["study"] != evidence["study"] or old["study"] != st["study"]:
            raise ValueError("CROSS_SOURCE_STUDY_CONFLICT")
        if old["FINAL_CLINICAL_APPROVAL"] != "NO_UNSIGNED" or evidence["FINAL_APPROVAL"] != "NO_UNSIGNED":
            raise ValueError("UNEXPECTED_HISTORICAL_CLINICAL_RELEASE")
        if st["technical_hold"] != "NO":
            raise ValueError("V5_PRIMARY_WITH_STAGE41_TECHNICAL_HOLD")
        for name in ("cohort_approved", "primary_tumour_certified", "pretreatment_certified", "patient_independence_certified", "atlas_independence_certified"):
            if st[name] != "NO":
                raise ValueError("UNEXPECTED_STAGE41_CERTIFICATION " + name)
        bucket = evidence["unsigned_evidence_bucket"]
        excluded = bucket == "EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE"
        row = {
            "study": old["study"],
            "experiment_acc": acc,
            "sra_sample_acc_LOCAL_ONLY": old.get("source_SRA_sample_acc", ""),
            "biosample_LOCAL_ONLY": old["biosample"],
            "administrative_patient_key_LOCAL_ONLY": old["patient_key"],
            "technical_screen": "STAGE41_NONHELD_NOT_FULL_ASSAY_CERTIFICATION",
            "specimen_evidence_bucket": bucket,
            "historical_tissue_grade": old["tissue_evidence_grade"],
            "historical_treatment_grade": old["treatment_evidence_grade"],
            "historical_patient_grade": old["patient_evidence_grade"],
            "sampled_site_gate": "SOURCE_NON_TUMOUR_EXCLUDE" if excluded else (
                "EXPLICIT_PRIMARY_SPECIMEN_OTHER_GATES_OPEN" if bucket == "PRIMARY_TUMOUR_SPECIMEN_ONLY_OTHER_GATES_OPEN" else "NOT_CERTIFIED"),
            "pretreatment_gate": "NOT_CERTIFIED",
            "patient_independence_gate": "NOT_CERTIFIED",
            "atlas_independence_gate": "NOT_CERTIFIED",
            "expression_assay_compatibility_gate": "NOT_TESTED_IN_THIS_PREFLIGHT",
            "model_transfer_gate": "NOT_TESTED_IN_THIS_PREFLIGHT",
            "disposition": "EXCLUDE_SOURCE_NON_TUMOUR" if excluded else "HOLD_PENDING_GATE_EVIDENCE",
            "final_clinical_approval": "NO_UNSIGNED",
            "source_snapshot": "FROZEN_V5_V6_STAGE41_EXACT_JOIN",
        }
        output.append(row)
        bystudy[old["study"]].append(row)
        k = old["patient_key"].strip()
        if k:
            patient_keys[old["study"]].append(k)

    assert len(output) == 2055 and len(bystudy) == 77
    assert sum(r["disposition"] == "EXCLUDE_SOURCE_NON_TUMOUR" for r in output) == 22
    assert all(r["final_clinical_approval"] == "NO_UNSIGNED" for r in output)

    study_rows = []
    for study, items in bystudy.items():
        buckets = Counter(r["specimen_evidence_bucket"] for r in items)
        grades = Counter(r["historical_treatment_grade"] for r in items)
        keys = patient_keys[study]
        study_rows.append({
            "study": study,
            "n_screen": len(items),
            "n_exclude_source_non_tumour": buckets["EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE"],
            "n_primary_specimen_source_supported": buckets["PRIMARY_TUMOUR_SPECIMEN_ONLY_OTHER_GATES_OPEN"],
            "n_sarcoma_sampled_site_hold": buckets["HOLD_SAMPLED_SITE_METASTASIS_DISEASE_FLAG"],
            "n_context_evidence_other_gates_open": buckets["V5_TUMOUR_OR_CONTEXT_SPECIMEN_GATES_OPEN"],
            "n_historical_only": buckets["HISTORICAL_EVIDENCE_REVIEW_PENDING"],
            "n_historical_cohort_pretreatment_support": grades["COHORT_PRETREATMENT_SOURCE_SUPPORTED"],
            "n_historical_geo_baseline_support": grades["GEO_BASELINE_TEXT_SOURCE_SUPPORTED"],
            "n_administrative_patient_keys": len(keys),
            "n_distinct_administrative_patient_keys": len(set(keys)),
            "n_final_approved": 0,
            "state": "UNSIGNED",
        })
    study_rows.sort(key=lambda r: (-r["n_historical_cohort_pretreatment_support"], -r["n_primary_specimen_source_supported"], -r["n_screen"], r["study"]))

    out = args.output.resolve()
    try:
        out.relative_to(repo)
    except ValueError:
        raise ValueError("OUTPUT_MUST_BE_UNDER_REPO_FOR_PROVENANCE")
    if out.exists():
        raise FileExistsError("REFUSING_TO_OVERWRITE_OUTPUT " + str(out))
    out.mkdir(parents=True)
    os.chmod(out, 0o700)
    individual = out / "2055_gate_matrix_LOCAL_ONLY.tsv"
    studies = out / "77_study_triage_LOCAL_ONLY.tsv"
    write_tsv(individual, list(output[0]), output)
    write_tsv(studies, list(study_rows[0]), study_rows)
    summary = {
        "version": "v7_integrated_candidate_gates_READONLY_v1",
        "input_sha256": EXPECTED,
        "input_paths": source_names,
        "n_original_v5_primary_screen": len(output),
        "n_studies": len(bystudy),
        "evidence_buckets": dict(sorted(actual_buckets.items())),
        "n_excluded_source_non_tumour": 22,
        "n_held_pending_clinical_gates": 2033,
        "n_clinically_certified": 0,
        "n_eligible_for_definitive_spatialecotyper_validation": 0,
        "individual_matrix_sha256": sha_bytes(read_bytes(individual)),
        "study_table_sha256": sha_bytes(read_bytes(studies)),
        "disclaimer": "Historical clinical grades and unique administrative keys are not specimen/patient certification; do not release data or infer eligibility.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("V7_COHORT_WIDE_GATE_PREFLIGHT_PASS")
    print("V5_PRIMARY_EXPERIMENTS", len(output))
    print("V5_STUDIES", len(bystudy))
    for name, count in sorted(actual_buckets.items()):
        print("SOURCE_BUCKET", name, count)
    print("SOURCE_NON_TUMOUR_EXCLUSIONS 22")
    print("CLINICAL_GATES_OPEN 2033")
    print("CLINICAL_ADMISSIONS 0")
    print("TOP_STUDIES_BY_HISTORICAL_PRETREATMENT_SUPPORT_NOT_APPROVAL")
    for r in study_rows[:12]:
        print("STUDY", r["study"], "SCREEN", r["n_screen"], "COHORT_TREATMENT_SUPPORT", r["n_historical_cohort_pretreatment_support"], "EXPLICIT_PRIMARY", r["n_primary_specimen_source_supported"], "NON_TUMOUR", r["n_exclude_source_non_tumour"])
    print("LOCAL_MATRIX_SHA256", summary["individual_matrix_sha256"])
    print("LOCAL_STUDY_TABLE_SHA256", summary["study_table_sha256"])
    print("NO_SOURCE_FILES_MODIFIED")
    print("NO_MODEL_EXECUTION")
    print("V7_COHORT_WIDE_GATE_PREFLIGHT_COMPLETE")

if __name__ == "__main__":
    main()
