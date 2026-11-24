#!/usr/bin/env python3
"""Append-only, source-grounded clinical evidence addendum.
No clinical admissions or changes to frozen V5/V6 inputs.
"""
import csv, gzip, hashlib, json, re, sys
from pathlib import Path
from collections import Counter, defaultdict

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
folder = repo / "diagnostics/20261124_v6_priority_64_geo_exact"
out = repo / "diagnostics/20261124_v6_priority64_evidence_addendum_v1"

detail = folder / "64_geo_exact_source_evidence_LOCAL_ONLY.tsv"
queue = repo / (
    "diagnostics/20261124_v6_599_adjudication_queue/"
    "599_clinical_review_queue_LOCAL_ONLY.tsv"
)
v5 = root / (
    "26_recount3_external_primary_v1/"
    "clinical_specimen_patient_v5_unsigned_20261122/"
    "01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_CROSSWALK_UNSIGNED.tsv"
)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

expected = {
    detail: "71f788ee1fd34edb3adb8b7fe26dc0cd06067a3a87b31e377b9d5124748fdae5",
    queue: "b5fdf922e4791b4b30c2ffa798c6ffed99b43478a0681c02a1faa2aaa248f3bd",
    v5: "ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a",
}

series = {
    "SRP094795": ("GSE91065", 25, 23,
        "4dd7d20aabbaafb4c2d0b16d541da000787fa880d652e2db0b1aa1b52a7ac730"),
    "SRP072494": ("GSE79671", 36, 20,
        "25d6a138fa91e90b75bbefc375b06cdb34e01b25608e98a38fa6b89b6266eedc"),
    "SRP194626": ("GSE130469", 20, 4,
        "44317d378636086bfa304756520251776709efa1f1a36eacbd9d74d3effc9c58"),
    "SRP112751": ("GSE101588", 17, 17,
        "2974a161614d170a59d5d640192ce808ffb414c8f8ca6d2f25c577866cc1c407"),
}

def read(p):
    with p.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def parse_soft(p):
    samples, current = {}, None
    with gzip.open(p, "rt", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith("^"):
                if current is not None:
                    assert current["gsm"] not in samples
                    samples[current["gsm"]] = current
                current = None
                if line.startswith("^SAMPLE = "):
                    current = {
                        "gsm": line.split(" = ", 1)[1],
                        "title": "",
                        "traits": defaultdict(list),
                        "protocol": [],
                    }
                continue
            if current is None or " = " not in line:
                continue
            key, val = line.split(" = ", 1)
            val = val.strip()
            if key == "!Sample_title":
                current["title"] = val
            elif key == "!Sample_treatment_protocol_ch1":
                current["protocol"].append(val)
            elif key == "!Sample_characteristics_ch1" and ":" in val:
                k, v = val.split(":", 1)
                current["traits"][k.strip().lower()].append(v.strip())
    if current is not None:
        assert current["gsm"] not in samples
        samples[current["gsm"]] = current
    return samples

def trait(s, key):
    vals = s["traits"].get(key, [])
    return vals[0] if len(vals) == 1 else "UNRESOLVED"

def prefix(s):
    m = re.match(r"^P(\d+)[-_]", s["title"])
    return m.group(1) if m else None

def group_hash(study, token):
    if token is None or token == "UNRESOLVED":
        return ""
    return hashlib.sha256(
        (study + "\0" + token).encode()
    ).hexdigest()[:20]

for p, digest in expected.items():
    assert p.is_file() and sha(p) == digest, p

geo = {}
for study, (gse, n_full, _, digest) in series.items():
    p = folder / "geo_soft_LOCAL_ONLY" / f"{gse}_family.soft.gz"
    assert p.is_file() and sha(p) == digest, p
    geo[study] = parse_soft(p)
    assert len(geo[study]) == n_full

selected = read(detail)
canonical = read(queue)
full = read(v5)

assert len(selected) == 64
assert len(canonical) == 599
assert len(full) == 5241

q_by_acc = {r["experiment_acc"]: r for r in canonical}
assert len(q_by_acc) == 599

# The additional cervical tumour remains a separate technical QC hold.
qc = [
    r for r in full if r["library_id"] == "SRX2407571"
]
assert len(qc) == 1
assert qc[0]["study"] == "SRP094795"
assert qc[0]["v5_original_sra_sample_name"] == "GSM2420406"
assert qc[0]["stage41_partition_NOT_ADMISSION"] == "QC_HOLD_NOT_ADMITTED"
assert qc[0]["recommended_provisional_disposition"] == "HOLD_TECHNICAL_QC"
assert qc[0]["FINAL_CLINICAL_APPROVAL"] == "NO_UNSIGNED"
assert "SRX2407571" not in q_by_acc

gbm_subjects = defaultdict(set)
for gsm, s in geo["SRP072494"].items():
    sid = trait(s, "subject id")
    if sid != "UNRESOLVED":
        gbm_subjects[sid].add(gsm)

breast_groups = defaultdict(set)
for gsm, s in geo["SRP194626"].items():
    p = prefix(s)
    if p is not None:
        breast_groups[p].add(gsm)

rows = []
metrics = defaultdict(Counter)

for r in sorted(selected, key=lambda x: (x["study"], x["SRX"])):
    study, gsm, srx = r["study"], r["GSM"], r["SRX"]
    gse = series[study][0]
    sample = geo[study][gsm]
    old = q_by_acc[srx]

    assert old["study"] == study
    assert r["GEO_relation_check"] == "EXACT_GSM_SRX_RELATION"
    assert old["final_clinical_approval"] == "NO_UNSIGNED"
    assert old["primary_sampled_site_review"] == "NOT_ADJUDICATED"
    assert old["collection_vs_treatment_review"] == "NOT_ADJUDICATED"

    group = ""
    linked = 0
    diagnosis_context = "NOT_APPLICABLE"

    if study == "SRP094795":
        assert trait(sample, "tissue") == "Uterine cervix tumor"
        assert "before patient treatment" in (
            " ".join(sample["protocol"]).lower()
        )
        specimen = "GEO_CERVICAL_TUMOUR_SOURCE_LABEL"
        timing = "GEO_PROTOCOL_BEFORE_PATIENT_TREATMENT"
        patient = "GEO_PATIENT_ID_NOT_ESTABLISHED"

    elif study == "SRP072494":
        assert trait(sample, "time point") == (
            "before bevacizumab combination treatment"
        )
        assert trait(sample, "first-line therapy") in (
            "RT + TMZ", "RT/TMZ+TMZ"
        )
        diagnosis_context = trait(
            sample, "sample from primary gbm diagnosis"
        )
        assert diagnosis_context in ("yes", "no")
        sid = trait(sample, "subject id")
        assert sid != "UNRESOLVED"
        group = group_hash(study, sid)
        partners = gbm_subjects[sid] - {gsm}
        assert all(
            trait(geo[study][x], "time point") ==
            "after bevacizumab combination treatment"
            for x in partners
        )
        linked = len(partners)
        specimen = "GEO_GBM_SOURCE_LABEL"
        timing = "BEFORE_BEVACIZUMAB_PRIOR_TREATMENT_TIMELINE_OPEN"
        patient = (
            "GEO_SUBJECT_TOKEN_WITH_LONGITUDINAL_PARTNER"
            if linked else "GEO_SUBJECT_TOKEN_WITHOUT_SERIES_PARTNER"
        )

    elif study == "SRP194626":
        assert trait(sample, "time point") == "pre-treatment biopsy"
        assert trait(sample, "sequencing platform") == "Lexogen QuantSeq"
        p = prefix(sample)
        assert p is not None
        group = group_hash(study, p)
        linked = len(breast_groups[p] - {gsm})
        assert linked >= 1
        specimen = "GEO_BREAST_CANCER_BIOPSY_SOURCE_LABEL"
        timing = "GEO_PRETREATMENT_BIOPSY_THERAPY_SCOPE_OPEN"
        patient = "GEO_PATIENT_TITLE_PREFIX_LINKED_NOT_CERTIFIED"

    else:
        assert study == "SRP112751"
        assert trait(sample, "cohort") == "Singapore"
        assert trait(sample, "tissue") == "Primary CRC total RNA"
        assert trait(sample, "site") != "UNRESOLVED"
        specimen = "GEO_PRIMARY_CRC_SOURCE_LABEL"
        timing = "PER_SPECIMEN_PRETREATMENT_NOT_ESTABLISHED"
        patient = "GEO_PATIENT_ID_NOT_ESTABLISHED"

    metrics[study]["n"] += 1
    metrics[study]["linked_other_GEO_records"] += linked
    metrics[study]["with_partner"] += int(linked > 0)
    if study == "SRP072494":
        metrics[study]["primary_diagnosis_" + diagnosis_context] += 1

    rows.append({
        "study": study,
        "GSE": gse,
        "GSM": gsm,
        "SRX": srx,
        "source_snapshot_sha256": series[study][3],
        "specimen_source_evidence": specimen,
        "treatment_source_evidence": timing,
        "within_study_group_evidence": patient,
        "source_primary_diagnosis_label": diagnosis_context,
        "source_group_token_HASH_LOCAL_ONLY": group,
        "n_linked_other_GEO_records": linked,
        "historical_treatment_grade":
            old["historical_treatment_grade"],
        "historical_patient_grade":
            old["historical_patient_grade"],
        "technical_release": "NOT_CERTIFIED",
        "patient_independence": "NOT_CERTIFIED",
        "treatment_naive": "NOT_CERTIFIED",
        "cross_study_independence": "NOT_CERTIFIED",
        "clinical_approval": "NO_UNSIGNED",
    })

assert len(rows) == 64
assert {s: c["n"] for s, c in metrics.items()} == {
    s: info[2] for s, info in series.items()
}
assert metrics["SRP072494"]["with_partner"] == 16
assert metrics["SRP072494"]["primary_diagnosis_yes"] == 10
assert metrics["SRP072494"]["primary_diagnosis_no"] == 10
assert metrics["SRP194626"]["with_partner"] == 4
assert metrics["SRP194626"]["linked_other_GEO_records"] == 5

# New outputs only. Never overwrite an earlier evidence addendum.
out.mkdir(parents=True, exist_ok=True)
local = out / "64_source_evidence_addendum_LOCAL_ONLY.tsv"
summary_path = out / "summary.json"
assert not local.exists() and not summary_path.exists(), (
    "Evidence addendum already exists; inspect it before rerunning"
)

with local.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(rows[0]),
        delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)

summary = {
    "status": "SOURCE_EVIDENCE_CLASSIFIED_NOT_CERTIFIED",
    "n_records": 64,
    "n_studies": 4,
    "n_GEO_SRX_relations_already_verified": 64,
    "study_metrics": {
        s: dict(c) for s, c in sorted(metrics.items())
    },
    "cervical_additional_v5_experiment": {
        "experiment": "SRX2407571",
        "partition": "QC_HOLD_NOT_ADMITTED",
        "clinical_approval": "NO_UNSIGNED",
        "included_in_64": False,
    },
    "input_hashes": {
        "geo_exact_64": sha(detail),
        "canonical_599_queue": sha(queue),
        "frozen_v5_5241": sha(v5),
    },
    "local_addendum_sha256": sha(local),
    "clinical_approvals": 0,
}
summary_path.write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)

print("V6_PRIORITY64_SOURCE_EVIDENCE_ADDENDUM")
for study, c in sorted(metrics.items()):
    print(study, dict(c))
print("CERVICAL_24TH_TUMOUR: QC_HOLD_RETAINED")
print("Local addendum SHA256:", sha(local))
print("Summary JSON SHA256:", sha(summary_path))
print("Clinical approvals: 0")
print("SOURCE_EVIDENCE_CLASSIFIED_NOT_CERTIFIED")
