#!/usr/bin/env python3
"""V6 NSCLC82: exact TCIA clinical linkage and treatment-field audit.

No frozen metadata changes. Patient-linked output stays private.
"""

import csv, hashlib, io, json, re, sys, urllib.request
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"

source = (
    repo / "private_data/v6_geo240_clinical_semantics_v1/"
    "240_exact_geo_sra_clinical_evidence_LOCAL_ONLY.tsv"
)
private = repo / "private_data/v6_nsclc82_tcia_clinical_v1"
out = repo / "diagnostics/20261124_v6_nsclc82_tcia_clinical_v1"

SOURCE_SHA = (
    "0740b18fb374e4f51134469aecda7a01aeb601c78eacefff3aec84c6a1cec292"
)
TCIA_SHA = (
    "0a1f7cdf02cef58450db17f116a3fb0ecfd4d737ca6843bb1939dbe421071476"
)

def sha(data):
    return hashlib.sha256(data).hexdigest()

assert sha(source.read_bytes()) == SOURCE_SHA

with source.open(newline="", encoding="utf-8-sig") as f:
    lung = [
        r for r in csv.DictReader(f, delimiter="\t")
        if r["study"] == "SRP117020"
    ]

assert len(lung) == 82
assert len({r["srx"] for r in lung}) == 82

url = (
    "https://wiki.cancerimagingarchive.net/download/"
    "attachments/28672347/"
    "NSCLCR01Radiogenomic_DATA_LABELS_2018-05-22_1500-shifted.csv"
    "?api=v2&modificationDate=1531967714295&version=1"
)

req = urllib.request.Request(
    url, headers={"User-Agent": "SpatialEcoTyper-V6-NSCLC82/1.0"}
)
with urllib.request.urlopen(req, timeout=45) as response:
    data = response.read(2_000_001)

assert len(data) <= 2_000_000
assert sha(data) == TCIA_SHA, "TCIA source version changed"

reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
fields = reader.fieldnames
clinical = list(reader)

assert fields and len(fields) == 40
assert len(clinical) == 211

fieldmap = {c.strip().casefold(): c for c in fields}
assert len(fieldmap) == len(fields)
assert "case id" in fieldmap
case_field = fieldmap["case id"]

def norm(x):
    return " ".join(str(x).strip().casefold().split())

cases = {}
for r in clinical:
    key = norm(r[case_field])
    assert key and key not in cases
    cases[key] = r

assert len(cases) == 211

selected = {}
for r in lung:
    title = norm(r["geo_title"])
    assert title and title in cases
    assert title not in selected
    selected[title] = (r, cases[title])

assert len(selected) == 82

def histology(x):
    x = norm(x)
    if "adenocarcinoma" in x:
        return "ADENOCARCINOMA"
    if "squamous" in x:
        return "SQUAMOUS_CELL_CARCINOMA"
    return "OTHER_OR_UNCLASSIFIED"

assert "histology" in fieldmap
hist_field = fieldmap["histology"]

agreement = Counter()
for geo, tcia in selected.values():
    chars = json.loads(geo["geo_characteristics_json"])
    labels = chars.get("histology", [])
    assert len(labels) == 1
    a = histology(labels[0])
    b = histology(tcia[hist_field])
    agreement[(a, b)] += 1

audit_fields = [
    "Histology",
    "Histopathological Grade",
    "Pathological T stage",
    "Pathological N stage",
    "Pathological M stage",
    "Adjuvant Treatment",
    "Chemotherapy",
    "Radiation",
]

missing_tokens = {
    "", "na", "n/a", "nan", "none", "null",
    "unknown", "not available", "-", "--"
}

def aggregate(field):
    c = Counter(
        " ".join(str(t[field]).strip().split())
        for _, t in selected.values()
    )
    missing = sum(
        n for value, n in c.items()
        if norm(value) in missing_tokens
    )
    # Only display short categorical values, never accession-like IDs.
    safe = (
        len(c) <= 15
        and all(len(v) <= 55 for v in c)
        and not any(
            re.search(r"\b(?:R01|AMC|GSM|SRX)-?\d+\b", v, re.I)
            for v in c
        )
    )
    return {
        "available_records": len(selected) - missing,
        "missing_or_unknown": missing,
        "n_distinct_literal_values": len(c),
        "literal_value_counts": dict(sorted(c.items())) if safe else None,
    }

summaries = {}
for name in audit_fields:
    if name.casefold() in fieldmap:
        original = fieldmap[name.casefold()]
        summaries[name] = aggregate(original)

timing_field_names = [
    f for f in fields
    if any(k in f.casefold() for k in (
        "neoadjuvant", "neoadj", "preoperative",
        "pre-operative", "induction", "prior therapy",
        "prior treatment"
    ))
]

# Complete validation before writing any output.
assert not private.exists(), "Prior private audit exists"
assert not (out / "summary.json").exists(), "Prior summary exists"

private.mkdir(parents=True)
out.mkdir(parents=True, exist_ok=True)

(private / "TCIA_NSCLC_RADIOGENOMICS_211_SOURCE.csv").write_bytes(data)

joinfields = [
    "study", "srx", "gsm", "geo_title", "v5_patient_key"
] + ["TCIA__" + x for x in fields]

joined = []
for title in sorted(selected):
    geo, tcia = selected[title]
    record = {
        "study": geo["study"],
        "srx": geo["srx"],
        "gsm": geo["gsm"],
        "geo_title": geo["geo_title"],
        "v5_patient_key": geo["v5_patient_key"],
    }
    record.update({"TCIA__" + k: tcia[k] for k in fields})
    joined.append(record)

detail = private / "82_exact_tcia_clinical_crosswalk_LOCAL_ONLY.tsv"
with detail.open("x", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(
        f, fieldnames=joinfields,
        delimiter="\t", lineterminator="\n"
    )
    w.writeheader()
    w.writerows(joined)

histology_match = (
    sum(n for (a, b), n in agreement.items() if a == b)
)

status = (
    "NSCLC82_TCIA_LINK_AND_HISTOLOGY_CONCORDANT"
    if histology_match == 82
    else "NSCLC82_TCIA_LINK_PASS_HISTOLOGY_DISCORDANCE_REVIEW"
)

summary = {
    "audit": "V6_NSCLC82_TCIA_CLINICAL_V1",
    "study": "SRP117020",
    "geo": "GSE103584",
    "selected_experiments": 82,
    "clinical_source_records": 211,
    "exact_unique_subject_links": 82,
    "histology_agreement": {
        "matched": histology_match,
        "cross_tab": [
            {"geo": a, "tcia": b, "n": n}
            for (a, b), n in sorted(agreement.items())
        ]
    },
    "clinical_fields": summaries,
    "explicit_preoperative_treatment_fields": timing_field_names,
    "published_treatment_naive_statement": (
        "COHORT_SOURCE_SUPPORTED_NOT_INDIVIDUALLY_ATTESTED"
    ),
    "clinical_approval": "NO_UNSIGNED",
    "source_sha256": sha(data),
    "private_crosswalk_sha256": sha(detail.read_bytes()),
    "status": status,
}

summary_path = out / "summary.json"
with summary_path.open("x", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, sort_keys=True)
    f.write("\n")

print("V6_NSCLC82_TCIA_CLINICAL_RECONCILIATION")
print("Clinical source records:", len(clinical))
print("Selected SRX-GEO-TCIA links:", len(selected))
print("TCIA source SHA256:", sha(data))
print("Histology agreement:", histology_match, "/ 82")
print("Histology cross-tab:", summary["histology_agreement"]["cross_tab"])

for name, values in summaries.items():
    print("\nCLINICAL_FIELD", repr(name))
    print(json.dumps(values, sort_keys=True))

print("\nEXPLICIT_PREOPERATIVE_TREATMENT_FIELDS",
      timing_field_names)
print("PRIVATE_CROSSWALK_SHA256", summary["private_crosswalk_sha256"])
print("SUMMARY_SHA256", sha(summary_path.read_bytes()))
print("STATUS:", status)
print("No patient identifiers printed")
print("No frozen metadata changed")
print("NO CLINICAL ADMISSIONS")
