#!/usr/bin/env python3
"""Exact-source 240-record GEO clinical semantics audit.

Append-only diagnostics; frozen V5 and clinical grades untouched.
Patient-linked source values are never printed or committed.
"""

import csv
import gzip
import hashlib
import io
import json
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"

private = repo / "private_data/v6_geo240_clinical_semantics_v1"
report = repo / "diagnostics/20261124_v6_geo240_clinical_semantics_v1"

studies = {
    "SRP113755": (
        "GSE102055", 75, 77,
        "dff7f303f4d8d8e581de88aa093bdc462e72011532cc29fd137483f6be210b65"
    ),
    "SRP036848": (
        "GSE54460", 83, 106,
        "ec27d7c8a5871dbede0cfda2fc1b8643b2fccbf30fce850444ddbc3cd84e8ddd"
    ),
    "SRP117020": (
        "GSE103584", 82, 130,
        "26e514f9f57e56909f906df304adfe189f6ec5165134746c1cfafca25c911924"
    ),
}

inputs = {
    "ledger": (
        repo / "diagnostics/20261124_v6_historical1696_integrated_v1/"
        "1696_clinical_source_evidence_LOCAL_ONLY.tsv",
        "d2cecf185f15dcd8c72ad8a2c4edc7e6481b5e0767c75245401f7db1df7762dc"
    ),
    "599": (
        repo / "diagnostics/20261123_v6_primary_unattested_review_v1_1/"
        "source_evidence_599_LOCAL_ONLY.tsv",
        "8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb"
    ),
    "1097": (
        repo / "diagnostics/20261124_v6_remaining1097_source_v1/"
        "1097_original_source_evidence_LOCAL_ONLY.tsv",
        "194a214e48b0105e9e6a40b1b72a8d6c5c7cf08abf4a728d274f80575f04bf1c"
    ),
    "v5": (
        root / "26_recount3_external_primary_v1/"
        "clinical_specimen_patient_v5_unsigned_20261122/"
        "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv",
        "06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49"
    ),
}

def digest(b):
    return hashlib.sha256(b).hexdigest()

def load(path, expected):
    data = path.read_bytes()
    assert digest(data) == expected, path
    return list(csv.DictReader(
        io.StringIO(data.decode("utf-8-sig")), delimiter="\t"
    ))

tables = {name: load(*spec) for name, spec in inputs.items()}

assert len(tables["ledger"]) == 1696
assert len(tables["599"]) == 599
assert len(tables["1097"]) == 1097
assert len(tables["v5"]) == 2055

ledger = {
    (r["study"], r["experiment_acc"]): r
    for r in tables["ledger"]
}
v5 = {r["library_id"]: r for r in tables["v5"]}

sra = {}
for origin in ("599", "1097"):
    for row in tables[origin]:
        key = (row["study"], row["experiment_acc"])
        assert key not in sra
        sra[key] = (origin, row)

assert len(ledger) == len(sra) == 1696
assert set(ledger) == set(sra)

def parse_sra(text):
    out = defaultdict(set)
    for segment in text.split("|"):
        if ";;" in segment:
            k, v = segment.split(";;", 1)
            if v.strip():
                out[k.strip().lower()].add(v.strip())
    return out

def fetch_geo(geo, expected_sha):
    group = re.sub(r"\d{3}$", "nnn", geo)
    url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{group}/{geo}/soft/{geo}_family.soft.gz"
    )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "V6-clinical-evidence-archive/1.0"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        data = response.read(20 * 1024 * 1024 + 1)
    assert len(data) <= 20 * 1024 * 1024
    assert digest(data) == expected_sha, (
        f"{geo}: downloaded source differs from verified original"
    )
    return data

def parse_soft(data):
    samples = {}
    current = None
    table_mode = False

    with gzip.open(io.BytesIO(data), "rt",
                   encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")

            if line.startswith("^SAMPLE = "):
                gsm = line.split("=", 1)[1].strip()
                assert re.fullmatch(r"GSM\d+", gsm)
                assert gsm not in samples
                current = {
                    "gsm": gsm,
                    "title": "",
                    "source_name": "",
                    "srx": set(),
                    "chars": defaultdict(set),
                }
                samples[gsm] = current
                table_mode = False
                continue

            if line.startswith("^"):
                current = None
                table_mode = False
                continue

            if current is None:
                continue

            if line.lower().startswith("!sample_table_begin"):
                table_mode = True
                continue
            if line.lower().startswith("!sample_table_end"):
                table_mode = False
                continue
            if table_mode or not line.startswith("!Sample_"):
                continue

            key, sep, value = line.partition(" = ")
            if not sep:
                continue

            key = key[len("!Sample_"):].lower()
            value = value.strip()

            if key == "relation":
                current["srx"].update(
                    x.upper() for x in
                    re.findall(r"\bSRX\d+\b", value, re.I)
                    if "SRA:" in value.upper()
                )
            elif key == "title":
                current["title"] = value
            elif key.startswith("source_name"):
                current["source_name"] = value
            elif key.startswith("characteristics") and ":" in value:
                k, v = value.split(":", 1)
                if v.strip():
                    current["chars"][k.strip().lower()].add(v.strip())

    return samples

def literal(v):
    return " ".join(v.strip().casefold().split())

def grouping(rows, field):
    values = defaultdict(set)
    for r in rows:
        val = r[field]
        if val:
            values[literal(val)].add(r["srx"])
    duplicates = [v for v in values.values() if len(v) > 1]
    return {
        "populated": sum(bool(r[field]) for r in rows),
        "distinct": len(values),
        "repeated_groups": len(duplicates),
        "experiments_in_repeated_groups": sum(map(len, duplicates)),
        "largest_group": max(map(len, values.values()), default=0),
    }

safe_fields = {
    "SRP113755": ("tissue type", "tissue", "disease", "disease_stage"),
    "SRP036848": ("tissue", "collection site", "pstage", "bcr"),
    "SRP117020": (
        "histology", "disease state", "histopathological grade"
    ),
}

def categorical_counts(rows, source_type, key):
    counts = Counter()
    for row in rows:
        values = row["geo_chars"].get(key, set()) \
            if source_type == "GEO" else row["sra_attrs"].get(key, set())
        for v in values:
            counts[v] += 1

    # Only predetermined clinical category fields, not patient-linked IDs.
    safe = {
        v: n for v, n in counts.items()
        if len(v) <= 90 and not re.search(
            r"\b(?:GSM|SRX|SRR|SAMN)\d+\b", v, re.I
        )
    }
    return dict(sorted(safe.items(), key=lambda x: (-x[1], x[0])))

# Verify all three sources and exact links before creating any output.
snapshots = {}
for study, (geo, expected_n, full_n, sha) in studies.items():
    data = fetch_geo(geo, sha)
    samples = parse_soft(data)
    assert len(samples) == full_n

    direct = {}
    for g in samples.values():
        assert len(g["srx"]) == 1
        accession = next(iter(g["srx"]))
        assert accession not in direct
        direct[accession] = g

    selected = {
        acc for (st, acc) in ledger if st == study
    }
    assert len(selected) == expected_n
    assert selected <= set(direct)
    snapshots[study] = (data, samples, direct)

assert not private.exists(), "Local snapshot exists; inspect before rerunning"
private.mkdir(parents=True)
report.mkdir(parents=True, exist_ok=True)

print("V6_GEO240_CLINICAL_SEMANTICS_AND_PATIENT_GROUPING")
summary = {
    "status": "EXACT_GEO240_LINK_PASS_CLINICAL_SEMANTICS_REVIEW",
    "clinical_approvals": 0,
    "studies": {},
}
detail_rows = []

for study, (geo, expected_n, full_n, sha) in studies.items():
    data, samples, direct = snapshots[study]
    (private / f"{geo}_family.soft.gz").write_bytes(data)

    selected = {
        acc for (st, acc) in ledger if st == study
    }

    all_rows = []
    for acc, g in direct.items():
        sra_attrs = {}
        if acc in selected:
            origin, source = sra[(study, acc)]
            raw = (
                source["original_sample_attributes"]
                if origin == "599"
                else source["source_sample_attributes_LOCAL_ONLY"]
            )
            sra_attrs = parse_sra(raw)

        def one(values):
            return next(iter(values)) if len(values) == 1 else ""

        chars = g["chars"]
        row = {
            "srx": acc,
            "gsm": g["gsm"],
            "geo_title": g["title"],
            "geo_source_name": g["source_name"],
            "external_id": one(chars.get("external id", set())),
            "cgsr_id": one(chars.get("cgsr id", set())),
            "sra_isolate": one(sra_attrs.get("isolate", set())),
            "v5_patient_key": (
                v5[acc]["patient_key"] if acc in selected else ""
            ),
            "geo_chars": chars,
            "sra_attrs": sra_attrs,
            "is_selected": acc in selected,
        }
        all_rows.append(row)

        if acc in selected:
            detail_rows.append({
                "study": study,
                "srx": acc,
                "gsm": g["gsm"],
                "geo_title": g["title"],
                "geo_source_name": g["source_name"],
                "external_id": row["external_id"],
                "cgsr_id": row["cgsr_id"],
                "sra_isolate": row["sra_isolate"],
                "v5_patient_key": row["v5_patient_key"],
                "geo_characteristics_json": json.dumps(
                    {k: sorted(v) for k, v in chars.items()},
                    sort_keys=True
                ),
                "sra_attributes_json": json.dumps(
                    {k: sorted(v) for k, v in sra_attrs.items()},
                    sort_keys=True
                ),
            })

    selected_rows = [r for r in all_rows if r["is_selected"]]
    assert len(selected_rows) == expected_n

    print("\nSTUDY", study, "GEO", geo)
    print("EXACT_SELECTED", len(selected_rows))
    print("FULL_GEO", len(all_rows))
    print("GEO_SOFT_SHA256", sha)

    groups = {}
    for field in (
        "geo_title", "external_id", "cgsr_id",
        "sra_isolate", "v5_patient_key"
    ):
        groups[field] = {
            "full_GEO": grouping(all_rows, field),
            "selected": grouping(selected_rows, field)
        }
        print("IDENTIFIER_STRUCTURE", field,
              json.dumps(groups[field], sort_keys=True))

    attr_agree = {}
    for key in ("external id", "cgsr id", "source_name"):
        comparable = match = 0
        for r in selected_rows:
            a = r["sra_attrs"].get(key, set())
            b = r["geo_chars"].get(key, set())
            if key == "source_name":
                target = r["geo_source_name"]
                if len(a) == 1 and target:
                    comparable += 1
                    match += literal(next(iter(a))) == literal(target)
            elif len(a) == len(b) == 1:
                comparable += 1
                match += literal(next(iter(a))) == literal(next(iter(b)))
        attr_agree[key] = {"comparable": comparable, "agree": match}
        print("SRA_VS_GEO_FIELD", key, attr_agree[key])

    categorical = {}
    for key in safe_fields[study]:
        categorical[key] = {
            "GEO": categorical_counts(selected_rows, "GEO", key),
            "SRA": categorical_counts(selected_rows, "SRA", key),
        }
        print("CLINICAL_CATEGORY", key,
              json.dumps(categorical[key], sort_keys=True))

    summary["studies"][study] = {
        "GEO": geo,
        "selected_exact": len(selected_rows),
        "full_GEO": len(all_rows),
        "GEO_SOFT_SHA256": sha,
        "identifier_structure": groups,
        "source_field_agreement": attr_agree,
        "categorical_evidence": categorical,
        "patient_identity_certified": False,
        "clinical_approval": "NO_UNSIGNED",
    }

assert len(detail_rows) == 240

detail = private / "240_exact_geo_sra_clinical_evidence_LOCAL_ONLY.tsv"
with detail.open("x", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(
        f, fieldnames=list(detail_rows[0]),
        delimiter="\t", lineterminator="\n"
    )
    w.writeheader()
    w.writerows(detail_rows)

summary["local_evidence_sha256"] = hashlib.sha256(
    detail.read_bytes()
).hexdigest()

summary_path = report / "summary.json"
with summary_path.open("x", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, sort_keys=True)
    f.write("\n")

print("\nLOCAL_EVIDENCE_SHA256", summary["local_evidence_sha256"])
print("SUMMARY_SHA256", hashlib.sha256(
    summary_path.read_bytes()
).hexdigest())
print("LOCAL_EVIDENCE_ROWS", len(detail_rows))
print("STATUS EXACT_GEO240_LINK_PASS_CLINICAL_SEMANTICS_REVIEW")
print("NO patient identifiers printed")
print("NO frozen data modified")
print("NO clinical approvals")
