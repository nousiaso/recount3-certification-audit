#!/usr/bin/env python3
"""Verify original GEO GSM-to-SRX relations for 64 selected experiments.
Uses frozen priority linkage and downloaded public GEO SOFT snapshots.
No clinical decisions, frozen-input modifications or model execution.
"""
import csv, gzip, hashlib, json, re, sys
from collections import Counter, defaultdict
from pathlib import Path

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"
prior = repo / (
    "diagnostics/20261124_v6_priority_64_source_linkage/"
    "64_source_linkage_LOCAL_ONLY.tsv"
)
out = repo / "diagnostics/20261124_v6_priority_64_geo_exact"
snap = out / "geo_soft_LOCAL_ONLY"

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

assert sha(prior) == (
    "eb9ccdab319e9531f35b9c5ca2133feb9952fb8f56979b86ae426aea6a244e6b"
)

with prior.open(newline="", encoding="utf-8-sig") as f:
    selected = list(csv.DictReader(f, delimiter="\t"))

series = {
    "SRP094795": ("GSE91065", 25, 23),
    "SRP072494": ("GSE79671", 36, 20),
    "SRP194626": ("GSE130469", 20, 4),
    "SRP112751": ("GSE101588", 17, 17),
}

assert len(selected) == 64
assert len({r["experiment_acc"] for r in selected}) == 64
assert Counter(r["study"] for r in selected) == {
    s: n for s, (_, _, n) in series.items()
}
assert all(r["final_approval"] == "NO_UNSIGNED" for r in selected)

def parse_soft(path):
    samples = {}
    current = None

    with gzip.open(
        path, "rt", encoding="utf-8", errors="replace"
    ) as f:
        for line in f:
            line = line.rstrip("\r\n")

            if line.startswith("^"):
                if current is not None:
                    assert current["gsm"] not in samples
                    samples[current["gsm"]] = current
                current = None

                if line.startswith("^SAMPLE = "):
                    gsm = line.split(" = ", 1)[1].strip()
                    current = {
                        "gsm": gsm, "srx": set(),
                        "title": "", "source": "",
                        "traits": defaultdict(list),
                        "protocol": [],
                    }
                continue

            if (
                current is None
                or not line.startswith("!Sample_")
                or " = " not in line
            ):
                continue

            key, val = line.split(" = ", 1)
            val = val.strip()

            if key == "!Sample_relation":
                current["srx"].update(
                    re.findall(r"\bSRX\d+\b", val)
                )
            elif key == "!Sample_title":
                current["title"] = val
            elif key == "!Sample_source_name_ch1":
                current["source"] = val
            elif key == "!Sample_characteristics_ch1" and ":" in val:
                k, v = val.split(":", 1)
                current["traits"][k.strip().lower()].append(
                    v.strip()
                )
            elif key == "!Sample_treatment_protocol_ch1":
                current["protocol"].append(val)

    if current is not None:
        assert current["gsm"] not in samples
        samples[current["gsm"]] = current

    return samples

geo = {}
source_hashes = {}

for study, (gse, expected_count, _) in series.items():
    path = snap / f"{gse}_family.soft.gz"
    assert path.is_file(), f"GEO_SNAPSHOT_MISSING: {path}"

    samples = parse_soft(path)
    assert len(samples) == expected_count, (
        gse, len(samples), expected_count
    )
    geo[study] = samples
    source_hashes[gse] = sha(path)

stats = defaultdict(Counter)
local = []

for row in sorted(
    selected, key=lambda r: (r["study"], r["experiment_acc"])
):
    study = row["study"]
    srx = row["experiment_acc"]
    gsm = row["v5_GSM_alias_LOCAL_ONLY"].strip()
    sample = geo[study].get(gsm)

    if sample is None:
        status = "HOLD_GSM_ABSENT_FROM_GEO"
    elif not sample["srx"]:
        status = "HOLD_GEO_SRX_RELATION_MISSING"
    elif len(sample["srx"]) != 1:
        status = "HOLD_MULTIPLE_GEO_SRX_RELATIONS"
    elif sample["srx"] != {srx}:
        status = "HOLD_SRX_RELATION_MISMATCH"
    else:
        status = "EXACT_GSM_SRX_RELATION"

    stats[study][status] += 1

    if sample and (
        re.search(
            r"(?i)\bcell[ -]?line\b|caski",
            sample["source"]
        )
        or sample["traits"].get("cell line")
    ):
        stats[study]["selected_possible_cell_line"] += 1

    local.append({
        "study": study,
        "GSE": series[study][0],
        "GSM": gsm,
        "SRX": srx,
        "GEO_relation_check": status,
        "GEO_title_LOCAL_ONLY":
            sample["title"] if sample else "",
        "GEO_source_LOCAL_ONLY":
            sample["source"] if sample else "",
        "GEO_SRX_relations_LOCAL_ONLY":
            "|".join(sorted(sample["srx"])) if sample else "",
        "GEO_traits_json_LOCAL_ONLY":
            json.dumps(sample["traits"]) if sample else "",
        "GEO_treatment_protocol_json_LOCAL_ONLY":
            json.dumps(sample["protocol"]) if sample else "",
        "final_clinical_approval": "NO_UNSIGNED",
    })

# Original GEO study-level grouping checks.
gbm = geo["SRP072494"]
gbm_subjects = Counter(
    x["traits"]["subject id"][0]
    for x in gbm.values()
    if len(x["traits"].get("subject id", [])) == 1
)

# Breast titles contain candidate patient prefixes, e.g. P1-2.
# These are not automatically certified biological patient IDs.
breast_prefixes = Counter()
for r in local:
    if r["study"] == "SRP194626":
        m = re.match(r"^(P\d+)-", r["GEO_title_LOCAL_ONLY"])
        if m:
            breast_prefixes[m.group(1)] += 1

out.mkdir(parents=True, exist_ok=True)
details = out / "64_geo_exact_source_evidence_LOCAL_ONLY.tsv"

with details.open("w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=list(local[0]),
        delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(local)

n_exact = sum(
    c["EXACT_GSM_SRX_RELATION"] for c in stats.values()
)
n_possible_cell_line = sum(
    c["selected_possible_cell_line"] for c in stats.values()
)

status = (
    "EXACT_GEO_RELATIONS_PASS_NOT_CERTIFIED"
    if n_exact == 64 and n_possible_cell_line == 0
    else "GEO_EXACT_RELATION_EXCEPTIONS_HOLD"
)

summary = {
    "status": status,
    "prior_64_sha256": sha(prior),
    "geo_family_soft_sha256": source_hashes,
    "n_selected": 64,
    "n_exact_gsm_srx": n_exact,
    "n_possible_cell_line_selected": n_possible_cell_line,
    "studies": {
        study: {
            "n_selected": n,
            "n_original_GEO_samples": total,
            "relation_status_counts": dict(sorted(stats[study].items())),
        }
        for study, (_, total, n) in series.items()
    },
    "gbm_full_series_distinct_submitted_subject_tokens":
        len(gbm_subjects),
    "gbm_full_series_subject_token_multiplicity":
        dict(sorted(Counter(gbm_subjects.values()).items())),
    "breast_selected_distinct_title_patient_prefixes_UNVERIFIED":
        len(breast_prefixes),
    "breast_selected_reused_title_patient_prefix_groups_UNVERIFIED":
        sum(v > 1 for v in breast_prefixes.values()),
    "local_detail_sha256": sha(details),
    "clinical_approvals": 0,
}

report = out / "summary.json"
report.write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n"
)

print("V6_PRIORITY_64_GEO_EXACT_RELATION_AUDIT")
for study, d in summary["studies"].items():
    print(
        study,
        "selected", d["n_selected"],
        "GEO_samples", d["n_original_GEO_samples"],
        "relations", d["relation_status_counts"]
    )

print(
    "GBM full-series subject-token multiplicity:",
    summary["gbm_full_series_subject_token_multiplicity"]
)
print(
    "Breast reused title-prefix groups (UNVERIFIED):",
    summary[
        "breast_selected_reused_title_patient_prefix_groups_UNVERIFIED"
    ]
)
print("Exact GSM-to-SRX relations:", n_exact, "of 64")
print("Possible selected cell lines:", n_possible_cell_line)
print("GEO snapshot hashes:", json.dumps(source_hashes, sort_keys=True))
print("Local detail SHA256:", sha(details))
print("Aggregate SHA256:", sha(report))
print("STATUS:", status)
print("NO CLINICAL ADMISSIONS")
