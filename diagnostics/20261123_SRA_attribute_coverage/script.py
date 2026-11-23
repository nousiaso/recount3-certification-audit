import csv, gzip, re, sys
from pathlib import Path
from collections import Counter, defaultdict

rec = Path(sys.argv[1])
v5 = rec / "clinical_specimen_patient_v5_unsigned_20261122"
screen = v5 / "03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv"

csv.field_size_limit(50_000_000)

with screen.open(newline="") as f:
    candidates = list(csv.DictReader(f, delimiter="\t"))

by_study = defaultdict(dict)
for r in candidates:
    by_study[r["study"]][r["library_id"]] = r

KEYWORDS = (
    "tissue", "specimen", "organ", "anatomic",
    "tumor", "tumour", "metasta",
    "sample type", "disease", "treatment", "therapy"
)

def extract_attributes(text):
    values = []
    for item in (text or "").split("|"):
        key, sep, value = item.partition(";;")
        if sep and any(k in key.lower() for k in KEYWORDS):
            values.append((key.strip(), value.strip()))
    return values

coverage = Counter()
tissue_values = Counter()
study_summary = []
potential_nonprimary = []
missing_files = []
errors = []

for study, wanted in sorted(by_study.items()):
    path = (
        rec / "raw_studies" / study / "metadata"
        / f"sra.sra.{study}.MD.gz"
    )

    if not path.exists():
        missing_files.append((study, len(wanted)))
        continue

    found = defaultdict(list)

    try:
        with gzip.open(path, "rt", errors="replace",
                       newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                exp = (row.get("experiment_acc") or "").strip()
                if exp in wanted:
                    found[exp].append(row)
    except Exception as e:
        errors.append((study, str(e)))
        continue

    n_attr = 0
    n_informative_old = 0
    examples = []

    for exp, original in wanted.items():
        rows = found.get(exp, [])

        if not rows:
            coverage["missing_experiment"] += 1
            continue

        coverage["matched_experiment"] += 1

        attributes = set()
        for row in rows:
            attributes.update(extract_attributes(
                row.get("sample_attributes", "")
            ))

        if attributes:
            coverage["has_specific_attributes"] += 1
            n_attr += 1

            grade = original.get("tissue_evidence_grade", "")

            if grade in (
                "PRIMARY_TUMOUR_NOT_ATTESTED",
                "PRIMARY_SURGERY_CONTEXT_NOT_PATHOLOGY_CERTIFIED",
                "TISSUE_CONFLICT_UNRESOLVED"
            ):
                n_informative_old += 1

            for key, val in attributes:
                if "tissue" in key.lower():
                    tissue_values[val] += 1

                    if re.search(
                        r"\b(normal|healthy|benign|adjacent)\b"
                        r"|metasta",
                        val, re.I
                    ):
                        potential_nonprimary.append(
                            (study, exp, key, val,
                             original.get(
                                 "v5_specimen_evidence_class", ""
                             ))
                        )

            if len(examples) < 2:
                examples.append(
                    (exp, sorted(attributes)[:5])
                )

    study_summary.append((
        study, len(wanted), len(found),
        n_attr, n_informative_old, examples
    ))

print("=== GLOBAL RAW SRA COVERAGE ===")
print("Primary-screen records:", len(candidates))
print("Primary-screen studies:", len(by_study))
print("Studies with raw source:", len(study_summary))
print("Missing raw source files:", len(missing_files))
print("Coverage:", dict(coverage))

print("\n=== HIGHEST POTENTIAL EVIDENCE RECOVERY ===")
for s, total, matched, n, old, examples in sorted(
    study_summary, key=lambda x: (-x[4], -x[3])
)[:25]:
    print(
        f"{s}: screen={total}, matched={matched}, "
        f"attribute_evidence={n}, "
        f"previously_unattested_or_ambiguous={old}"
    )
    for exp, attrs in examples:
        print(" ", exp, str(attrs)[:220])

print("\n=== FREQUENT EXACT TISSUE VALUES ===")
for val, n in tissue_values.most_common(30):
    print(n, repr(val[:150]))

print("\n=== NONPRIMARY TISSUE LABELS TO REVIEW ===")
print("Total review flags:", len(potential_nonprimary))
for item in potential_nonprimary[:35]:
    print(item)

print("\n=== MISSING SOURCE FILES ===")
for item in missing_files[:40]:
    print(item)

print("\n=== SOURCE READ ERRORS ===")
for item in errors[:20]:
    print(item)

print("\nREAD-ONLY COVERAGE AUDIT COMPLETE")
