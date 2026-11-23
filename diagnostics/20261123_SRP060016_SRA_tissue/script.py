from pathlib import Path
import sys, csv, gzip

root = Path(sys.argv[1])
name = "sra.sra.SRP060016.MD.gz"

print("=== SRP060016 RAW SOURCE AUDIT ===")

matches = sorted(root.rglob(name))
print("Matching source files:", len(matches))

if not matches:
    print("SOURCE FILE NOT FOUND IN PROJECT ROOT")
    sys.exit(0)

p = matches[0]
print("Using:", p)

with gzip.open(p, "rt", errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")
    cols = reader.fieldnames or []
    print("All columns:", cols)

    rows = list(reader)
    print("Total source rows:", len(rows))

    wanted = [
        c for c in cols
        if any(term in c.lower() for term in (
            "experiment", "sample", "biosample",
            "alias", "title", "tissue",
            "specimen", "organ", "source"
        ))
    ]

    print("\n=== EXPERIMENT AND SPECIMEN FIELDS ===")
    print("Selected columns:", wanted)

    for row in rows:
        accession = " ".join(
            row.get(c, "") for c in cols
            if "experiment" in c.lower()
        )

        if not accession.strip():
            continue

        print("\nEXPERIMENT:", accession)

        for c in wanted:
            val = row.get(c, "").strip()
            if val:
                print(f"  {c}: {val[:350]}")

print("\nREAD-ONLY AUDIT COMPLETE")
