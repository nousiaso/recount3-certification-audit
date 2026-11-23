#!/usr/bin/env python3
"""Archive hash-verified historical files. Source files are read-only."""
from pathlib import Path
import sys, hashlib, shutil, csv

root = Path(sys.argv[1])
repo = root / "RECOUNT3_CERTIFICATION_GIT"

frozen = (
    "RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/"
    "PROJECT_TREE/"
)
rec = "26_recount3_external_primary_v1/"

# Each source hash was observed in the cluster inventory.
items = [
    ("frozen_source",
     frozen + "06_audit_ARCHS4_primary_expression_compatibility_v1.py",
     "01191dfcd29ad729391ca14b837605009ccbdb7f5509bc8b4183236ff91eef32"),

    ("frozen_source",
     frozen + "09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py",
     "46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500"),

    ("frozen_source",
     frozen + "10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R",
     "16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c"),

    ("frozen_source",
     frozen + "11_build_ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1.R",
     "529551d7ecce253b1ab8626956afd076690631bf7d8b51cfcf0c613440a603f5"),

    ("frozen_source",
     frozen + "12_run_ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3.R",
     "a63750e5bfae63383a397aefb8af94a3e29b17bf3952b2cf47d8c3fc9b858f69"),

    ("frozen_source",
     frozen + "13_build_ARCHS4_projection_use_analysis_ready_atlas_v1_0_3.R",
     "18ec4761bd1374e9ca6fc63d52048f43e7680cc418d0f61d1fb013026c561907"),

    ("contract_script",
     "31_frozen_execution_contract_audit_v1_1_20261119/"
     "audit_frozen_outputs_READONLY_v1_1.py",
     "a6c1005a64483aba8da6afa3c5e0c3dfe282bb08b994c35a2d798ae96df4af05"),

    ("contract_result",
     "31_frozen_execution_contract_audit_v1_1_20261119/"
     "frozen_output_contract_v1_1.json",
     "855f6f7f06f5c3a6117410bd82e32717b004f041597c5ca218ac547d9187f3b8"),

    ("source_review",
     rec + "recount3_accession_source_evidence_v1_0_0.py",
     "ce0ad20af5dbe711aa28eca1cf15da363eb47e90d8612e5d4d98c2fcd04920f0"),

    ("v5_builder",
     rec + "RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122/"
     "scripts/build_source_crosswalk_v5.py",
     "74bf23c25c55febe99a42f51ed61e92d14c2e8a93db07b5aff712b5ecd7b634e"),

    ("w499_reference",
     frozen +
     "11_pan_cancer_fixed_anchor_projection_v1_0_4/"
     "16A_frozen_W_gene_symbol_harmonization.tsv",
     "fee127d95fcd94db2fd8d110c222428cf7298f838b8e5e3ce239de7415e71706"),
]

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()

def destination(group, src):
    if group in ("contract_result", "w499_reference"):
        return repo / "evidence" / "history" / group / src.name
    return repo / "scripts" / "originals" / group / src.name

# Verify everything BEFORE creating any copies.
verified = []
for group, relative, expected in items:
    src = root / relative
    if not src.is_file():
        raise SystemExit("MISSING SOURCE: " + str(src))
    actual = sha256(src)
    if actual != expected:
        raise SystemExit(
            f"HASH MISMATCH: {src}\n"
            f"expected={expected}\nactual={actual}"
        )
    dst = destination(group, src)
    if dst.exists() and sha256(dst) != expected:
        raise SystemExit("CONFLICTING EXISTING COPY: " + str(dst))
    verified.append((group, src, dst, expected))

print("PREFLIGHT PASS:", len(verified), "source files verified")

rows = []
for group, src, dst, expected in verified:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if not dst.exists():
        shutil.copyfile(src, dst)
    if sha256(dst) != expected:
        raise SystemExit("COPY VERIFICATION FAILED: " + str(dst))
    print("ARCHIVED:", dst.relative_to(repo))
    rows.append((
        group, str(src), str(dst.relative_to(repo)),
        expected, "RECOVERED_NOT_REEXECUTED"
    ))

manifest = repo / "docs" / "HISTORICAL_SOURCE_RECOVERY_MANIFEST.tsv"
with manifest.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "group", "original_cluster_path", "repository_path",
        "sha256", "execution_status"
    ])
    writer.writerows(rows)

print("MANIFEST:", manifest)
print("RECOVERY COMPLETE — ORIGINAL SOURCES UNCHANGED")
