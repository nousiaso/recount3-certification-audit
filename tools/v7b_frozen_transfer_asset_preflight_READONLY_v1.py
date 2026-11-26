#!/usr/bin/env python3
"""Targeted READ-ONLY asset preflight before the frozen 499-gene assay-transfer experiment.

No expression values loaded; no specimen identities printed; no clinical admissions.
The 143 SRP074349 same-specimen pairs are calibration, NOT independent validation.
"""
import argparse
import csv
import hashlib
import importlib.util
import os
import shutil
from pathlib import Path

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def describe_tsv(label, path):
    print("ASSET", label, "PRESENT", path.is_file())
    if not path.is_file():
        return False
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames or []
        n = sum(1 for _ in reader)
    print("ASSET", label, "ROWS", n, "FIELDS", len(headers))
    print("ASSET", label, "FIELD_NAMES", ";".join(headers[:60]))
    print("ASSET", label, "SHA256", sha(path))
    return True

def describe_hdf5(label, path, expected_runs):
    print("H5", label, "PRESENT", path.is_file())
    if not path.is_file():
        return False
    if importlib.util.find_spec("h5py") is None:
        print("H5", label, "H5PY_MISSING")
        return False
    import h5py
    with h5py.File(path, "r") as f:
        print("H5", label, "DATASETS", ";".join(sorted(f.keys())))
        for key in ("coverage_bp", "gene_id", "library_id", "run_id"):
            if key not in f:
                print("H5", label, "MISSING_DATASET", key)
                return False
            ds = f[key]
            print("H5", label, key, "SHAPE", "x".join(map(str, ds.shape)), "DTYPE", str(ds.dtype))
        if f["coverage_bp"].shape != (63856, expected_runs):
            raise ValueError("H5_AXIS_CONTRACT_MISMATCH " + label)
    return True

def scoped_find(base, wanted, max_paths=15000):
    """Search only the authoritative frozen release project-tree for two known files."""
    found = {n: [] for n in wanted}
    if not base.is_dir():
        return found
    visited = 0
    for directory, folders, files in os.walk(base):
        folders[:] = sorted(d for d in folders if d not in (".git", "__pycache__"))
        for filename in files:
            visited += 1
            if filename in found:
                found[filename].append(Path(directory) / filename)
        if visited > max_paths:
            print("FROZEN_TREE_SEARCH_LIMIT_REACHED", max_paths)
            break
    return found

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    rec = root / "26_recount3_external_primary_v1"
    repo = root / "RECOUNT3_CERTIFICATION_GIT"
    print("V7B_FROZEN_TRANSFER_ASSET_PREFLIGHT_V1")
    print("PYTHON_H5PY_AVAILABLE", importlib.util.find_spec("h5py") is not None)

    summary = repo / "private_data/v7_integrated_candidate_gates_v1/summary.json"
    print("V7_BASELINE_PRESENT", summary.is_file())
    if summary.is_file():
        import json
        with summary.open(encoding="utf-8") as f:
            d = json.load(f)
        print("V7_BASELINE_EXPERIMENTS", d.get("n_original_v5_primary_screen"))
        print("V7_BASELINE_STUDIES", d.get("n_studies"))
        print("V7_CLINICAL_RELEASES", d.get("n_clinically_certified"))
        if (d.get("n_original_v5_primary_screen"), d.get("n_studies"), d.get("n_clinically_certified")) != (2055, 77, 0):
            raise ValueError("V7_PRECONDITION_MISMATCH")

    stage5 = rec / "stage5_output"
    for label, name in (
        ("MATCHED_BRIDGE", "03_bridge_165_availability.tsv"),
        ("GENE_ROW_MAPPING", "00_gene_row_mapping.tsv"),
        ("LIBRARY_PREANALYSIS", "04_library_preanalysis_decision_3639.tsv"),
    ):
        describe_tsv(label, stage5 / name)

    h5_results = []
    for label, study, expected_runs in (
        ("MATCHED_CALIBRATION_SRP074349", "SRP074349", 144),
        ("FUSCC_TNBC_SRP157974", "SRP157974", 360),
    ):
        path = rec / "stage2_output/studies" / study / "selected_run_coverage_G026.h5"
        h5_results.append(describe_hdf5(label, path, expected_runs))

    w_candidates = [
        root.parent / "02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds",
        Path("/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds"),
    ]
    w_found = set()
    for path in w_candidates:
        if path.is_file() and str(path.resolve()) not in w_found:
            w_found.add(str(path.resolve()))
            print("FROZEN_W_PRESENT", True, "BYTES", path.stat().st_size, "SHA256", sha(path))
    print("FROZEN_W_ASSET_COUNT", len(w_found))

    release = root / "RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
    wanted = (
        "09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py",
        "10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R",
    )
    print("FROZEN_RELEASE_PROJECT_TREE_PRESENT", release.is_dir())
    located = scoped_find(release, wanted)
    for name in wanted:
        print("FROZEN_SCRIPT", name, "COUNT", len(located[name]))
        for p in located[name][:2]:
            print("FROZEN_SCRIPT_LOCATION", name, str(p.relative_to(root)))
            print("FROZEN_SCRIPT_SHA256", name, sha(p))

    rscript = shutil.which("Rscript")
    print("RSCRIPT_ON_PATH", bool(rscript))
    print("KNOWN_R_ENV_RSCRIPTPATH_PRESENT",
          Path("/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript").is_file())

    all_h5 = all(h5_results)
    print("CALIBRATION_H5_AXIS_CONTRACTS_PASS", all_h5)
    print("STAGE5_BRIDGE_IS_TECHNICAL_NOT_INDEPENDENT_VALIDATION")
    print("RECOUNT3_MEASURE_IS_GENE_BODY_COVERAGE_NOT_COUNTS")
    print("NEXT_EXPERIMENT_499_GENE_MATCHED_H_AND_FRACTION_COMPARISON")
    print("NO_MODEL_REFIT_OR_FROZEN_SOURCE_MODIFICATION")
    print("CLINICAL_ADMISSIONS 0")
    print("V7B_FROZEN_TRANSFER_ASSET_PREFLIGHT_V1_COMPLETE")

if __name__ == "__main__":
    main()
