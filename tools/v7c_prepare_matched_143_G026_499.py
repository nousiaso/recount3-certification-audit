#!/usr/bin/env python3
"""V7C: locked SRP074349 matched-assay 499-gene input preparation.

Reads the exact 143-pair Stage5 bridge, original Stage09 processing helpers,
and recount3 G026 coverage HDF5. Reuses frozen Stage09 gene collapse, TMM
and UQ implementations; NEVER treats base-pair coverage as RNA read counts.
No original files, W, frozen anchors or historical cohorts are changed.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

STUDY = "SRP074349"
BRANCHES = ("CPM", "TMM_CPM", "UQ_CPM")

def sha_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()

def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        return rd.fieldnames or [], list(rd)

def decode(x):
    if isinstance(x, (bytes, np.bytes_)):
        return bytes(x).decode("utf-8").strip()
    return str(x).strip()

def import_original(path):
    spec = importlib.util.spec_from_file_location("frozen_stage09_unchanged", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot import frozen Stage09 source")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def save_matrix(path, genes, ids, values):
    if values.shape != (len(genes), len(ids)):
        raise ValueError("Matrix serialization axis mismatch")
    with gzip.open(path, "wt", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["model_gene", *ids])
        for gene, vec in zip(genes, values):
            w.writerow([gene, *(format(float(x), ".17g") for x in vec)])

def main():
    p = argparse.ArgumentParser()
    p.add_argument("root", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    root = args.root.resolve()
    out = args.output.resolve()
    if out.exists():
        raise FileExistsError("OUTPUT_EXISTS_NO_OVERWRITE " + str(out))
    repo = root / "RECOUNT3_CERTIFICATION_GIT"
    if out.parent != (repo / "private_data").resolve():
        raise ValueError("Output must be directly under repo/private_data")
    out.mkdir(parents=True, mode=0o700)

    rec = root / "26_recount3_external_primary_v1"
    stage5 = rec / "stage5_output"
    rel = root / "RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
    original_py = rel / "09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py"
    if sha_file(original_py) != "46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500":
        raise ValueError("Original Stage09 source changed")
    old = import_original(original_py)

    bridge_path = stage5 / "03_bridge_165_availability.tsv"
    _, b = read_tsv(bridge_path)
    if len(b) != 165:
        raise ValueError("Frozen bridge has changed")
    pairs = [r for r in b if r["study"] == STUDY]
    if len(pairs) != 143:
        raise ValueError("Expected exactly 143 same-specimen bridge rows")
    pairs.sort(key=lambda r: r["ARCHS4_GSM"])
    gsm = [r["ARCHS4_GSM"].strip() for r in pairs]
    libs = [r["library_id"].strip() for r in pairs]
    indices = [int(r["ARCHS4_expression_column_index"]) for r in pairs]
    if not all(gsm) or len(set(gsm)) != 143 or len(set(libs)) != 143:
        raise ValueError("Duplicate/blank bridge GSM or library")
    if len(set(indices)) != 143 or min(indices) < 0:
        raise ValueError("Invalid frozen ARCHS4 HDF5 column indices")
    if not all(r["is_study_gene_file_present"].strip().upper() in ("TRUE", "YES", "1") for r in pairs):
        raise ValueError("Bridge includes source file not recorded as present")
    run_ids = []
    for r in pairs:
        runs = re.findall(r"(?:SRR|ERR|DRR)\d+", r["recount3_run_ids"], re.I)
        if len(runs) != 1 or len(set(runs)) != 1 or int(r["n_runs"]) != 1:
            raise ValueError("Expected one unambiguous recount3 run per bridge pair")
        run_ids.append(runs[0])
    if len(set(run_ids)) != 143:
        raise ValueError("Duplicated run in matched bridge")

    arch_h5 = root / "00_download/human_gene_v2.latest.h5"
    with h5py.File(arch_h5, "r") as f:
        access = f["meta/samples/geo_accession"]
        if max(indices) >= len(access):
            raise ValueError("Bridge ARCHS4 column outside accession axis")
        observed = [decode(access[i]) for i in indices]
        if observed != gsm:
            raise ValueError("ARCHS4 bridge GSM/index mismatch: no automatic index shift")
    print("EXACT_BRIDGE_ARCHS4_ACCESSION_INDEX_PASS", len(gsm))
    print("BRIDGE_ONE_TO_ONE_RUNS_PASS", len(run_ids))

    h5path = rec / "stage2_output/studies/SRP074349/selected_run_coverage_G026.h5"
    mapping_path = stage5 / "00_gene_row_mapping.tsv"
    cm, mapped = read_tsv(mapping_path)
    need = {"row_zero_based", "source_gene_id", "gene_base", "gene_symbol"}
    if not need <= set(cm) or len(mapped) != 63856:
        raise ValueError("Stage5 source gene map changed")
    mapped.sort(key=lambda r: int(r["row_zero_based"]))
    if [int(x["row_zero_based"]) for x in mapped] != list(range(63856)):
        raise ValueError("Gene row indices are not exact 0-based axis")
    symbols = [r["gene_symbol"] for r in mapped]

    with h5py.File(h5path, "r") as f:
        raw_genes = [decode(x) for x in f["gene_id"][:]]
        names = [decode(x) for x in f["run_id"][:]]
        library_axis = [decode(x) for x in f["library_id"][:]]
        if len(raw_genes) != 63856 or len(names) != 144 or len(library_axis) != 144:
            raise ValueError("Unexpected expression HDF5 axes")
        if len(set(names)) != 144 or len(set(library_axis)) != 144:
            raise ValueError("Repeated source run or library IDs")
        mismatches = sum(
            not (name == m["source_gene_id"] or name.split(".", 1)[0] == m["gene_base"])
            for name, m in zip(raw_genes, mapped)
        )
        if mismatches:
            raise ValueError("H5 gene axis disagrees with Stage5 gene map: " + str(mismatches))
        run_to_idx = {x: i for i, x in enumerate(names)}
        if not set(run_ids) <= set(run_to_idx):
            raise ValueError("A frozen matched-bridge run is absent from source HDF5")
        positions = [run_to_idx[x] for x in run_ids]
        for lib, position in zip(libs, positions):
            if library_axis[position] != lib:
                raise ValueError("Run-to-library membership mismatch")
        order = np.argsort(positions)
        sorted_pos = np.asarray(positions, dtype=int)[order]
        selected = np.asarray(f["coverage_bp"][:, sorted_pos], dtype=np.float64)
        restore = np.argsort(order)
        selected = selected[:, restore]
    if selected.shape != (63856, 143):
        raise ValueError("Selected coverage axis mismatch")
    if not np.isfinite(selected).all() or (selected < 0).any():
        raise ValueError("Negative or nonfinite recount3 coverage is not eligible")
    print("RECOUNT3_G026_SOURCE_ROWS", len(raw_genes))
    print("RECOUNT3_MATCHED_RUNS", selected.shape[1])
    print("RECOUNT3_G026_NOT_RNA_COUNTS")

    stage09_dir = root / "10_preprojection_normalization_v1_0_1"
    audit = stage09_dir / "08_model_gene_availability.tsv"
    if not audit.is_file():
        alt = rel / "10_preprojection_normalization_v1_0_1/08_model_gene_availability.tsv"
        audit = alt
    cols, availability = read_tsv(audit)
    if "model_gene" not in cols or len(availability) != 499:
        raise ValueError("Frozen Stage09 model gene audit not available")
    model_genes = [r["model_gene"].strip() for r in availability]
    if len(set(model_genes)) != 499:
        raise ValueError("Duplicated model gene symbol")
    gene_maps = old.build_gene_maps(pd.Index(symbols), model_genes)
    if len(gene_maps["available_model_genes"]) != 499:
        raise ValueError("Full 499-gene W input not reconstructed")
    print("MODEL_GENES_RESOLVED", len(gene_maps["available_model_genes"]))
    print("ALIAS_RESOLVED", int(gene_maps["model_gene_resolution"]["resolution"].eq("alias_resolved").sum()))

    collapsed = old.collapse_batch(
        selected,
        gene_maps["valid_positions"],
        gene_maps["row_to_unique"],
        len(gene_maps["unique_genes"]),
    )
    L = collapsed.sum(axis=0)
    if (L <= 0).any() or not np.isfinite(L).all():
        raise ValueError("Invalid all-valid-symbol normalization denominators")

    model_data = np.zeros((499, 143), dtype=float)
    np.add.at(
        model_data,
        gene_maps["model_row_indices"],
        selected[gene_maps["model_row_positions"], :],
    )
    if (model_data < 0).any() or not np.isfinite(model_data).all():
        raise ValueError("Invalid model-gene matrix")

    q75 = old.quantile75_fraction(collapsed, L)
    ref_idx = old.select_tmm_reference(gsm, q75, L)
    factors = np.ones(143, dtype=float)
    retained = np.zeros(143, dtype=int)
    for j in range(143):
        factors[j], retained[j] = old.tmm_factor(
            collapsed[:, j], collapsed[:, ref_idx],
            float(L[j]), float(L[ref_idx]), .30, .05
        )
    factors /= old.geometric_mean_positive(factors)
    if (factors <= 0).any() or not np.isfinite(factors).all():
        raise ValueError("Invalid Stage09 TMM library-size factors")

    uq_evaluable = np.isfinite(q75) & (q75 > 0)
    uq = q75.copy()
    if not uq_evaluable.any():
        raise ValueError("UQ has no evaluable samples")
    uq[~uq_evaluable] = np.median(uq[uq_evaluable])
    uq /= old.geometric_mean_positive(uq)

    norm = {
        "CPM": model_data / L[None, :] * 1e6,
        "TMM_CPM": model_data / (L * factors)[None, :] * 1e6,
        "UQ_CPM": model_data / (L * uq)[None, :] * 1e6,
    }
    for branch in BRANCHES:
        if not np.isfinite(norm[branch]).all():
            raise ValueError("Nonfinite normalized gene values in " + branch)
        save_matrix(out / ("499_" + branch + "_LOCAL_ONLY.tsv.gz"), model_genes, gsm, norm[branch])
    with (out / "143_pair_manifest_LOCAL_ONLY.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["ARCHS4_GSM", "ARCHS4_expression_column_index", "library_id", "run_id", "UQ_evaluable", "TMM_factor", "UQ_factor", "valid_symbol_library_mass"])
        for i, r in enumerate(pairs):
            w.writerow([gsm[i], indices[i], libs[i], run_ids[i], str(bool(uq_evaluable[i])).upper(), format(float(factors[i]), ".17g"), format(float(uq[i]), ".17g"), format(float(L[i]), ".17g")])
    provenance = {
        "matched_pairs": len(pairs),
        "study": STUDY,
        "assay": "recount3_G026_gene_body_base_pair_coverage_NOT_READ_COUNTS",
        "normalization": "frozen_Stage09_methods_recomputed_within_143_bridge_samples",
        "tmm_reference_sample_position_0based_LOCAL_ONLY": ref_idx,
        "uq_evaluable": int(uq_evaluable.sum()),
        "full_symbol_denominator": True,
        "model_genes": 499,
        "source_hashes": {
            "bridge": sha_file(bridge_path),
            "Stage5_gene_row_map": sha_file(mapping_path),
            "frozen_Stage09_source": sha_file(original_py),
            "Stage09_499_gene_audit": sha_file(audit),
        },
        "output_sha256": {
            name: sha_file(out / name)
            for name in [
                "499_CPM_LOCAL_ONLY.tsv.gz",
                "499_TMM_CPM_LOCAL_ONLY.tsv.gz",
                "499_UQ_CPM_LOCAL_ONLY.tsv.gz",
                "143_pair_manifest_LOCAL_ONLY.tsv",
            ]
        },
    }
    (out / "preparation_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print("FULL_499_GENE_NORM_BRANCHES_PASS", len(BRANCHES))
    print("UPPER_QUARTILE_EVALUABLE", int(uq_evaluable.sum()))
    print("MATCHED_MANIFEST_SHA256", provenance["output_sha256"]["143_pair_manifest_LOCAL_ONLY.tsv"])
    print("NO_FROZEN_OR_SOURCE_FILES_MODIFIED")
    print("V7C_MATCHED_143_NORMALIZATION_INPUT_PREPARATION_COMPLETE")

if __name__ == "__main__":
    main()
