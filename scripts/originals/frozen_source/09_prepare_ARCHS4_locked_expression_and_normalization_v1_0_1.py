#!/usr/bin/env python3
"""
ARCHS4 Stage 09: locked expression extraction and preprojection normalization.

This stage consumes the authoritative Stage 08 v1.0.1 biological lock and
extracts expression directly from the current ARCHS4 human HDF5 release using
the frozen zero-based expression-column indices.

It does not run SpatialEcoTyper.

Primary outputs
---------------
1. Ordered canonical and cross-sectional anchor manifests.
2. Direct-HDF5 model-gene ARCHS4 pseudocount matrix.
3. CPM, TMM-CPM and upper-quartile CPM matrices.
4. Study- and cancer-balanced anchor weights.
5. Complete HDF5 alignment, model-gene and normalization QC contracts.

The primary projection branch remains library-size CPM. TMM-CPM and UQ-CPM
remain sensitivity branches. These are normalized ARCHS4 Kallisto
pseudocount-abundance matrices, not literal transcript-level TPM.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import socket
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np
import pandas as pd

try:
    import pyarrow.parquet as pq
except ImportError:
    pq = None


SCRIPT_VERSION = "1.0.1"
EXPECTED_STAGE08_SCRIPT_VERSION = "1.0.1"

PROJECT_ROOT = Path(
    "/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6"
)
DEFAULT_LOCKED_DIR = (
    PROJECT_ROOT / "09_patient_canonical_tier_resolution_v1_0_1"
)
DEFAULT_H5 = PROJECT_ROOT / "00_download/human_gene_v2.latest.h5"
DEFAULT_MODEL_WEIGHTS = Path(
    "/home/tpn890/orestis/ST_ECOTYPES/"
    "00_preflight_ARCHS4_solid_tumors/official_Bulk_SE_Recovery_W.tsv"
)
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "10_preprojection_normalization_v1_0_1"

EXPECTED_COMBINED = 9536
EXPECTED_CANONICAL = 9465
EXPECTED_CROSS_SECTIONAL = 8962
EXPECTED_GOLD = 7821
EXPECTED_SILVER = 1715
EXPECTED_CANCERS = 18
EXPECTED_H5_GENES = 67186
EXPECTED_H5_SAMPLES = 1093742
EXPECTED_MODEL_GENES = 499
EXPECTED_AVAILABLE_MODEL_GENES = 499
EXPECTED_ALIAS_RESOLVED_MODEL_GENES = 7

GOLD = "GOLD_CONFIRMED_WHOLE_TUMOUR_WHOLE_TRANSCRIPTOME"
SILVER = "SILVER_CONFIRMED_WHOLE_TUMOUR_ACCEPTABLE_COVERAGE"

MODEL_GENE_ALIASES = {
    "CXORF36": "DIPK2B",
    "FAM71A": "GARIN4",
    "GPR1": "CMKLR2",
    "H2AFY": "MACROH2A1",
    "MRVI1": "IRAG1",
    "RARRES3": "PLAAT4",
    "WARS": "WARS1",
}

COMPONENTS = ["NonSE", *[f"SE{i}" for i in range(1, 10)]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract the Stage 08 locked canonical ARCHS4 matrix directly "
            "from HDF5 and create preprojection normalization branches."
        )
    )
    parser.add_argument("--locked-dir", type=Path, default=DEFAULT_LOCKED_DIR)
    parser.add_argument("--h5", type=Path, default=DEFAULT_H5)
    parser.add_argument("--model-weights", type=Path, default=DEFAULT_MODEL_WEIGHTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--sample-batch-size", type=int, default=16)
    parser.add_argument("--tmm-logratio-trim", type=float, default=0.30)
    parser.add_argument("--tmm-sum-trim", type=float, default=0.05)
    parser.add_argument("--expected-combined", type=int, default=EXPECTED_COMBINED)
    parser.add_argument("--expected-canonical", type=int, default=EXPECTED_CANONICAL)
    parser.add_argument(
        "--expected-cross-sectional",
        type=int,
        default=EXPECTED_CROSS_SECTIONAL,
    )
    parser.add_argument("--expected-cancers", type=int, default=EXPECTED_CANCERS)
    parser.add_argument("--expected-h5-genes", type=int, default=EXPECTED_H5_GENES)
    parser.add_argument("--expected-h5-samples", type=int, default=EXPECTED_H5_SAMPLES)
    parser.add_argument(
        "--expected-available-model-genes",
        type=int,
        default=EXPECTED_AVAILABLE_MODEL_GENES,
    )
    parser.add_argument(
        "--expected-alias-resolved-model-genes",
        type=int,
        default=EXPECTED_ALIAS_RESOLVED_MODEL_GENES,
    )
    parser.add_argument("--self-test", "--run-self-tests", action="store_true")
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(message: str) -> None:
    print(f"[{utc_now()}] {message}", flush=True)


def decode_string_array(values: np.ndarray) -> np.ndarray:
    decoded: list[str] = []
    for value in values:
        if isinstance(value, (bytes, np.bytes_)):
            decoded.append(bytes(value).decode("utf-8", errors="replace"))
        else:
            decoded.append(str(value).strip())
    return np.asarray(decoded, dtype=object)


def h5_take_columns_preserve_order(
    dataset: h5py.Dataset,
    column_indices: np.ndarray,
) -> np.ndarray:
    """Read arbitrary unique HDF5 columns and restore the requested order.

    h5py requires fancy indices to be strictly increasing. Stage 08 manifests
    are ordered biologically, not by HDF5 column. This helper sorts indices only
    for the physical read and then returns columns in the exact locked order.
    """
    columns = np.asarray(column_indices, dtype=np.int64)
    if columns.ndim != 1:
        raise ValueError("column_indices must be one-dimensional.")
    if len(columns) == 0:
        return np.empty((dataset.shape[0], 0), dtype=dataset.dtype)
    if len(np.unique(columns)) != len(columns):
        raise RuntimeError("HDF5 column request contains duplicated indices.")

    read_order = np.argsort(columns, kind="mergesort")
    sorted_columns = columns[read_order]
    if np.any(np.diff(sorted_columns) <= 0):
        raise RuntimeError("Sorted HDF5 column indices are not strictly increasing.")

    sorted_values = np.asarray(dataset[:, sorted_columns])
    restore_order = np.empty_like(read_order)
    restore_order[read_order] = np.arange(len(read_order), dtype=read_order.dtype)
    return sorted_values[:, restore_order]


def h5_take_vector_preserve_order(
    dataset: h5py.Dataset,
    indices: np.ndarray,
) -> np.ndarray:
    """Read arbitrary unique positions from a one-dimensional HDF5 dataset."""
    positions = np.asarray(indices, dtype=np.int64)
    if positions.ndim != 1:
        raise ValueError("indices must be one-dimensional.")
    if len(positions) == 0:
        return np.asarray([], dtype=dataset.dtype)
    if len(np.unique(positions)) != len(positions):
        raise RuntimeError("HDF5 vector request contains duplicated indices.")

    read_order = np.argsort(positions, kind="mergesort")
    sorted_positions = positions[read_order]
    if np.any(np.diff(sorted_positions) <= 0):
        raise RuntimeError("Sorted HDF5 vector indices are not strictly increasing.")

    sorted_values = np.asarray(dataset[sorted_positions])
    restore_order = np.empty_like(read_order)
    restore_order[read_order] = np.arange(len(read_order), dtype=read_order.dtype)
    return sorted_values[restore_order]


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            hasher.update(block)
    return hasher.hexdigest()


def require_nonempty(path: Path, label: str) -> None:
    if not path.exists() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing or empty {label}: {path}")


def bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        raise RuntimeError(f"Required boolean column is missing: {column}")
    return frame[column].fillna(False).astype(bool)


def prepare_downstream_manifest(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["sample_id"] = result["sample_id"].fillna("").astype(str).str.strip()
    result["primary_gse_id"] = (
        result["primary_gse_id"].fillna("").astype(str).str.strip()
    )
    result["analysis_cancer_class"] = (
        result["stage05_resolved_cancer_class"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    tier_map = {GOLD: "GOLD", SILVER: "SILVER"}
    result["provisional_cohort"] = result["stage08_tier"].map(tier_map)
    if result["provisional_cohort"].isna().any():
        bad = sorted(result.loc[result["provisional_cohort"].isna(), "stage08_tier"].astype(str).unique())
        raise RuntimeError(f"Unexpected Stage 08 tier values: {bad}")
    return result


def validate_stage08_inputs(
    locked_dir: Path,
    expected_combined: int,
    expected_canonical: int,
    expected_anchor: int,
    expected_cancers: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], pd.DataFrame]:
    status_path = locked_dir / "16_patient_canonical_resolution_status.json"
    qc_path = locked_dir / "15_patient_canonical_resolution_QC_contract.tsv"
    combined_path = locked_dir / "22_FINAL_LOCKED_combined_primary.tsv.gz"
    canonical_path = locked_dir / "23_FINAL_LOCKED_canonical_library_atlas.tsv.gz"
    anchor_path = locked_dir / "24_FINAL_LOCKED_cross_sectional_anchor.tsv.gz"
    normalization_ids_path = (
        locked_dir / "25_FINAL_LOCKED_normalization_ready_sample_ids.tsv"
    )

    for path, label in [
        (status_path, "Stage 08 status"),
        (qc_path, "Stage 08 QC contract"),
        (combined_path, "Stage 08 combined lock"),
        (canonical_path, "Stage 08 canonical lock"),
        (anchor_path, "Stage 08 cross-sectional anchor"),
        (normalization_ids_path, "Stage 08 normalization-ready IDs"),
    ]:
        require_nonempty(path, label)

    status = json.loads(status_path.read_text())
    required_status = {
        "status": "PASS",
        "script_version": EXPECTED_STAGE08_SCRIPT_VERSION,
        "mode": "lock",
        "audit_state": "FINAL_LOCK_VALIDATED",
        "final_lock_created": True,
        "normalization_ready_manifest_created": True,
        "expression_matrix_exported": False,
        "manual_review_condition_groups": 0,
        "qc_failures": 0,
        "input_expression_eligible_samples": expected_combined,
        "proposed_canonical_libraries": expected_canonical,
        "proposed_cross_sectional_anchor_libraries": expected_anchor,
    }
    for key, expected in required_status.items():
        observed = status.get(key)
        if observed != expected:
            raise RuntimeError(
                f"Stage 08 status mismatch for {key}: "
                f"observed={observed!r};expected={expected!r}"
            )

    stage08_qc = pd.read_csv(qc_path, sep="\t", low_memory=False)
    if stage08_qc.empty or set(stage08_qc["status"].astype(str)) != {"PASS"}:
        raise RuntimeError("Stage 08 QC contract is not an all-PASS contract.")

    combined = prepare_downstream_manifest(
        pd.read_csv(combined_path, sep="\t", compression="gzip", low_memory=False)
    )
    canonical = prepare_downstream_manifest(
        pd.read_csv(canonical_path, sep="\t", compression="gzip", low_memory=False)
    )
    anchor = prepare_downstream_manifest(
        pd.read_csv(anchor_path, sep="\t", compression="gzip", low_memory=False)
    )
    normalization_ids = pd.read_csv(
        normalization_ids_path, sep="\t", low_memory=False
    )

    required_columns = {
        "sample_id",
        "primary_gse_id",
        "expression_column_index",
        "expression_column_number",
        "stage05_resolved_cancer_class",
        "stage08_tier",
        "stage08_canonical_library",
        "stage08_cross_sectional_anchor",
        "expression_profile_sha256",
        "library_sum_archs4_pseudocounts",
        "stage07_expression_quality_tier",
    }
    for label, frame in [
        ("combined", combined),
        ("canonical", canonical),
        ("anchor", anchor),
    ]:
        missing = sorted(required_columns - set(frame.columns))
        if missing:
            raise RuntimeError(f"Stage 08 {label} manifest lacks columns: {missing}")
        validate_unique_ids(frame, f"Stage 08 {label}")
        if frame["primary_gse_id"].eq("").any():
            raise RuntimeError(f"Stage 08 {label} has blank primary_gse_id.")
        if frame["analysis_cancer_class"].eq("").any():
            raise RuntimeError(
                f"Stage 08 {label} has blank analysis_cancer_class."
            )
        frame["expression_column_index"] = pd.to_numeric(
            frame["expression_column_index"], errors="raise"
        ).astype(np.int64)
        frame["expression_column_number"] = pd.to_numeric(
            frame["expression_column_number"], errors="raise"
        ).astype(np.int64)
        if not np.array_equal(
            frame["expression_column_number"].to_numpy(),
            frame["expression_column_index"].to_numpy() + 1,
        ):
            raise RuntimeError(
                f"Stage 08 {label} has invalid one-based expression columns."
            )
        if frame["expression_column_index"].duplicated().any():
            raise RuntimeError(
                f"Stage 08 {label} has duplicated HDF5 expression columns."
            )

    observed_sizes = (len(combined), len(canonical), len(anchor))
    expected_sizes = (expected_combined, expected_canonical, expected_anchor)
    if observed_sizes != expected_sizes:
        raise RuntimeError(
            f"Stage 08 lock sizes differ: observed={observed_sizes};"
            f"expected={expected_sizes}"
        )

    combined_ids = set(combined["sample_id"])
    canonical_ids = set(canonical["sample_id"])
    anchor_ids = set(anchor["sample_id"])
    if not canonical_ids.issubset(combined_ids):
        raise RuntimeError("Canonical lock is not a subset of the combined lock.")
    if not anchor_ids.issubset(canonical_ids):
        raise RuntimeError("Cross-sectional anchor is not a subset of canonical.")
    normalization_sample_ids = (
        normalization_ids["sample_id"].fillna("").astype(str).str.strip().tolist()
    )
    canonical_sample_ids = canonical["sample_id"].astype(str).tolist()
    if normalization_sample_ids != canonical_sample_ids:
        if set(normalization_sample_ids) != canonical_ids:
            missing = sorted(canonical_ids - set(normalization_sample_ids))
            extra = sorted(set(normalization_sample_ids) - canonical_ids)
            raise RuntimeError(
                "Normalization-ready Stage 08 IDs do not equal the canonical lock. "
                f"Missing={missing[:20]};extra={extra[:20]}"
            )
        raise RuntimeError(
            "Normalization-ready Stage 08 IDs contain the canonical set but not "
            "the exact locked canonical order."
        )
    if not bool_series(canonical, "stage08_canonical_library").all():
        raise RuntimeError("Canonical lock contains a noncanonical row.")
    if not bool_series(anchor, "stage08_cross_sectional_anchor").all():
        raise RuntimeError("Anchor lock contains a non-anchor row.")

    tier_counts = combined["stage08_tier"].value_counts().to_dict()
    if tier_counts.get(GOLD, 0) != EXPECTED_GOLD:
        raise RuntimeError(f"Unexpected Gold count: {tier_counts}")
    if tier_counts.get(SILVER, 0) != EXPECTED_SILVER:
        raise RuntimeError(f"Unexpected Silver count: {tier_counts}")

    n_cancers = canonical["analysis_cancer_class"].nunique()
    if n_cancers != expected_cancers:
        raise RuntimeError(
            f"Expected {expected_cancers} canonical cancer classes, "
            f"observed {n_cancers}."
        )

    duplicate_hash = (
        canonical["expression_profile_sha256"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    duplicate_hash = duplicate_hash[duplicate_hash.ne("")]
    if duplicate_hash.duplicated().any():
        examples = duplicate_hash[duplicate_hash.duplicated(False)].head(20).tolist()
        raise RuntimeError(
            "Exact expression-profile duplicates remain in canonical lock: "
            f"{examples}"
        )

    # Preserve the exact Stage 08 biological lock order. HDF5 indices are
    # sorted only inside the low-level read helper and restored immediately.
    canonical = canonical.reset_index(drop=True)
    canonical["stage09_matrix_column_index"] = np.arange(
        len(canonical), dtype=np.int64
    )
    canonical["stage09_matrix_column_number"] = (
        canonical["stage09_matrix_column_index"] + 1
    )

    canonical_position = canonical.set_index("sample_id")[
        "stage09_matrix_column_index"
    ]
    anchor = anchor.copy().reset_index(drop=True)
    anchor["stage09_matrix_column_index"] = anchor["sample_id"].map(
        canonical_position
    )
    if anchor["stage09_matrix_column_index"].isna().any():
        raise RuntimeError("At least one anchor is absent from canonical order.")
    anchor["stage09_matrix_column_index"] = anchor[
        "stage09_matrix_column_index"
    ].astype(np.int64)
    if len(anchor) > 1 and np.any(
        np.diff(anchor["stage09_matrix_column_index"].to_numpy()) <= 0
    ):
        raise RuntimeError(
            "Stage 08 anchor order is not the exact canonical subsequence order."
        )
    anchor["stage09_matrix_column_number"] = (
        anchor["stage09_matrix_column_index"] + 1
    )

    return combined, canonical, anchor, status, stage08_qc


def validate_h5(
    h5: h5py.File,
    expected_genes: int,
    expected_samples: int,
) -> tuple[h5py.Dataset, h5py.Dataset, np.ndarray, np.ndarray]:
    required = [
        "data/expression",
        "meta/samples/geo_accession",
        "meta/genes/symbol",
        "meta/genes/ensembl_gene",
    ]
    missing = [name for name in required if name not in h5]
    if missing:
        raise RuntimeError(f"ARCHS4 HDF5 lacks datasets: {missing}")
    expression = h5["data/expression"]
    if expression.ndim != 2:
        raise RuntimeError(
            f"ARCHS4 data/expression must be 2D, observed {expression.shape}."
        )
    observed_shape = tuple(map(int, expression.shape))
    expected_shape = (expected_genes, expected_samples)
    if observed_shape != expected_shape:
        raise RuntimeError(
            f"ARCHS4 HDF5 shape differs: observed={observed_shape};"
            f"expected={expected_shape}"
        )
    accessions = h5["meta/samples/geo_accession"]
    symbols = decode_string_array(h5["meta/genes/symbol"][:])
    ensembl = decode_string_array(h5["meta/genes/ensembl_gene"][:])
    if len(accessions) != expected_samples:
        raise RuntimeError("ARCHS4 accession axis is not sample-aligned.")
    if len(symbols) != expected_genes or len(ensembl) != expected_genes:
        raise RuntimeError("ARCHS4 gene metadata axes are not expression-aligned.")
    return expression, accessions, symbols, ensembl


def clean_text(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return str(value).strip()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def write_tsv(
    frame: pd.DataFrame,
    path: Path,
    compression: str | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        compression=compression,
        na_rep="NA",
    )


def matrix_to_tsv_gz(
    matrix: np.ndarray,
    row_names: list[str],
    col_names: list[str],
    path: Path,
) -> None:
    if matrix.shape != (len(row_names), len(col_names)):
        raise RuntimeError(
            f"Matrix shape mismatch for {path.name}: "
            f"{matrix.shape} versus {(len(row_names), len(col_names))}"
        )
    frame = pd.DataFrame(matrix, index=row_names, columns=col_names)
    frame.index.name = "gene"
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", compression="gzip", float_format="%.10g")


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(
        path,
        sep="\t",
        compression="gzip" if path.suffix == ".gz" else None,
        low_memory=False,
    )


def parquet_columns(pf: pq.ParquetFile) -> tuple[list[str], list[str]]:
    names = list(pf.schema_arrow.names)
    index_columns: list[str] = []
    raw = (pf.schema_arrow.metadata or {}).get(b"pandas")
    if raw:
        try:
            metadata = json.loads(raw.decode())
            for item in metadata.get("index_columns", []):
                if isinstance(item, str):
                    index_columns.append(item)
        except Exception:
            pass
    for name in names:
        if name.startswith("__index_level_") and name not in index_columns:
            index_columns.append(name)
    data_columns = [name for name in names if name not in index_columns]
    return data_columns, index_columns


def validate_unique_ids(frame: pd.DataFrame, name: str) -> None:
    if "sample_id" not in frame.columns:
        raise RuntimeError(f"{name} lacks sample_id.")
    if frame["sample_id"].isna().any():
        raise RuntimeError(f"{name} contains missing sample_id.")
    duplicated = frame["sample_id"].astype(str).duplicated(keep=False)
    if duplicated.any():
        ids = sorted(frame.loc[duplicated, "sample_id"].astype(str).unique())
        raise RuntimeError(f"{name} contains duplicated sample IDs: {ids[:20]}")


def validate_manifests(
    combined: pd.DataFrame,
    canonical: pd.DataFrame,
    anchor: pd.DataFrame,
    exclusions: pd.DataFrame,
    expected_combined: int,
    expected_canonical: int,
    expected_anchor: int,
) -> None:
    required = {
        "sample_id",
        "primary_gse_id",
        "provisional_cohort",
        "analysis_cancer_class",
    }
    for label, frame in [
        ("combined", combined),
        ("canonical", canonical),
        ("cross-sectional", anchor),
    ]:
        missing = required - set(frame.columns)
        if missing:
            raise RuntimeError(f"{label} manifest missing columns: {sorted(missing)}")
        validate_unique_ids(frame, label)

    observed = {
        "combined": len(combined),
        "canonical": len(canonical),
        "anchor": len(anchor),
    }
    expected = {
        "combined": expected_combined,
        "canonical": expected_canonical,
        "anchor": expected_anchor,
    }
    if observed != expected:
        raise RuntimeError(
            "Final-lock manifest sizes differ from the locked contract.\n"
            f"Expected: {expected}\nObserved: {observed}"
        )

    combined_ids = set(combined["sample_id"].astype(str))
    canonical_ids = set(canonical["sample_id"].astype(str))
    anchor_ids = set(anchor["sample_id"].astype(str))

    if not canonical_ids.issubset(combined_ids):
        raise RuntimeError("Canonical manifest is not a subset of combined primary.")
    if not anchor_ids.issubset(canonical_ids):
        raise RuntimeError("Cross-sectional anchor is not a subset of canonical.")

    if not exclusions.empty and "sample_id" in exclusions.columns:
        excluded_ids = set(exclusions["sample_id"].dropna().astype(str))
        overlap = sorted(combined_ids & excluded_ids)
        if overlap:
            raise RuntimeError(
                f"Final primary manifest overlaps the exclusion ledger: {overlap[:20]}"
            )

    for label, frame in [("canonical", canonical), ("anchor", anchor)]:
        for column in [
            "primary_gse_id",
            "provisional_cohort",
            "analysis_cancer_class",
        ]:
            blank = frame[column].fillna("").astype(str).str.strip().eq("")
            if blank.any():
                raise RuntimeError(
                    f"{label} has {int(blank.sum())} blank values in {column}."
                )


def load_model_genes(path: Path) -> tuple[list[str], pd.DataFrame]:
    if not path.exists():
        raise FileNotFoundError(f"Missing official model-weight table: {path}")
    frame = pd.read_csv(path, sep="\t", low_memory=False)

    if "model_gene" in frame.columns:
        gene_column = "model_gene"
    elif "gene" in frame.columns:
        gene_column = "gene"
    else:
        raise RuntimeError("Model-weight table lacks gene/model_gene column.")

    missing_components = set(COMPONENTS) - set(frame.columns)
    if missing_components:
        raise RuntimeError(
            f"Model-weight table missing components: {sorted(missing_components)}"
        )

    frame = frame.copy()
    frame["model_gene"] = (
        frame[gene_column].fillna("").astype(str).str.strip().str.upper()
    )
    if frame["model_gene"].eq("").any():
        raise RuntimeError("Model-weight table contains blank gene names.")

    for component in COMPONENTS:
        frame[component] = pd.to_numeric(
            frame[component],
            errors="coerce",
        ).fillna(0.0)

    collapsed = (
        frame.groupby("model_gene", sort=True, as_index=False)[COMPONENTS]
        .sum()
    )
    genes = collapsed["model_gene"].tolist()
    if len(genes) != EXPECTED_MODEL_GENES:
        raise RuntimeError(
            f"Expected {EXPECTED_MODEL_GENES} unique model genes, observed {len(genes)}."
        )
    return genes, collapsed


def gene_axis_from_parquet(path: Path) -> tuple[pd.Index, list[str], list[str]]:
    pf = pq.ParquetFile(path)
    data_columns, index_columns = parquet_columns(pf)
    if not data_columns:
        raise RuntimeError(f"No expression columns in {path}")
    probe = pd.read_parquet(path, columns=[data_columns[0]])
    axis = pd.Index(probe.index.astype(str).str.strip())
    return axis, data_columns, index_columns


def axis_digest(axis: Iterable[str]) -> str:
    hasher = hashlib.sha256()
    for value in axis:
        hasher.update(str(value).encode("utf-8"))
        hasher.update(b"\0")
    return hasher.hexdigest()


def build_gene_maps(
    axis: pd.Index,
    model_genes: list[str],
) -> dict[str, Any]:
    normalized = axis.astype(str).str.strip().str.upper()
    valid_mask = ~normalized.isin(["", "NA", "NAN", "NONE", "NULL"])
    valid_positions = np.flatnonzero(valid_mask)
    valid_names = normalized[valid_mask]

    unique_genes = sorted(set(valid_names))
    unique_index = {gene: i for i, gene in enumerate(unique_genes)}
    row_to_unique = np.asarray(
        [unique_index[gene] for gene in valid_names],
        dtype=np.int64,
    )

    if len(model_genes) != len(set(model_genes)):
        raise RuntimeError("Official model-gene list contains duplicated symbols.")

    axis_gene_set = set(valid_names)
    source_gene_to_model_index: dict[str, int] = {}
    resolution_rows: list[dict[str, Any]] = []
    available_model: list[str] = []

    for model_index, model_gene in enumerate(model_genes):
        alias_gene = MODEL_GENE_ALIASES.get(model_gene, "")
        exact_present = model_gene in axis_gene_set
        alias_present = bool(alias_gene) and alias_gene in axis_gene_set

        # Never silently choose between two distinct rows carrying the old and
        # current symbols. Such a matrix requires explicit manual resolution.
        if exact_present and alias_present and alias_gene != model_gene:
            raise RuntimeError(
                "Ambiguous model-gene resolution: both the frozen symbol "
                f"{model_gene} and its current alias {alias_gene} occur in "
                "the ARCHS4 axis."
            )

        if exact_present:
            source_gene = model_gene
            resolution = "exact"
        elif alias_present:
            source_gene = alias_gene
            resolution = "alias_resolved"
        else:
            source_gene = ""
            resolution = "unresolved"

        if source_gene:
            previous = source_gene_to_model_index.get(source_gene)
            if previous is not None and previous != model_index:
                raise RuntimeError(
                    f"ARCHS4 source symbol {source_gene} resolves to more than "
                    "one frozen model row."
                )
            source_gene_to_model_index[source_gene] = model_index
            available_model.append(model_gene)

        resolution_rows.append(
            {
                "model_gene": model_gene,
                "ARCHS4_source_gene": source_gene or "NA",
                "resolution": resolution,
                "legacy_alias": alias_gene or "NA",
                "frozen_symbol_present": exact_present,
                "alias_symbol_present": alias_present,
                "available_in_ARCHS4_gene_axis": bool(source_gene),
            }
        )

    model_row_positions: list[int] = []
    model_row_indices: list[int] = []
    for position, gene in zip(valid_positions, valid_names):
        model_index = source_gene_to_model_index.get(gene)
        if model_index is not None:
            model_row_positions.append(int(position))
            model_row_indices.append(model_index)

    duplicate_summary = (
        pd.Series(valid_names)
        .value_counts()
        .rename_axis("gene")
        .reset_index(name="n_rows")
    )
    duplicate_summary = duplicate_summary[
        duplicate_summary["n_rows"] > 1
    ].sort_values(["n_rows", "gene"], ascending=[False, True])

    resolution_table = pd.DataFrame(resolution_rows)

    return {
        "normalized_axis": normalized,
        "valid_positions": valid_positions,
        "valid_names": list(valid_names),
        "unique_genes": unique_genes,
        "row_to_unique": row_to_unique,
        "model_row_positions": np.asarray(model_row_positions, dtype=np.int64),
        "model_row_indices": np.asarray(model_row_indices, dtype=np.int64),
        "available_model_genes": available_model,
        "model_gene_resolution": resolution_table,
        "duplicate_summary": duplicate_summary,
    }


def collapse_batch(
    values: np.ndarray,
    valid_positions: np.ndarray,
    row_to_unique: np.ndarray,
    n_unique: int,
) -> np.ndarray:
    selected = values[valid_positions, :]
    collapsed = np.zeros((n_unique, values.shape[1]), dtype=np.float64)
    np.add.at(collapsed, row_to_unique, selected)
    return collapsed


def clean_expression(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    finite = np.isfinite(values)
    negative = (finite & (values < 0)).sum(axis=0)
    nonfinite = (~finite).sum(axis=0)
    clean = np.where(finite & (values > 0), values, 0.0)
    return clean, negative, nonfinite


def quantile75_fraction(counts: np.ndarray, library_size: np.ndarray) -> np.ndarray:
    q75 = np.quantile(counts, 0.75, axis=0)
    return np.divide(
        q75,
        library_size,
        out=np.zeros_like(q75, dtype=np.float64),
        where=library_size > 0,
    )


def select_tmm_reference(
    sample_ids: list[str],
    f75: np.ndarray,
    library_sizes: np.ndarray,
) -> int:
    positive = np.isfinite(f75) & (f75 > 0)
    if positive.any() and np.median(f75[positive]) >= 1e-20:
        target = float(np.mean(f75[positive]))
        distance = np.where(positive, np.abs(f75 - target), np.inf)
        return int(np.argmin(distance))

    positive_lib = library_sizes[library_sizes > 0]
    if len(positive_lib) == 0:
        raise RuntimeError("No positive library sizes available for TMM reference.")
    target = math.sqrt(float(positive_lib.min()) * float(positive_lib.max()))
    return int(np.argmin(np.abs(library_sizes - target)))


def rank_keep_mask(values: np.ndarray, trim: float) -> np.ndarray:
    n = len(values)
    if n == 0:
        return np.zeros(0, dtype=bool)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(n, dtype=np.int64)
    ranks[order] = np.arange(1, n + 1)
    lo = math.floor(n * trim) + 1
    hi = n + 1 - lo
    return (ranks >= lo) & (ranks <= hi)


def tmm_factor(
    sample: np.ndarray,
    reference: np.ndarray,
    sample_library: float,
    reference_library: float,
    logratio_trim: float,
    sum_trim: float,
) -> tuple[float, int]:
    if sample_library <= 0 or reference_library <= 0:
        return 1.0, 0

    keep = (sample > 0) & (reference > 0)
    if int(keep.sum()) < 10:
        return 1.0, int(keep.sum())

    y = sample[keep].astype(np.float64)
    r = reference[keep].astype(np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        log_r = np.log2((y / sample_library) / (r / reference_library))
        abs_e = 0.5 * (
            np.log2(y / sample_library)
            + np.log2(r / reference_library)
        )
        variance = (
            (sample_library - y) / (sample_library * y)
            + (reference_library - r) / (reference_library * r)
        )

    finite = (
        np.isfinite(log_r)
        & np.isfinite(abs_e)
        & np.isfinite(variance)
        & (variance > 0)
    )
    log_r = log_r[finite]
    abs_e = abs_e[finite]
    variance = variance[finite]

    if len(log_r) < 10:
        return 1.0, len(log_r)

    retained = (
        rank_keep_mask(log_r, logratio_trim)
        & rank_keep_mask(abs_e, sum_trim)
    )
    if int(retained.sum()) < 10:
        return 1.0, int(retained.sum())

    weights = 1.0 / variance[retained]
    mean_log_r = np.sum(log_r[retained] * weights) / np.sum(weights)
    return float(2.0 ** mean_log_r), int(retained.sum())


def geometric_mean_positive(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    positive = values[np.isfinite(values) & (values > 0)]
    if len(positive) == 0:
        return 1.0
    return float(np.exp(np.mean(np.log(positive))))


def make_balanced_anchor_weights(anchor: pd.DataFrame) -> pd.DataFrame:
    frame = anchor[
        ["sample_id", "primary_gse_id", "analysis_cancer_class"]
    ].copy()
    frame["sample_id"] = frame["sample_id"].astype(str)
    frame["primary_gse_id"] = frame["primary_gse_id"].astype(str)
    frame["analysis_cancer_class"] = frame["analysis_cancer_class"].astype(str)

    n_cancers = frame["analysis_cancer_class"].nunique()
    study_counts = (
        frame.groupby("analysis_cancer_class")["primary_gse_id"]
        .nunique()
        .rename("n_studies_in_cancer")
    )
    sample_counts = (
        frame.groupby(["analysis_cancer_class", "primary_gse_id"])
        .size()
        .rename("n_samples_in_cancer_study")
    )

    frame = frame.join(study_counts, on="analysis_cancer_class")
    frame = frame.join(
        sample_counts,
        on=["analysis_cancer_class", "primary_gse_id"],
    )

    frame["pan_cancer_study_balanced_weight"] = (
        1.0
        / n_cancers
        / frame["n_studies_in_cancer"]
        / frame["n_samples_in_cancer_study"]
    )
    frame["within_cancer_study_balanced_weight"] = (
        1.0
        / frame["n_studies_in_cancer"]
        / frame["n_samples_in_cancer_study"]
    )

    total = frame["pan_cancer_study_balanced_weight"].sum()
    if not np.isclose(total, 1.0, atol=1e-12):
        raise RuntimeError(f"Pan-cancer anchor weights sum to {total}, not 1.")

    within = (
        frame.groupby("analysis_cancer_class")[
            "within_cancer_study_balanced_weight"
        ]
        .sum()
    )
    if not np.allclose(within.to_numpy(), 1.0, atol=1e-12):
        raise RuntimeError("Within-cancer study-balanced weights do not sum to 1.")

    return frame.sort_values(
        ["analysis_cancer_class", "primary_gse_id", "sample_id"]
    )


def run_self_tests() -> None:
    # Gene-axis normalization, invalid-name removal and duplicate mapping.
    maps = build_gene_maps(
        pd.Index(["TP53", "egfr ", "NA", "", "TP53"]),
        ["EGFR", "TP53"],
    )
    assert maps["valid_positions"].tolist() == [0, 1, 4]
    assert maps["valid_names"] == ["TP53", "EGFR", "TP53"]
    assert maps["available_model_genes"] == ["EGFR", "TP53"]
    assert (
        maps["duplicate_summary"]
        .set_index("gene")
        .loc["TP53", "n_rows"]
        == 2
    )

    # Legacy frozen-model symbols resolve from current ARCHS4 symbols while
    # retaining the original frozen model-row name.
    alias_maps = build_gene_maps(
        pd.Index(["DIPK2B", "TP53"]),
        ["CXORF36", "TP53"],
    )
    assert alias_maps["available_model_genes"] == ["CXORF36", "TP53"]
    alias_resolution = alias_maps["model_gene_resolution"].set_index("model_gene")
    assert alias_resolution.loc["CXORF36", "ARCHS4_source_gene"] == "DIPK2B"
    assert alias_resolution.loc["CXORF36", "resolution"] == "alias_resolved"
    assert alias_maps["model_row_positions"].tolist() == [0, 1]
    assert alias_maps["model_row_indices"].tolist() == [0, 1]

    # The extractor must not silently choose when both old and current symbols
    # occur in the same axis.
    ambiguity_failed = False
    try:
        build_gene_maps(
            pd.Index(["CXORF36", "DIPK2B"]),
            ["CXORF36"],
        )
    except RuntimeError as exc:
        ambiguity_failed = "Ambiguous model-gene resolution" in str(exc)
    assert ambiguity_failed

    # Duplicate-gene collapse.
    values = np.array(
        [
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ]
    )
    collapsed = collapse_batch(
        values,
        np.array([0, 1, 2]),
        np.array([0, 0, 1]),
        2,
    )
    assert np.allclose(collapsed, np.array([[4.0, 6.0], [5.0, 6.0]]))

    # TMM identity.
    x = np.arange(1, 101, dtype=float)
    factor, retained = tmm_factor(
        x,
        x.copy(),
        float(x.sum()),
        float(x.sum()),
        0.30,
        0.05,
    )
    assert np.isclose(factor, 1.0, atol=1e-12)
    assert retained >= 10

    # TMM global scaling is removed by library-size normalization.
    factor2, _ = tmm_factor(
        x * 2,
        x,
        float((x * 2).sum()),
        float(x.sum()),
        0.30,
        0.05,
    )
    assert np.isclose(factor2, 1.0, atol=1e-12)

    # Balanced weights.
    anchor = pd.DataFrame(
        {
            "sample_id": ["A", "B", "C", "D"],
            "primary_gse_id": ["S1", "S1", "S2", "S3"],
            "analysis_cancer_class": ["C1", "C1", "C1", "C2"],
        }
    )
    weights = make_balanced_anchor_weights(anchor)
    assert np.isclose(weights["pan_cancer_study_balanced_weight"].sum(), 1.0)
    by_cancer = weights.groupby("analysis_cancer_class")[
        "pan_cancer_study_balanced_weight"
    ].sum()
    assert np.allclose(by_cancer.to_numpy(), 0.5)

    # HDF5 fancy-index helpers must preserve arbitrary requested order.
    import tempfile
    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_h5 = Path(temporary_directory) / "order_test.h5"
        expected_matrix = np.arange(24, dtype=np.int64).reshape(4, 6)
        with h5py.File(temporary_h5, "w") as handle:
            matrix_dataset = handle.create_dataset("matrix", data=expected_matrix)
            vector_dataset = handle.create_dataset(
                "vector",
                data=np.asarray([b"A", b"B", b"C", b"D", b"E", b"F"]),
            )
            requested = np.asarray([4, 1, 5, 0], dtype=np.int64)
            observed_matrix = h5_take_columns_preserve_order(
                matrix_dataset, requested
            )
            observed_vector = h5_take_vector_preserve_order(
                vector_dataset, requested
            )
            assert np.array_equal(
                observed_matrix,
                expected_matrix[:, requested],
            )
            assert decode_string_array(observed_vector).tolist() == [
                "E", "B", "F", "A"
            ]

    print("ALL PREPROJECTION NORMALIZATION V1.0.1 SELF-TESTS PASSED")





def main() -> int:
    args = parse_args()
    if args.self_test:
        run_self_tests()
        # HDF5 decoding and index conventions specific to this revision.
        assert decode_string_array(np.array([b"GSM1", np.bytes_(b"GSM2")])).tolist() == [
            "GSM1",
            "GSM2",
        ]
        print("ALL ARCHS4 STAGE 09 V1.0.1 DIRECT-HDF5 TESTS PASSED")
        return 0

    if args.sample_batch_size < 1 or args.sample_batch_size > 128:
        raise ValueError("--sample-batch-size must be between 1 and 128.")
    if not (0 <= args.tmm_logratio_trim < 0.5):
        raise ValueError("--tmm-logratio-trim must be in [0, 0.5).")
    if not (0 <= args.tmm_sum_trim < 0.5):
        raise ValueError("--tmm-sum-trim must be in [0, 0.5).")

    started = utc_now()
    started_epoch = time.time()
    locked_dir = args.locked_dir.resolve()
    h5_path = args.h5.resolve()
    model_weights_path = args.model_weights.resolve()
    output_dir = args.output_dir.resolve()

    if not locked_dir.is_dir():
        raise FileNotFoundError(f"Stage 08 lock directory is missing: {locked_dir}")
    require_nonempty(h5_path, "ARCHS4 HDF5")
    require_nonempty(model_weights_path, "official model-weight table")

    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            f"Output directory is non-empty: {output_dir}. "
            "Stage outputs are immutable; use a new versioned directory."
        )
    staging = output_dir.parent / f"{output_dir.name}.staging.{os.getpid()}"
    if staging.exists():
        raise FileExistsError(f"Staging directory already exists: {staging}")
    staging.mkdir(parents=True, exist_ok=False)

    try:
        log("Validating the authoritative Stage 08 v1.0.1 lock.")
        combined, canonical, anchor, stage08_status, stage08_qc = (
            validate_stage08_inputs(
                locked_dir,
                args.expected_combined,
                args.expected_canonical,
                args.expected_cross_sectional,
                args.expected_cancers,
            )
        )
        model_genes, model_weight_table = load_model_genes(model_weights_path)

        sample_ids = canonical["sample_id"].astype(str).tolist()
        column_indices = canonical["expression_column_index"].to_numpy(
            dtype=np.int64
        )
        sample_to_position = {
            sample_id: position for position, sample_id in enumerate(sample_ids)
        }
        n_samples = len(canonical)
        n_model = len(model_genes)

        log(f"Opening ARCHS4 HDF5: {h5_path}")
        with h5py.File(h5_path, "r") as h5:
            expression, accessions, symbols, ensembl = validate_h5(
                h5,
                args.expected_h5_genes,
                args.expected_h5_samples,
            )
            if column_indices.min() < 0 or column_indices.max() >= expression.shape[1]:
                raise RuntimeError("At least one locked HDF5 column is out of bounds.")
            if len(np.unique(column_indices)) != len(column_indices):
                raise RuntimeError("Canonical HDF5 indices are not unique.")

            axis = pd.Index(symbols.astype(str))
            gene_maps = build_gene_maps(axis, model_genes)
            available_model_genes = gene_maps["available_model_genes"]
            model_gene_resolution = gene_maps["model_gene_resolution"].copy()
            alias_resolved_model_genes = int(
                model_gene_resolution["resolution"].eq("alias_resolved").sum()
            )
            exact_model_genes = int(
                model_gene_resolution["resolution"].eq("exact").sum()
            )
            unresolved_model_genes = int(
                model_gene_resolution["resolution"].eq("unresolved").sum()
            )
            if len(available_model_genes) != args.expected_available_model_genes:
                raise RuntimeError(
                    "Resolved model-gene count differs from contract: "
                    f"observed={len(available_model_genes)};"
                    f"expected={args.expected_available_model_genes}"
                )
            if alias_resolved_model_genes != args.expected_alias_resolved_model_genes:
                raise RuntimeError(
                    "Alias-resolved model-gene count differs from contract: "
                    f"observed={alias_resolved_model_genes};"
                    f"expected={args.expected_alias_resolved_model_genes}"
                )
            if unresolved_model_genes != 0:
                unresolved = model_gene_resolution.loc[
                    model_gene_resolution["resolution"].eq("unresolved"),
                    "model_gene",
                ].tolist()
                raise RuntimeError(f"Unresolved official model genes: {unresolved}")

            n_unique = len(gene_maps["unique_genes"])
            model_counts = np.zeros((n_model, n_samples), dtype=np.float64)
            library_all_rows = np.zeros(n_samples, dtype=np.float64)
            library_valid_symbols = np.zeros(n_samples, dtype=np.float64)
            detected_unique_genes = np.zeros(n_samples, dtype=np.int64)
            negative_counts = np.zeros(n_samples, dtype=np.int64)
            nonfinite_counts = np.zeros(n_samples, dtype=np.int64)
            f75 = np.zeros(n_samples, dtype=np.float64)
            observed_accessions = np.empty(n_samples, dtype=object)

            model_row_positions = gene_maps["model_row_positions"]
            model_row_indices = gene_maps["model_row_indices"]

            log(
                "PASS 1: direct-HDF5 extraction of canonical model counts and "
                "normalization metrics."
            )
            for start in range(0, n_samples, args.sample_batch_size):
                end = min(start + args.sample_batch_size, n_samples)
                cols = column_indices[start:end]
                raw = h5_take_columns_preserve_order(expression, cols)
                observed = decode_string_array(h5_take_vector_preserve_order(accessions, cols))
                expected = np.asarray(sample_ids[start:end], dtype=object)
                observed_accessions[start:end] = observed
                if not np.array_equal(observed, expected):
                    mismatch = np.flatnonzero(observed != expected)
                    examples = [
                        {
                            "position": int(start + i),
                            "expected": str(expected[i]),
                            "observed": str(observed[i]),
                            "h5_column": int(cols[i]),
                        }
                        for i in mismatch[:20]
                    ]
                    raise RuntimeError(
                        f"HDF5 accession alignment failed: {examples}"
                    )

                values, negative, nonfinite = clean_expression(
                    raw.astype(np.float64, copy=False)
                )
                collapsed = collapse_batch(
                    values,
                    gene_maps["valid_positions"],
                    gene_maps["row_to_unique"],
                    n_unique,
                )
                library_all_rows[start:end] = values.sum(axis=0)
                valid_library = collapsed.sum(axis=0)
                library_valid_symbols[start:end] = valid_library
                detected_unique_genes[start:end] = (collapsed > 0).sum(axis=0)
                negative_counts[start:end] = negative
                nonfinite_counts[start:end] = nonfinite
                f75[start:end] = quantile75_fraction(
                    collapsed, valid_library
                )

                batch_model = np.zeros((n_model, end - start), dtype=np.float64)
                selected = values[model_row_positions, :]
                np.add.at(batch_model, model_row_indices, selected)
                model_counts[:, start:end] = batch_model

                log(f"PASS 1: {end:,}/{n_samples:,} canonical samples")

            if np.any(library_valid_symbols <= 0):
                bad = [
                    sample_ids[i]
                    for i in np.flatnonzero(library_valid_symbols <= 0)[:20]
                ]
                raise RuntimeError(
                    f"Non-positive valid-symbol library denominators: {bad}"
                )
            if int(negative_counts.sum()) != 0:
                raise RuntimeError(
                    f"Canonical HDF5 columns contain "
                    f"{int(negative_counts.sum())} negative values."
                )
            if int(nonfinite_counts.sum()) != 0:
                raise RuntimeError(
                    f"Canonical HDF5 columns contain "
                    f"{int(nonfinite_counts.sum())} non-finite values."
                )

            frozen_library = pd.to_numeric(
                canonical["library_sum_archs4_pseudocounts"],
                errors="coerce",
            ).to_numpy(dtype=np.float64)
            library_replay_difference = library_all_rows - frozen_library
            max_abs_library_replay_difference = float(
                np.nanmax(np.abs(library_replay_difference))
            )
            if (
                not np.isfinite(max_abs_library_replay_difference)
                or max_abs_library_replay_difference > 0.5
            ):
                raise RuntimeError(
                    "Direct-HDF5 library sums do not replay frozen Stage 08 "
                    f"metrics; max absolute difference="
                    f"{max_abs_library_replay_difference}"
                )

            reference_position = select_tmm_reference(
                sample_ids, f75, library_valid_symbols
            )
            reference_sample = sample_ids[reference_position]
            reference_col = np.asarray(
                [column_indices[reference_position]], dtype=np.int64
            )
            reference_raw = h5_take_columns_preserve_order(expression, reference_col)
            reference_values, ref_negative, ref_nonfinite = clean_expression(
                reference_raw.astype(np.float64, copy=False)
            )
            if int(ref_negative.sum()) or int(ref_nonfinite.sum()):
                raise RuntimeError("TMM reference contains invalid values.")
            reference_collapsed = collapse_batch(
                reference_values,
                gene_maps["valid_positions"],
                gene_maps["row_to_unique"],
                n_unique,
            )[:, 0]
            reference_library = float(reference_collapsed.sum())

            tmm_factors = np.ones(n_samples, dtype=np.float64)
            tmm_retained_genes = np.zeros(n_samples, dtype=np.int64)
            log(
                f"PASS 2: estimating TMM factors against {reference_sample}."
            )
            for start in range(0, n_samples, args.sample_batch_size):
                end = min(start + args.sample_batch_size, n_samples)
                cols = column_indices[start:end]
                raw = h5_take_columns_preserve_order(expression, cols)
                values, negative, nonfinite = clean_expression(
                    raw.astype(np.float64, copy=False)
                )
                if int(negative.sum()) or int(nonfinite.sum()):
                    raise RuntimeError(
                        "Invalid expression appeared during TMM replay."
                    )
                collapsed = collapse_batch(
                    values,
                    gene_maps["valid_positions"],
                    gene_maps["row_to_unique"],
                    n_unique,
                )
                for local_index in range(end - start):
                    global_index = start + local_index
                    factor, retained = tmm_factor(
                        collapsed[:, local_index],
                        reference_collapsed,
                        library_valid_symbols[global_index],
                        reference_library,
                        args.tmm_logratio_trim,
                        args.tmm_sum_trim,
                    )
                    tmm_factors[global_index] = factor
                    tmm_retained_genes[global_index] = retained
                log(f"PASS 2: {end:,}/{n_samples:,} canonical samples")

            h5_inventory = {
                "h5_path": str(h5_path),
                "h5_basename": h5_path.name,
                "h5_size_bytes": int(h5_path.stat().st_size),
                "expression_rows": int(expression.shape[0]),
                "expression_columns": int(expression.shape[1]),
                "expression_dtype": str(expression.dtype),
                "expression_chunks": str(expression.chunks),
                "expression_compression": expression.compression or "",
                "gene_symbol_axis_sha256": axis_digest(symbols),
                "ensembl_axis_sha256": axis_digest(ensembl),
                "nonblank_unique_gene_symbols": n_unique,
                "blank_or_missing_gene_symbol_rows": int(
                    len(symbols) - len(gene_maps["valid_positions"])
                ),
                "duplicated_gene_symbol_names": int(
                    len(gene_maps["duplicate_summary"])
                ),
            }

        # edgeR-style geometric centering.
        tmm_factors /= geometric_mean_positive(tmm_factors)

        uq_factors = f75.copy()
        positive_uq = np.isfinite(uq_factors) & (uq_factors > 0)
        if not positive_uq.any():
            raise RuntimeError("No positive upper-quartile fractions.")
        if not positive_uq.all():
            uq_factors[~positive_uq] = float(np.median(uq_factors[positive_uq]))
        uq_factors /= geometric_mean_positive(uq_factors)

        effective_library_tmm = library_valid_symbols * tmm_factors
        effective_library_uq = library_valid_symbols * uq_factors

        cpm = model_counts / library_valid_symbols[None, :] * 1e6
        tmm_cpm = model_counts / effective_library_tmm[None, :] * 1e6
        uq_cpm = model_counts / effective_library_uq[None, :] * 1e6

        model_detected = (model_counts > 0).sum(axis=0)
        excluded_symbol_mass = library_all_rows - library_valid_symbols
        excluded_symbol_mass_fraction = np.divide(
            excluded_symbol_mass,
            library_all_rows,
            out=np.zeros_like(excluded_symbol_mass),
            where=library_all_rows > 0,
        )

        sample_qc = canonical[
            [
                "sample_id",
                "primary_gse_id",
                "provisional_cohort",
                "analysis_cancer_class",
                "stage08_tier",
                "stage07_expression_quality_tier",
                "expression_column_index",
                "stage09_matrix_column_index",
            ]
        ].copy()
        sample_qc["h5_geo_accession_observed"] = observed_accessions
        sample_qc["h5_accession_match"] = (
            sample_qc["sample_id"].to_numpy(dtype=object)
            == observed_accessions
        )
        sample_qc["library_sum_all_HDF5_gene_rows"] = library_all_rows
        sample_qc["frozen_stage08_library_sum"] = frozen_library
        sample_qc["library_sum_replay_difference"] = library_replay_difference
        sample_qc["normalization_library_sum_valid_gene_symbols"] = (
            library_valid_symbols
        )
        sample_qc["excluded_blank_symbol_mass"] = excluded_symbol_mass
        sample_qc["excluded_blank_symbol_mass_fraction"] = (
            excluded_symbol_mass_fraction
        )
        sample_qc["detected_unique_gene_symbols_recomputed"] = (
            detected_unique_genes
        )
        sample_qc["model_genes_detected_recomputed"] = model_detected
        sample_qc["model_genes_available"] = len(available_model_genes)
        sample_qc["model_genes_official"] = len(model_genes)
        sample_qc["upper_quartile_fraction"] = f75
        sample_qc["tmm_reference_sample"] = reference_sample
        sample_qc["tmm_factor"] = tmm_factors
        sample_qc["tmm_genes_retained"] = tmm_retained_genes
        sample_qc["uq_factor"] = uq_factors
        sample_qc["effective_library_tmm"] = effective_library_tmm
        sample_qc["effective_library_uq"] = effective_library_uq
        sample_qc["negative_values"] = negative_counts
        sample_qc["nonfinite_values"] = nonfinite_counts

        method_rows: list[dict[str, Any]] = []
        for comparison, alternative in [
            ("CPM_vs_TMM_CPM", tmm_cpm),
            ("CPM_vs_UQ_CPM", uq_cpm),
        ]:
            correlations: list[float] = []
            for index in range(n_samples):
                x = np.log2(cpm[:, index] + 1)
                y = np.log2(alternative[:, index] + 1)
                if np.std(x) > 0 and np.std(y) > 0:
                    correlations.append(float(np.corrcoef(x, y)[0, 1]))
            if not correlations:
                raise RuntimeError(
                    f"No valid sample correlations for {comparison}."
                )
            method_rows.append(
                {
                    "comparison": comparison,
                    "n_samples_evaluated": len(correlations),
                    "minimum_sample_pearson_log2": float(np.min(correlations)),
                    "median_sample_pearson_log2": float(np.median(correlations)),
                    "mean_sample_pearson_log2": float(np.mean(correlations)),
                    "maximum_sample_pearson_log2": float(np.max(correlations)),
                }
            )
        method_agreement = pd.DataFrame(method_rows)

        anchor_weights = make_balanced_anchor_weights(anchor)
        anchor_positions = anchor["stage09_matrix_column_index"].to_numpy(
            dtype=np.int64
        )
        anchor_ids = anchor["sample_id"].astype(str).tolist()

        provenance_rows = [
            {
                "input": "stage08_status",
                "path": str(
                    locked_dir / "16_patient_canonical_resolution_status.json"
                ),
                "sha256": sha256_file(
                    locked_dir / "16_patient_canonical_resolution_status.json"
                ),
            },
            {
                "input": "stage08_qc_contract",
                "path": str(
                    locked_dir
                    / "15_patient_canonical_resolution_QC_contract.tsv"
                ),
                "sha256": sha256_file(
                    locked_dir
                    / "15_patient_canonical_resolution_QC_contract.tsv"
                ),
            },
            {
                "input": "stage08_canonical_manifest",
                "path": str(
                    locked_dir / "23_FINAL_LOCKED_canonical_library_atlas.tsv.gz"
                ),
                "sha256": sha256_file(
                    locked_dir / "23_FINAL_LOCKED_canonical_library_atlas.tsv.gz"
                ),
            },
            {
                "input": "stage08_anchor_manifest",
                "path": str(
                    locked_dir / "24_FINAL_LOCKED_cross_sectional_anchor.tsv.gz"
                ),
                "sha256": sha256_file(
                    locked_dir / "24_FINAL_LOCKED_cross_sectional_anchor.tsv.gz"
                ),
            },
            {
                "input": "official_model_weights",
                "path": str(model_weights_path),
                "sha256": sha256_file(model_weights_path),
            },
            {
                "input": "ARCHS4_HDF5",
                "path": str(h5_path),
                "sha256": "NOT_REHASHED_LARGE_FROZEN_SOURCE",
            },
        ]
        write_tsv(pd.DataFrame(provenance_rows), staging / "00_input_provenance.tsv")

        write_tsv(
            canonical,
            staging / "01_canonical_library_manifest.tsv.gz",
            compression="gzip",
        )
        write_tsv(
            anchor,
            staging / "02_cross_sectional_anchor_manifest.tsv.gz",
            compression="gzip",
        )
        write_tsv(
            anchor_weights,
            staging / "03_cross_sectional_balanced_anchor_weights.tsv",
        )
        write_tsv(
            sample_qc,
            staging / "04_sample_normalization_QC.tsv.gz",
            compression="gzip",
        )

        matrix_manifest = pd.DataFrame(
            [
                {
                    **h5_inventory,
                    "canonical_samples_requested": n_samples,
                    "canonical_samples_aligned": int(
                        sample_qc["h5_accession_match"].sum()
                    ),
                    "cross_sectional_anchor_samples": len(anchor),
                    "matrix_column_order": "exact_Stage08_locked_canonical_order",
                    "normalization_denominator": (
                        "sum_of_nonnegative_expression_over_nonblank_"
                        "case-normalized_gene_symbols"
                    ),
                }
            ]
        )
        write_tsv(matrix_manifest, staging / "05_matrix_extraction_manifest.tsv")
        write_tsv(
            gene_maps["duplicate_summary"],
            staging / "06_duplicate_gene_symbol_rows.tsv",
        )
        write_tsv(
            model_weight_table,
            staging / "07_official_model_gene_weights.tsv",
        )
        write_tsv(
            model_gene_resolution,
            staging / "08_model_gene_availability.tsv",
        )

        matrix_to_tsv_gz(
            model_counts,
            model_genes,
            sample_ids,
            staging / "09_model_gene_ARCHS4_pseudocounts_canonical.tsv.gz",
        )
        matrix_to_tsv_gz(
            cpm,
            model_genes,
            sample_ids,
            staging / "10_model_gene_CPM_canonical.tsv.gz",
        )
        matrix_to_tsv_gz(
            tmm_cpm,
            model_genes,
            sample_ids,
            staging / "11_model_gene_TMM_CPM_canonical.tsv.gz",
        )
        matrix_to_tsv_gz(
            uq_cpm,
            model_genes,
            sample_ids,
            staging / "12_model_gene_UQ_CPM_canonical.tsv.gz",
        )
        matrix_to_tsv_gz(
            cpm[:, anchor_positions],
            model_genes,
            anchor_ids,
            staging / "13_model_gene_CPM_cross_sectional_anchor.tsv.gz",
        )
        matrix_to_tsv_gz(
            tmm_cpm[:, anchor_positions],
            model_genes,
            anchor_ids,
            staging / "14_model_gene_TMM_CPM_cross_sectional_anchor.tsv.gz",
        )
        matrix_to_tsv_gz(
            uq_cpm[:, anchor_positions],
            model_genes,
            anchor_ids,
            staging / "15_model_gene_UQ_CPM_cross_sectional_anchor.tsv.gz",
        )
        write_tsv(
            method_agreement,
            staging / "16_normalization_method_agreement.tsv",
        )

        checks = [
            (
                "stage08_final_lock_validated",
                stage08_status.get("audit_state") == "FINAL_LOCK_VALIDATED",
                f"audit_state={stage08_status.get('audit_state')}",
            ),
            (
                "stage08_qc_all_pass",
                set(stage08_qc["status"].astype(str)) == {"PASS"},
                f"rows={len(stage08_qc)}",
            ),
            (
                "combined_manifest_size",
                len(combined) == args.expected_combined,
                f"observed={len(combined)};expected={args.expected_combined}",
            ),
            (
                "canonical_manifest_size",
                len(canonical) == args.expected_canonical,
                f"observed={len(canonical)};expected={args.expected_canonical}",
            ),
            (
                "cross_sectional_manifest_size",
                len(anchor) == args.expected_cross_sectional,
                f"observed={len(anchor)};"
                f"expected={args.expected_cross_sectional}",
            ),
            (
                "cancer_class_count",
                canonical["analysis_cancer_class"].nunique()
                == args.expected_cancers,
                f"observed={canonical['analysis_cancer_class'].nunique()};"
                f"expected={args.expected_cancers}",
            ),
            (
                "canonical_HDF5_indices_unique",
                bool(
                    not canonical["expression_column_index"].duplicated().any()
                    and len(np.unique(column_indices)) == len(column_indices)
                ),
                (
                    f"min={column_indices.min()};max={column_indices.max()};"
                    f"strictly_increasing_in_locked_order="
                    f"{bool(np.all(np.diff(column_indices) > 0))}"
                ),
            ),
            (
                "all_HDF5_accessions_aligned",
                bool(sample_qc["h5_accession_match"].all()),
                f"aligned={int(sample_qc['h5_accession_match'].sum())};"
                f"expected={n_samples}",
            ),
            (
                "frozen_library_sums_exactly_replayed",
                max_abs_library_replay_difference <= 0.5,
                f"max_abs_difference={max_abs_library_replay_difference:.12g}",
            ),
            (
                "official_model_gene_count",
                len(model_genes) == EXPECTED_MODEL_GENES,
                f"observed={len(model_genes)};expected={EXPECTED_MODEL_GENES}",
            ),
            (
                "resolved_model_gene_count",
                len(available_model_genes)
                == args.expected_available_model_genes,
                f"observed={len(available_model_genes)};"
                f"expected={args.expected_available_model_genes}",
            ),
            (
                "alias_resolved_model_gene_count",
                alias_resolved_model_genes
                == args.expected_alias_resolved_model_genes,
                f"observed={alias_resolved_model_genes};"
                f"expected={args.expected_alias_resolved_model_genes}",
            ),
            (
                "zero_unresolved_model_genes",
                unresolved_model_genes == 0,
                f"unresolved={unresolved_model_genes}",
            ),
            (
                "positive_normalization_libraries",
                bool(np.all(library_valid_symbols > 0)),
                f"minimum={library_valid_symbols.min():.12g}",
            ),
            (
                "zero_negative_expression_values",
                int(negative_counts.sum()) == 0,
                f"negative_values={int(negative_counts.sum())}",
            ),
            (
                "zero_nonfinite_expression_values",
                int(nonfinite_counts.sum()) == 0,
                f"nonfinite_values={int(nonfinite_counts.sum())}",
            ),
            (
                "TMM_factors_finite_positive",
                bool(
                    np.all(np.isfinite(tmm_factors))
                    and np.all(tmm_factors > 0)
                ),
                f"min={tmm_factors.min():.12g};max={tmm_factors.max():.12g}",
            ),
            (
                "TMM_geometric_mean_one",
                bool(
                    np.isclose(
                        geometric_mean_positive(tmm_factors),
                        1.0,
                        atol=1e-12,
                    )
                ),
                "geometric_mean="
                f"{geometric_mean_positive(tmm_factors):.12g}",
            ),
            (
                "pan_cancer_anchor_weights_sum_one",
                bool(
                    np.isclose(
                        anchor_weights[
                            "pan_cancer_study_balanced_weight"
                        ].sum(),
                        1.0,
                        atol=1e-12,
                    )
                ),
                "sum="
                f"{anchor_weights['pan_cancer_study_balanced_weight'].sum():.12g}",
            ),
            (
                "all_output_matrices_finite_nonnegative",
                bool(
                    np.isfinite(cpm).all()
                    and np.isfinite(tmm_cpm).all()
                    and np.isfinite(uq_cpm).all()
                    and (cpm >= 0).all()
                    and (tmm_cpm >= 0).all()
                    and (uq_cpm >= 0).all()
                ),
                "CPM,TMM_CPM,UQ_CPM",
            ),
            (
                "anchor_is_exact_canonical_matrix_subset",
                bool(
                    np.array_equal(
                        np.asarray(anchor_ids, dtype=object),
                        np.asarray(sample_ids, dtype=object)[anchor_positions],
                    )
                ),
                f"anchor_samples={len(anchor_ids)}",
            ),
        ]
        qc = pd.DataFrame(
            [
                {
                    "check": check,
                    "status": "PASS" if passed else "FAIL",
                    "detail": detail,
                }
                for check, passed, detail in checks
            ]
        )
        write_tsv(
            qc, staging / "17_preprojection_normalization_QC_contract.tsv"
        )
        failed = qc.loc[qc["status"].eq("FAIL")]

        status = {
            "started_at_utc": started,
            "completed_at_utc": utc_now(),
            "runtime_seconds": round(time.time() - started_epoch, 3),
            "host": socket.gethostname(),
            "script": Path(__file__).name,
            "script_version": SCRIPT_VERSION,
            "source_stage08_dir": str(locked_dir),
            "source_stage08_script_version": stage08_status.get(
                "script_version"
            ),
            "source_stage08_status_sha256": sha256_file(
                locked_dir / "16_patient_canonical_resolution_status.json"
            ),
            "h5_path": str(h5_path),
            "h5_genes": h5_inventory["expression_rows"],
            "h5_samples": h5_inventory["expression_columns"],
            "combined_locked_samples": len(combined),
            "canonical_samples": len(canonical),
            "cross_sectional_samples": len(anchor),
            "gold_samples": int(
                (combined["stage08_tier"] == GOLD).sum()
            ),
            "silver_samples": int(
                (combined["stage08_tier"] == SILVER).sum()
            ),
            "cancer_classes": int(
                canonical["analysis_cancer_class"].nunique()
            ),
            "model_genes_official": len(model_genes),
            "model_genes_available": len(available_model_genes),
            "model_genes_exact": exact_model_genes,
            "model_genes_alias_resolved": alias_resolved_model_genes,
            "model_genes_unresolved": unresolved_model_genes,
            "unique_ARCHS4_gene_symbols": n_unique,
            "duplicated_gene_symbol_names": len(
                gene_maps["duplicate_summary"]
            ),
            "maximum_absolute_frozen_library_sum_replay_difference": (
                max_abs_library_replay_difference
            ),
            "maximum_excluded_blank_symbol_mass_fraction": float(
                np.max(excluded_symbol_mass_fraction)
            ),
            "median_excluded_blank_symbol_mass_fraction": float(
                np.median(excluded_symbol_mass_fraction)
            ),
            "tmm_reference_sample": reference_sample,
            "tmm_reference_cancer": canonical.loc[
                reference_position, "analysis_cancer_class"
            ],
            "tmm_reference_tier": canonical.loc[
                reference_position, "provisional_cohort"
            ],
            "tmm_factor_min": float(tmm_factors.min()),
            "tmm_factor_median": float(np.median(tmm_factors)),
            "tmm_factor_max": float(tmm_factors.max()),
            "expression_matrix_exported": failed.empty,
            "projection_run": False,
            "qc_failures": int(len(failed)),
            "status": "PASS" if failed.empty else "FAIL",
            "output_dir": str(output_dir),
        }
        atomic_json(
            staging / "18_preprojection_normalization_status.json", status
        )

        report_lines = [
            "ARCHS4 STAGE 09 LOCKED EXPRESSION AND NORMALIZATION V1.0.1",
            "=" * 88,
            "",
            f"Status: {status['status']}",
            f"Source Stage 08: {locked_dir}",
            f"ARCHS4 HDF5: {h5_path}",
            f"Combined locked samples: {len(combined):,}",
            f"Canonical libraries exported: {len(canonical):,}",
            f"Cross-sectional anchors exported: {len(anchor):,}",
            f"Cancer classes: {status['cancer_classes']:,}",
            f"Gold / Silver: {status['gold_samples']:,} / "
            f"{status['silver_samples']:,}",
            f"Model genes: {len(available_model_genes):,}/"
            f"{len(model_genes):,} resolved",
            f"Exact / alias / unresolved: {exact_model_genes:,} / "
            f"{alias_resolved_model_genes:,} / {unresolved_model_genes:,}",
            f"TMM reference: {reference_sample}",
            f"TMM factor range: {tmm_factors.min():.5g}–"
            f"{tmm_factors.max():.5g}",
            "",
            "NORMALIZATION DENOMINATOR",
            "-------------------------",
            (
                "CPM denominators use the sum of nonnegative ARCHS4 "
                "pseudocounts after excluding rows with blank/missing gene "
                "symbols and collapsing duplicated case-normalized symbols."
            ),
            (
                "The complete all-row HDF5 library sum was independently "
                "recomputed and required to reproduce the frozen Stage 08 "
                "library-size metric."
            ),
            "",
            "PRIMARY PROJECTION INPUT",
            "10_model_gene_CPM_canonical.tsv.gz",
            "",
            "NORMALIZATION SENSITIVITY INPUTS",
            "11_model_gene_TMM_CPM_canonical.tsv.gz",
            "12_model_gene_UQ_CPM_canonical.tsv.gz",
            "",
            "FIXED-ANCHOR INPUTS",
            "13_model_gene_CPM_cross_sectional_anchor.tsv.gz",
            "14_model_gene_TMM_CPM_cross_sectional_anchor.tsv.gz",
            "15_model_gene_UQ_CPM_cross_sectional_anchor.tsv.gz",
            "",
            "INTERPRETATION",
            (
                "The output matrices are normalized ARCHS4 Kallisto "
                "pseudocount-abundance matrices, not literal TPM."
            ),
            (
                "No SpatialEcoTyper projection was run in this stage. "
                "Pan-cancer and cancer-specific fixed-anchor projections "
                "must remain separate downstream score spaces."
            ),
            "",
            "QC CONTRACT",
            qc.to_string(index=False),
        ]
        (
            staging
            / "ARCHS4_preprojection_normalization_v1_0_1_report.txt"
        ).write_text("\n".join(report_lines) + "\n")

        if not failed.empty:
            print(failed.to_string(index=False), file=sys.stderr)
            raise RuntimeError(
                f"Stage 09 completed with {len(failed)} failed QC checks. "
                f"Inspect staging directory: {staging}"
            )

        if output_dir.exists():
            # An existing empty directory is permitted, but publication is
            # atomic and therefore requires removing it first.
            output_dir.rmdir()
        staging.replace(output_dir)

        print("=" * 96)
        print("ARCHS4 STAGE 09 LOCKED EXPRESSION/NORMALIZATION: PASS")
        print("=" * 96)
        print(f"Canonical libraries:       {len(canonical):,}")
        print(f"Cross-sectional anchors:   {len(anchor):,}")
        print(f"Cancer classes:            {status['cancer_classes']:,}")
        print(f"Resolved model genes:      {len(available_model_genes):,}")
        print(f"TMM reference sample:      {reference_sample}")
        print(f"Output:                    {output_dir}")
        return 0

    except Exception:
        print(
            f"ERROR: Stage 09 failed. Staging output retained at: {staging}",
            file=sys.stderr,
            flush=True,
        )
        raise


if __name__ == "__main__":
    raise SystemExit(main())
