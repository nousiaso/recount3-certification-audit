#!/usr/bin/env python3
"""
Audit direct ARCHS4 expression compatibility for the frozen Stage 05 primary set.

This stage reads the current ARCHS4 HDF5 expression matrix directly for every
Stage 05 primary-pending-expression-QC sample. It does not publish a final
Gold/Silver lock and does not export a final expression matrix.

The stage verifies exact GSM-to-HDF5-column alignment, expression validity,
gene-detection breadth, library complexity, SpatialEcoTyper model-feature
coverage, exact expression-profile duplicates, and robust technical outliers.

Output statuses are intentionally conservative:
  COMPATIBLE
  COMPATIBLE_WITH_WARNING
  REVIEW_TECHNICAL_OUTLIER
  INCOMPATIBLE_HARD_FAIL

Only objective structural/invalid-expression conditions are hard failures.
Distributional outliers are review flags, not automatic exclusions.
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import os
import shutil
import socket
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import h5py
import numpy as np
import pandas as pd


SCRIPT_VERSION = "1.0.0"

DEFAULT_PROJECT = Path(
    "/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6"
)
DEFAULT_H5 = DEFAULT_PROJECT / "00_download/human_gene_v2.latest.h5"
DEFAULT_STAGE05 = DEFAULT_PROJECT / "06_manual_resolution_v1_0_1"
DEFAULT_MANIFEST = DEFAULT_STAGE05 / "03_manual_resolution_manifest.parquet"
DEFAULT_STAGE05_STATUS = DEFAULT_STAGE05 / "12_manual_resolution_status.json"
DEFAULT_MODEL_WEIGHTS = Path(
    "/home/tpn890/orestis/ST_ECOTYPES/00_preflight_ARCHS4_solid_tumors/"
    "official_Bulk_SE_Recovery_W.tsv"
)
DEFAULT_OUTPUT = DEFAULT_PROJECT / "07_expression_compatibility_qc_v1"

EXPECTED_H5_GENES = 67186
EXPECTED_H5_SAMPLES = 1093742
EXPECTED_CANDIDATES = 220052
EXPECTED_PRIMARY = 9954
EXPECTED_PREVIOUS_PRIMARY = 4209
EXPECTED_NEW_PRIMARY = 5745
EXPECTED_MODEL_GENES = 499
EXPECTED_ALIAS_RESOLVED = 7

COMPONENTS = ["NonSE", *[f"SE{i}" for i in range(1, 10)]]

# Frozen model names that were renamed on the current ARCHS4 symbol axis.
MODEL_GENE_ALIASES = {
    "CXORF36": "DIPK2B",
    "FAM71A": "GARIN4",
    "GPR1": "CMKLR2",
    "H2AFY": "MACROH2A1",
    "MRVI1": "IRAG1",
    "RARRES3": "PLAAT4",
    "WARS": "WARS1",
}

MISSING_TEXT = {"", "na", "n/a", "nan", "none", "null", "unknown", "not available"}

PRIMARY_BUCKET = "PRIMARY_PENDING_EXPRESSION_QC"

REQUIRED_MANIFEST_COLUMNS = {
    "sample_id",
    "expression_column_index",
    "expression_column_number",
    "primary_gse_id",
    "previous_primary_status",
    "previous_lock_tier",
    "stage05_resolution_bucket",
    "stage05_resolved_cancer_class",
    "stage05_independence_group",
    "stage05_independence_review_required",
    "stage05_canonical_priority_hint",
    "stage05_primary_pending_expression_qc",
}


@dataclass(frozen=True)
class ModelMapping:
    model_genes: list[str]
    weights: np.ndarray
    component_names: list[str]
    model_row_positions: np.ndarray
    model_row_indices: np.ndarray
    resolution: pd.DataFrame
    total_abs_by_component: np.ndarray
    total_abs_overall: float


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def log(message: str) -> None:
    print(f"[{utc_now()}] {message}", flush=True)


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    return str(value).strip()


def normalize_symbol(value: Any) -> str:
    text = clean_text(value).upper()
    return "" if text.casefold() in MISSING_TEXT else text


def decode_string_array(values: np.ndarray) -> np.ndarray:
    result: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            result.append(value.decode("utf-8", errors="replace"))
        elif isinstance(value, np.bytes_):
            result.append(bytes(value).decode("utf-8", errors="replace"))
        else:
            result.append(clean_text(value))
    return np.asarray(result, dtype=object)


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def write_tsv(frame: pd.DataFrame, path: Path, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, compression=compression, na_rep="NA")


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            hasher.update(block)
    return hasher.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit ARCHS4 expression compatibility for Stage 05 primary samples."
    )
    parser.add_argument("--h5", type=Path, default=DEFAULT_H5)
    parser.add_argument("--stage05-manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--stage05-status", type=Path, default=DEFAULT_STAGE05_STATUS)
    parser.add_argument("--model-weights", type=Path, default=DEFAULT_MODEL_WEIGHTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--sample-batch-size", type=int, default=16)
    parser.add_argument("--expected-h5-genes", type=int, default=EXPECTED_H5_GENES)
    parser.add_argument("--expected-h5-samples", type=int, default=EXPECTED_H5_SAMPLES)
    parser.add_argument("--expected-candidates", type=int, default=EXPECTED_CANDIDATES)
    parser.add_argument("--expected-primary", type=int, default=EXPECTED_PRIMARY)
    parser.add_argument("--expected-previous-primary", type=int, default=EXPECTED_PREVIOUS_PRIMARY)
    parser.add_argument("--expected-new-primary", type=int, default=EXPECTED_NEW_PRIMARY)
    parser.add_argument("--expected-model-genes", type=int, default=EXPECTED_MODEL_GENES)
    parser.add_argument("--expected-alias-resolved", type=int, default=EXPECTED_ALIAS_RESOLVED)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def load_stage05_inputs(
    manifest_path: Path,
    status_path: Path,
    expected_candidates: int,
    expected_primary: int,
    expected_previous: int,
    expected_new: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not manifest_path.exists() or manifest_path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing Stage 05 manifest: {manifest_path}")
    if not status_path.exists() or status_path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing Stage 05 status: {status_path}")

    status = json.loads(status_path.read_text())
    if status.get("status") != "PASS":
        raise RuntimeError(f"Stage 05 status is not PASS: {status}")
    if int(status.get("qc_failures", -1)) != 0:
        raise RuntimeError(f"Stage 05 has QC failures: {status}")
    if bool(status.get("expression_extracted", True)):
        raise RuntimeError("Stage 05 unexpectedly reports expression_extracted=true.")
    if bool(status.get("final_lock_created", True)):
        raise RuntimeError("Stage 05 unexpectedly reports final_lock_created=true.")

    full = pd.read_parquet(manifest_path)
    missing = sorted(REQUIRED_MANIFEST_COLUMNS - set(full.columns))
    if missing:
        raise RuntimeError(f"Stage 05 manifest lacks required columns: {missing}")
    if len(full) != expected_candidates:
        raise RuntimeError(
            f"Expected {expected_candidates:,} Stage 05 rows, observed {len(full):,}."
        )

    primary_mask = (
        full["stage05_resolution_bucket"].astype(str).eq(PRIMARY_BUCKET)
        & full["stage05_primary_pending_expression_qc"].fillna(False).astype(bool)
    )
    primary = full.loc[primary_mask].copy()

    if len(primary) != expected_primary:
        raise RuntimeError(
            f"Expected {expected_primary:,} primary rows, observed {len(primary):,}."
        )

    primary["sample_id"] = primary["sample_id"].fillna("").astype(str).str.strip()
    primary["primary_gse_id"] = primary["primary_gse_id"].fillna("").astype(str).str.strip()
    primary["expression_column_index"] = pd.to_numeric(
        primary["expression_column_index"], errors="raise"
    ).astype(np.int64)
    primary["expression_column_number"] = pd.to_numeric(
        primary["expression_column_number"], errors="raise"
    ).astype(np.int64)

    if primary["sample_id"].eq("").any():
        raise RuntimeError("Primary manifest contains blank sample IDs.")
    if primary["sample_id"].duplicated().any():
        dup = primary.loc[primary["sample_id"].duplicated(False), "sample_id"].tolist()
        raise RuntimeError(f"Duplicated primary sample IDs: {dup[:20]}")
    if primary["expression_column_index"].duplicated().any():
        dup = primary.loc[
            primary["expression_column_index"].duplicated(False),
            ["sample_id", "expression_column_index"],
        ]
        raise RuntimeError(f"Duplicated HDF5 expression indices:\n{dup.head(20)}")
    if not (
        primary["expression_column_number"].to_numpy()
        == primary["expression_column_index"].to_numpy() + 1
    ).all():
        raise RuntimeError("expression_column_number is not expression_column_index + 1.")

    previous = primary["previous_primary_status"].fillna(False).astype(bool)
    if int(previous.sum()) != expected_previous:
        raise RuntimeError(
            f"Expected {expected_previous:,} previous primary rows, observed {int(previous.sum()):,}."
        )
    if int((~previous).sum()) != expected_new:
        raise RuntimeError(
            f"Expected {expected_new:,} new primary rows, observed {int((~previous).sum()):,}."
        )

    primary = primary.sort_values("expression_column_index", kind="mergesort").reset_index(drop=True)
    return primary, status


def load_model_weights(path: Path, expected_genes: int) -> tuple[list[str], np.ndarray, pd.DataFrame]:
    if not path.exists() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing model-weight table: {path}")
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    gene_column = "model_gene" if "model_gene" in frame.columns else "gene" if "gene" in frame.columns else None
    if gene_column is None:
        raise RuntimeError("Model-weight table lacks model_gene/gene column.")
    missing_components = sorted(set(COMPONENTS) - set(frame.columns))
    if missing_components:
        raise RuntimeError(f"Model-weight table lacks components: {missing_components}")

    work = frame[[gene_column, *COMPONENTS]].copy()
    work["model_gene"] = work[gene_column].map(normalize_symbol)
    if work["model_gene"].eq("").any():
        raise RuntimeError("Model-weight table contains blank gene names.")
    for component in COMPONENTS:
        work[component] = pd.to_numeric(work[component], errors="coerce").fillna(0.0)
    collapsed = work.groupby("model_gene", sort=True, as_index=False)[COMPONENTS].sum()
    if len(collapsed) != expected_genes:
        raise RuntimeError(
            f"Expected {expected_genes} unique model genes, observed {len(collapsed)}."
        )
    genes = collapsed["model_gene"].tolist()
    weights = collapsed[COMPONENTS].to_numpy(dtype=np.float64, copy=True)
    return genes, weights, collapsed


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
    missing = [path for path in required if path not in h5]
    if missing:
        raise RuntimeError(f"HDF5 is missing required datasets: {missing}")

    expression = h5["data/expression"]
    if expression.ndim != 2:
        raise RuntimeError(f"data/expression must be 2D, observed {expression.shape}")
    n_genes, n_samples = map(int, expression.shape)
    if n_genes != expected_genes or n_samples != expected_samples:
        raise RuntimeError(
            f"Unexpected HDF5 dimensions {n_genes:,} x {n_samples:,}; "
            f"expected {expected_genes:,} x {expected_samples:,}."
        )

    accession = h5["meta/samples/geo_accession"]
    symbols = decode_string_array(h5["meta/genes/symbol"][:])
    ensembl = decode_string_array(h5["meta/genes/ensembl_gene"][:])
    if len(accession) != n_samples or len(symbols) != n_genes or len(ensembl) != n_genes:
        raise RuntimeError("HDF5 metadata axes are not aligned to expression dimensions.")
    return expression, accession, symbols, ensembl


def build_model_mapping(
    symbols: np.ndarray,
    model_genes: list[str],
    weights: np.ndarray,
    expected_alias_resolved: int,
) -> ModelMapping:
    symbol_upper = np.asarray([normalize_symbol(value) for value in symbols], dtype=object)
    positions: dict[str, list[int]] = collections.defaultdict(list)
    for row, symbol in enumerate(symbol_upper):
        if symbol:
            positions[symbol].append(row)

    resolution_rows: list[dict[str, Any]] = []
    model_row_positions: list[int] = []
    model_row_indices: list[int] = []

    for model_index, gene in enumerate(model_genes):
        exact_rows = positions.get(gene, [])
        alias = MODEL_GENE_ALIASES.get(gene, "")
        alias_rows = positions.get(alias, []) if alias else []
        if exact_rows:
            resolved_symbol = gene
            rows = exact_rows
            resolution = "exact"
        elif alias_rows:
            resolved_symbol = alias
            rows = alias_rows
            resolution = "alias_resolved"
        else:
            resolved_symbol = ""
            rows = []
            resolution = "unresolved"

        for row in rows:
            model_row_positions.append(row)
            model_row_indices.append(model_index)

        resolution_rows.append(
            {
                "model_gene": gene,
                "resolved_archs4_symbol": resolved_symbol,
                "resolution": resolution,
                "archs4_row_count": len(rows),
                "archs4_row_indices_zero_based": ";".join(str(x) for x in rows),
                "alias_target": alias,
                "model_abs_weight_total": float(np.abs(weights[model_index, :]).sum()),
                **{
                    f"abs_weight_{component}": float(abs(weights[model_index, j]))
                    for j, component in enumerate(COMPONENTS)
                },
            }
        )

    resolution = pd.DataFrame(resolution_rows)
    unresolved = resolution.loc[resolution["resolution"].eq("unresolved"), "model_gene"].tolist()
    if unresolved:
        raise RuntimeError(f"Unresolved SpatialEcoTyper model genes: {unresolved}")
    alias_count = int(resolution["resolution"].eq("alias_resolved").sum())
    if alias_count != expected_alias_resolved:
        raise RuntimeError(
            f"Expected {expected_alias_resolved} alias-resolved model genes, observed {alias_count}."
        )

    total_abs_by_component = np.abs(weights).sum(axis=0)
    if np.any(total_abs_by_component <= 0):
        raise RuntimeError("At least one model component has zero total absolute weight.")
    total_abs_overall = float(np.abs(weights).sum())
    return ModelMapping(
        model_genes=model_genes,
        weights=weights,
        component_names=COMPONENTS.copy(),
        model_row_positions=np.asarray(model_row_positions, dtype=np.int64),
        model_row_indices=np.asarray(model_row_indices, dtype=np.int64),
        resolution=resolution,
        total_abs_by_component=total_abs_by_component,
        total_abs_overall=total_abs_overall,
    )


def make_unique_symbol_codes(symbols: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    normalized = np.asarray([normalize_symbol(value) for value in symbols], dtype=object)
    valid = normalized != ""
    codes = np.full(len(normalized), -1, dtype=np.int64)
    unique_symbols = pd.unique(normalized[valid])
    lookup = {symbol: index for index, symbol in enumerate(unique_symbols)}
    codes[valid] = np.asarray([lookup[x] for x in normalized[valid]], dtype=np.int64)
    return normalized, codes, len(unique_symbols)


def exact_profile_sha256(values: np.ndarray, original_dtype: np.dtype) -> str:
    hasher = hashlib.sha256()
    hasher.update(str(np.dtype(original_dtype)).encode("utf-8"))
    hasher.update(b"\0")
    hasher.update(np.asarray(values).shape.__repr__().encode("utf-8"))
    hasher.update(b"\0")
    hasher.update(np.ascontiguousarray(values).tobytes(order="C"))
    return hasher.hexdigest()


def safe_quantile(values: np.ndarray, q: float) -> float:
    valid = np.asarray(values, dtype=np.float64)
    valid = valid[np.isfinite(valid)]
    return float(np.quantile(valid, q)) if len(valid) else math.nan


def robust_mad(values: np.ndarray) -> float:
    valid = np.asarray(values, dtype=np.float64)
    valid = valid[np.isfinite(valid)]
    if not len(valid):
        return math.nan
    median = np.median(valid)
    return float(np.median(np.abs(valid - median)) * 1.4826)


def lower_warning_threshold(values: np.ndarray, objective_floor: float) -> tuple[float, dict[str, float]]:
    valid = np.asarray(values, dtype=np.float64)
    valid = valid[np.isfinite(valid)]
    median = float(np.median(valid))
    mad = robust_mad(valid)
    q005 = safe_quantile(valid, 0.005)
    q01 = safe_quantile(valid, 0.01)
    robust = median - 4.5 * mad if np.isfinite(mad) else q005
    candidate = min(q005, robust) if np.isfinite(robust) else q005
    threshold = max(float(objective_floor), float(candidate))
    return threshold, {
        "reference_median": median,
        "reference_mad_scaled": mad,
        "reference_q005": q005,
        "reference_q01": q01,
        "robust_candidate": float(robust),
    }


def upper_warning_threshold(values: np.ndarray, objective_cap: float) -> tuple[float, dict[str, float]]:
    valid = np.asarray(values, dtype=np.float64)
    valid = valid[np.isfinite(valid)]
    median = float(np.median(valid))
    mad = robust_mad(valid)
    q99 = safe_quantile(valid, 0.99)
    q995 = safe_quantile(valid, 0.995)
    robust = median + 4.5 * mad if np.isfinite(mad) else q995
    candidate = max(q995, robust) if np.isfinite(robust) else q995
    threshold = min(float(objective_cap), float(candidate))
    return threshold, {
        "reference_median": median,
        "reference_mad_scaled": mad,
        "reference_q99": q99,
        "reference_q995": q995,
        "robust_candidate": float(robust),
    }


def append_code(existing: str, code: str) -> str:
    values = [x for x in clean_text(existing).split(";") if x]
    if code and code not in values:
        values.append(code)
    return ";".join(values)


def summarize_by(frame: pd.DataFrame, group_column: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for group, sub in frame.groupby(group_column, dropna=False, sort=True):
        counts = sub["expression_qc_status"].value_counts().to_dict()
        rows.append(
            {
                group_column: clean_text(group),
                "samples": len(sub),
                "compatible": int(counts.get("COMPATIBLE", 0)),
                "compatible_with_warning": int(counts.get("COMPATIBLE_WITH_WARNING", 0)),
                "review_technical_outlier": int(counts.get("REVIEW_TECHNICAL_OUTLIER", 0)),
                "incompatible_hard_fail": int(counts.get("INCOMPATIBLE_HARD_FAIL", 0)),
                "median_library_sum": float(sub["library_sum_archs4_pseudocounts"].median()),
                "median_detected_unique_symbols": float(sub["detected_unique_gene_symbols"].median()),
                "median_model_genes_detected": float(sub["model_genes_detected"].median()),
                "median_min_component_weight_coverage": float(sub["min_component_weight_coverage"].median()),
            }
        )
    return pd.DataFrame(rows)


def run_self_tests() -> None:
    assert normalize_symbol(" wars ") == "WARS"
    assert normalize_symbol(None) == ""
    x = np.array([1, 2, 3, 4, 100], dtype=float)
    assert robust_mad(x) > 0
    lo, _ = lower_warning_threshold(np.arange(1, 101), 1)
    hi, _ = upper_warning_threshold(np.arange(1, 101), 1000)
    assert lo >= 1
    assert hi <= 1000
    assert append_code("A;B", "B") == "A;B"
    assert append_code("A", "C") == "A;C"

    # Exact/alias model mapping synthetic test.
    symbols = np.array(["WARS1", "TP53", "TP53", "GAPDH"], dtype=object)
    genes = ["WARS", "TP53"]
    weights = np.ones((2, len(COMPONENTS)), dtype=float)
    mapping = build_model_mapping(symbols, genes, weights, expected_alias_resolved=1)
    assert mapping.resolution["resolution"].tolist() == ["exact", "alias_resolved"] or mapping.resolution["resolution"].tolist() == ["alias_resolved", "exact"]
    # Order is genes as supplied: WARS alias, TP53 exact.
    row = mapping.resolution.set_index("model_gene")
    assert row.loc["WARS", "resolution"] == "alias_resolved"
    assert row.loc["TP53", "archs4_row_count"] == 2

    # Synthetic HDF5 expression access and fingerprint stability.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "test.h5"
        with h5py.File(path, "w") as h5:
            h5.create_dataset("data/expression", data=np.arange(20, dtype=np.int32).reshape(4, 5))
            h5.create_dataset("meta/samples/geo_accession", data=np.array([f"GSM{i}".encode() for i in range(5)]))
            h5.create_dataset("meta/genes/symbol", data=np.array([b"A", b"B", b"C", b"D"]))
            h5.create_dataset("meta/genes/ensembl_gene", data=np.array([b"E1", b"E2", b"E3", b"E4"]))
        with h5py.File(path, "r") as h5:
            expression, accession, syms, ens = validate_h5(h5, 4, 5)
            assert expression.shape == (4, 5)
            assert decode_string_array(accession[[1, 3]]).tolist() == ["GSM1", "GSM3"]
            a = expression[:, 1]
            assert exact_profile_sha256(a, expression.dtype) == exact_profile_sha256(a.copy(), expression.dtype)
            assert syms.tolist() == ["A", "B", "C", "D"]
            assert ens.tolist() == ["E1", "E2", "E3", "E4"]

    print("ALL ARCHS4 PRIMARY EXPRESSION COMPATIBILITY V1 SELF-TESTS PASSED")


def main(args: argparse.Namespace) -> None:
    if args.sample_batch_size < 1 or args.sample_batch_size > 256:
        raise ValueError("--sample-batch-size must be between 1 and 256.")

    started = utc_now()
    started_epoch = time.time()
    h5_path = args.h5.resolve()
    manifest_path = args.stage05_manifest.resolve()
    status_path = args.stage05_status.resolve()
    weights_path = args.model_weights.resolve()
    output_dir = args.output_dir.resolve()

    for path in [h5_path, manifest_path, status_path, weights_path]:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Required input is missing or empty: {path}")

    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output directory is non-empty: {output_dir}")
    staging = output_dir.parent / f"{output_dir.name}.staging.{os.getpid()}"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True, exist_ok=False)

    try:
        log("Loading frozen Stage 05 v1.0.1 manifest.")
        primary, stage05_status = load_stage05_inputs(
            manifest_path,
            status_path,
            args.expected_candidates,
            args.expected_primary,
            args.expected_previous_primary,
            args.expected_new_primary,
        )
        model_genes, model_weights, _ = load_model_weights(weights_path, args.expected_model_genes)

        n_samples = len(primary)
        sample_ids = primary["sample_id"].tolist()
        column_indices = primary["expression_column_index"].to_numpy(dtype=np.int64)
        previous_mask = primary["previous_primary_status"].fillna(False).astype(bool).to_numpy()

        log(f"Opening ARCHS4 HDF5: {h5_path}")
        with h5py.File(h5_path, "r") as h5:
            expression, accession, symbols, ensembl = validate_h5(
                h5, args.expected_h5_genes, args.expected_h5_samples
            )
            expression_dtype_text = str(expression.dtype)
            expression_chunks_text = str(expression.chunks)
            expression_compression_text = expression.compression or ""
            if column_indices.min() < 0 or column_indices.max() >= expression.shape[1]:
                raise RuntimeError("At least one Stage 05 expression column index is out of bounds.")

            mapping = build_model_mapping(
                symbols,
                model_genes,
                model_weights,
                args.expected_alias_resolved,
            )
            normalized_symbols, symbol_codes, n_unique_symbols = make_unique_symbol_codes(symbols)
            valid_symbol_rows = symbol_codes >= 0

            inventory_rows = [
                {"metric": "h5_path", "value": str(h5_path)},
                {"metric": "h5_basename", "value": h5_path.name},
                {"metric": "h5_size_bytes", "value": h5_path.stat().st_size},
                {"metric": "expression_shape", "value": f"{expression.shape[0]}x{expression.shape[1]}"},
                {"metric": "expression_dtype", "value": expression_dtype_text},
                {"metric": "expression_chunks", "value": expression_chunks_text},
                {"metric": "expression_compression", "value": expression_compression_text},
                {"metric": "gene_rows", "value": expression.shape[0]},
                {"metric": "sample_columns", "value": expression.shape[1]},
                {"metric": "nonblank_unique_gene_symbols", "value": n_unique_symbols},
                {"metric": "blank_gene_symbol_rows", "value": int((normalized_symbols == "").sum())},
                {"metric": "duplicated_symbol_rows", "value": int(valid_symbol_rows.sum() - n_unique_symbols)},
                {"metric": "model_weight_path", "value": str(weights_path)},
                {"metric": "model_weight_sha256", "value": sha256_file(weights_path)},
                {"metric": "model_genes", "value": len(model_genes)},
                {"metric": "model_genes_exact", "value": int(mapping.resolution["resolution"].eq("exact").sum())},
                {"metric": "model_genes_alias_resolved", "value": int(mapping.resolution["resolution"].eq("alias_resolved").sum())},
                {"metric": "model_genes_unresolved", "value": int(mapping.resolution["resolution"].eq("unresolved").sum())},
            ]
            write_tsv(pd.DataFrame(inventory_rows), staging / "01_expression_source_and_gene_axis_inventory.tsv")
            write_tsv(mapping.resolution, staging / "02_model_gene_resolution.tsv")

            metrics: dict[str, np.ndarray] = {
                "h5_accession_match": np.zeros(n_samples, dtype=bool),
                "library_sum_archs4_pseudocounts": np.zeros(n_samples, dtype=np.float64),
                "detected_gene_rows": np.zeros(n_samples, dtype=np.int32),
                "detected_unique_gene_symbols": np.zeros(n_samples, dtype=np.int32),
                "zero_gene_row_fraction": np.zeros(n_samples, dtype=np.float64),
                "positive_expression_q75": np.zeros(n_samples, dtype=np.float64),
                "top1_gene_fraction": np.zeros(n_samples, dtype=np.float64),
                "top10_gene_fraction": np.zeros(n_samples, dtype=np.float64),
                "model_genes_detected": np.zeros(n_samples, dtype=np.int16),
                "model_gene_detection_fraction": np.zeros(n_samples, dtype=np.float64),
                "overall_model_weight_coverage": np.zeros(n_samples, dtype=np.float64),
                "min_component_weight_coverage": np.zeros(n_samples, dtype=np.float64),
                "median_component_weight_coverage": np.zeros(n_samples, dtype=np.float64),
                "components_weight_coverage_ge_0_80": np.zeros(n_samples, dtype=np.int16),
                "negative_expression_values": np.zeros(n_samples, dtype=np.int32),
                "nonfinite_expression_values": np.zeros(n_samples, dtype=np.int32),
            }
            top_gene_symbol = np.empty(n_samples, dtype=object)
            profile_hashes = np.empty(n_samples, dtype=object)
            accession_observed = np.empty(n_samples, dtype=object)
            component_coverage = np.zeros((n_samples, len(COMPONENTS)), dtype=np.float32)

            log(
                f"Scanning {n_samples:,} primary samples in batches of {args.sample_batch_size}; "
                f"HDF5 chunks={expression.chunks}."
            )
            for start in range(0, n_samples, args.sample_batch_size):
                end = min(start + args.sample_batch_size, n_samples)
                cols = column_indices[start:end]
                # Primary rows are sorted by expression_column_index, satisfying h5py's
                # strictly increasing fancy-index requirement.
                raw = expression[:, cols]
                observed = decode_string_array(accession[cols])
                expected = np.asarray(sample_ids[start:end], dtype=object)
                accession_observed[start:end] = observed
                metrics["h5_accession_match"][start:end] = observed == expected

                values = np.asarray(raw, dtype=np.float64)
                nonfinite = ~np.isfinite(values)
                negative = values < 0
                metrics["nonfinite_expression_values"][start:end] = nonfinite.sum(axis=0)
                metrics["negative_expression_values"][start:end] = negative.sum(axis=0)
                clean = values.copy()
                clean[nonfinite | negative] = 0.0

                library = clean.sum(axis=0)
                positive = clean > 0
                detected_rows = positive.sum(axis=0)
                metrics["library_sum_archs4_pseudocounts"][start:end] = library
                metrics["detected_gene_rows"][start:end] = detected_rows
                metrics["zero_gene_row_fraction"][start:end] = 1.0 - detected_rows / clean.shape[0]

                # Top-gene dominance and positive-value quartile.
                if clean.shape[0] >= 10:
                    top10 = np.partition(clean, clean.shape[0] - 10, axis=0)[-10:, :].sum(axis=0)
                else:
                    top10 = clean.sum(axis=0)
                top_indices = np.argmax(clean, axis=0)
                top_values = clean[top_indices, np.arange(clean.shape[1])]
                safe_library = np.where(library > 0, library, np.nan)
                metrics["top1_gene_fraction"][start:end] = np.nan_to_num(top_values / safe_library, nan=0.0)
                metrics["top10_gene_fraction"][start:end] = np.nan_to_num(top10 / safe_library, nan=0.0)
                top_gene_symbol[start:end] = [clean_text(symbols[i]) for i in top_indices]

                for local in range(clean.shape[1]):
                    global_index = start + local
                    positive_rows = np.flatnonzero(positive[:, local] & valid_symbol_rows)
                    metrics["detected_unique_gene_symbols"][global_index] = (
                        np.unique(symbol_codes[positive_rows]).size if len(positive_rows) else 0
                    )
                    positive_values = clean[positive[:, local], local]
                    metrics["positive_expression_q75"][global_index] = (
                        float(np.quantile(positive_values, 0.75)) if len(positive_values) else 0.0
                    )
                    profile_hashes[global_index] = exact_profile_sha256(
                        raw[:, local], expression.dtype
                    )

                model_counts = np.zeros((len(model_genes), clean.shape[1]), dtype=np.float64)
                selected = clean[mapping.model_row_positions, :]
                np.add.at(model_counts, mapping.model_row_indices, selected)
                model_detected_mask = model_counts > 0
                model_detected = model_detected_mask.sum(axis=0)
                metrics["model_genes_detected"][start:end] = model_detected
                metrics["model_gene_detection_fraction"][start:end] = model_detected / len(model_genes)

                abs_weights = np.abs(mapping.weights)
                coverage = model_detected_mask.T.astype(np.float64) @ abs_weights
                coverage /= mapping.total_abs_by_component[None, :]
                component_coverage[start:end, :] = coverage.astype(np.float32)
                metrics["min_component_weight_coverage"][start:end] = coverage.min(axis=1)
                metrics["median_component_weight_coverage"][start:end] = np.median(coverage, axis=1)
                metrics["components_weight_coverage_ge_0_80"][start:end] = (coverage >= 0.80).sum(axis=1)
                overall = (
                    model_detected_mask.T.astype(np.float64) @ abs_weights.sum(axis=1)
                ) / mapping.total_abs_overall
                metrics["overall_model_weight_coverage"][start:end] = overall

                log(f"Expression QC: {end:,}/{n_samples:,} samples")

        # Construct the per-sample QC table after closing HDF5.
        qc = primary[
            [
                "sample_id",
                "expression_column_index",
                "expression_column_number",
                "primary_gse_id",
                "previous_lock_tier",
                "previous_primary_status",
                "stage05_resolved_cancer_class",
                "stage05_independence_group",
                "stage05_independence_group_source",
                "stage05_independence_review_required",
                "stage05_canonical_priority_hint",
                "stage05_canonical_priority_reason",
            ]
        ].copy()
        qc["h5_geo_accession_observed"] = accession_observed
        for name, values in metrics.items():
            qc[name] = values
        qc["top_expressed_gene_symbol"] = top_gene_symbol
        qc["expression_profile_sha256"] = profile_hashes
        for component_index, component in enumerate(COMPONENTS):
            qc[f"weight_coverage_{component}"] = component_coverage[:, component_index]

        # Exact expression-profile duplicate groups.
        hash_counts = qc["expression_profile_sha256"].value_counts()
        duplicate_hashes = set(hash_counts.loc[hash_counts > 1].index)
        qc["exact_profile_duplicate_group_size"] = qc["expression_profile_sha256"].map(hash_counts).astype(np.int32)
        qc["exact_profile_duplicate_flag"] = qc["expression_profile_sha256"].isin(duplicate_hashes)
        qc["exact_profile_duplicate_scope"] = ""
        duplicate_rows: list[dict[str, Any]] = []
        duplicate_group_number = 0
        for profile_hash, sub in qc.loc[qc["exact_profile_duplicate_flag"]].groupby(
            "expression_profile_sha256", sort=True
        ):
            duplicate_group_number += 1
            gses = sorted(sub["primary_gse_id"].astype(str).unique())
            groups = sorted(sub["stage05_independence_group"].astype(str).unique())
            cancers = sorted(sub["stage05_resolved_cancer_class"].astype(str).unique())
            if len(gses) > 1:
                scope = "CROSS_STUDY"
            elif len(groups) > 1:
                scope = "WITHIN_STUDY_ACROSS_INDEPENDENCE_GROUPS"
            else:
                scope = "WITHIN_INDEPENDENCE_GROUP"
            qc.loc[sub.index, "exact_profile_duplicate_scope"] = scope
            group_id = f"EXACT_EXPR_DUP_{duplicate_group_number:05d}"
            for row in sub.itertuples(index=False):
                duplicate_rows.append(
                    {
                        "duplicate_group_id": group_id,
                        "expression_profile_sha256": profile_hash,
                        "duplicate_group_size": len(sub),
                        "duplicate_scope": scope,
                        "group_GSE_count": len(gses),
                        "group_independence_count": len(groups),
                        "group_cancer_class_count": len(cancers),
                        "sample_id": row.sample_id,
                        "primary_gse_id": row.primary_gse_id,
                        "stage05_independence_group": row.stage05_independence_group,
                        "stage05_resolved_cancer_class": row.stage05_resolved_cancer_class,
                    }
                )
        duplicate_table = pd.DataFrame(
            duplicate_rows,
            columns=[
                "duplicate_group_id",
                "expression_profile_sha256",
                "duplicate_group_size",
                "duplicate_scope",
                "group_GSE_count",
                "group_independence_count",
                "group_cancer_class_count",
                "sample_id",
                "primary_gse_id",
                "stage05_independence_group",
                "stage05_resolved_cancer_class",
            ],
        )

        # Thresholds are calibrated from the frozen historical primary cohort.
        reference = qc.loc[qc["previous_primary_status"].fillna(False).astype(bool)].copy()
        if len(reference) != args.expected_previous_primary:
            raise RuntimeError("Previous-primary reference cohort size changed unexpectedly.")

        threshold_rows: list[dict[str, Any]] = []
        lower_specs = {
            "library_sum_archs4_pseudocounts": 1.0,
            "detected_unique_gene_symbols": 500.0,
            "model_genes_detected": 25.0,
            "overall_model_weight_coverage": 0.10,
            "min_component_weight_coverage": 0.05,
        }
        upper_specs = {
            "top1_gene_fraction": 0.95,
            "top10_gene_fraction": 0.995,
            "zero_gene_row_fraction": 0.995,
        }
        thresholds: dict[str, float] = {}
        for metric, floor in lower_specs.items():
            threshold, detail = lower_warning_threshold(reference[metric].to_numpy(), floor)
            thresholds[f"lower::{metric}"] = threshold
            threshold_rows.append(
                {
                    "metric": metric,
                    "direction": "lower",
                    "warning_threshold": threshold,
                    "objective_bound": floor,
                    "reference_samples": len(reference),
                    **detail,
                }
            )
        for metric, cap in upper_specs.items():
            threshold, detail = upper_warning_threshold(reference[metric].to_numpy(), cap)
            thresholds[f"upper::{metric}"] = threshold
            threshold_rows.append(
                {
                    "metric": metric,
                    "direction": "upper",
                    "warning_threshold": threshold,
                    "objective_bound": cap,
                    "reference_samples": len(reference),
                    **detail,
                }
            )
        threshold_table = pd.DataFrame(threshold_rows)

        # Objective failures and conservative review flags.
        hard_codes = np.full(len(qc), "", dtype=object)
        warning_codes = np.full(len(qc), "", dtype=object)
        severe_review = np.zeros(len(qc), dtype=bool)
        warning_count = np.zeros(len(qc), dtype=np.int16)

        def add_hard(mask: np.ndarray | pd.Series, code: str) -> None:
            array = np.asarray(mask, dtype=bool)
            for i in np.flatnonzero(array):
                hard_codes[i] = append_code(hard_codes[i], code)

        def add_warning(mask: np.ndarray | pd.Series, code: str, severe: bool = False) -> None:
            array = np.asarray(mask, dtype=bool)
            warning_count[array] += 1
            if severe:
                severe_review[array] = True
            for i in np.flatnonzero(array):
                warning_codes[i] = append_code(warning_codes[i], code)

        add_hard(~qc["h5_accession_match"].to_numpy(), "H5_ACCESSION_MISMATCH")
        add_hard(qc["nonfinite_expression_values"].to_numpy() > 0, "NONFINITE_EXPRESSION")
        add_hard(qc["negative_expression_values"].to_numpy() > 0, "NEGATIVE_EXPRESSION")
        add_hard(qc["library_sum_archs4_pseudocounts"].to_numpy() <= 0, "NONPOSITIVE_LIBRARY")
        add_hard(qc["detected_unique_gene_symbols"].to_numpy() < 200, "GROSSLY_LOW_GENE_DETECTION")
        add_hard(qc["model_genes_detected"].to_numpy() < 10, "GROSSLY_LOW_MODEL_FEATURE_DETECTION")
        add_hard(qc["overall_model_weight_coverage"].to_numpy() < 0.02, "GROSSLY_LOW_MODEL_WEIGHT_COVERAGE")

        for metric in lower_specs:
            add_warning(
                qc[metric].to_numpy() < thresholds[f"lower::{metric}"],
                f"LOW_{metric.upper()}",
            )
        for metric in upper_specs:
            add_warning(
                qc[metric].to_numpy() > thresholds[f"upper::{metric}"],
                f"HIGH_{metric.upper()}",
            )

        add_warning(qc["detected_unique_gene_symbols"].to_numpy() < 1000, "SEVERELY_LOW_GENE_DETECTION", severe=True)
        add_warning(qc["model_genes_detected"].to_numpy() < 50, "SEVERELY_LOW_MODEL_FEATURE_DETECTION", severe=True)
        add_warning(qc["overall_model_weight_coverage"].to_numpy() < 0.25, "SEVERELY_LOW_MODEL_WEIGHT_COVERAGE", severe=True)
        add_warning(qc["min_component_weight_coverage"].to_numpy() < 0.10, "SEVERELY_LOW_COMPONENT_COVERAGE", severe=True)
        add_warning(qc["top1_gene_fraction"].to_numpy() > 0.80, "EXTREME_SINGLE_GENE_DOMINANCE", severe=True)

        cross_group_duplicate = qc["exact_profile_duplicate_scope"].isin(
            ["CROSS_STUDY", "WITHIN_STUDY_ACROSS_INDEPENDENCE_GROUPS"]
        ).to_numpy()
        within_group_duplicate = qc["exact_profile_duplicate_scope"].eq(
            "WITHIN_INDEPENDENCE_GROUP"
        ).to_numpy()
        add_warning(cross_group_duplicate, "EXACT_PROFILE_DUPLICATE_ACROSS_INDEPENDENCE_GROUPS", severe=True)
        add_warning(within_group_duplicate, "EXACT_PROFILE_DUPLICATE_WITHIN_INDEPENDENCE_GROUP")

        qc["expression_qc_hard_failure_codes"] = hard_codes
        qc["expression_qc_warning_codes"] = warning_codes
        qc["expression_qc_warning_count"] = warning_count
        hard_fail = pd.Series(hard_codes).astype(str).ne("").to_numpy()
        review = (~hard_fail) & (severe_review | (warning_count >= 2))
        minor = (~hard_fail) & (~review) & (warning_count == 1)
        status = np.full(len(qc), "COMPATIBLE", dtype=object)
        status[minor] = "COMPATIBLE_WITH_WARNING"
        status[review] = "REVIEW_TECHNICAL_OUTLIER"
        status[hard_fail] = "INCOMPATIBLE_HARD_FAIL"
        qc["expression_qc_status"] = status
        qc["expression_compatible_for_next_stage"] = np.isin(
            status, ["COMPATIBLE", "COMPATIBLE_WITH_WARNING"]
        )
        qc["expression_qc_review_required"] = status == "REVIEW_TECHNICAL_OUTLIER"
        qc["expression_qc_hard_incompatible"] = status == "INCOMPATIBLE_HARD_FAIL"

        # Write primary QC products.
        qc.to_parquet(staging / "03_primary_expression_qc_metrics.parquet", index=False)
        write_tsv(
            qc.loc[qc["expression_compatible_for_next_stage"]],
            staging / "04_expression_compatible_primary_candidates.tsv.gz",
            compression="gzip",
        )
        write_tsv(
            qc.loc[qc["expression_qc_review_required"]],
            staging / "05_expression_qc_review_candidates.tsv.gz",
            compression="gzip",
        )
        write_tsv(
            qc.loc[qc["expression_qc_hard_incompatible"]],
            staging / "06_expression_incompatible_candidates.tsv.gz",
            compression="gzip",
        )
        write_tsv(threshold_table, staging / "07_expression_qc_thresholds.tsv")
        write_tsv(
            duplicate_table,
            staging / "08_exact_expression_duplicate_groups.tsv.gz",
            compression="gzip",
        )
        write_tsv(
            summarize_by(qc, "stage05_resolved_cancer_class"),
            staging / "09_expression_qc_summary_by_cancer.tsv",
        )
        write_tsv(
            summarize_by(qc, "primary_gse_id"),
            staging / "10_expression_qc_summary_by_study.tsv",
        )

        previous_summary = (
            qc.assign(
                cohort=np.where(
                    qc["previous_primary_status"].fillna(False).astype(bool),
                    "previous_primary",
                    "new_primary",
                )
            )
            .groupby(["cohort", "expression_qc_status"], as_index=False)
            .size()
            .rename(columns={"size": "samples"})
        )
        write_tsv(previous_summary, staging / "11_previous_primary_expression_qc_summary.tsv")

        fingerprint_manifest = qc[
            [
                "sample_id",
                "primary_gse_id",
                "stage05_independence_group",
                "expression_column_index",
                "expression_profile_sha256",
                "exact_profile_duplicate_group_size",
                "exact_profile_duplicate_scope",
            ]
        ].copy()
        write_tsv(
            fingerprint_manifest,
            staging / "14_expression_profile_fingerprint_manifest.tsv.gz",
            compression="gzip",
        )
        alignment = qc[
            [
                "sample_id",
                "expression_column_index",
                "expression_column_number",
                "h5_geo_accession_observed",
                "h5_accession_match",
            ]
        ].copy()
        write_tsv(
            alignment,
            staging / "15_sample_index_alignment_audit.tsv.gz",
            compression="gzip",
        )

        coverage_summary_rows: list[dict[str, Any]] = []
        for component in COMPONENTS:
            column = f"weight_coverage_{component}"
            coverage_summary_rows.append(
                {
                    "component": component,
                    "samples": len(qc),
                    "minimum": float(qc[column].min()),
                    "q01": float(qc[column].quantile(0.01)),
                    "median": float(qc[column].median()),
                    "q99": float(qc[column].quantile(0.99)),
                    "maximum": float(qc[column].max()),
                    "samples_below_0_50": int((qc[column] < 0.50).sum()),
                    "samples_below_0_80": int((qc[column] < 0.80).sum()),
                }
            )
        write_tsv(
            pd.DataFrame(coverage_summary_rows),
            staging / "16_model_component_weight_coverage_summary.tsv",
        )

        enriched = primary.merge(
            qc,
            on=[
                "sample_id",
                "expression_column_index",
                "expression_column_number",
                "primary_gse_id",
                "previous_lock_tier",
                "previous_primary_status",
                "stage05_resolved_cancer_class",
                "stage05_independence_group",
                "stage05_independence_group_source",
                "stage05_independence_review_required",
                "stage05_canonical_priority_hint",
                "stage05_canonical_priority_reason",
            ],
            how="left",
            validate="one_to_one",
        )
        enriched.to_parquet(staging / "17_primary_manifest_with_expression_qc.parquet", index=False)

        counts = qc["expression_qc_status"].value_counts().to_dict()
        previous_hard = int(
            (
                qc["previous_primary_status"].fillna(False).astype(bool)
                & qc["expression_qc_hard_incompatible"]
            ).sum()
        )
        previous_review = int(
            (
                qc["previous_primary_status"].fillna(False).astype(bool)
                & qc["expression_qc_review_required"]
            ).sum()
        )
        new_hard = int(
            (
                ~qc["previous_primary_status"].fillna(False).astype(bool)
                & qc["expression_qc_hard_incompatible"]
            ).sum()
        )
        new_review = int(
            (
                ~qc["previous_primary_status"].fillna(False).astype(bool)
                & qc["expression_qc_review_required"]
            ).sum()
        )

        qc_checks: list[dict[str, Any]] = []
        def check(name: str, condition: bool, detail: str, warning_only: bool = False) -> None:
            qc_checks.append(
                {
                    "check": name,
                    "status": "PASS" if condition else ("WARN" if warning_only else "FAIL"),
                    "detail": detail,
                }
            )

        check("stage05_status_pass", stage05_status.get("status") == "PASS", str(stage05_status.get("status")))
        check("input_candidate_rows_exact", int(stage05_status.get("input_candidate_samples", -1)) == args.expected_candidates, f"rows={stage05_status.get('input_candidate_samples')}")
        check("primary_rows_exact", len(qc) == args.expected_primary, f"rows={len(qc)}")
        check("previous_primary_rows_exact", int(qc["previous_primary_status"].fillna(False).astype(bool).sum()) == args.expected_previous_primary, f"rows={int(qc['previous_primary_status'].fillna(False).astype(bool).sum())}")
        check("new_primary_rows_exact", int((~qc["previous_primary_status"].fillna(False).astype(bool)).sum()) == args.expected_new_primary, f"rows={int((~qc['previous_primary_status'].fillna(False).astype(bool)).sum())}")
        check("unique_sample_ids", not qc["sample_id"].duplicated().any(), f"duplicates={int(qc['sample_id'].duplicated().sum())}")
        check("unique_expression_indices", not qc["expression_column_index"].duplicated().any(), f"duplicates={int(qc['expression_column_index'].duplicated().sum())}")
        check("all_h5_accessions_match", bool(qc["h5_accession_match"].all()), f"mismatches={int((~qc['h5_accession_match']).sum())}")
        check("all_model_genes_resolved", int(mapping.resolution["resolution"].eq("unresolved").sum()) == 0, f"unresolved={int(mapping.resolution['resolution'].eq('unresolved').sum())}")
        check("model_gene_count_exact", len(mapping.model_genes) == args.expected_model_genes, f"genes={len(mapping.model_genes)}")
        check("alias_resolved_count_exact", int(mapping.resolution["resolution"].eq("alias_resolved").sum()) == args.expected_alias_resolved, f"aliases={int(mapping.resolution['resolution'].eq('alias_resolved').sum())}")
        check("expression_status_partition_complete", sum(counts.values()) == len(qc), f"status_total={sum(counts.values())};primary={len(qc)}")
        check("negative_values_zero", int(qc["negative_expression_values"].sum()) == 0, f"negative_values={int(qc['negative_expression_values'].sum())}")
        check("nonfinite_values_zero", int(qc["nonfinite_expression_values"].sum()) == 0, f"nonfinite_values={int(qc['nonfinite_expression_values'].sum())}")
        check("previous_primary_hard_failures_zero", previous_hard == 0, f"hard_failures={previous_hard}", warning_only=True)
        check("previous_primary_review_rows_reported", previous_review == previous_review, f"review_rows={previous_review}")
        check("final_lock_not_created", True, "final_lock_created=false")
        check("final_expression_matrix_not_exported", True, "expression_matrix_exported=false")

        contract = pd.DataFrame(qc_checks)
        write_tsv(contract, staging / "12_expression_QC_contract.tsv")
        failures = int(contract["status"].eq("FAIL").sum())
        warnings = int(contract["status"].eq("WARN").sum())

        status_payload = {
            "script": Path(__file__).name,
            "script_version": SCRIPT_VERSION,
            "status": "PASS" if failures == 0 else "FAIL",
            "started_at_utc": started,
            "completed_at_utc": utc_now(),
            "runtime_seconds": round(time.time() - started_epoch, 3),
            "host": socket.gethostname(),
            "h5_path": str(h5_path),
            "h5_expression_genes": args.expected_h5_genes,
            "h5_expression_samples": args.expected_h5_samples,
            "h5_expression_dtype": expression_dtype_text,
            "h5_expression_chunks": expression_chunks_text,
            "stage05_manifest": str(manifest_path),
            "stage05_status": str(status_path),
            "model_weights": str(weights_path),
            "model_weights_sha256": sha256_file(weights_path),
            "output_dir": str(output_dir),
            "primary_samples_scanned": len(qc),
            "previous_primary_samples": int(previous_mask.sum()),
            "new_primary_samples": int((~previous_mask).sum()),
            "compatible_samples": int(counts.get("COMPATIBLE", 0)),
            "compatible_with_warning_samples": int(counts.get("COMPATIBLE_WITH_WARNING", 0)),
            "technical_review_samples": int(counts.get("REVIEW_TECHNICAL_OUTLIER", 0)),
            "hard_incompatible_samples": int(counts.get("INCOMPATIBLE_HARD_FAIL", 0)),
            "previous_primary_technical_review_samples": previous_review,
            "previous_primary_hard_incompatible_samples": previous_hard,
            "new_primary_technical_review_samples": new_review,
            "new_primary_hard_incompatible_samples": new_hard,
            "exact_duplicate_profiles": int(len(duplicate_hashes)),
            "samples_in_exact_duplicate_profiles": int(qc["exact_profile_duplicate_flag"].sum()),
            "cross_study_exact_duplicate_samples": int(qc["exact_profile_duplicate_scope"].eq("CROSS_STUDY").sum()),
            "model_genes": len(mapping.model_genes),
            "model_genes_exact": int(mapping.resolution["resolution"].eq("exact").sum()),
            "model_genes_alias_resolved": int(mapping.resolution["resolution"].eq("alias_resolved").sum()),
            "model_genes_unresolved": int(mapping.resolution["resolution"].eq("unresolved").sum()),
            "qc_failures": failures,
            "qc_warnings": warnings,
            "expression_qc_completed": failures == 0,
            "expression_matrix_exported": False,
            "final_lock_created": False,
        }
        atomic_json(staging / "13_expression_QC_status.json", status_payload)

        if failures:
            raise RuntimeError("Expression compatibility QC contract contains failures.")

        output_dir.parent.mkdir(parents=True, exist_ok=True)
        staging.replace(output_dir)

        log("=" * 96)
        log("ARCHS4 PRIMARY EXPRESSION COMPATIBILITY V1: PASS")
        log("=" * 96)
        log(f"Primary samples scanned:             {len(qc):,}")
        log(f"Compatible:                          {int(counts.get('COMPATIBLE', 0)):,}")
        log(f"Compatible with warning:             {int(counts.get('COMPATIBLE_WITH_WARNING', 0)):,}")
        log(f"Technical review:                    {int(counts.get('REVIEW_TECHNICAL_OUTLIER', 0)):,}")
        log(f"Hard incompatible:                   {int(counts.get('INCOMPATIBLE_HARD_FAIL', 0)):,}")
        log(f"Previous-primary technical review:   {previous_review:,}")
        log(f"Previous-primary hard incompatible:  {previous_hard:,}")
        log(f"Exact duplicate profile groups:      {len(duplicate_hashes):,}")
        log(f"Output:                               {output_dir}")

    except Exception:
        # Preserve staging outputs for diagnosis only when files have been written.
        raise


if __name__ == "__main__":
    args = parse_args()
    if args.self_test:
        run_self_tests()
        sys.exit(0)
    main(args)
