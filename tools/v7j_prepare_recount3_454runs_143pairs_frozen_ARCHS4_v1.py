#!/usr/bin/env python3
"""V7J: original-ARCHS4-frozen-baseline, recount3-only 499-gene preparation.

Reads original gzip gene-body coverage from all 454 SRR runs; aggregates within
143 exact SRX/GSM pairs; calls the authentic Stage09 symbol mapping and
normalization methods for recount3 ONLY. The ARCHS4 baseline must be read from
original immutable 9465-sample Stage09/10 outputs by the R projection runner.
TMM and UQ on recount3 are explicit 143-cohort sensitivity analyses.

No cohort admissions; sample-linked outputs only under ignored private_data.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import re
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

RUN_RE = re.compile(r"^(?:SRR|ERR|DRR)\d+$", re.I)
RUN_IN_TEXT = re.compile(r"(?:SRR|ERR|DRR)\d+", re.I)
EXPECTED_BRIDGE_SHA = "c7d6d7c754e3ca77d62424f26627a40ceaec10fac98a0741c945773f3d3b1b67"
EXPECTED_MAP_SHA = "68be266dd6bd30bfa79a2bd8ac79f0e49c6f5d0f19d2d083e2ca57850ef17391"
EXPECTED_RAW_SHA = "b056331eb16d6ec612409fd37a266abfe10c4c15329b486dbb9df213612797c5"
EXPECTED_STAGE09_SHA = "46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500"

def digest(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()

def read_tsv(p):
    with p.open(encoding="utf-8-sig", newline="") as f:
        d = csv.DictReader(f, delimiter="\t")
        return d.fieldnames or [], list(d)

def decode(v):
    return bytes(v).decode("utf-8").strip() if isinstance(v, (bytes, np.bytes_)) else str(v).strip()

def original_module(path):
    spec = importlib.util.spec_from_file_location("frozen_stage09_v101", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot import authoritative frozen Stage09 source")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def write_matrix(path, genes, ids, value):
    assert value.shape == (len(genes), len(ids))
    with gzip.open(path, "wt", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["model_gene", *ids])
        for g, r in zip(genes, value):
            w.writerow([g, *(format(float(v), ".17g") for v in r)])

def record_matrix(raw, gene_symbols, genes, old, assay, gsm, out):
    """Original Stage09 gene-collapse and factors, on recount3 (NOT ARCHS4)."""
    if raw.ndim != 2 or raw.shape[1] != 143:
        raise ValueError(assay + " expected matrix of 143 paired specimens")
    if not np.isfinite(raw).all() or (raw < 0).any():
        raise ValueError(assay + " negative or non-finite raw signal")
    maps = old.build_gene_maps(pd.Index(gene_symbols), genes)
    if len(maps["available_model_genes"]) != 499:
        raise ValueError(assay + " not all frozen W model genes resolve")
    unresolved = maps["model_gene_resolution"]["resolution"].eq("unresolved").sum()
    if unresolved:
        raise ValueError(assay + " unresolved model genes")

    collapsed = old.collapse_batch(
        raw, maps["valid_positions"], maps["row_to_unique"],
        len(maps["unique_genes"])
    )
    libsum = collapsed.sum(axis=0)
    if not np.isfinite(libsum).all() or (libsum <= 0).any():
        raise ValueError(assay + " invalid all-symbol normalization mass")

    model = np.zeros((499, 143), dtype=float)
    np.add.at(model, maps["model_row_indices"], raw[maps["model_row_positions"], :])
    uq75 = old.quantile75_fraction(collapsed, libsum)
    tmm_ref = old.select_tmm_reference(gsm, uq75, libsum)
    tmm = np.ones(143, dtype=float)
    retained = np.zeros(143, dtype=int)
    for i in range(143):
        tmm[i], retained[i] = old.tmm_factor(
            collapsed[:, i], collapsed[:, tmm_ref],
            float(libsum[i]), float(libsum[tmm_ref]), 0.30, 0.05
        )
    tmm /= old.geometric_mean_positive(tmm)
    if not np.isfinite(tmm).all() or (tmm <= 0).any():
        raise ValueError(assay + " invalid TMM")
    uq_valid = np.isfinite(uq75) & (uq75 > 0)
    if not uq_valid.any():
        raise ValueError(assay + " no evaluable UQ samples")
    uq = uq75.copy()
    uq[~uq_valid] = np.median(uq75[uq_valid])
    uq /= old.geometric_mean_positive(uq)
    if not np.isfinite(uq).all() or (uq <= 0).any():
        raise ValueError(assay + " invalid UQ factors")
    norms = {
        "CPM": libsum,
        "TMM_CPM": libsum*tmm,
        "UQ_CPM": libsum*uq,
    }
    for branch,denom in norms.items():
        values = model / denom[None, :] * 1e6
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(assay + " invalid normalized " + branch)
        write_matrix(out / ("499_" + assay + "_" + branch + "_LOCAL_ONLY.tsv.gz"),
                     genes, gsm, values)
    alias_n = int(maps["model_gene_resolution"]["resolution"].eq("alias_resolved").sum())
    print("RECOUNT3_ASSAY",assay,"RAW_ROWS",raw.shape[0],"VALID_SYMBOL_GROUPS",len(maps["unique_genes"]),
          "MODEL_GENES",len(maps["available_model_genes"]),"ALIAS_RESOLVED",alias_n,
          "UQ_EVALUABLE",int(uq_valid.sum()),"TMM_REF_CHOSEN",True,flush=True)
    return {"uq_valid":uq_valid, "tmm":tmm, "uq":uq, "library_mass":libsum,
            "tmm_ref_index_local_only":int(tmm_ref),
            "source_alias_resolved":alias_n}

def main():
    a = argparse.ArgumentParser()
    a.add_argument("root",type=Path)
    a.add_argument("output",type=Path)
    args = a.parse_args()
    root=args.root.resolve()
    out=args.output.resolve()
    repo=root/"RECOUNT3_CERTIFICATION_GIT"
    if out.parent != (repo/"private_data").resolve():
        raise ValueError("Only local ignored private_data output is allowed")
    if out.exists():
        raise FileExistsError("REFUSE_OVERWRITE "+str(out))

    rec=root/"26_recount3_external_primary_v1"
    stage5=rec/"stage5_output"
    rawpath=rec/"raw_studies/SRP074349/gene_sums/sra.gene_sums.SRP074349.G026.gz"
    bridgepath=stage5/"03_bridge_165_availability.tsv"
    mappath=stage5/"00_gene_row_mapping.tsv"
    rel=root/"RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
    frozen=rel/"09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py"
    expected={bridgepath:EXPECTED_BRIDGE_SHA,mappath:EXPECTED_MAP_SHA,
              rawpath:EXPECTED_RAW_SHA,frozen:EXPECTED_STAGE09_SHA}
    for p,s in expected.items():
        if digest(p) != s:
            raise ValueError("FROZEN_INPUT_SHA_MISMATCH "+str(p))
    old=original_module(frozen)
    bcols,bridge=read_tsv(bridgepath)
    assert len(bridge)==165
    assert {"study","library_id","ARCHS4_GSM","ARCHS4_expression_column_index","recount3_run_ids","n_runs"} <= set(bcols)
    pairs=sorted((r for r in bridge if r["study"]=="SRP074349"),key=lambda r:r["ARCHS4_GSM"])
    if len(pairs)!=143:
        raise ValueError("Expected 143 calibration-only bridge experiments")
    if set(r["bridge_class"] for r in pairs)!={"CALIBRATION_ONLY_NOT_IN_PRIMARY"}:
        raise ValueError("Bridge contains primary-cohort experiments")

    gsm=[r["ARCHS4_GSM"].strip() for r in pairs]
    libs=[r["library_id"].strip() for r in pairs]
    indices=[int(r["ARCHS4_expression_column_index"]) for r in pairs]
    if len(set(gsm))!=143 or len(set(libs))!=143 or len(set(indices))!=143 or min(indices)<0:
        raise ValueError("Duplicate GSM/SRX/index in bridge")
    run_to_experiment={}
    nrun=[]
    for exp,r in enumerate(pairs):
        runs=[x.upper() for x in RUN_IN_TEXT.findall(r["recount3_run_ids"] or "")]
        if len(runs)!=int(r["n_runs"]) or len(runs)!=len(set(runs)):
            raise ValueError("Incorrect declared run multiplicity")
        for run in runs:
            if run in run_to_experiment:
                raise ValueError("Cross-experiment duplicated run")
            run_to_experiment[run]=exp
        nrun.append(len(runs))
    if len(run_to_experiment)!=454 or Counter(nrun)!=Counter({1:84,2:4,3:22,8:1,9:32}):
        raise ValueError("Unexpected bridge run cardinality")
    print("BRIDGE_143_EXPERIMENTS_454_RUNS_PASS",flush=True)

    _,mapping=read_tsv(mappath)
    if len(mapping)!=63856 or [int(x["row_zero_based"]) for x in mapping]!=list(range(63856)):
        raise ValueError("G026 row-map axis mismatch")
    symbols=[r["gene_symbol"] for r in mapping]
    # In Stage09, valid genes are cleaned/uppercased and duplicated rows summed.
    audit_candidates=[
        root/"10_preprojection_normalization_v1_0_1/08_model_gene_availability.tsv",
        rel/"10_preprojection_normalization_v1_0_1/08_model_gene_availability.tsv",
    ]
    audit=next((p for p in audit_candidates if p.is_file()),None)
    if audit is None:
        raise FileNotFoundError("Frozen 499 gene availability ledger is missing")
    cols,availability=read_tsv(audit)
    if len(availability)!=499 or "model_gene" not in cols:
        raise ValueError("Frozen 499 gene ledger changed")
    genes=[r["model_gene"].strip() for r in availability]
    if len(set(genes))!=499:
        raise ValueError("Duplicated model genes")

    # Identity check ONLY; archived Stage09/10 matrices remain the reference.
    archpath=root/"00_download/human_gene_v2.latest.h5"
    with h5py.File(archpath,"r") as h:
        access=h["meta/samples/geo_accession"]
        if max(indices)>=len(access):
            raise ValueError("ARCHS4 H5 bridge index out of bounds")
        observed=[decode(access[i]) for i in indices]
        if observed!=gsm:
            raise ValueError("ARCHS4 exact GSM/index disagreement")
    print("ARCHS4_H5_GSM_INDEX_LOCK_PASS",len(gsm),flush=True)
    print("ARCHS4_BASELINE_ARCHIVED_FROZEN_NO_RENORMALIZATION",flush=True)

    # The raw G026 gzip header is line 3 (index 2) in V7F; discover by run IDs,
    # rather than assuming any fixed header offset or column order.
    g026=np.zeros((63856,143),dtype=np.float64)
    with gzip.open(rawpath,"rt",encoding="utf-8-sig",errors="strict",newline="") as f:
        header=None
        for i in range(30):
            line=f.readline()
            if not line:
                raise ValueError("Raw G026 header not found")
            vals=line.rstrip("\r\n").split("\t")
            if len(vals)==669 and all(RUN_RE.fullmatch(x.strip().strip('"')) for x in vals[1:]):
                header=vals
                header_line=i
                break
        if header is None:
            raise ValueError("Raw G026 668-run header not identified")
        headers=[x.strip().strip('"').upper() for x in header[1:]]
        if len(headers)!=668 or len(set(headers))!=668:
            raise ValueError("Raw G026 run axis is not 668 unique SRRs")
        if not set(run_to_experiment)<=set(headers):
            raise ValueError("Raw G026 missing matched runs")
        run_cols=np.asarray([j for j,x in enumerate(headers) if x in run_to_experiment],dtype=np.int64)
        groups=np.asarray([run_to_experiment[headers[j]] for j in run_cols],dtype=np.int64)
        if len(run_cols)!=454 or len(set(run_cols))!=454:
            raise ValueError("Run header projection axis mismatch")
        print("RAW_G026_HEADER_LINE_0BASED",header_line,
              "RAW_RUN_COLUMNS",len(headers),
              "MATCHED_RUN_COLUMNS",len(run_cols),flush=True)
        n=0
        max_relative_aggregation_mass_error=0.0
        for line in f:
            if not line.strip():
                raise ValueError("Unexpected blank row in raw G026")
            gene, sep, values = line.rstrip("\r\n").partition("\t")
            if not sep:
                raise ValueError("Expected tab-delimited G026 gene row")
            if n >= 63856:
                raise ValueError("Extra gene row in original G026 source")
            name=gene.strip().strip('"')
            source=mapping[n]
            if name!=source["source_gene_id"] and name.split(".",1)[0]!=source["gene_base"]:
                raise ValueError("Raw G026 gene order differs from Stage5 locked map; row="+str(n))
            numeric=np.fromstring(values,sep="\t",dtype=np.float64)
            if len(numeric)!=668:
                raise ValueError("Wrong number of G026 columns at row "+str(n))
            if not np.isfinite(numeric).all() or (numeric<0).any():
                raise ValueError("Invalid G026 gene-body coverage at row "+str(n))
            g026[n,:]=np.bincount(groups,weights=numeric[run_cols],minlength=143)
            source_total=float(np.sum(numeric[run_cols]))
            group_total=float(np.sum(g026[n,:]))
            max_relative_aggregation_mass_error=max(
                max_relative_aggregation_mass_error,
                abs(source_total-group_total)/max(1.0,source_total))
            n+=1
            if n%12000==0:
                print("G026_SOURCE_GENE_ROWS_PROCESSED",n,flush=True)
        if n!=63856:
            raise ValueError("Raw G026 gene rows incomplete: "+str(n))
    if max_relative_aggregation_mass_error > 1e-12:
        raise ValueError("Raw technical-run coverage mass was not conserved")
    print("RAW_G026_SOURCE_GENE_ROWS_PASS",n,flush=True)
    print("TECHNICAL_RUNS_SUMMED_WITHIN_EXPERIMENTS",454,143,flush=True)
    print("MAX_RELATIVE_RAW_COVERAGE_AGGREGATION_MASS_ERROR",
          format(max_relative_aggregation_mass_error,".12g"),flush=True)

    out.mkdir(parents=True,mode=0o700)
    recount_stats=record_matrix(g026,symbols,genes,old,"recount3",gsm,out)
    pair_file=out/"143_pair_manifest_LOCAL_ONLY.tsv"
    with pair_file.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["ARCHS4_GSM","ARCHS4_expression_column_index","library_id",
                    "n_technical_runs","UQ_evaluable_recount3","source_data_type"])
        for i in range(143):
            w.writerow([gsm[i],indices[i],libs[i],nrun[i],
                        str(bool(recount_stats["uq_valid"][i])).upper(),
                        "G026_gene_body_bp_coverage_sum_over_all_SRRs"])
    run_file=out/"454_run_to_experiment_LOCAL_ONLY.tsv"
    with run_file.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["run_id","library_id","ARCHS4_GSM"])
        for run,idx in sorted(run_to_experiment.items()):
            w.writerow([run,libs[idx],gsm[idx]])
    tracked=sorted(out.glob("*.tsv*"))
    manifest={
        "version":"v7j_frozen_archs4_vs_recount3_full499_v1",
        "study":"SRP074349",
        "n_paired_specimens":143,
        "n_matched_technical_runs":454,
        "n_total_raw_g026_run_columns":668,
        "n_raw_g026_gene_rows":63856,
        "n_model_genes":499,
        "reference_ARCHS4":"Original Stage09/10 immutable 9465-cohort per-GSM normalized matrices, scaled Z, raw H and fraction scores",
        "new_recount3_source":"G026 gene-body base-pair coverage summed across all technical SRRs within each SRX",
        "new_recount3_normalization":"Original Stage09 methods CPM/TMM/UQ; recount3 TMM and UQ factors estimated on 143 matched SRXs, not original 9465 ARCHS4 reference",
        "primary_branch":"CPM",
        "normalization_sensitivity_branches":["TMM_CPM","UQ_CPM"],
        "no_frozen_anchor_or_W_refit":True,
        "clinical_admissions":0,
        "input_sha256":{str(p):v for p,v in expected.items()},
        "frozen_model_gene_audit_sha256":digest(audit),
        "recount3_factors_LOCAL_ONLY":{
            "uq_evaluable":int(recount_stats["uq_valid"].sum()),
            "alias_resolved_in_recount3_source":recount_stats["source_alias_resolved"],
            "tmm_reference_index":recount_stats["tmm_ref_index_local_only"]},
        "outputs_sha256":{p.name:digest(p) for p in tracked},
    }
    (out/"input_provenance.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print("RECOUNT3_ONLY_FULL_499_INPUTS_PREPARED",flush=True)
    print("RECOUNT3_UQ_EVALUABLE",int(recount_stats["uq_valid"].sum()),flush=True)
    print("ARCHS4_ONLY_ORIGINAL_IMMUTABLE_STAGE09_STAGE10_REFERENCE",flush=True)
    print("NO_SOURCE_OR_FROZEN_FILE_MODIFICATION",flush=True)
    print("V7J_RECOUNT3_ONLY_143_MATCHED_FULL499_PREPARATION_PASS",flush=True)

if __name__=="__main__":
    main()
