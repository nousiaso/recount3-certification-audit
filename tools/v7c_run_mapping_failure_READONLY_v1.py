#!/usr/bin/env python3
"""Read-only diagnosis of failed V7C 143-pair run mapping and frozen input contracts.

Reports only aggregate counts and categories. Does NOT print study accessions,
patient identifiers, or source values; does NOT load expression matrix values.
Preserves failed V7C artifacts. No repairs or reruns.
"""
import csv
import importlib.util
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(sys.argv[1]).resolve()
REC = ROOT / "26_recount3_external_primary_v1"
STAGE5 = REC / "stage5_output"
RELEASE = ROOT / "RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
RX = re.compile(r"(?:SRR|ERR|DRR)\d+", re.I)

def tab(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        return r.fieldnames or [], list(r)

def decode(v):
    return bytes(v).decode("utf-8").strip() if isinstance(v, (bytes, np.bytes_)) else str(v).strip()

def report(label, counts):
    print(label, " ".join(f"{k}={v}" for k,v in sorted(counts.items(), key=lambda kv: str(kv[0]))))

print("V7C_RUN_MAPPING_DIAGNOSTIC_V1")
_, bridge = tab(STAGE5 / "03_bridge_165_availability.tsv")
pairs = [r for r in bridge if r["study"] == "SRP074349"]
assert len(bridge) == 165 and len(pairs) == 143
print("BRIDGE_ROWS", len(pairs))

bylib = {r["library_id"].strip(): r for r in pairs}
assert len(bylib) == 143
observed_n = Counter()
token_n = Counter()
n_vs_tokens = Counter()
presence = Counter()
reported_coverage = Counter()
bridge_class = Counter()
bridge_runs = {}
for lib,r in bylib.items():
    n = str(r.get("n_runs", "")).strip()
    tokens = RX.findall(r.get("recount3_run_ids", "") or "")
    unique = set(tokens)
    observed_n[n] += 1
    token_n[len(tokens)] += 1
    n_vs_tokens[(n,len(tokens),len(unique))] += 1
    presence[str(r.get("is_study_gene_file_present", "")).strip()] += 1
    reported_coverage[str(r.get("n_recount3_runs_present_in_downloaded_raw", "")).strip()] += 1
    bridge_class[str(r.get("bridge_class", "")).strip()] += 1
    bridge_runs[lib] = set(tokens)
report("BRIDGE_REPORTED_N_RUNS", observed_n)
report("BRIDGE_RUN_TOKEN_COUNT", token_n)
report("BRIDGE_NRUNS_TOKEN_UNIQUE_CROSS", n_vs_tokens)
report("BRIDGE_RAW_PRESENT_COUNTS", reported_coverage)
report("BRIDGE_CLASS_COUNTS", bridge_class)
report("BRIDGE_FILE_PRESENT", presence)
print("BRIDGE_NRUN_ONE_TOKEN_ONE", sum(k[0]=="1" and k[1]==1 and k[2]==1 for k,v in n_vs_tokens.items() for _ in range(v)))
print("BRIDGE_MULTI_RUN_ROWS", sum(v for (n,t,u),v in n_vs_tokens.items() if (n.isdigit() and int(n)>1) or t>1))
print("BRIDGE_NO_EXTRACTABLE_RUN_TOKENS", token_n.get(0,0))

h5path = REC / "stage2_output/studies/SRP074349/selected_run_coverage_G026.h5"
with h5py.File(h5path) as h:
    rids=[decode(v) for v in h["run_id"][:]]
    lids=[decode(v) for v in h["library_id"][:]]
    geneids=[decode(v) for v in h["gene_id"][:]]
    print("H5_AXES", len(geneids), len(rids), len(lids))
    assert len(rids)==len(lids)==144 and len(geneids)==63856
    by_h5 = defaultdict(set)
    for lib,run in zip(lids,rids):
        by_h5[lib].add(run)
    print("H5_DISTINCT_LIBRARIES", len(by_h5))
    report("H5_RUNS_PER_LIBRARY", Counter(len(x) for x in by_h5.values()))
    print("H5_RUN_DUPLICATES", len(rids)-len(set(rids)))
    print("BRIDGE_LIBRARIES_IN_H5", len(set(bylib) & set(by_h5)))
    print("H5_LIBRARIES_NOT_IN_BRIDGE", len(set(by_h5) - set(bylib)))
    print("BRIDGE_LIBRARIES_NOT_IN_H5", len(set(bylib) - set(by_h5)))
    complete = sum(bridge_runs[lib]==by_h5.get(lib,set()) and bool(bridge_runs[lib]) for lib in bylib)
    missing = sum(bool(by_h5.get(lib,set())-bridge_runs[lib]) for lib in bylib)
    extra = sum(bool(bridge_runs[lib]-by_h5.get(lib,set())) for lib in bylib)
    print("EXACT_BRIDGE_VS_H5_RUN_SET_MATCH", complete)
    print("LIBRARIES_WITH_H5_RUNS_ABSENT_FROM_BRIDGE_TOKENS", missing)
    print("LIBRARIES_WITH_BRIDGE_TOKENS_ABSENT_FROM_H5", extra)
    report("BRIDGE_DECLARED_VS_H5_RUN_COUNT",
           Counter((r.get("n_runs","").strip(),len(by_h5.get(lib,set()))) for lib,r in bylib.items()))

    _,mapping=tab(STAGE5 / "00_gene_row_mapping.tsv")
    print("GENE_MAP_ROWS",len(mapping))
    if len(mapping)==63856:
        matches=sum(g==m["source_gene_id"] or g.split(".",1)[0]==m["gene_base"]
                    for g,m in zip(geneids,mapping))
        print("H5_GENE_MAP_ORDER_MATCH_ROWS",matches)
        audit = ROOT/"10_preprojection_normalization_v1_0_1/08_model_gene_availability.tsv"
        if not audit.is_file():
            audit = RELEASE/"10_preprojection_normalization_v1_0_1/08_model_gene_availability.tsv"
        print("STAGE09_499_AUDIT_PRESENT",audit.is_file())
        if audit.is_file():
            fields,rows=tab(audit)
            print("STAGE09_AUDIT_ROWS",len(rows),"HAS_MODEL_GENE",("model_gene" in fields))
            if len(rows)==499 and "model_gene" in fields:
                source=RELEASE/"09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py"
                spec=importlib.util.spec_from_file_location("v7c_diagnostic_stage09",source)
                mod=importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                try:
                    result=mod.build_gene_maps(pd.Index([m["gene_symbol"] for m in mapping]),
                                               [r["model_gene"] for r in rows])
                    print("STAGE09_GENE_MODEL_RESOLVED",len(result["available_model_genes"]))
                    print("STAGE09_GENE_ALIAS_RESOLVED",int(result["model_gene_resolution"]["resolution"].eq("alias_resolved").sum()))
                    print("STAGE09_GENE_MODEL_UNRESOLVED",int(result["model_gene_resolution"]["resolution"].eq("unresolved").sum()))
                except Exception as e:
                    print("STAGE09_GENE_MODEL_AUDIT_ERROR",type(e).__name__,str(e)[:200])

_,lqc=tab(REC/"stage2_output/library_qc_all.tsv")
selected=[r for r in lqc if r["study"]=="SRP074349" and r["library_id"] in bylib]
print("STAGE2_LIBRARY_QC_MATCHED",len(selected))
for f in ("n_runs_expected","n_runs_expression_verified","run_count_status","not_collapsed"):
    report("STAGE2_"+f,Counter(r.get(f,"") for r in selected))

# Verify frozen archival presence, without loading or claiming recalibration.
for fn in ("03_pan_cancer_anchor_scaling_CPM.tsv.gz",
           "04_pan_cancer_anchor_scaling_TMM_CPM.tsv.gz",
           "05_pan_cancer_anchor_scaling_UQ_CPM.tsv.gz",
           "06_fixed_anchor_scaled_input_CPM.rds",
           "07_fixed_anchor_scaled_input_TMM_CPM.rds",
           "08_fixed_anchor_scaled_input_UQ_CPM.rds",
           "09_raw_H_CPM.tsv.gz","11_raw_H_TMM_CPM.tsv.gz","13_raw_H_UQ_CPM.tsv.gz"):
    yes=any((base/fn).is_file() for base in (
        ROOT/"11_pan_cancer_fixed_anchor_projection_v1_0_4",
        RELEASE/"11_pan_cancer_fixed_anchor_projection_v1_0_4"))
    print("FROZEN_STAGE10_ASSET",fn,"PRESENT",yes)

# No actual specimen/run identifiers are printed.
print("FAILURE_REPRODUCED_AT_BRIDGE_CARDINALITY_ONLY")
print("NO_COHORT_OR_FROZEN_FILES_MODIFIED")
print("NO_MODEL_PROJECTION")
print("CLINICAL_ADMISSIONS 0")
print("V7C_RUN_MAPPING_DIAGNOSTIC_V1_COMPLETE")
