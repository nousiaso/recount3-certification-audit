#!/usr/bin/env python3
"""V7I: read-only post-V7H2 preflight for exact 143-pair transfer inputs.

No abundance values or run IDs printed. Confirms source SHA256, 143/454 identity,
W-gene resolution in BOTH assays, and 143-versus-frozen atlas membership.
Never edits V7G, V7C or frozen Stage09-12 outputs. No full matrix extraction.
"""
from collections import Counter
from pathlib import Path
import csv, gzip, hashlib, importlib.util, re, sys
import h5py
import numpy as np
import pandas as pd

ROOT=Path(sys.argv[1]).resolve()
REC=ROOT/"26_recount3_external_primary_v1"
REL=ROOT/"RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
STAGE5=REC/"stage5_output"
S9="10_preprojection_normalization_v1_0_1"
S10="11_pan_cancer_fixed_anchor_projection_v1_0_4"
SCRIPT=ROOT/"RECOUNT3_CERTIFICATION_GIT/tools/v7g_prepare_143_matched_dual_assay_G026_ARCHS4_v1.py"
FROZEN=REL/"09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py"
BRIDGE=STAGE5/"03_bridge_165_availability.tsv"
GENEMAP=STAGE5/"00_gene_row_mapping.tsv"
G026=REC/"raw_studies/SRP074349/gene_sums/sra.gene_sums.SRP074349.G026.gz"
ARCH=ROOT/"00_download/human_gene_v2.latest.h5"
EXPECTED_RAW="b056331eb16d6ec612409fd37a266abfe10c4c15329b486dbb9df213612797c5"
EXPECTED_FROZEN="46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500"
RUNS=re.compile(r"(?:SRR|ERR|DRR)\d+",re.I)
def read_tsv(p):
    with p.open("r",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        return reader.fieldnames or [],list(reader)
def hash_file(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()
def decode(x):
    return bytes(x).decode("utf-8").strip() if isinstance(x,(bytes,np.bytes_)) else str(x).strip()
def frozen_file(stage,name):
    for base in (ROOT,REL):
        p=base/stage/name
        if p.is_file(): return p
    raise FileNotFoundError(stage+"/"+name)

print("V7I_POST_V7H2_TRANSFER_INPUT_CONTRACT_V1")
for name,p in (
    ("BRIDGE",BRIDGE),("GENE_MAP",GENEMAP),
    ("RAW_G026",G026),("FROZEN_STAGE09_CODE",FROZEN),
    ("V7G_ORIGINAL_PREPARATION_CODE",SCRIPT)):
    assert p.is_file(),f"MISSING_ASSET {name}"
    h=hash_file(p)
    print("SHA256",name,h)
    if name=="RAW_G026": assert h==EXPECTED_RAW,"Unexpected original G026 checksum"
    if name=="FROZEN_STAGE09_CODE": assert h==EXPECTED_FROZEN,"Original Stage09 source changed"

# The previous V7G script has unverified fixed expectations. Report exact mismatch.
import ast
tree=ast.parse(SCRIPT.read_text())
hardcoded={}
for node in tree.body:
    if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
        if node.targets[0].id in ("EXPECTED_BRIDGE_SHA","EXPECTED_MAP_SHA") and isinstance(node.value,ast.Constant):
            hardcoded[node.targets[0].id]=node.value.value
assert set(hardcoded)=={"EXPECTED_BRIDGE_SHA","EXPECTED_MAP_SHA"}
print("V7G_V1_BRIDGE_SHA_MATCH",hardcoded["EXPECTED_BRIDGE_SHA"]==hash_file(BRIDGE))
print("V7G_V1_GENEMAP_SHA_MATCH",hardcoded["EXPECTED_MAP_SHA"]==hash_file(GENEMAP))

cols,b=read_tsv(BRIDGE)
pairs=[r for r in b if r["study"]=="SRP074349"]
assert len(b)==165 and len(pairs)==143
assert set(r["bridge_class"] for r in pairs)=={"CALIBRATION_ONLY_NOT_IN_PRIMARY"}
gsms={r["ARCHS4_GSM"].strip() for r in pairs}
srx={r["library_id"].strip() for r in pairs}
bios={r["biosample"].strip() for r in pairs}
all_runs=[r for p in pairs for r in RUNS.findall(p["recount3_run_ids"] or "")]
assert len(gsms)==len(srx)==len(bios)==143
assert len(all_runs)==len(set(all_runs))==454
assert Counter(int(r["n_runs"]) for r in pairs)==Counter({1:84,2:4,3:22,8:1,9:32})
print("EXACT_CALIBRATION_IDENTITY_PASS","SRX=143","GSM=143","SRR=454","BIOSAMPLE=143")

mapcols,mapped=read_tsv(GENEMAP)
assert len(mapped)==63856
assert [int(r["row_zero_based"]) for r in mapped]==list(range(63856))
print("G026_SOURCE_GENE_MAP_63856_PASS")

_,geneaudit=read_tsv(frozen_file(S9,"08_model_gene_availability.tsv"))
assert len(geneaudit)==499
model_genes=[r["model_gene"].strip() for r in geneaudit]
assert len(set(model_genes))==499

spec=importlib.util.spec_from_file_location("frozen_original_stage09",FROZEN)
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
raw_gene_symbols=[r["gene_symbol"] for r in mapped]
rec_map=old.build_gene_maps(pd.Index(raw_gene_symbols),model_genes)
assert len(rec_map["available_model_genes"])==499
print("RECOUNT3_W_GENE_MAPPING_PASS",
      "RESOLVED",len(rec_map["available_model_genes"]),
      "ALIASES",int(rec_map["model_gene_resolution"]["resolution"].eq("alias_resolved").sum()))

index_by_gsm={r["ARCHS4_GSM"].strip():int(r["ARCHS4_expression_column_index"]) for r in pairs}
with h5py.File(ARCH,"r") as h:
    accessor=h["meta/samples/geo_accession"]
    assert all(0<=ix<len(accessor) and decode(accessor[ix])==gsm
               for gsm,ix in index_by_gsm.items()),"ARCHS4 GSM index mismatch"
    source_symbols=old.decode_string_array(h["meta/genes/symbol"][:])
    assert len(source_symbols)==67186
    arch_map=old.build_gene_maps(pd.Index(source_symbols),model_genes)
    assert len(arch_map["available_model_genes"])==499
print("ARCHS4_GSM_AND_GENE_MAPPING_PASS",
      "GSM",143,
      "MODEL_GENES",len(arch_map["available_model_genes"]),
      "ALIASES",int(arch_map["model_gene_resolution"]["resolution"].eq("alias_resolved").sum()))

canonical=frozen_file(S9,"01_canonical_library_manifest.tsv.gz")
with gzip.open(canonical,"rt",encoding="utf-8-sig",newline="") as f:
    reader=csv.DictReader(f,delimiter="\t")
    original_ids={r["sample_id"] for r in reader}
assert len(original_ids)==9465
print("MATCHED_GSMS_IN_ORIGINAL_9465_ATLAS",len(gsms & original_ids))
print("MATCHED_GSMS_OUTSIDE_ORIGINAL_9465_ATLAS",len(gsms-original_ids))

# Check original Stage10 frozen W mapping audit and direct model replay status.
mcols,ledger=read_tsv(frozen_file(S10,"16A_frozen_W_gene_symbol_harmonization.tsv"))
assert len(ledger)==1497 and Counter(r["branch"] for r in ledger)==Counter({"CPM":499,"TMM_CPM":499,"UQ_CPM":499})
print("FROZEN_STAGE10_W_MAPPING_LEDGER_PASS","499_PER_BRANCH")

# Explicitly locate the replay output without conflating presence with execution.
replay=ROOT/"RECOUNT3_CERTIFICATION_GIT/diagnostics/20261127_v7h2_ARCHS4_9465_original_numeric_replay_v2/03_ARCHS4_9465_numerical_replay_stdout.txt"
print("V7H2_REPLAY_LOG_PRESENT",replay.is_file())
if replay.is_file():
    t=replay.read_text(encoding="utf-8",errors="replace")
    checks=("V7H2_ORIGINAL_ARCHS4_FULL_NUMERICAL_REPLAY_COMPLETE",
            "V7H2_FROZEN_STAGE09_NORMALIZATION_AND_STAGE10_Z_REPLAY_PASS",
            "STAGE10_ORIGINAL_H_REPLAY CPM",
            "STAGE10_ORIGINAL_H_REPLAY TMM_CPM",
            "STAGE10_ORIGINAL_H_REPLAY UQ_CPM")
    print("V7H2_COMPLETION_AND_ALL_BRANCH_LABELS_PRESENT",all(x in t for x in checks))

print("NO_EXPRESSION_COUNTS_LOADED_NO_EXTRACTION")
print("NO_UNAPPROVED_V7G_RUN_NO_SOURCE_MODIFICATIONS")
print("CLINICAL_ADMISSIONS 0")
print("V7I_POST_V7H2_TRANSFER_INPUT_CONTRACT_V1_COMPLETE")
