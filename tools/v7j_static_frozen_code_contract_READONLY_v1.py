#!/usr/bin/env python3
"""V7J static fidelity preflight: READ ONLY. No expression values processed.

Confirms original Stage09/Stage10 immutable scripts, exact transcribed Stage10
R helper bodies, new V7J Python AST, source hash locks, model code structure.
"""
import ast
import hashlib
import sys
from pathlib import Path

ROOT=Path(sys.argv[1]).resolve()
REPO=ROOT/"RECOUNT3_CERTIFICATION_GIT"
REL=ROOT/"RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE"
S9=REL/"09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py"
S10=REL/"10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R"
PREP=REPO/"tools/v7j_prepare_recount3_454runs_143pairs_frozen_ARCHS4_v1.py"
PROJ=REPO/"tools/v7j_project_recount3_vs_original_frozen_ARCHS4_143pairs_v1.R"
FILES={
    "STAGE09_SOURCE":(S9,"46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500"),
    "STAGE10_SOURCE":(S10,"16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c"),
    "BRIDGE":(ROOT/"26_recount3_external_primary_v1/stage5_output/03_bridge_165_availability.tsv",
              "c7d6d7c754e3ca77d62424f26627a40ceaec10fac98a0741c945773f3d3b1b67"),
    "GENE_MAP":(ROOT/"26_recount3_external_primary_v1/stage5_output/00_gene_row_mapping.tsv",
                "68be266dd6bd30bfa79a2bd8ac79f0e49c6f5d0f19d2d083e2ca57850ef17391"),
    "ORIGINAL_G026":(ROOT/"26_recount3_external_primary_v1/raw_studies/SRP074349/gene_sums/sra.gene_sums.SRP074349.G026.gz",
                     "b056331eb16d6ec612409fd37a266abfe10c4c15329b486dbb9df213612797c5"),
}
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):h.update(data)
    return h.hexdigest()

def rfunction(source,name):
    begin=source.find(name+" <- function(")
    assert begin!=-1, f"Cannot find original function {name}"
    end=source.find("\n}\n",begin)
    assert end!=-1, f"Cannot delimit original function {name}"
    return source[begin:end+2]

print("V7J_ORIGINAL_STAGE09_STAGE10_FUNCTION_FIDELITY_PREFLIGHT")
for label,(path,sha) in FILES.items():
    assert path.is_file(),f"Missing expected source {label}"
    actual=digest(path)
    assert actual==sha,f"Frozen/source checksum changed: {label}"
    print("HASH_MATCH",label,actual)

prep=PREP.read_text(encoding="utf-8")
ast.parse(prep)
required_py=(
  "old.build_gene_maps(","old.collapse_batch(","old.quantile75_fraction(",
  "old.select_tmm_reference(","old.tmm_factor(","old.geometric_mean_positive(",
  "EXPECTED_RAW_SHA","G026_SOURCE_GENE_ROWS_PASS",
  "ARCHS4_BASELINE_ARCHIVED_FROZEN_NO_RENORMALIZATION",
)
assert all(x in prep for x in required_py),"New recount3 prep missing original Stage09 stages"
assert "arch_raw=" not in prep and "arch_stats=" not in prep, "ARCHS4 recalculation still present"
print("V7J_PREPARATION_PYTHON_AST_PASS")
print("V7J_PREPARATION_ARCHS4_UNCHANGED_REFERENCE_PASS")
print("V7J_PREPARATION_SHA256",digest(PREP))

original=S10.read_text(encoding="utf-8")
new=PROJ.read_text(encoding="utf-8")
function_names=(
  "score_fractions","clr_matrix","dominant_component",
  "harmonize_expression_to_W","make_directional_bundle",
  "project_raw_h","calculate_reconstruction_adequacy",
  "normalization_sensitivity",
)
for name in function_names:
    f=rfunction(original,name)
    assert new.count(f)==1,f"Original Stage10 helper not byte-identical: {name}"
    print("ORIGINAL_STAGE10_FUNCTION_EXACT",name,len(f))
assert all(x in new for x in (
    "V7J_FROZEN_ARCHS4_BASELINE_143_IN_9465_PASS",
    "V7J_ORIGINAL_ARCHS4_FROZEN_Z_REPLAY",
    "V7J_ARCHIVED_ARCHS4_H_FRACTION_PASS",
    "V7J_ARCHS4_FROZEN_PRIMARY_ADEQUACY_REPLAY_MAXABS",
    "V7J_RECOUNT3_DIRECTIONAL_FEATURES",
    "V7J_RECOUNT3_ORIGINAL_PROJECTION_PASS",
    "V7J_RECOUNT3_ORDER_SENTINEL",
    "V7J_FROZEN_ARCHS4_VS_RECOUNT3_143_PAIRS_499_GENES_COMPLETED_PENDING_REVIEW",
)), "New R projection missing original-baseline safety gates"
assert "boot_median" not in new, "Specimen-level bootstrap unsupported for unknown shared patient identities"
print("V7J_PROJECTION_SOURCE_SHA256",digest(PROJ))
print("STAGE10_AUTHENTIC_FUNCTIONS_IDENTICAL",len(function_names))
print("STATIC_PARSER_R_WILL_BE_RUN_BY_WRAPPER")
print("NO_DATA_FILES_MODIFIED")
print("NO_MATRIX_VALUES_LOADED_OR_MODEL_EXECUTED")
print("NO_CLINICAL_ADMISSIONS")
print("V7J_STATIC_FIDELITY_CONTRACT_PASS")
