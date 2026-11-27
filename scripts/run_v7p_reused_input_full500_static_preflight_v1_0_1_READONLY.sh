#!/usr/bin/env bash
# V7P static preflight only: source parse, original V7J input integrity,
# exact V7O numerical decision, no R model fitting nor raw recount3 extraction.
set -euo pipefail
umask 077
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
RSCRIPT=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/Rscript
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
SOURCE="$REPO/private_data/v7j_143_frozen_archs4_vs_recount3_full499_v1"
V7O="$REPO/diagnostics/20261127_v7o_143_frozen_ARCHS4_500_vs_496_KL_v1/01_500_vs_496_KL_stdout.txt"
SCRIPT="$REPO/tools/v7p_frozen_ARCHS4_vs_recount3_143pairs_full500_KL_v1_0_1.R"
PRE="$REPO/diagnostics/20261127_v7p_static_preflight_v1_0_1"
OUTPUT="$REPO/private_data/v7p_143_frozen_ARCHS4_vs_recount3_full500_KL_v1_0_1"
AUDIT="$REPO/diagnostics/20261127_v7p_143_frozen_ARCHS4_vs_recount3_full500_KL_v1_0_1"

test -x "$PYTHON" && test -x "$RSCRIPT" && test -f "$SCRIPT"
test -f "$V7O" && test -f "$SOURCE/input_provenance.json"
test ! -e "$OUTPUT" && test ! -e "$AUDIT" && test ! -e "$PRE"
grep -q '^V7O_CONVEX_KL_OPTIMUM_STABILITY_PASS$' "$V7O"
grep -q '^V7O_FROZEN_143_FIXED_W_KL_FEATURE_AND_OPTIMIZER_DIAGNOSTIC_COMPLETE$' "$V7O"
mkdir -p "$PRE"

"$RSCRIPT" -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("V7P_R_SYNTAX_PASS\n")' "$SCRIPT" 2>&1 | tee "$PRE/01_R_syntax_stdout.txt"
"$RSCRIPT" -e 'suppressPackageStartupMessages({library(data.table);library(NMF);library(SpatialEcoTyper);library(jsonlite)});stopifnot(as.character(packageVersion("SpatialEcoTyper"))=="1.0.3");cat("V7P_R_PACKAGES_PASS\n")' 2>&1 | tee "$PRE/02_R_packages_stdout.txt"
"$PYTHON" - "$SOURCE" <<'PY' 2>&1 | tee "$PRE/03_V7J_reused_inputs_stdout.txt"
import hashlib,json,sys
from pathlib import Path
src=Path(sys.argv[1])
prov=json.loads((src/"input_provenance.json").read_text())
if prov.get("version")!="v7j_frozen_archs4_vs_recount3_full499_v1":raise RuntimeError("Unexpected source version")
if (prov.get("n_paired_specimens"),prov.get("n_matched_technical_runs"),prov.get("n_model_genes"))!=(143,454,499):
    raise RuntimeError("Prepared source provenance dimensions changed")
files=prov.get("outputs_sha256")
if not isinstance(files,dict) or len(files)<4:raise RuntimeError("Missing prepared source checksums")
for filename,digest in sorted(files.items()):
    p=src/filename
    if not p.is_file():raise RuntimeError(f"Prepared output missing: {filename}")
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    if h.hexdigest()!=digest:raise RuntimeError(f"Checksum mismatch {filename}")
for b in ("CPM","TMM_CPM","UQ_CPM"):
    if not (src/f"499_recount3_{b}_LOCAL_ONLY.tsv.gz").is_file():
        raise RuntimeError(f"Missing reuse-ready {b} matrix")
if not (src/"143_pair_manifest_LOCAL_ONLY.tsv").is_file():
    raise RuntimeError("Missing matching GSM manifest")
print("V7P_V7J_PREPARED_INPUT_REUSE_SHA256_PASS",len(files),"files")
print("V7P_NO_RECOUNT3_G026_REEXTRACTION")
print("V7P_MATCHED_SPECIMENS_143_RUNS_454_GENES_499")
PY
"$RSCRIPT" --vanilla - "$SCRIPT" <<'RS' 2>&1 | tee "$PRE/04_assay_key_regression_stdout.txt"
suppressPackageStartupMessages(library(data.table))
exprs <- parse(file=commandArgs(trailingOnly=TRUE)[[1L]])
env <- new.env()
for (nm in c("COMPONENTS","CLR_FLOOR","clr_matrix","dominant_component",
             "normalization_sensitivity","assay_storage_key")) {
  hit <- Filter(function(e) is.call(e) && identical(e[[1L]],as.name("<-")) &&
                length(e)==3L && identical(e[[2L]],as.name(nm)),as.list(exprs))
  stopifnot(length(hit)==1L)
  eval(hit[[1L]],envir=env)
}
with(env,{
  f <- outer(seq_len(143L),seq_len(10L),function(i,j) 1+((i*j+j*j)%%149))
  f <- f/rowSums(f)
  dimnames(f) <- list(paste0("synthetic",seq_len(143L)),COMPONENTS)
  alt <- f[,c(2:10,1)]
  colnames(alt) <- COMPONENTS
  x <- list(CPM=list(arch=list(f=f),recount3=list(f=f)),
            TMM_CPM=list(arch=list(f=alt),recount3=list(f=alt)),
            UQ_CPM=list(arch=list(f=alt),recount3=list(f=alt)))
  stopifnot(is.null(x[["CPM"]][["ARCHS4"]]),
            identical(assay_storage_key("ARCHS4"),"arch"),
            identical(assay_storage_key("recount3"),"recount3"))
  for (assay in c("ARCHS4","recount3")) {
    k <- assay_storage_key(assay)
    for (branch in c("TMM_CPM","UQ_CPM")) {
      mask <- rep(TRUE,143L)
      if (branch=="UQ_CPM") mask[143L] <- FALSE
      result <- normalization_sensitivity(x[["CPM"]][[k]]$f,
                 x[[branch]][[k]]$f,paste0("CPM_vs_",branch),evaluable_mask=mask)
      stopifnot(nrow(result$sample)==143L,nrow(result$component)==10L,
                sum(result$sample$comparison_evaluable)==sum(mask),
                all(is.finite(result$sample$total_variation_distance[mask])),
                all(is.na(result$sample$total_variation_distance[!mask])))
    }
  }
})
cat("V7P_V1_0_1_ASSAY_KEY_REGRESSION_PASS\n")
RS
(
 cd "$PRE"
 sha256sum "$SCRIPT" ./*.txt > SHA256SUMS.txt
 sha256sum -c SHA256SUMS.txt
)
echo "V7P_V1_0_1_STATIC_INPUT_PREFLIGHT_CAPTURED_PENDING_REVIEW"
echo "NO_PROJECTION_EXECUTED_NO_CLINICAL_ADMISSIONS"
