#!/usr/bin/env bash
set -euo pipefail
umask 077
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
ENV=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human
PYTHON="$ENV/bin/python"
RSCRIPT="$ENV/bin/Rscript"
STUDY="${1:?Supply configured study accession}"
case "$STUDY" in SRP074349|SRP157215) ;; *) echo "Unknown study" >&2; exit 1;; esac
CONFIG="$REPO/configs/v7r_matched_calibration_studies_v1.json"
PRE="$REPO/diagnostics/20261127_v7r_${STUDY}_preflight_v1"
AUDIT="$REPO/diagnostics/20261127_v7r_${STUDY}_calibration_v1"
SOURCE="$REPO/private_data/v7r_${STUDY}_prepared_v1"
OUT="$REPO/private_data/v7r_${STUDY}_projection_v1"
P="$REPO/tools/v7r_prepare_configured_matched_G026_full499_v1.py"
R="$REPO/tools/v7r_project_configured_matched_full500_KL_v1.R"
RAW="$ROOT/26_recount3_external_primary_v1/raw_studies/$STUDY/gene_sums/sra.gene_sums.$STUDY.G026.gz"
cd "$REPO"

for d in "$PRE" "$AUDIT" "$SOURCE" "$OUT"; do test ! -e "$d"; done
for f in "$REPO/scripts/run_v7r_configured_static_preflight_v1.sh" "$REPO/scripts/run_v7r_configured_matched_calibration_v1.sh" "$REPO/scripts/slurm_v7r_configured_matched_calibration_v1.sbatch"; do bash -n "$f"; done
mkdir "$PRE"
"$PYTHON" - "$P" "$CONFIG" "$STUDY" "$RAW" <<'PY' | tee "$PRE/01_python_contract.txt"
import ast,hashlib,json,sys
from pathlib import Path
p,c,s,raw=sys.argv[1:]
ast.parse(Path(p).read_text())
cfg=json.loads(Path(c).read_text())
assert cfg["version"]=="v7r_study_contract_v1"
v=cfg["studies"][s]
assert sum(int(k)*n for k,n in v["multiplicity"].items())==v["n_runs"]
assert sum(v["multiplicity"].values())==v["n_pairs"]
h=hashlib.sha256()
with open(raw,"rb") as f:
    for b in iter(lambda:f.read(1048576),b""):h.update(b)
assert h.hexdigest()==v["raw_sha256"]
import h5py,numpy,pandas
print("V7R_PYTHON_AND_RAW_CONTRACT_PASS",s,v["n_pairs"],v["n_runs"],v["n_raw_runs"])
PY
"$RSCRIPT" --vanilla -e 'invisible(parse(file=commandArgs(trailingOnly=TRUE)[[1]]));cat("V7R_R_SYNTAX_PASS\n")' "$R" | tee "$PRE/02_R_syntax.txt"
"$RSCRIPT" --vanilla -e 'suppressPackageStartupMessages({library(data.table);library(NMF);library(SpatialEcoTyper);library(jsonlite)});stopifnot(as.character(packageVersion("SpatialEcoTyper"))=="1.0.3");cat("V7R_R_PACKAGES_PASS\n")' | tee "$PRE/03_R_packages.txt"
"$RSCRIPT" --vanilla - "$R" <<'RS' 2>&1 | tee "$PRE/04_assay_key_regression_stdout.txt"
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
  f <- outer(seq_len(22L),seq_len(10L),function(i,j) 1+((i*j+j*j)%%149))
  f <- f/rowSums(f)
  dimnames(f) <- list(paste0("synthetic",seq_len(22L)),COMPONENTS)
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
      mask <- rep(TRUE,22L)
      if (branch=="UQ_CPM") mask[22L] <- FALSE
      result <- normalization_sensitivity(x[["CPM"]][[k]]$f,
                 x[[branch]][[k]]$f,paste0("CPM_vs_",branch),evaluable_mask=mask)
      stopifnot(nrow(result$sample)==22L,nrow(result$component)==10L,
                sum(result$sample$comparison_evaluable)==sum(mask),
                all(is.finite(result$sample$total_variation_distance[mask])),
                all(is.na(result$sample$total_variation_distance[!mask])))
    }
  }
})
cat("V7R_ASSAY_KEY_REGRESSION_PASS\n")
RS
(
 cd "$PRE"
 sha256sum "$P" "$R" "$CONFIG" "$REPO/scripts/run_v7r_configured_static_preflight_v1.sh" "$REPO/scripts/run_v7r_configured_matched_calibration_v1.sh" "$REPO/scripts/slurm_v7r_configured_matched_calibration_v1.sbatch" "$RAW" ./*.txt > SHA256SUMS.txt
 sha256sum -c SHA256SUMS.txt
)
echo "V7R_STATIC_PREFLIGHT_PASS_PENDING_REVIEW"
echo "NO_EXTRACTION_NO_MODEL_EXECUTION_CLINICAL_ADMISSIONS_0"
