#!/usr/bin/env Rscript
# V7N v1 — read-only four-feature impact diagnostic on original ARCHS4 CPM.
# Observe the exact features lost from the 143-sample directional subset and
# their contribution to the original frozen H's KL reconstruction objective.
# This is an EVALUATION of archived H, not optimization nor causal attribution.
suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
  library(SpatialEcoTyper)
})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: Rscript v7n.R ROOT",call.=FALSE)
ROOT <- normalizePath(args[[1L]],mustWork=TRUE)
REPO <- file.path(ROOT,"RECOUNT3_CERTIFICATION_GIT")
REL <- file.path(ROOT,"RELEASES","ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0","PROJECT_TREE")
frozen <- function(stage,name) {
  p <- c(file.path(ROOT,stage,name),file.path(REL,stage,name))
  x <- p[file.exists(p)]
  if (!length(x)) stop("FROZEN_ASSET_MISSING ",stage,"/",name,call.=FALSE)
  x[[1L]]
}
assert <- function(ok,msg) if (!isTRUE(ok)) stop(msg,call.=FALSE)
assert(as.character(packageVersion("SpatialEcoTyper"))=="1.0.3",
       "Unexpected SpatialEcoTyper version")
W <- readRDS("/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds")
assert(is.matrix(W)&&identical(dim(W),c(500L,10L)),
       "Frozen W dimensional contract mismatch")
w_na <- sum(is.na(W))
W[is.na(W)] <- 0
assert(all(is.finite(W))&&all(W>=0),
       "Frozen W cannot be interpreted as nonnegative component weights")
Wgenes <- unique(sub("__.*$","",rownames(W)))
assert(length(Wgenes)==499L, "Unexpected W gene count")
Z <- readRDS(frozen("11_pan_cancer_fixed_anchor_projection_v1_0_4",
                    "06_fixed_anchor_scaled_input_CPM.rds"))
assert(is.matrix(Z)&&identical(dim(Z),c(499L,9465L)) &&
       identical(rownames(Z),Wgenes) &&
       !anyDuplicated(colnames(Z)), "Frozen Z axes invalid")
manifest <- fread(file.path(REPO,"private_data",
    "v7j_143_frozen_archs4_vs_recount3_full499_v1",
    "143_pair_manifest_LOCAL_ONLY.tsv"),
    sep="\t",quote="",showProgress=FALSE)
assert(nrow(manifest)==143L && "ARCHS4_GSM" %in% names(manifest),
       "V7J frozen 143 GSM manifest missing")
ids <- as.character(manifest$ARCHS4_GSM)
assert(!anyDuplicated(ids)&&all(ids %in% colnames(Z)),
       "V7J 143 GSM manifest not aligned with frozen Z")
Z143 <- Z[,ids,drop=FALSE]
rm(Z); invisible(gc(verbose=FALSE))

# Match the original Stage10/.nmf.predict pre-zero-variance-filter
# directional expansion and naming conventions; keep all rows for audit.
v <- as.matrix(NMF::posneg(Z143[Wgenes,,drop=FALSE]))
dup <- duplicated(rownames(v))
rownames(v)[!dup] <- paste0(rownames(v)[!dup],"__pos")
rownames(v)[dup] <- paste0(rownames(v)[dup],"__neg")
varrow <- apply(v,1L,var)
eligible <- is.finite(varrow)&varrow>0
features <- intersect(rownames(v)[eligible],rownames(W))
missing <- setdiff(rownames(W),features)
expected_missing <- c("F12__pos","JUND__pos","PSMA3__neg","SURF4__neg")
assert(setequal(missing,expected_missing) &&
       length(missing)==4L && length(features)==496L,
       "Expected original V7M four-feature pattern changed; STOP")
for (k in missing)
  assert(!is.na(match(k,rownames(v))),
         paste("Expected directional row missing",k))

Harchive <- fread(frozen("11_pan_cancer_fixed_anchor_projection_v1_0_4",
                        "09_raw_H_CPM.tsv.gz"),sep="\t",
                  check.names=FALSE,showProgress=FALSE)
comps <- colnames(W)
assert(nrow(Harchive)==9465L &&
       all(c("sample_id",comps) %in% names(Harchive)) &&
       !anyDuplicated(Harchive$sample_id),
       "Original archived H matrix contract invalid")
hi <- match(ids,as.character(Harchive$sample_id))
assert(!anyNA(hi),"143 original GSM IDs absent in frozen H")
H <- as.matrix(Harchive[hi,..comps])
storage.mode(H) <- "double"
rownames(H) <- ids
colnames(H) <- comps
assert(identical(dim(H),c(143L,10L))&&
       all(is.finite(H))&&all(H>=0),
       "Archived H values invalid")
V <- v[rownames(W),,drop=FALSE]
assert(identical(dim(V),c(500L,143L)) &&
       all(is.finite(V))&&all(V>=0),
       "Frozen directional V invalid")
Vhat <- W %*% t(H)
dimnames(Vhat) <- dimnames(V)
assert(all(is.finite(Vhat))&&all(Vhat>=0),
       "Frozen-H reconstructed V invalid")
eps <- 1e-12
K <- Vhat-V
positive <- V>0
K[positive] <- V[positive]*log(V[positive]/pmax(Vhat[positive],eps)) -
               V[positive]+Vhat[positive]
assert(all(is.finite(K)) && min(K)>=-1e-8,
       "Invalid frozen KL objective calculation")
K <- pmax(K,0)
Kall <- colSums(K)
Kmissing <- colSums(K[missing,,drop=FALSE])
pct <- ifelse(Kall>0,100*Kmissing/Kall,NA_real_)
fmt <- function(x) format(x,digits=12,scientific=TRUE)
q95 <- function(x) unname(quantile(x,.95,na.rm=TRUE))
cat("V7N_FROZEN_143_ORIGINAL_W_FEATURE_IMPACT_V1\n")
cat("MODEL_W_DIM",nrow(W),ncol(W),"ORIGINAL_W_NA_TO_ZERO",w_na,"\n")
cat("ARCHIVED_H_SAMPLES",nrow(H),"COMPONENTS",ncol(H),"\n")
cat("DIRECTIONAL_ROWS_ALL",nrow(v),
    "NONVARIABLE_ALL_DIRECTIONAL",sum(!eligible),
    "ORIGINAL_W_FEATURES_RETAINED",length(features),
    "ORIGINAL_W_FEATURES_MISSING",length(missing),"\n")
cat("MISSING_W_FEATURE_IDS",paste(missing,collapse=";"),"\n")
for (f in missing) {
  observed <- as.numeric(V[f,])
  predicted <- as.numeric(Vhat[f,])
  weights <- as.numeric(W[f,])
  kk <- as.numeric(K[f,])
  cat("FEATURE",f,
      "MEASURED_MIN",fmt(min(observed)),
      "MEASURED_MAX",fmt(max(observed)),
      "MEASURED_MEAN",fmt(mean(observed)),
      "MEASURED_NONZERO",sum(observed>0),
      "W_NONZERO_COMPONENTS",sum(weights>0),
      "W_L1",fmt(sum(weights)),
      "W_MAX",fmt(max(weights)),
      "FROZEN_H_PREDICTED_MEDIAN",fmt(median(predicted)),
      "FROZEN_H_KL_MEDIAN",fmt(median(kk)),
      "FROZEN_H_KL_MAX",fmt(max(kk)),"\n")
}
cat("FROZEN_H_FULL500_KL_MEDIAN",fmt(median(Kall)),
    "P95",fmt(q95(Kall)),"\n")
cat("FROZEN_H_MISSING4_KL_MEDIAN",fmt(median(Kmissing)),
    "P95",fmt(q95(Kmissing)),
    "MAX",fmt(max(Kmissing)),"\n")
cat("MISSING4_PERCENT_OF_FULL500_KL_MEDIAN",fmt(median(pct,na.rm=TRUE)),
    "P95",fmt(q95(pct)),
    "MAX",fmt(max(pct,na.rm=TRUE)),"\n")
cat("IMPORTANT_MISSING4_KL_PENALTY_IS_EVALUATED_AT_ARCHIVED_H_NOT_CAUSAL_REESTIMATION\n")
cat("NO_NMF_OPTIMIZATION_OR_RECOUNT3_PROJECTION\n")
cat("NO_FROZEN_FILES_MODIFIED\nCLINICAL_ADMISSIONS 0\n")
cat("V7N_FROZEN_143_ORIGINAL_W_FEATURE_IMPACT_V1_COMPLETE\n")
