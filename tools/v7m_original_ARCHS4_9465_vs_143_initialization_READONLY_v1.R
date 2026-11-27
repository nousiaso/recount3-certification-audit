#!/usr/bin/env Rscript
# V7M: frozen ARCHS4 9465-vs-143 cohort-dependence diagnostic, READ ONLY.
# Measures DECONVOLUTESE auto-log predicate, variable directional W features,
# and NMF::rnmf(seed=39) initialization of identical 143 GSMs in 143 vs 9465.
# Does NOT run .nmf.predict, DeconvoluteSE, or overwrite any source artifacts.
suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
  library(SpatialEcoTyper)
})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
ROOT <- normalizePath(args[[1L]],mustWork=TRUE)
REPO <- file.path(ROOT,"RECOUNT3_CERTIFICATION_GIT")
REL <- file.path(ROOT,"RELEASES","ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0","PROJECT_TREE")
file_frozen <- function(stage,name) {
  p <- c(file.path(ROOT,stage,name),file.path(REL,stage,name))
  f <- p[file.exists(p)]
  if (!length(f)) stop("FROZEN_FILE_NOT_FOUND ",stage,"/",name,call.=FALSE)
  f[[1L]]
}
assert <- function(test,msg) if (!isTRUE(test)) stop(msg,call.=FALSE)
assert(identical(as.character(packageVersion("SpatialEcoTyper")),"1.0.3"),
       "Unexpected installed SpatialEcoTyper version")
W <- readRDS("/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds")
assert(is.matrix(W)&&identical(dim(W),c(500L,10L)),"Frozen W axes differ")
Wgenes <- unique(sub("__.*$","",rownames(W)))
assert(length(Wgenes)==499L,"Frozen W model-gene identity unexpected")
Z <- readRDS(file_frozen("11_pan_cancer_fixed_anchor_projection_v1_0_4",
                          "06_fixed_anchor_scaled_input_CPM.rds"))
assert(is.matrix(Z) && identical(dim(Z),c(499L,9465L)) &&
       identical(rownames(Z),Wgenes) &&
       !anyDuplicated(colnames(Z)),"Original full-9465 scaled Z invalid")
manifest <- fread(file.path(REPO,"private_data",
                "v7j_143_frozen_archs4_vs_recount3_full499_v1",
                "143_pair_manifest_LOCAL_ONLY.tsv"),
                sep="\t",quote="",showProgress=FALSE)
assert(nrow(manifest)==143L && "ARCHS4_GSM" %in% names(manifest),
       "V7J prepared 143-pair manifest missing/invalid")
ids <- as.character(manifest$ARCHS4_GSM)
assert(!anyDuplicated(ids) && all(ids %in% colnames(Z)),
       "Matched GSM IDs not unique or not in archived 9465 reference")
Z143 <- Z[,ids,drop=FALSE]

cat("V7M_FROZEN_ARCHS4_9465_VS_143_NMF_CONTEXT_DIAGNOSTIC\n")
cat("SPATIALECOTYPER_VERSION 1.0.3\n")
cat("W_DIMENSIONS",nrow(W),ncol(W),"\n")
cat("Z_MATRIX_DIMENSIONS",nrow(Z),ncol(Z),"\n")
cat("MATCHED_SAMPLES",ncol(Z143),"\n")
for (nm in c("FULL_9465","MATCHED_143")) {
  m <- if (nm=="FULL_9465") Z else Z143
  mn <- min(m);mx <- max(m)
  auto <- all(m>=0) && mx>80
  cat("WRAPPER_AUTOLOG",nm,
      "MIN",format(mn,digits=9,scientific=TRUE),
      "MAX",format(mx,digits=9,scientific=TRUE),
      "ALL_NONNEGATIVE",all(m>=0),
      "WOULD_AUTOLOG",auto,"\n")
}
directional <- function(m) {
  original_gene_order <- Wgenes
  use <- m[original_gene_order,,drop=FALSE]
  x <- as.matrix(NMF::posneg(use))
  duplicate <- duplicated(rownames(x))
  rownames(x)[!duplicate] <- paste0(rownames(x)[!duplicate],"__pos")
  rownames(x)[duplicate] <- paste0(rownames(x)[duplicate],"__neg")
  variances <- apply(x,1L,var)
  valid <- is.finite(variances) & variances>0
  features <- intersect(rownames(x)[valid],rownames(W))
  list(
    V=x[match(features,rownames(x)),,drop=FALSE],
    features=features,
    n_zero_or_nonfinite=sum(!valid),
    n_missing_W=length(setdiff(rownames(W),features)),
    missing_W=setdiff(rownames(W),features),
    var=variances
  )
}
a <- directional(Z)
b <- directional(Z143)
cat("FULL_9465_DIRECTIONAL_RETAINED",length(a$features),
    "MISSING_W",a$n_missing_W,
    "NONVARIABLE_OR_INVALID_ALL",a$n_zero_or_nonfinite,"\n")
cat("MATCHED_143_DIRECTIONAL_RETAINED",length(b$features),
    "MISSING_W",b$n_missing_W,
    "NONVARIABLE_OR_INVALID_ALL",b$n_zero_or_nonfinite,"\n")
cat("DIRECTIONAL_FEATURE_SETS_IDENTICAL",
    identical(a$features,b$features),"\n")
cat("MATCHED_143_W_FEATURES_ABSENT",
    if(length(b$missing_W))paste(b$missing_W,collapse=";") else "none","\n")
cat("FULL_9465_W_FEATURES_ABSENT",
    if(length(a$missing_W))paste(a$missing_W,collapse=";") else "none","\n")
assert(length(a$features)>=5L&&length(b$features)>=5L,
       "Insufficient directional features for initialization inspection")

# Faithfully duplicate only ORIGINAL .nmf.predict random H initialization.
# Important: both NMF::rnmf calls must reset the seed, as .nmf.predict does.
# No KL iterations are performed here.
set.seed(39)
first <- NMF::rnmf(ncol(W),a$V)
set.seed(39)
second <- NMF::rnmf(ncol(W),b$V)
Hfull <- as.matrix(first@H)
H143 <- as.matrix(second@H)
assert(identical(dim(Hfull),c(10L,9465L)) &&
       identical(dim(H143),c(10L,143L)),
       "rnmf initial H dimensions unexpected")
idx <- match(ids,colnames(a$V))
assert(!anyNA(idx),"Manifest IDs not in full directional input")
Hfull_matched <- Hfull[,idx,drop=FALSE]
assert(all(is.finite(Hfull_matched))&&all(is.finite(H143)),
       "Random initialization contains nonfinite values")
d <- abs(Hfull_matched-H143)
a_frac <- sweep(Hfull_matched,2L,colSums(Hfull_matched),"/")
b_frac <- sweep(H143,2L,colSums(H143),"/")
tvd <- colSums(abs(a_frac-b_frac))/2
cat("RANDOM_INITIALIZER_FIXED_SEED 39\n")
cat("MATCHED_143_INITIAL_H_MAX_ABS_DIFFERENCE",
    format(max(d),digits=12,scientific=TRUE),"\n")
cat("MATCHED_143_INITIAL_H_MEDIAN_ABS_DIFFERENCE",
    format(median(d),digits=12,scientific=TRUE),"\n")
cat("MATCHED_143_INITIAL_H_COLUMNS_DIFFERING_GT_1E_12",
    sum(colSums(d)>1e-12),"OF",length(ids),"\n")
cat("MATCHED_143_INITIAL_FRACTION_TVD_MEDIAN",
    format(median(tvd),digits=12,scientific=TRUE),"\n")
cat("MATCHED_143_INITIAL_FRACTION_TVD_MAX",
    format(max(tvd),digits=12,scientific=TRUE),"\n")

cat("STOPPING_RULE_IN_INSTALLED_SOURCE connectivity\n")
cat("STOPPING_RULE_DEPENDENCE_NOT_QUANTIFIED_BY_THIS_DIAGNOSTIC\n")
cat("NO_DECONVOLUTESE_OR_NMF_OPTIMIZER_EXECUTED\n")
cat("NO_RECOUNT3_EXPRESSION_LOADED\n")
cat("NO_FROZEN_FILES_MODIFIED\nCLINICAL_ADMISSIONS 0\n")
cat("V7M_FROZEN_ARCHS4_9465_VS_143_NMF_CONTEXT_DIAGNOSTIC_COMPLETE\n")
