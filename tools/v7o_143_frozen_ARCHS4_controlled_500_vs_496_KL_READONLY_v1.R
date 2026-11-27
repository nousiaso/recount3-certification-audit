#!/usr/bin/env Rscript
# V7O: original ARCHS4 143 matched GSMs. Controlled fixed-W KL projection
# with 500 versus 496 original directional features; original package 143-H
# is computed as a reference only. Optimizer here is INDEPENDENT stats::optim,
# not claimed to be the original NMF solver. No recount3 expression is loaded.
# Four variants per sample: W500 and W496 x archived-H start / uniform start.
# These isolate changes to the KL objective feature universe and measure
# sensitivity of a converged solution to initialization.
suppressPackageStartupMessages({
  library(data.table)
  library(NMF)
  library(SpatialEcoTyper)
})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: Rscript v7o.R ROOT",call.=FALSE)
ROOT <- normalizePath(args[[1L]],mustWork=TRUE)
REPO <- file.path(ROOT,"RECOUNT3_CERTIFICATION_GIT")
REL <- file.path(ROOT,"RELEASES","ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0","PROJECT_TREE")
frozen <- function(stage,name) {
  p <- c(file.path(ROOT,stage,name),file.path(REL,stage,name))
  f <- p[file.exists(p)]
  if (!length(f)) stop("FROZEN_ASSET_MISSING ",stage,"/",name,call.=FALSE)
  f[[1L]]
}
assert <- function(p,s) if (!isTRUE(p)) stop(s,call.=FALSE)
assert(as.character(packageVersion("SpatialEcoTyper"))=="1.0.3",
       "Original SpatialEcoTyper 1.0.3 required")
W <- readRDS("/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds")
assert(is.matrix(W)&&identical(dim(W),c(500L,10L)),"Frozen W shape invalid")
assert(all(is.finite(W))&&all(W>=0),"Frozen W negative or nonfinite")
component_names <- colnames(W)
Wgenes <- unique(sub("__.*$","",rownames(W)))
assert(length(Wgenes)==499L,"Frozen W base genes expected 499")

full <- readRDS(frozen("11_pan_cancer_fixed_anchor_projection_v1_0_4",
                       "06_fixed_anchor_scaled_input_CPM.rds"))
assert(is.matrix(full)&&identical(dim(full),c(499L,9465L)) &&
       identical(rownames(full),Wgenes) &&
       !anyDuplicated(colnames(full)),"Frozen 9465 Z invalid")
manifest <- fread(file.path(REPO,"private_data",
             "v7j_143_frozen_archs4_vs_recount3_full499_v1",
             "143_pair_manifest_LOCAL_ONLY.tsv"),
             sep="\t",quote="",showProgress=FALSE)
ids <- as.character(manifest$ARCHS4_GSM)
assert(length(ids)==143L&&!anyDuplicated(ids)&&all(ids %in% colnames(full)),
       "Exact 143 reference GSM identities unavailable")
Z <- full[,ids,drop=FALSE]
rm(full); invisible(gc(verbose=FALSE))

# Stage10/installed .nmf.predict deterministic directional transform.
Vall <- as.matrix(NMF::posneg(Z[Wgenes,,drop=FALSE]))
d <- duplicated(rownames(Vall))
rownames(Vall)[!d] <- paste0(rownames(Vall)[!d],"__pos")
rownames(Vall)[d] <- paste0(rownames(Vall)[d],"__neg")
assert(all(rownames(W) %in% rownames(Vall)),"Frozen W directional features missing")
V500 <- Vall[rownames(W),,drop=FALSE]
assert(identical(dim(V500),c(500L,143L)), "Full directional matrix invalid")
rowvar <- apply(Vall,1L,stats::var)
retained <- intersect(rownames(Vall)[is.finite(rowvar)&rowvar>0],rownames(W))
lost <- setdiff(rownames(W),retained)
assert(length(retained)==496L && setequal(lost,
       c("F12__pos","JUND__pos","PSMA3__neg","SURF4__neg")),
       "Original matched 143 W-feature retention changed")
V496 <- Vall[retained,,drop=FALSE]
W496 <- W[retained,,drop=FALSE]
assert(all(V500>=0)&&all(V496>=0),"Invalid nonnegative V")
rm(Vall)
cat("V7O_FROZEN_143_CONTROLLED_FEATURE_OBJECTIVE_DIAGNOSTIC\n")
cat("ARCHIVED_REFERENCE_SAMPLES 143\n")
cat("FROZEN_W_BASE_GENES 499\n")
cat("FEATURE_UNIVERSE_500 500\nFEATURE_UNIVERSE_496 496\n")
cat("DROPPED_W_FEATURES",paste(lost,collapse=";"),"\n")

h_dt <- fread(frozen("11_pan_cancer_fixed_anchor_projection_v1_0_4",
                     "09_raw_H_CPM.tsv.gz"),sep="\t",quote="",check.names=FALSE,
               showProgress=FALSE)
assert(nrow(h_dt)==9465L && all(c("sample_id",component_names) %in% names(h_dt)),
       "Original H archive invalid")
idx <- match(ids,as.character(h_dt$sample_id))
assert(!anyNA(idx),"Original H GSM index mismatch")
Harchive <- as.matrix(h_dt[idx,..component_names])
storage.mode(Harchive) <- "double"
rownames(Harchive) <- ids
colnames(Harchive) <- component_names
assert(identical(dim(Harchive),c(143L,10L)) &&
       all(is.finite(Harchive))&&all(Harchive>=0)&&all(rowSums(Harchive)>0),
       "Archived original H matrix invalid")

# Original official 143-context H, independently re-executed using its fixed W.
# Expect same mismatch as V7J; no accepted assay calibration results produced.
H143 <- SpatialEcoTyper::DeconvoluteSE(
  dat=Z, scale=FALSE,W=W,nsample.per.run=10000L,sum2one=FALSE,ncores=1L)
H143 <- as.matrix(H143)
assert(all(ids %in% rownames(H143)) &&
       identical(colnames(H143),component_names),"Original 143 H axes changed")
H143 <- H143[ids,component_names,drop=FALSE]
assert(all(is.finite(H143))&&all(H143>=0)&&all(rowSums(H143)>0),
       "Original solver 143 H invalid")
fractions <- function(h) h/rowSums(h)
tvd <- function(a,b) rowSums(abs(fractions(a)-fractions(b)))/2
same_original <- tvd(Harchive,H143)
cat("ORIGINAL_143_VS_9465_H_MAX_ABS",
    format(max(abs(H143-Harchive)),digits=12,scientific=TRUE),"\n")
cat("ORIGINAL_143_VS_9465_TVD_MEDIAN",
    format(median(same_original),digits=12,scientific=TRUE),"\n")
cat("ORIGINAL_143_VS_9465_TVD_MAX",
    format(max(same_original),digits=12,scientific=TRUE),"\n")

# Full W and reduced W share the exact same original 143 sample list and Z.
# For fixed W, each sample's KL objective is convex in its nonnegative h.
# Use stats::optim L-BFGS-B to find a controlled stationary point for each
# sample, with TWO different starts to test solution stability.
# IMPORTANT: this is a mathematical diagnostic; it is NOT SpatialEcoTyper NMF.
kl_fit <- function(Wsub, v, h0) {
  v <- as.numeric(v)
  Ws <- as.matrix(Wsub)
  eps <- 1e-12
  nz <- v>0
  objective <- function(h) {
    p <- pmax(as.numeric(Ws %*% h),eps)
    sum(p-v) + sum(v[nz]*log(v[nz]/p[nz]))
  }
  gradient <- function(h) {
    p <- pmax(as.numeric(Ws %*% h),eps)
    as.numeric(crossprod(Ws,1-v/p))
  }
  sol <- stats::optim(
    par=as.numeric(pmax(h0,1e-7)),
    fn=objective,gr=gradient,method="L-BFGS-B",
    lower=rep(0,ncol(Ws)),
    control=list(maxit=2000L,factr=1e7,pgtol=1e-8)
  )
  h <- pmax(sol$par,0)
  c(h,objective=objective(h),convergence=sol$convergence,
    iterations=if (!is.null(sol$counts[["function"]]))
      as.numeric(sol$counts[["function"]]) else NA_real_)
}
opt500_archive <- matrix(NA_real_,143,10,dimnames=list(ids,component_names))
opt500_uniform <- opt500_archive
opt496_archive <- opt500_archive
opt496_uniform <- opt500_archive
objective_deltas <- rep(NA_real_,143)
solver_codes <- matrix(NA_integer_,143,4)
max_relative_objective_start_diff <- 0
for (i in seq_along(ids)) {
  init_archive <- pmax(as.numeric(Harchive[i,]),1e-7)
  init_flat <- rep(max(mean(init_archive),1e-3),10L)
  f5a <- kl_fit(W,V500[,i],init_archive)
  f5u <- kl_fit(W,V500[,i],init_flat)
  f4a <- kl_fit(W496,V496[,i],init_archive)
  f4u <- kl_fit(W496,V496[,i],init_flat)
  opt500_archive[i,] <- f5a[1:10]
  opt500_uniform[i,] <- f5u[1:10]
  opt496_archive[i,] <- f4a[1:10]
  opt496_uniform[i,] <- f4u[1:10]
  solver_codes[i,] <- as.integer(c(f5a["convergence"],f5u["convergence"],
                                   f4a["convergence"],f4u["convergence"]))
  if (any(!is.finite(c(f5a,f5u,f4a,f4u))))
    stop("NONFINITE_KL_OPTIMIZATION_RESULT",call.=FALSE)
  relative5 <- abs(f5a["objective"]-f5u["objective"])/max(1,abs(f5a["objective"]))
  relative4 <- abs(f4a["objective"]-f4u["objective"])/max(1,abs(f4a["objective"]))
  max_relative_objective_start_diff <- max(
    max_relative_objective_start_diff,relative5,relative4)
  objective_deltas[i] <- f5a["objective"]-f4a["objective"]
  if (i%%25L==0L) cat("KL_CONVEX_SOLVER_SAMPLES_COMPLETED",i,"\n")
}
assert(all(rowSums(opt500_archive)>0)&&all(rowSums(opt496_archive)>0),
       "KL optimized H has empty vectors")
t_start500 <- tvd(opt500_archive,opt500_uniform)
t_start496 <- tvd(opt496_archive,opt496_uniform)
feature_tvd <- tvd(opt500_archive,opt496_archive)
arch_to_opt500 <- tvd(Harchive,opt500_archive)
official143_to_opt496 <- tvd(H143,opt496_archive)
fmt <- function(x) format(x,digits=12,scientific=TRUE)
q95 <- function(x) unname(quantile(x,.95))
cat("OPTIM_METHOD LBFGSB_PER_SPECIMEN_CONVEX_FIXED_W_KL\n")
cat("OPTIM_IS_DIAGNOSTIC_NOT_ORIGINAL_SPATIALECOTYPER\n")
cat("OPTIM_RUNS",length(ids)*4L,"\n")
cat("OPTIM_NONZERO_CONVERGENCE_CODES",sum(solver_codes!=0L),"\n")
cat("OPTIM_STARTING_POINT_OBJECTIVE_MAX_RELATIVE_DIFFERENCE",
    fmt(max_relative_objective_start_diff),"\n")
cat("OPTIM_STARTING_POINT_TVD_MAX_500",fmt(max(t_start500)),"\n")
cat("OPTIM_STARTING_POINT_TVD_MAX_496",fmt(max(t_start496)),"\n")
cat("CONTROLLED_500_VS_496_FRACTION_TVD_MEDIAN",fmt(median(feature_tvd)),
    "P95",fmt(q95(feature_tvd)),
    "MAX",fmt(max(feature_tvd)),"\n")
cat("ARCHIVED_ORIGINAL_VS_OPTIMAL_500_FRACTION_TVD_MEDIAN",
    fmt(median(arch_to_opt500)),
    "P95",fmt(q95(arch_to_opt500)),
    "MAX",fmt(max(arch_to_opt500)),"\n")
cat("OFFICIAL_143_VS_OPTIMAL_496_FRACTION_TVD_MEDIAN",
    fmt(median(official143_to_opt496)),
    "P95",fmt(q95(official143_to_opt496)),
    "MAX",fmt(max(official143_to_opt496)),"\n")
cat("FULL500_VS_REDUCED496_KL_OBJECTIVE_DELTA_MEDIAN",
    fmt(median(objective_deltas)),
    "MAX",fmt(max(objective_deltas)),"\n")
cat("INTERPRETATION_TVD_DOES_NOT_DECOMPOSE_ADDITIVELY\n")
cat("CONVERGENCE_IS_ASSESSED_BY_OPTIM_CODES_AND_TWO_START_AGREEMENT\n")
if (any(solver_codes!=0L) ||
    max_relative_objective_start_diff>1e-5 ||
    max(t_start500,t_start496)>1e-3) {
  cat("V7O_CONVEX_KL_OPTIMUM_STABILITY_HOLD\n")
} else {
  cat("V7O_CONVEX_KL_OPTIMUM_STABILITY_PASS\n")
}
cat("NO_RECOUNT3_EXPRESSION_OR_MODEL_PROJECTION\n")
cat("NO_FROZEN_FILES_MODIFIED\nCLINICAL_ADMISSIONS 0\n")
cat("V7O_FROZEN_143_FIXED_W_KL_FEATURE_AND_OPTIMIZER_DIAGNOSTIC_COMPLETE\n")
