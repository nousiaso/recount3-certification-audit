#!/usr/bin/env Rscript
# V7G same-specimen dual-assay full-499 frozen SpatialEcoTyper transfer.
# Two 143-specimen input matrices already normalized by original frozen Stage09.
# Model W and pan-cancer anchor mean/SD are immutable; no fitting/recalibration.
suppressPackageStartupMessages({
  library(data.table)
  library(SpatialEcoTyper)
  library(NMF)
})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: Rscript V7G.R ROOT PRIVATE_OUTPUT", call.=FALSE)
ROOT <- normalizePath(args[[1L]],mustWork=TRUE)
OUT <- normalizePath(args[[2L]],mustWork=TRUE)
REL <- file.path(ROOT,"RELEASES","ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0","PROJECT_TREE")
BRANCHES <- c("CPM","TMM_CPM","UQ_CPM")
ASSAYS <- c("ARCHS4","recount3")
COMPONENTS <- c("NonSE",paste0("SE",1:9))
assert <- function(p,s) if (!isTRUE(p)) stop(s,call.=FALSE)
ref_file <- function(name,required=TRUE) {
  p <- c(file.path(ROOT,"11_pan_cancer_fixed_anchor_projection_v1_0_4",name),
         file.path(REL,"11_pan_cancer_fixed_anchor_projection_v1_0_4",name))
  f <- p[file.exists(p)]
  if (!length(f)) {
    if (required) stop("Frozen Stage10 reference missing: ",name,call.=FALSE)
    return(NA_character_)
  }
  f[[1L]]
}

W_FILE <- "/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds"
assert(file.exists(W_FILE),"Original frozen W RDS missing")
W <- readRDS(W_FILE)
assert(is.matrix(W) && identical(dim(W),c(500L,10L)),
       "Frozen W must be 500 directional features by 10 components")
assert(identical(colnames(W),COMPONENTS),"Frozen W component order differs")
installed <- file.path(system.file("extdata",package="SpatialEcoTyper"),
                       "Bulk_SE_Recovery_W.rds")
assert(file.exists(installed) && identical(W,readRDS(installed)),
       "Frozen W is not identical to installed SpatialEcoTyper W")
assert(as.character(packageVersion("SpatialEcoTyper"))=="1.0.3",
       "Expected original SpatialEcoTyper 1.0.3")

manifest <- fread(file.path(OUT,"143_pair_manifest_LOCAL_ONLY.tsv"),
                  sep="\t",check.names=FALSE)
assert(nrow(manifest)==143L && uniqueN(manifest$ARCHS4_GSM)==143L,
       "143 exact unique bridge GSMs required")
gsm <- as.character(manifest$ARCHS4_GSM)
uq_arch <- toupper(as.character(manifest$UQ_evaluable_ARCHS4))=="TRUE"
uq_rec <- toupper(as.character(manifest$UQ_evaluable_recount3))=="TRUE"
assert(!anyNA(uq_arch) && !anyNA(uq_rec),"UQ evaluability incomplete")
W_genes <- unique(sub("__.*$","",rownames(W)))
k <- function(x) toupper(trimws(as.character(x)))
assert(length(W_genes)==499L && !anyDuplicated(k(W_genes)),
       "Frozen W does not contain 499 distinct symbols")
nine <- c("FAM20C","FCAR","ISLR","MAL2","MATR3","PODNL1",
          "PSME2","SLC39A4","TNNT3")

read_expr <- function(assay,branch) {
  path <- file.path(OUT,paste0("499_",assay,"_",branch,"_LOCAL_ONLY.tsv.gz"))
  x <- fread(path,sep="\t",check.names=FALSE)
  assert("model_gene" %in% names(x) && nrow(x)==499L,
         paste("Incorrect expression matrix",assay,branch))
  assert(!anyDuplicated(k(x$model_gene)),
         paste("Duplicate gene identifiers",assay,branch))
  assert(identical(names(x)[-1L],gsm),
         paste("Matched sample order differs",assay,branch))
  ix <- match(k(W_genes),k(x$model_gene))
  assert(!anyNA(ix),paste("Missing W genes",assay,branch))
  v <- as.matrix(x[ix, ..gsm])
  storage.mode(v) <- "double"
  rownames(v) <- W_genes
  colnames(v) <- gsm
  assert(identical(dim(v),c(499L,143L)) &&
         all(is.finite(v)) && all(v>=0),
         paste("Invalid expression values",assay,branch))
  v
}
read_params <- function(branch) {
  name <- switch(branch,
       CPM="03_pan_cancer_anchor_scaling_CPM.tsv.gz",
       TMM_CPM="04_pan_cancer_anchor_scaling_TMM_CPM.tsv.gz",
       UQ_CPM="05_pan_cancer_anchor_scaling_UQ_CPM.tsv.gz")
  x <- fread(ref_file(name),sep="\t")
  cols <- c("model_gene","weighted_anchor_mean_log2","weighted_anchor_sd_log2")
  assert(all(cols %in% names(x)) && nrow(x)==499L,
         paste("Frozen anchor scaling contract failed",branch))
  assert(!anyDuplicated(k(x$model_gene)),"Duplicate frozen scaling genes")
  i <- match(k(W_genes),k(x$model_gene))
  assert(!anyNA(i),"Frozen scaling missing W genes")
  mu <- as.numeric(x$weighted_anchor_mean_log2[i])
  sd <- as.numeric(x$weighted_anchor_sd_log2[i])
  assert(all(is.finite(mu)) && all(is.finite(sd)) && all(sd>0),
         "Invalid frozen means or SDs")
  list(mu=mu,sd=sd)
}
to_directional <- function(Z) {
  all <- as.matrix(NMF::posneg(Z))
  duplicated_rows <- duplicated(rownames(all))
  rownames(all)[!duplicated_rows] <- paste0(rownames(all)[!duplicated_rows],"__pos")
  rownames(all)[duplicated_rows] <- paste0(rownames(all)[duplicated_rows],"__neg")
  i <- match(rownames(W),rownames(all))
  assert(!anyNA(i),"Directional NMF::posneg failed to recover all 500 frozen W rows")
  V <- all[i,,drop=FALSE]
  assert(identical(rownames(V),rownames(W)) &&
         all(is.finite(V)) && all(V>=0),
         "Invalid W-aligned nonnegative directional matrix")
  V
}
project <- function(Z) {
  h <- SpatialEcoTyper::DeconvoluteSE(
    dat=Z,scale=FALSE,W=W,nsample.per.run=10000L,
    sum2one=FALSE,ncores=1L
  )
  h <- as.matrix(h)
  storage.mode(h) <- "double"
  assert(setequal(rownames(h),gsm) &&
         identical(colnames(h),COMPONENTS),
         "DeconvoluteSE returned unexpected sample/component membership")
  h <- h[gsm,COMPONENTS,drop=FALSE]
  assert(all(is.finite(h)) && all(h>=0) && all(rowSums(h)>0),
         "Invalid frozen model raw H")
  h
}
metrics <- function(Z,H) {
  V <- to_directional(Z)
  Vhat <- W %*% t(H)
  r <- V - Vhat
  centered <- sweep(V,2L,colMeans(V),"-")
  sst <- colSums(centered^2)
  list(
    reconstruction_rmse=sqrt(colMeans(r^2)),
    reconstruction_mae=colMeans(abs(r)),
    reconstruction_r2=1-colSums(r^2)/pmax(sst,1e-12),
    directional_input_mass=colSums(V),
    raw_H_total=rowSums(H)
  )
}
describe <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(c(median=NA_real_,p90=NA_real_,p95=NA_real_,max=NA_real_))
  c(median=median(x),p90=unname(quantile(x,.90)),
    p95=unname(quantile(x,.95)),max=max(x))
}
boot_median <- function(x,B=500L) {
  if (length(x)<2L) return(c(NA_real_,NA_real_))
  set.seed(20261127L)
  d <- replicate(B,median(sample(x,length(x),replace=TRUE)))
  unname(quantile(d,c(.025,.975)))
}

by_sample <- list()
by_gene <- list()
summary <- list()
projections <- list()

for (branch in BRANCHES) {
  message("V7G_FROZEN_DUAL_ASSAY_BRANCH_BEGIN ",branch)
  par <- read_params(branch)
  obj <- list()
  for (assay in ASSAYS) {
    X <- read_expr(assay,branch)
    Z <- sweep(sweep(log2(X+1),1L,par$mu,"-"),1L,par$sd,"/")
    rownames(Z) <- W_genes
    colnames(Z) <- gsm
    assert(all(is.finite(Z)),paste("Bad frozen-Z",assay,branch))
    H <- project(Z)
    obj[[assay]] <- list(Z=Z,H=H,f=H/rowSums(H),qc=metrics(Z,H))
    message("V7G_PROJECTED ",assay," ",branch," SAMPLES ",nrow(H))
  }
  arch <- obj[["ARCHS4"]]
  rec <- obj[["recount3"]]
  use <- if (branch=="UQ_CPM") uq_arch & uq_rec else rep(TRUE,143L)
  assert(sum(use)>=2L,"Insufficient jointly evaluable matched UQ specimens")
  f_diff <- rec$f-arch$f
  tvd <- rowSums(abs(f_diff))/2
  dom_arch <- COMPONENTS[max.col(arch$f,ties.method="first")]
  dom_rec <- COMPONENTS[max.col(rec$f,ties.method="first")]
  same <- dom_arch==dom_rec
  H_diff <- rec$H-arch$H
  mean_abs_H_log2 <- rowMeans(abs(log2(rec$H+1e-8)-log2(arch$H+1e-8)))
  gene_Z_diff <- rec$Z-arch$Z
  Z_mae <- colMeans(abs(gene_Z_diff))
  primary_sample <- data.table(
      branch=branch, sample_id_LOCAL_ONLY=gsm,
      n_technical_runs=as.integer(manifest$n_technical_runs),
      UQ_joint_evaluable=uq_arch & uq_rec,
      fraction_TVD=tvd,
      same_dominant_component=same,
      ARCHS4_dominant=dom_arch,recount3_dominant=dom_rec,
      ARCHS4_H_total=arch$qc$raw_H_total,
      recount3_H_total=rec$qc$raw_H_total,
      H_mean_abs_log2_difference=mean_abs_H_log2,
      frozen_Z_MAE=Z_mae,
      ARCHS4_reconstruction_RMSE=arch$qc$reconstruction_rmse,
      recount3_reconstruction_RMSE=rec$qc$reconstruction_rmse,
      ARCHS4_reconstruction_R2=arch$qc$reconstruction_r2,
      recount3_reconstruction_R2=rec$qc$reconstruction_r2
  )
  for (comp in COMPONENTS) {
    primary_sample[[paste0("fraction_difference_",comp)]] <- f_diff[,comp]
    primary_sample[[paste0("raw_H_difference_",comp)]] <- H_diff[,comp]
  }
  by_sample[[branch]] <- primary_sample
  by_gene[[branch]] <- data.table(
      branch=branch,model_gene=W_genes,
      median_absolute_frozen_Z_difference=apply(abs(gene_Z_diff),1L,median),
      mean_frozen_Z_bias=rowMeans(gene_Z_diff),
      nine_different_source_Ensembl_membership=W_genes %in% nine
  )
  ci <- boot_median(tvd[use])
  s <- describe(tvd[use])
  summary[[branch]] <- data.table(
      branch=branch,n_pairs_total=143L,n_jointly_evaluable=sum(use),
      n_UQ_fallback=143L-sum(use),
      median_fraction_TVD=s[["median"]],
      p90_fraction_TVD=s[["p90"]],
      p95_fraction_TVD=s[["p95"]],
      max_fraction_TVD=s[["max"]],
      median_fraction_TVD_bootstrap_low_descriptive=ci[[1L]],
      median_fraction_TVD_bootstrap_high_descriptive=ci[[2L]],
      dominant_component_concordance=mean(same[use]),
      median_H_mean_abs_log2_difference=median(mean_abs_H_log2[use]),
      median_frozen_Z_mean_absolute_difference=median(Z_mae[use]),
      median_ARCHS4_reconstruction_RMSE=median(arch$qc$reconstruction_rmse[use]),
      median_recount3_reconstruction_RMSE=median(rec$qc$reconstruction_rmse[use]),
      n_distinct_ARCHS4_dominant_components=length(unique(dom_arch[use])),
      n_distinct_recount3_dominant_components=length(unique(dom_rec[use]))
  )
  message("V7G_BRANCH ",branch," N=",sum(use),
          " MEDIAN_FRACTION_TVD=",format(s[["median"]],digits=6),
          " P95_TVD=",format(s[["p95"]],digits=6),
          " DOMINANT_CONCORDANCE=",format(mean(same[use]),digits=6))
  projections[[branch]] <- obj
}

# Within-assay normalization robustness, not conflated with cross-assay transfer.
norm_summary <- list()
for (assay in ASSAYS) {
  primary <- projections[["CPM"]][[assay]]$f
  primary_dom <- COMPONENTS[max.col(primary,ties.method="first")]
  for (branch in c("TMM_CPM","UQ_CPM")) {
    alternative <- projections[[branch]][[assay]]$f
    use <- if (branch=="UQ_CPM") {
      if (assay=="ARCHS4") uq_arch else uq_rec
    } else {
      rep(TRUE,143L)
    }
    if (!any(use)) next
    distance <- rowSums(abs(primary-alternative))/2
    other_dom <- COMPONENTS[max.col(alternative,ties.method="first")]
    stats <- describe(distance[use])
    norm_summary[[paste(assay,branch,sep="_")]] <- data.table(
      assay=assay,comparison=paste0("CPM_vs_",branch),
      n_evaluable=sum(use),
      median_fraction_TVD=stats[["median"]],
      p95_fraction_TVD=stats[["p95"]],
      dominant_component_agreement=mean(primary_dom[use]==other_dom[use])
    )
  }
}
norm_aggregate <- rbindlist(norm_summary)
fwrite(norm_aggregate,
       file.path(OUT,"aggregate_normalization_sensitivity.tsv"),sep="\t")
message("V7G_WITHIN_ASSAY_NORMALIZATION_ROBUSTNESS_COMPLETE")

# Prespecified nine-annotation-difference sensitivity:
# neutralize exactly these nine scaled inputs (frozen Z = 0) in BOTH assays,
# without modifying W, anchors, or the primary full-499 results.
assert(sum(W_genes %in% nine)==9L,
       "The nine known source annotation difference genes are not all in W")
neutral <- list()
neutral_per_pair <- data.table(sample_id_LOCAL_ONLY=gsm)
for (assay in ASSAYS) {
  ref <- projections[["CPM"]][[assay]]
  z <- ref$Z
  z[W_genes %in% nine,] <- 0
  h <- project(z)
  neutral[[assay]] <- h/rowSums(h)
  neutral_per_pair[[paste0("within_",assay,"_fraction_TVD_full_vs_nine_neutralized")]] <-
    rowSums(abs(neutral[[assay]]-ref$f))/2
}
base_tvd <- rowSums(abs(
  projections[["CPM"]][["ARCHS4"]]$f -
    projections[["CPM"]][["recount3"]]$f))/2
neutral_tvd <- rowSums(abs(neutral[["ARCHS4"]]-neutral[["recount3"]]))/2
neutral_per_pair[["cross_assay_CPM_TVD_full_499"]] <- base_tvd
neutral_per_pair[["cross_assay_CPM_TVD_nine_neutralized"]] <- neutral_tvd
neutral_per_pair[["cross_assay_CPM_TVD_neutral_minus_full"]] <- neutral_tvd-base_tvd
fwrite(neutral_per_pair,
       file.path(OUT,"143_nine_gene_neutralization_LOCAL_ONLY.tsv"),sep="\t")
nine_agg <- data.table(
  sensitivity="set_frozen_Z_to_zero_for_nine_known_annotation_difference_genes",
  frozen_model_gene_count=499L,neutralized_genes=9L,matched_pairs=143L,
  median_CPM_TVD_full=median(base_tvd),
  median_CPM_TVD_neutralized=median(neutral_tvd),
  median_per_pair_delta_TVD=median(neutral_tvd-base_tvd),
  p95_absolute_per_pair_delta_TVD=unname(quantile(abs(neutral_tvd-base_tvd),.95))
)
fwrite(nine_agg,
       file.path(OUT,"aggregate_nine_gene_sensitivity.tsv"),sep="\t")
message("V7G_NINE_ANNOTATION_GENE_NEUTRALIZATION_SENSITIVITY_COMPLETE")

# Optional historical CPM H replay only for GSMs actually present in the
# frozen 9,465-sample atlas. Absent specimens are NOT silently substituted.
path_h <- ref_file("09_raw_H_CPM.tsv.gz",required=FALSE)
reference_replay_count <- 0L
if (!is.na(path_h)) {
  original <- fread(path_h,sep="\t",check.names=FALSE)
  assert("sample_id" %in% names(original) && all(COMPONENTS %in% names(original)),
         "Historical frozen H contract invalid")
  common <- gsm[gsm %in% original$sample_id]
  reference_replay_count <- length(common)
  if (length(common)>0L) {
    ix <- match(common,original$sample_id)
    orig_h <- as.matrix(original[ix, ..COMPONENTS])
    storage.mode(orig_h) <- "double"
    rownames(orig_h) <- common
    new_h <- projections[["CPM"]][["ARCHS4"]][["H"]][common,,drop=FALSE]
    replay <- data.table(
       sample_id_LOCAL_ONLY=common,
       fraction_TVD_to_historical_CPM_H=
         rowSums(abs(new_h/rowSums(new_h)-orig_h/rowSums(orig_h)))/2,
       raw_H_mean_absolute_difference=rowMeans(abs(new_h-orig_h))
    )
    fwrite(replay,
           file.path(OUT,"historical_ARCHS4_CPM_replay_LOCAL_ONLY.tsv"),
           sep="\t")
    message("V7G_HISTORICAL_CPM_REFERENCE_OVERLAP ",length(common),
            " MEDIAN_TVD ",format(median(replay$fraction_TVD_to_historical_CPM_H),
                                   digits=7))
  } else {
    message("V7G_HISTORICAL_CPM_REFERENCE_OVERLAP 0; dual assay projection self-contained")
  }
}
per_sample <- rbindlist(by_sample)
per_gene <- rbindlist(by_gene)
aggregate <- rbindlist(summary)
fwrite(per_sample,file.path(OUT,"143_per_pair_calibration_LOCAL_ONLY.tsv"),sep="\t")
fwrite(per_gene,file.path(OUT,"499_gene_calibration_LOCAL_ONLY.tsv"),sep="\t")
fwrite(aggregate,file.path(OUT,"aggregate_transfer_metrics.tsv"),sep="\t")
saveRDS(projections,file.path(OUT,"all_143_dual_assay_projections_LOCAL_ONLY.rds"),compress="gzip")
message("V7G_HISTORICAL_ARCHS4_CPM_H_OVERLAP ",reference_replay_count)
message("V7G_FULL_499_DUAL_ASSAY_COMPARISON_COMPLETE")
message("CALIBRATION_SPECIMENS 143")
message("SOURCE_TECHNICAL_RUNS 454")
message("CALIBRATION_ONLY_NOT_INDEPENDENT_BIOLOGICAL_VALIDATION")
message("DESCRIPTIVE_BOOTSTRAP_NOT_PATIENT_INDEPENDENCE_INFERENCE")
message("NO_W_OR_ANCHOR_REFIT_NO_FROZEN_FILES_MODIFIED")
message("CLINICAL_ADMISSIONS 0")
