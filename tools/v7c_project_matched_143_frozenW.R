#!/usr/bin/env Rscript
# V7C frozen-W matched assay-transfer calibration: real computation, NO refitting.
# Uses 143 exact G026/ARCHS4 pairs, frozen Stage10 scaling tables/H as reference.
# Writes all specimen-linked output only in ignored private_data.
suppressPackageStartupMessages({
  library(data.table)
  library(SpatialEcoTyper)
  library(NMF)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: Rscript v7c_project.R ROOT OUT", call.=FALSE)
ROOT <- normalizePath(args[[1L]], mustWork=TRUE)
OUT <- normalizePath(args[[2L]], mustWork=TRUE)
REL <- file.path(ROOT, "RELEASES", "ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0", "PROJECT_TREE")
BRANCHES <- c("CPM", "TMM_CPM", "UQ_CPM")
COMPONENTS <- c("NonSE", paste0("SE", 1:9))
N <- 143L

reference_file <- function(filename) {
  paths <- c(
    file.path(ROOT, "11_pan_cancer_fixed_anchor_projection_v1_0_4", filename),
    file.path(REL, "11_pan_cancer_fixed_anchor_projection_v1_0_4", filename)
  )
  found <- paths[file.exists(paths)]
  if (length(found) == 0L) stop("Missing frozen Stage10: ", filename, call.=FALSE)
  found[[1L]]
}
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call.=FALSE)

W_FILE <- "/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds"
assert(file.exists(W_FILE), "Frozen W missing")
W <- readRDS(W_FILE)
assert(is.matrix(W) && identical(dim(W), c(500L, 10L)), "W is not 500x10")
assert(identical(colnames(W), COMPONENTS), "W component order changed")
installed <- file.path(system.file("extdata", package="SpatialEcoTyper"), "Bulk_SE_Recovery_W.rds")
assert(file.exists(installed), "Installed W missing")
assert(identical(W, readRDS(installed)), "Frozen W differs from installed package W")
assert(as.character(packageVersion("SpatialEcoTyper")) == "1.0.3", "Expected SpatialEcoTyper 1.0.3")

pair <- fread(file.path(OUT, "143_pair_manifest_LOCAL_ONLY.tsv"),
              sep="\t", check.names=FALSE)
assert(nrow(pair) == N && uniqueN(pair$ARCHS4_GSM) == N, "Pair manifest count/duplicates")
gsm <- as.character(pair$ARCHS4_GSM)
uq_ok <- toupper(as.character(pair$UQ_evaluable)) == "TRUE"
assert(!anyNA(uq_ok), "UQ evaluability has NA")

W_gene <- unique(sub("__.*$", "", rownames(W)))
assert(length(W_gene) == 499L, "W does not contain exactly 499 base genes")
k <- function(x) toupper(trimws(as.character(x)))
assert(!anyDuplicated(k(W_gene)), "Ambiguous case-folded W gene")
gene_differences <- c("FAM20C", "FCAR", "ISLR", "MAL2",
                      "MATR3", "PODNL1", "PSME2", "SLC39A4", "TNNT3")

read_frozen_H <- function(branch) {
  fname <- switch(branch,
       CPM="09_raw_H_CPM.tsv.gz",
       TMM_CPM="11_raw_H_TMM_CPM.tsv.gz",
       UQ_CPM="13_raw_H_UQ_CPM.tsv.gz")
  x <- fread(reference_file(fname), sep="\t", check.names=FALSE)
  assert("sample_id" %in% names(x), "Frozen H lacks sample_id")
  assert(all(COMPONENTS %in% names(x)), "Frozen H lacks components")
  assert(!anyDuplicated(x$sample_id), "Frozen H duplicated sample_id")
  ix <- match(gsm, x$sample_id)
  assert(!anyNA(ix), paste("Not all 143 pairs in frozen", branch, "H"))
  m <- as.matrix(x[ix, ..COMPONENTS])
  storage.mode(m) <- "double"
  rownames(m) <- gsm; colnames(m) <- COMPONENTS
  assert(all(is.finite(m)) && all(m >= 0), "Invalid frozen H")
  m
}

read_scaling <- function(branch) {
  fname <- switch(branch,
      CPM="03_pan_cancer_anchor_scaling_CPM.tsv.gz",
      TMM_CPM="04_pan_cancer_anchor_scaling_TMM_CPM.tsv.gz",
      UQ_CPM="05_pan_cancer_anchor_scaling_UQ_CPM.tsv.gz")
  d <- fread(reference_file(fname), sep="\t")
  required <- c("model_gene", "weighted_anchor_mean_log2", "weighted_anchor_sd_log2")
  assert(all(required %in% names(d)), paste("Incomplete frozen scaling", branch))
  assert(nrow(d) == 499L, "Scaling table is not 499 genes")
  assert(!anyDuplicated(k(d$model_gene)), "Scaling table has duplicate genes")
  ix <- match(k(W_gene), k(d$model_gene))
  assert(!anyNA(ix), paste("W genes not in scaling", branch))
  a <- as.numeric(d$weighted_anchor_mean_log2[ix])
  s <- as.numeric(d$weighted_anchor_sd_log2[ix])
  assert(all(is.finite(a)) && all(is.finite(s) & s > 0),
         "Frozen scaling values not valid")
  list(mean=a, sd=s)
}

read_recount_matrix <- function(branch) {
  fname <- paste0("499_", branch, "_LOCAL_ONLY.tsv.gz")
  x <- fread(file.path(OUT, fname), sep="\t", check.names=FALSE)
  assert("model_gene" %in% names(x), "Recount matrix lacks model_gene")
  assert(nrow(x) == 499L, "Recount model matrix is not 499 genes")
  genes <- as.character(x$model_gene)
  assert(!anyDuplicated(k(genes)), "Ambiguous recount gene keys")
  assert(setequal(names(x)[-1], gsm), "Recount sample column membership disagrees with bridge")
  ix <- match(k(W_gene), k(genes))
  assert(!anyNA(ix), "Recount model gene mismatch")
  d <- as.matrix(x[ix, ..gsm])
  storage.mode(d) <- "double"
  rownames(d) <- W_gene; colnames(d) <- gsm
  assert(identical(dim(d), c(499L, N)), "Recount input dimensions invalid")
  assert(all(is.finite(d)) && all(d >= 0), "Invalid recount model gene expression")
  d
}

reference_z <- function(branch) {
  fname <- switch(branch,
      CPM="06_fixed_anchor_scaled_input_CPM.rds",
      TMM_CPM="07_fixed_anchor_scaled_input_TMM_CPM.rds",
      UQ_CPM="08_fixed_anchor_scaled_input_UQ_CPM.rds")
  x <- readRDS(reference_file(fname))
  assert(is.matrix(x) && nrow(x) == 499L && all(gsm %in% colnames(x)),
         "Frozen scaled RDS does not contain all 143 matched GSMs")
  ix <- match(k(W_gene), k(rownames(x)))
  assert(!anyNA(ix), "Frozen scaled matrix genes do not match W")
  x <- x[ix, gsm, drop=FALSE]
  rownames(x) <- W_gene
  x
}

directional <- function(Z) {
  V <- as.matrix(NMF::posneg(Z))
  dup <- duplicated(rownames(V))
  rownames(V)[!dup] <- paste0(rownames(V)[!dup], "__pos")
  rownames(V)[dup] <- paste0(rownames(V)[dup], "__neg")
  ix <- match(rownames(W), rownames(V))
  assert(!anyNA(ix), "NMF::posneg did not produce all W directional features")
  V <- V[ix, , drop=FALSE]
  assert(identical(rownames(V), rownames(W)), "Directional row order mismatch")
  assert(all(is.finite(V)) && all(V >= 0), "Invalid directional input")
  V
}

project <- function(Z) {
  raw <- SpatialEcoTyper::DeconvoluteSE(
    dat=Z, scale=FALSE, W=W,
    nsample.per.run=10000L, sum2one=FALSE, ncores=1L
  )
  raw <- as.matrix(raw)
  storage.mode(raw) <- "double"
  assert(setequal(rownames(raw), gsm), "Deconvolution sample IDs differ")
  assert(identical(colnames(raw), COMPONENTS), "Deconvolution component order differs")
  raw <- raw[gsm, COMPONENTS, drop=FALSE]
  assert(all(is.finite(raw)) && all(raw >= 0) && all(rowSums(raw) > 0),
         "Invalid projected raw H")
  raw
}

fractions <- function(H) {
  f <- H / rowSums(H)
  assert(all(abs(rowSums(f)-1) < 1e-8), "Fractions do not sum to one")
  f
}
qstats <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(rep(NA_real_, 5L))
  c(mean=mean(x), median=median(x),
    p90=unname(quantile(x, .9)), p95=unname(quantile(x, .95)), max=max(x))
}
bootstrap_ci <- function(x, B=1000L) {
  set.seed(20261127L)
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  v <- replicate(B, median(sample(x, length(x), replace=TRUE)))
  unname(quantile(v, c(.025,.975)))
}

all_sample <- list()
all_gene <- list()
aggregate <- list()

for (branch in BRANCHES) {
  message("V7C_TRANSFER_BRANCH_BEGIN ", branch)
  expr <- read_recount_matrix(branch)
  scale <- read_scaling(branch)
  Z <- sweep(sweep(log2(expr+1), 1L, scale$mean, "-"), 1L, scale$sd, "/")
  assert(all(is.finite(Z)), paste("Invalid frozen-Z branch", branch))
  rownames(Z) <- W_gene; colnames(Z) <- gsm
  refZ <- reference_z(branch)

  H <- project(Z)
  Href <- read_frozen_H(branch)
  f <- fractions(H)
  fref <- fractions(Href)
  V <- directional(Z)
  Vhat <- W %*% t(H)
  residual <- V-Vhat

  tvd <- rowSums(abs(f-fref))/2
  same <- apply(f,1L,which.max)==apply(fref,1L,which.max)
  raw_log2_difference <- rowMeans(abs(log2(H+1e-8)-log2(Href+1e-8)))
  model_z_mae <- colMeans(abs(Z-refZ))
  reconstruction_rmse <- sqrt(colMeans(residual^2))
  reconstruction_r2 <- 1 - colSums(residual^2) / pmax(colSums(sweep(V, 2L, colMeans(V), "-")^2), 1e-12)

  score <- data.table(
    branch=branch,
    sample_id_LOCAL_ONLY=gsm,
    UQ_formally_evaluable=uq_ok,
    H_reference_total=rowSums(Href),
    H_recount_total=rowSums(H),
    H_mean_absolute_log2_difference=raw_log2_difference,
    fraction_TVD=tvd,
    same_dominant_component=same,
    Z_mean_absolute_difference=model_z_mae,
    recount_reconstruction_RMSE=reconstruction_rmse,
    recount_reconstruction_R2=reconstruction_r2,
    reference_dominant=COMPONENTS[max.col(fref)],
    recount_dominant=COMPONENTS[max.col(f)]
  )
  all_sample[[branch]] <- score

  zdiff <- Z-refZ
  gene <- data.table(branch=branch, model_gene=W_gene,
     median_absolute_frozen_Z_difference=apply(abs(zdiff), 1L, median),
     mean_frozen_Z_bias=rowMeans(zdiff),
     nine_annotation_membership_difference=W_gene %in% gene_differences)
  all_gene[[branch]] <- gene

  # For UQ, fallback samples are reported but not counted as formal UQ-calibration.
  use <- if (branch=="UQ_CPM") uq_ok else rep(TRUE, N)
  assert(sum(use) >= 2L, paste("Insufficient comparable samples", branch))
  tv <- tvd[use]
  ci <- bootstrap_ci(tv)
  aggregate[[branch]] <- data.table(
    branch=branch, n_matched=N, n_formally_compared=sum(use),
    n_fallback=length(use)-sum(use),
    median_fraction_TVD=median(tv),
    fraction_TVD_p90=unname(quantile(tv,.9)),
    fraction_TVD_p95=unname(quantile(tv,.95)),
    fraction_TVD_max=max(tv),
    median_fraction_TVD_bootstrap95_low=ci[[1L]],
    median_fraction_TVD_bootstrap95_high=ci[[2L]],
    dominant_concordance=mean(same[use]),
    median_H_mean_abs_log2_difference=median(raw_log2_difference[use]),
    median_frozen_Z_MAE=median(model_z_mae[use]),
    median_reconstruction_RMSE=median(reconstruction_rmse[use]),
    median_reconstruction_R2=median(reconstruction_r2[use])
  )
  message("V7C_BRANCH ", branch,
          " N=", sum(use),
          " MEDIAN_TVD=", format(median(tv),digits=7),
          " P95_TVD=", format(unname(quantile(tv,.95)),digits=7),
          " DOMINANT_CONCORDANCE=", format(mean(same[use]),digits=7))
  saveRDS(list(H=H, Href=Href, fractions=f, fractions_ref=fref, Z=Z,
               sample=gsm, UQ_evaluable=uq_ok),
          file.path(OUT,paste0("143_",branch,"_matched_projection_LOCAL_ONLY.rds")),
          compress="gzip")
}
fwrite(rbindlist(all_sample), file.path(OUT, "143_per_pair_calibration_LOCAL_ONLY.tsv"),
       sep="\t")
fwrite(rbindlist(all_gene), file.path(OUT, "499_gene_calibration_LOCAL_ONLY.tsv"),
       sep="\t")
agg <- rbindlist(aggregate)
fwrite(agg, file.path(OUT, "aggregate_transfer_metrics.tsv"), sep="\t")
message("V7C_FULL499_FROZEN_TRANSFER_COMPUTATION_COMPLETE")
message("MATCHED_PAIRS ", N)
message("BRANCHES ", paste(BRANCHES, collapse=","))
message("NINE_DIFFERENTIAL_ANNOTATION_GENES_RETAINED ", sum(W_gene %in% gene_differences))
message("NO_W_REFIT_NO_ANCHOR_REFIT_NO_FROZEN_OUTPUT_MODIFICATION")
message("CLINICAL_ADMISSIONS 0")
message("EXTERNAL_BIOLOGICAL_VALIDATION_NOT_PERFORMED")
