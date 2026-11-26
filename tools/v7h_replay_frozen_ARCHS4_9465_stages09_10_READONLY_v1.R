#!/usr/bin/env Rscript
# V7H: full 9,465-sample ARCHS4 Stage09/10 numerical replay, READ ONLY.
# Does not execute the original scripts' output-writing main procedures.
# Re-derives all three Stage09 499-gene normalized matrices from frozen raw
# pseudocount outputs+frozen size factors, verifies frozen Stage10 scaled Z,
# then calls the ORIGINAL Stage10 SpatialEcoTyper deconvolution parameters on
# all 9,465 samples and compares against the archived raw H/fractions.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(SpatialEcoTyper)
  library(NMF)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: Rscript v7h ... ROOT",call.=FALSE)
ROOT <- normalizePath(args[[1L]],mustWork=TRUE)
REL <- file.path(ROOT,"RELEASES","ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0","PROJECT_TREE")
find_frozen <- function(stage,filename) {
  paths <- c(file.path(ROOT,stage,filename),file.path(REL,stage,filename))
  existing <- paths[file.exists(paths)]
  if (!length(existing)) stop("FROZEN_ASSET_MISSING ",stage,"/",filename,call.=FALSE)
  existing[[1L]]
}
assert <- function(x,s) if (!isTRUE(x)) stop(s,call.=FALSE)
components <- c("NonSE",paste0("SE",1:9))
B <- list(
  CPM=c(expression="10_model_gene_CPM_canonical.tsv.gz",
        scaling="03_pan_cancer_anchor_scaling_CPM.tsv.gz",
        z="06_fixed_anchor_scaled_input_CPM.rds",
        h="09_raw_H_CPM.tsv.gz",fractions="10_fraction_scores_CPM.tsv.gz"),
  TMM_CPM=c(expression="11_model_gene_TMM_CPM_canonical.tsv.gz",
            scaling="04_pan_cancer_anchor_scaling_TMM_CPM.tsv.gz",
            z="07_fixed_anchor_scaled_input_TMM_CPM.rds",
            h="11_raw_H_TMM_CPM.tsv.gz",fractions="12_fraction_scores_TMM_CPM.tsv.gz"),
  UQ_CPM=c(expression="12_model_gene_UQ_CPM_canonical.tsv.gz",
           scaling="05_pan_cancer_anchor_scaling_UQ_CPM.tsv.gz",
           z="08_fixed_anchor_scaled_input_UQ_CPM.rds",
           h="13_raw_H_UQ_CPM.tsv.gz",fractions="14_fraction_scores_UQ_CPM.tsv.gz")
)
S9="10_preprojection_normalization_v1_0_1"
S10="11_pan_cancer_fixed_anchor_projection_v1_0_4"
status9 <- fromJSON(find_frozen(S9,"18_preprojection_normalization_status.json"))
status10 <- fromJSON(find_frozen(S10,"28_pan_cancer_projection_status.json"))
assert(identical(status9$status,"PASS") && as.numeric(status9$qc_failures)==0,
       "Archived Stage09 normalization status not PASS")
assert(identical(status10$status,"PASS") && as.numeric(status10$qc_failures)==0,
       "Archived Stage10 model projection status not PASS")
assert(as.numeric(status10$canonical_samples)==9465,
       "Archived Stage10 canonical sample count mismatch")
assert(as.character(packageVersion("SpatialEcoTyper"))=="1.0.3",
       "Installed SpatialEcoTyper version differs from original 1.0.3")
Wpath="/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds"
assert(file.exists(Wpath),"Frozen W path missing")
W <- readRDS(Wpath)
installed <- file.path(system.file("extdata",package="SpatialEcoTyper"),
                       "Bulk_SE_Recovery_W.rds")
assert(file.exists(installed) && identical(W,readRDS(installed)),
       "Original installed W != archived W")
assert(is.matrix(W) && identical(dim(W),c(500L,10L)) &&
       identical(colnames(W),components),"W contract mismatch")
Wgenes <- unique(sub("__.*$","",rownames(W)))
assert(length(Wgenes)==499L,"Unexpected frozen W gene count")

read499 <- function(path) {
  x <- fread(path,sep="\t",header=TRUE,quote="",check.names=FALSE,
             showProgress=FALSE)
  assert(nrow(x)==499L && ncol(x)==9466L,
         paste("Expected 499 rows and 9465 samples",path))
  genes <- trimws(as.character(x[[1L]]))
  assert(all(nzchar(genes)) && !anyNA(genes) && !anyDuplicated(genes),
         paste("Ambiguous model genes",path))
  m <- as.matrix(x[,-1L,with=FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- genes
  assert(all(is.finite(m)) && all(m>=0),
         paste("Invalid count/CPM matrix",path))
  m
}
key <- function(x) toupper(trimws(as.character(x)))
readH <- function(path,ids) {
  dt <- fread(path,sep="\t",header=TRUE,quote="",check.names=FALSE,
              showProgress=FALSE)
  assert("sample_id" %in% names(dt) &&
         all(components %in% names(dt)),
         paste("Archived H/fractions column contract changed",path))
  assert(nrow(dt)==9465L && !anyDuplicated(dt$sample_id),
         paste("Archived H/fractions sample universe changed",path))
  ix <- match(ids,dt$sample_id)
  assert(!anyNA(ix),"Archived H/fractions sample mismatch")
  m <- as.matrix(dt[ix,..components])
  storage.mode(m) <- "double"
  rownames(m) <- ids
  colnames(m) <- components
  assert(all(is.finite(m)) && all(m>=0), "Archived H/fractions invalid")
  m
}
canonical <- fread(find_frozen(S9,"01_canonical_library_manifest.tsv.gz"),
                  sep="\t",showProgress=FALSE)
assert(nrow(canonical)==9465L && "sample_id" %in% names(canonical),
       "Original canonical 9465 manifest changed")
ids <- as.character(canonical$sample_id)
assert(!anyDuplicated(ids),"Canonical sample IDs not unique")
sample_qc <- fread(find_frozen(S9,"04_sample_normalization_QC.tsv.gz"),
                   sep="\t",showProgress=FALSE)
needed <- c("sample_id","normalization_library_sum_valid_gene_symbols",
            "effective_library_tmm","effective_library_uq")
assert(all(needed %in% names(sample_qc)) &&
       nrow(sample_qc)==9465L && !anyDuplicated(sample_qc$sample_id),
       "Stage09 sample-QC denominator contract changed")
q <- sample_qc[match(ids,sample_qc$sample_id)]
assert(identical(as.character(q$sample_id),ids),"Stage09 QC IDs not aligned")
denoms <- list(
  CPM=as.numeric(q$normalization_library_sum_valid_gene_symbols),
  TMM_CPM=as.numeric(q$effective_library_tmm),
  UQ_CPM=as.numeric(q$effective_library_uq)
)
for (b in names(denoms))
  assert(all(is.finite(denoms[[b]])) && all(denoms[[b]]>0),
         paste("Invalid frozen size factors",b))

raw <- read499(find_frozen(S9,"09_model_gene_ARCHS4_pseudocounts_canonical.tsv.gz"))
assert(identical(colnames(raw),ids),"Original 499-gene pseudocount columns not in canonical order")
print("V7H_ORIGINAL_ARCHS4_NUMERICAL_REPLAY")
cat("STAGE09_STATUS PASS\nSTAGE10_STATUS PASS\n")
cat("SAMPLES",length(ids),"\nMODEL_GENES",nrow(raw),
    "\nFROZEN_W_DIRECTIONAL_FEATURES",nrow(W),"\n")
cat("SPATIALECOTYPER_VERSION",as.character(packageVersion("SpatialEcoTyper")),"\n")
mapping <- fread(find_frozen(S10,"16A_frozen_W_gene_symbol_harmonization.tsv"),
                 sep="\t",showProgress=FALSE)
required_map <- c("branch","W_gene_index","W_model_gene","input_model_gene",
                  "harmonization")
assert(all(required_map %in% names(mapping)),
       "Original Stage10 W-harmonization ledger schema mismatch")

# Primary status is not changed by this diagnostic. Tight numerical checks
# cover deterministic archived expression transformations, while deconvolution
# equivalence is reported descriptively and reviewed separately.
for (branch in names(B)) {
  cat("V7H_BRANCH_BEGIN",branch,"\n")
  source <- read499(find_frozen(S9,B[[branch]][["expression"]]))
  assert(identical(rownames(source),rownames(raw)) &&
         identical(colnames(source),ids),
         paste("Original Stage09 input axis mismatch",branch))
  recomputed <- sweep(raw,2L,denoms[[branch]],"/")*1e6
  abs_err <- max(abs(recomputed-source))
  rel_err <- max(abs(recomputed-source)/pmax(1,abs(source)))
  cat("STAGE09_NORMALIZATION",branch,
      "MAX_ABSOLUTE_ERROR",format(abs_err,digits=15,scientific=TRUE),
      "MAX_RELATIVE_ERROR",format(rel_err,digits=15,scientific=TRUE),"\n")
  assert(rel_err<=1e-9,
         paste("Stage09 archived normalized expression not numerically reproduced",branch))

  params <- fread(find_frozen(S10,B[[branch]][["scaling"]]),sep="\t")
  required_params <- c("model_gene","weighted_anchor_mean_log2",
                       "weighted_anchor_sd_log2")
  assert(nrow(params)==499L && all(required_params %in% names(params)),
         "Frozen anchor parameters changed")
  ix <- match(rownames(source),params$model_gene)
  assert(!anyNA(ix) && !anyDuplicated(ix),"Stage10 anchor source-gene mapping ambiguous")
  mu <- as.numeric(params$weighted_anchor_mean_log2[ix])
  sd <- as.numeric(params$weighted_anchor_sd_log2[ix])
  assert(all(is.finite(mu)) && all(is.finite(sd)) && all(sd>0),
         "Frozen anchor mean/SD invalid")
  z_input <- sweep(sweep(log2(source+1),1L,mu,"-"),1L,sd,"/")
  assert(all(is.finite(z_input)),"Frozen standardized input invalid")
  map <- mapping[as.character(mapping[["branch"]])==branch]
  setorder(map,W_gene_index)
  assert(nrow(map)==499L &&
         identical(as.character(map$W_model_gene),Wgenes),
         "Original W-gene ledger order mismatch")
  match_idx <- match(as.character(map$input_model_gene),rownames(z_input))
  assert(!anyNA(match_idx) && !anyDuplicated(match_idx),
         "Stage10 frozen mapping cannot be reproduced")
  Z <- z_input[match_idx,,drop=FALSE]
  rownames(Z) <- Wgenes
  assert(identical(colnames(Z),ids),"Stage10 frozen sample order changed")
  frozen_Z <- readRDS(find_frozen(S10,B[[branch]][["z"]]))
  assert(is.matrix(frozen_Z) &&
         identical(dim(frozen_Z),c(499L,9465L)) &&
         identical(rownames(frozen_Z),Wgenes) &&
         identical(colnames(frozen_Z),ids),
         "Archived 499x9465 scaled matrix axis changed")
  z_abs <- max(abs(Z-frozen_Z))
  cat("STAGE10_FROZEN_Z",branch,
      "MAX_ABSOLUTE_ERROR",format(z_abs,digits=15,scientific=TRUE),"\n")
  assert(z_abs<=1e-8,
         paste("Stage10 archived frozen-Z not numerically reproduced",branch))

  oldH <- readH(find_frozen(S10,B[[branch]][["h"]]),ids)
  oldF <- readH(find_frozen(S10,B[[branch]][["fractions"]]),ids)
  derivedF <- oldH/rowSums(oldH)
  frac_err <- max(abs(derivedF-oldF))
  cat("STAGE10_ARCHIVED_H_FRACTION_CONSISTENCY",branch,
      "MAX_ABSOLUTE_ERROR",format(frac_err,digits=15,scientific=TRUE),"\n")
  assert(frac_err<=1e-8,"Archived frozen H and fractions not self-consistent")

  cat("STAGE10_ORIGINAL_DECONVOLUTESE_REPLAY_START",branch,"\n")
  # Exact original Stage10 call (same 9,465 sample batch; no refitting).
  replay <- SpatialEcoTyper::DeconvoluteSE(
    dat=Z,scale=FALSE,W=W,
    nsample.per.run=10000L,sum2one=FALSE,ncores=1L
  )
  replay <- as.matrix(replay)
  storage.mode(replay) <- "double"
  assert(setequal(rownames(replay),ids) &&
         identical(colnames(replay),components),
         "Original DeconvoluteSE output order/membership mismatch")
  replay <- replay[ids,components,drop=FALSE]
  assert(all(is.finite(replay)) &&
         all(replay>=0) && all(rowSums(replay)>0),
         "Original DeconvoluteSE invalid replayed H")
  replayF <- replay/rowSums(replay)
  raw_h_max <- max(abs(replay-oldH))
  H_rel <- max(abs(replay-oldH)/pmax(1,abs(oldH)))
  TVD <- rowSums(abs(replayF-oldF))/2
  argmax <- function(m) max.col(m,ties.method="first")
  concordance <- mean(argmax(replayF)==argmax(oldF))
  cat("STAGE10_ORIGINAL_H_REPLAY",branch,
      "MAX_ABSOLUTE_H_DIFFERENCE",format(raw_h_max,digits=12,scientific=TRUE),
      "MAX_RELATIVE_H_DIFFERENCE",format(H_rel,digits=12,scientific=TRUE),
      "MEDIAN_FRACTION_TVD",format(median(TVD),digits=12,scientific=TRUE),
      "P95_FRACTION_TVD",format(unname(quantile(TVD,.95)),digits=12,scientific=TRUE),
      "MAX_FRACTION_TVD",format(max(TVD),digits=12,scientific=TRUE),
      "DOMINANT_COMPONENT_CONCORDANCE",format(concordance,digits=12),"\n")
  # The original 1e-5 median/1e-3 max tolerances apply to reversed sample-order
  # sentinel testing, not necessarily to H replay; report rather than misuse.
  cat("V7H_BRANCH_NUMERICAL_AUDIT_RECORDED",branch,"\n")
  rm(source,recomputed,z_input,Z,frozen_Z,oldH,oldF,derivedF,replay,replayF,TVD)
  invisible(gc(verbose=FALSE))
}
cat("V7H_FROZEN_STAGE09_NORMALIZATION_AND_STAGE10_Z_REPLAY_PASS\n")
cat("V7H_STAGE10_H_REPLAY_METRICS_RECORDED_PENDING_REVIEW\n")
cat("NO_FROZEN_FILES_MODIFIED\nNO_RECOUNT3_PROJECTION\nCLINICAL_ADMISSIONS 0\n")
cat("V7H_ORIGINAL_ARCHS4_FULL_NUMERICAL_REPLAY_COMPLETE\n")
