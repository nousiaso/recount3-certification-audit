#!/usr/bin/env Rscript
# V7P v1: paired ARCHS4/recount3 technical calibration using independent full500 KL solver (not official .nmf.predict).
# Original Stage10 function bodies below are copied VERBATIM from verified
# 10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R
# SHA256 16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c.
# Frozen ARCHS4: archived 9465-cohort Stage09 normalization + Stage10 Z/H/f.
# Recount3: 454 SRRs -> 143 SRXs raw G026 coverage, original Stage09 formulae.
# NO anchor/W refit, NO patient admissions.
suppressPackageStartupMessages({
  library(data.table)
  library(SpatialEcoTyper)
  library(NMF)
  library(jsonlite)
})
argv <- commandArgs(trailingOnly=TRUE)
if (length(argv)!=3L) stop("Usage: Rscript V7P.R ROOT V7J_INPUT_DIR V7P_OUTPUT_DIR",call.=FALSE)
ROOT <- normalizePath(argv[[1L]],mustWork=TRUE)
OUT <- normalizePath(argv[[2L]],mustWork=TRUE)
FINAL <- normalizePath(argv[[3L]],mustWork=TRUE)
if (identical(OUT,FINAL)) stop("Refuse reuse V7J prepared directory for V7P output",call.=FALSE)
REL <- file.path(ROOT,"RELEASES","ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0","PROJECT_TREE")
S9 <- "10_preprojection_normalization_v1_0_1"
S10 <- "11_pan_cancer_fixed_anchor_projection_v1_0_4"
BRANCHES <- c("CPM","TMM_CPM","UQ_CPM")
COMPONENTS <- c("NonSE",paste0("SE",1:9))
EXPECTED_MODEL_GENES <- 499L
NSAMPLE_PER_RUN <- 10000L
NCORES <- 1L
CLR_FLOOR <- 1e-8
RECONSTRUCTION_EPS <- 1e-12
SUM_TO_ONE_TOLERANCE <- 1e-8
W_FILE <- "/home/tpn890/orestis/ST_ECOTYPES/02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds"
assert <- function(ok,msg) if (!isTRUE(ok)) stop(msg,call.=FALSE)
find_frozen <- function(stage,name) {
  p <- c(file.path(ROOT,stage,name),file.path(REL,stage,name))
  available <- p[file.exists(p)]
  if (!length(available)) stop("Required frozen asset missing: ",stage,"/",name,call.=FALSE)
  available[[1L]]
}
branch_assets <- list(
 CPM=list(expr="10_model_gene_CPM_canonical.tsv.gz",
          z="06_fixed_anchor_scaled_input_CPM.rds",
          h="09_raw_H_CPM.tsv.gz",
          f="10_fraction_scores_CPM.tsv.gz",
          params="03_pan_cancer_anchor_scaling_CPM.tsv.gz"),
 TMM_CPM=list(expr="11_model_gene_TMM_CPM_canonical.tsv.gz",
              z="07_fixed_anchor_scaled_input_TMM_CPM.rds",
              h="11_raw_H_TMM_CPM.tsv.gz",
              f="12_fraction_scores_TMM_CPM.tsv.gz",
              params="04_pan_cancer_anchor_scaling_TMM_CPM.tsv.gz"),
 UQ_CPM=list(expr="12_model_gene_UQ_CPM_canonical.tsv.gz",
             z="08_fixed_anchor_scaled_input_UQ_CPM.rds",
             h="13_raw_H_UQ_CPM.tsv.gz",
             f="14_fraction_scores_UQ_CPM.tsv.gz",
             params="05_pan_cancer_anchor_scaling_UQ_CPM.tsv.gz")
)
assert(file.exists(W_FILE),"Frozen W file missing")
W <- readRDS(W_FILE)
installed_w <- file.path(system.file("extdata",package="SpatialEcoTyper"),
                         "Bulk_SE_Recovery_W.rds")
assert(file.exists(installed_w) && identical(W,readRDS(installed_w)),
       "Frozen and installed W matrices are not identical")
assert(as.character(packageVersion("SpatialEcoTyper"))=="1.0.3",
       "SpatialEcoTyper package differs from original 1.0.3")
assert(is.matrix(W) && identical(dim(W),c(500L,10L)) &&
       identical(colnames(W),COMPONENTS),"W direction/components mismatch")
W_genes <- unique(sub("__.*$","",rownames(W)))
assert(length(W_genes)==499L,"Expected 499 original W base symbols")
manifest <- fread(file.path(OUT,"143_pair_manifest_LOCAL_ONLY.tsv"),sep="\t",
                  quote="",check.names=FALSE)
cols <- c("ARCHS4_GSM","library_id","n_technical_runs","UQ_evaluable_recount3")
assert(all(cols %in% names(manifest)) && nrow(manifest)==143L,
       "143-pair recount3 prepared manifest invalid")
gsm <- as.character(manifest$ARCHS4_GSM)
assert(uniqueN(gsm)==143L && uniqueN(manifest$library_id)==143L &&
       sum(manifest$n_technical_runs)==454L,"143 SRX/GSM/454 SRR contract mismatch")
uq_rec <- toupper(as.character(manifest$UQ_evaluable_recount3))=="TRUE"
assert(!anyNA(uq_rec),"Invalid recount3 upper-quartile evaluability mask")
reference_canonical <- fread(find_frozen(S9,"01_canonical_library_manifest.tsv.gz"),
                             sep="\t",quote="",showProgress=FALSE)
assert(nrow(reference_canonical)==9465L &&
       uniqueN(reference_canonical$sample_id)==9465L &&
       all(gsm %in% reference_canonical$sample_id),
       "All 143 GSMs must occur in frozen 9465 ARCHS4 samples")
gene_audit <- fread(find_frozen(S9,"08_model_gene_availability.tsv"),
                    sep="\t",quote="",showProgress=FALSE)
assert(nrow(gene_audit)==499L &&
       all(c("model_gene","ARCHS4_source_gene","resolution") %in% names(gene_audit)),
       "Frozen original Stage09 model-gene audit changed")
frozen_mapping <- fread(find_frozen(S10,"16A_frozen_W_gene_symbol_harmonization.tsv"),
                        sep="\t",quote="",showProgress=FALSE)
assert(nrow(frozen_mapping)==1497L &&
       all(vapply(BRANCHES,function(b)
           sum(frozen_mapping[["branch"]]==b)==499L,logical(1))),
       "Original Stage10 3-branch W harmonization ledger differs")
sample_qc <- fread(find_frozen(S9,"04_sample_normalization_QC.tsv.gz"),
                   sep="\t",quote="",showProgress=FALSE)
uq_qc <- sample_qc[match(gsm,as.character(sample_qc$sample_id))]
assert(nrow(uq_qc)==143L &&
       identical(as.character(uq_qc$sample_id),gsm) &&
       "upper_quartile_fraction" %in% names(uq_qc),
       "Frozen original ARCHS4 upper-quartile evaluability missing")
uq_arch <- as.numeric(uq_qc$upper_quartile_fraction)>0
assert(!anyNA(uq_arch),"Frozen ARCHS4 upper-quartile evaluability incomplete")
cat("V7J_FROZEN_ARCHS4_BASELINE_143_IN_9465_PASS\n")
cat("V7J_ORIGINAL_W_500x10_PACKAGE_1.0.3_PASS\n")

key <- function(x) toupper(trimws(as.character(x)))
read_499 <- function(path,expected_sample_ids,full_atlas=FALSE) {
  x <- fread(path,sep="\t",header=TRUE,quote="",check.names=FALSE,
             showProgress=FALSE)
  expected_n <- if (full_atlas) 9465L else 143L
  assert(nrow(x)==499L && ncol(x)==expected_n+1L,
         paste("Unexpected model matrix axes",path))
  genes <- trimws(as.character(x[[1L]]))
  assert(!anyNA(genes) && all(nzchar(genes)) && !anyDuplicated(key(genes)),
         paste("Duplicated model-gene symbols",path))
  have <- names(x)[-1L]
  assert(all(expected_sample_ids %in% have) &&
         !anyDuplicated(have),"Missing model input GSMs")
  m <- as.matrix(x[,..expected_sample_ids])
  storage.mode(m) <- "double"
  rownames(m) <- genes
  colnames(m) <- expected_sample_ids
  assert(all(is.finite(m)) && all(m>=0),"Nonfinite/negative model input")
  m
}
read_scores <- function(path) {
  x <- fread(path,sep="\t",header=TRUE,quote="",check.names=FALSE,
             showProgress=FALSE)
  assert(nrow(x)==9465L && "sample_id" %in% names(x) &&
         !anyDuplicated(x$sample_id) &&
         all(COMPONENTS %in% names(x)),"Original H/fraction archive invalid")
  ix <- match(gsm,as.character(x$sample_id))
  assert(!anyNA(ix),"Archived 143 GSMs not all found")
  z <- as.matrix(x[ix,..COMPONENTS])
  storage.mode(z) <- "double"
  rownames(z) <- gsm
  colnames(z) <- COMPONENTS
  assert(all(is.finite(z)) && all(z>=0),"Invalid original H/fractions")
  z
}
read_params <- function(branch) {
  p <- fread(find_frozen(S10,branch_assets[[branch]]$params),
             sep="\t",quote="",showProgress=FALSE)
  assert(nrow(p)==499L && !anyDuplicated(key(p$model_gene)) &&
         all(c("weighted_anchor_mean_log2","weighted_anchor_sd_log2") %in% names(p)),
         "Original pan-cancer anchor scaling metadata invalid")
  mu <- as.numeric(p$weighted_anchor_mean_log2)
  sd <- as.numeric(p$weighted_anchor_sd_log2)
  assert(all(is.finite(mu)) && all(is.finite(sd) & sd>0),
         "Original frozen anchor parameters invalid")
  names(mu) <- key(p$model_gene)
  names(sd) <- key(p$model_gene)
  list(mu=mu,sd=sd)
}
scaled_model <- function(X,params) {
  idx <- match(key(rownames(X)),names(params$mu))
  assert(!anyNA(idx) && !anyDuplicated(idx),
         "Frozen scaling gene identity not unique")
  Z <- sweep(sweep(log2(X+1),1L,params$mu[idx],"-"),1L,params$sd[idx],"/")
  assert(all(is.finite(Z)),"Invalid fixed-anchor scaled expression")
  Z
}
archiveZ <- function(branch) {
  full <- readRDS(find_frozen(S10,branch_assets[[branch]]$z))
  assert(is.matrix(full) &&
         identical(dim(full),c(499L,9465L)) &&
         identical(rownames(full),W_genes) &&
         all(gsm %in% colnames(full)),
         paste("Frozen Z full atlas axis mismatch",branch))
  z <- full[,gsm,drop=FALSE]
  rm(full)
  z
}
aligned_directional <- function(Z) {
  v <- as.matrix(NMF::posneg(Z[W_genes,,drop=FALSE]))
  d <- duplicated(rownames(v))
  rownames(v)[!d] <- paste0(rownames(v)[!d],"__pos")
  rownames(v)[d] <- paste0(rownames(v)[d],"__neg")
  idx <- match(rownames(W),rownames(v))
  assert(!anyNA(idx),"Original directional alignment missing W row")
  v <- v[idx,,drop=FALSE]
  assert(identical(rownames(v),rownames(W)),
         "Original directional output W-row order mismatch")
  v
}
summarize <- function(v) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (!length(v)) return(c(median=NA_real_,p90=NA_real_,p95=NA_real_,max=NA_real_))
  c(median=median(v),p90=unname(quantile(v,.9)),
    p95=unname(quantile(v,.95)),max=max(v))
}
# No specimen-level bootstrap: SRX is not established as independent patient.

# The following EXACT original Stage10 function bodies are not recreated from
# summaries. They were extracted from the source commit with verified SHA.
score_fractions <- function(raw_h) {
  totals <- rowSums(raw_h)

  if (
    any(!is.finite(totals)) ||
      any(totals <= 0)
  ) {
    stop("Raw H contains a non-positive sample total.", call. = FALSE)
  }

  fractions <- raw_h / totals

  if (
    any(!is.finite(fractions)) ||
      any(fractions < 0)
  ) {
    stop("Invalid compositional scores derived from raw H.", call. = FALSE)
  }

  fractions
}

clr_matrix <- function(score_matrix, floor_value = CLR_FLOOR) {
  x <- as.matrix(score_matrix)
  storage.mode(x) <- "double"

  if (
    any(!is.finite(x)) ||
      any(x < 0)
  ) {
    stop("Invalid values supplied to clr_matrix().", call. = FALSE)
  }

  x[x < floor_value] <- floor_value
  x <- x / rowSums(x)
  lx <- log(x)
  lx - rowMeans(lx)
}

dominant_component <- function(score_matrix) {
  COMPONENTS[
    max.col(
      score_matrix,
      ties.method = "first"
    )
  ]
}

harmonize_expression_to_W <- function(
    scaled_matrix,
    W,
    model_gene_audit,
    branch_name
) {
  required_audit_columns <- c(
    "model_gene",
    "ARCHS4_source_gene",
    "resolution"
  )

  missing_audit_columns <- setdiff(
    required_audit_columns,
    names(model_gene_audit)
  )

  if (length(missing_audit_columns) > 0L) {
    stop(
      branch_name,
      ": model-gene audit lacks columns: ",
      paste(missing_audit_columns, collapse = ", "),
      call. = FALSE
    )
  }

  W_model_genes <- unique(
    sub(
      "__.*$",
      "",
      rownames(W)
    )
  )

  input_genes <- rownames(
    scaled_matrix
  )

  normalize_symbol_key <- function(x) {
    toupper(
      trimws(
        as.character(x)
      )
    )
  }

  W_keys <- normalize_symbol_key(
    W_model_genes
  )

  input_keys <- normalize_symbol_key(
    input_genes
  )

  if (anyDuplicated(W_keys)) {
    duplicated_W <- unique(
      W_model_genes[
        duplicated(W_keys) |
          duplicated(
            W_keys,
            fromLast = TRUE
          )
      ]
    )

    stop(
      branch_name,
      ": frozen W contains case-insensitive duplicate base-gene symbols: ",
      paste(
        duplicated_W,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(input_keys)) {
    duplicated_input <- unique(
      input_genes[
        duplicated(input_keys) |
          duplicated(
            input_keys,
            fromLast = TRUE
          )
      ]
    )

    stop(
      branch_name,
      ": Stage 08 matrix contains case-insensitive duplicate row symbols: ",
      paste(
        duplicated_input,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  input_key_to_gene <- setNames(
    input_genes,
    input_keys
  )

  audit <- copy(
    as.data.table(
      model_gene_audit
    )
  )

  audit[
    ,
    model_gene := trimws(
      as.character(
        model_gene
      )
    )
  ]

  audit[
    ,
    ARCHS4_source_gene := trimws(
      as.character(
        ARCHS4_source_gene
      )
    )
  ]

  audit[
    is.na(ARCHS4_source_gene) |
      ARCHS4_source_gene %in% c(
        "",
        "NA",
        "NAN",
        "NONE",
        "NULL"
      ),
    ARCHS4_source_gene := NA_character_
  ]

  audit[
    ,
    model_gene_key := normalize_symbol_key(
      model_gene
    )
  ]

  audit[
    ,
    ARCHS4_source_gene_key := normalize_symbol_key(
      ARCHS4_source_gene
    )
  ]

  mapping_rows <- lapply(
    seq_along(
      W_model_genes
    ),
    function(i) {
      W_gene <- W_model_genes[[i]]
      W_key <- W_keys[[i]]

      direct_candidate <- unname(
        input_key_to_gene[
          W_key
        ]
      )

      if (
        length(direct_candidate) == 0L ||
          is.na(direct_candidate) ||
          direct_candidate == ""
      ) {
        direct_candidate <- character()
      }

      audit_candidates <- audit[
        !is.na(
          ARCHS4_source_gene
        ) &
          ARCHS4_source_gene_key ==
            W_key &
          model_gene_key %in%
            input_keys,
        unique(
          model_gene
        )
      ]

      candidates <- unique(
        c(
          direct_candidate,
          audit_candidates
        )
      )

      candidate_keys <- normalize_symbol_key(
        candidates
      )

      candidates <- candidates[
        !duplicated(
          candidate_keys
        )
      ]

      if (length(candidates) == 0L) {
        return(
          data.table(
            W_model_gene = W_gene,
            input_model_gene = NA_character_,
            harmonization = "unresolved"
          )
        )
      }

      if (length(candidates) > 1L) {
        stop(
          branch_name,
          ": frozen W gene ",
          W_gene,
          " maps ambiguously to Stage 08 rows: ",
          paste(
            candidates,
            collapse = ", "
          ),
          call. = FALSE
        )
      }

      input_gene <- candidates[[1L]]

      harmonization <- if (
        identical(
          W_gene,
          input_gene
        )
      ) {
        "exact_symbol"
      } else if (
        identical(
          normalize_symbol_key(
            W_gene
          ),
          normalize_symbol_key(
            input_gene
          )
        )
      ) {
        "case_harmonized_to_W"
      } else {
        "alias_harmonized_to_W"
      }

      data.table(
        W_model_gene = W_gene,
        input_model_gene = input_gene,
        harmonization = harmonization
      )
    }
  )

  mapping <- rbindlist(
    mapping_rows,
    use.names = TRUE
  )

  unresolved <- mapping[
    harmonization ==
      "unresolved"
  ]

  if (nrow(unresolved) > 0L) {
    stop(
      branch_name,
      ": ",
      nrow(unresolved),
      " frozen W genes cannot be resolved from the Stage 08 matrix: ",
      paste(
        unresolved$W_model_gene,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  duplicated_source_keys <- normalize_symbol_key(
    mapping$input_model_gene
  )

  duplicated_sources <- mapping[
    duplicated(
      duplicated_source_keys
    ) |
      duplicated(
        duplicated_source_keys,
        fromLast = TRUE
      )
  ]

  if (nrow(duplicated_sources) > 0L) {
    stop(
      branch_name,
      ": at least one Stage 08 row maps to more than one frozen W gene:\n",
      paste(
        capture.output(
          print(
            duplicated_sources,
            row.names = FALSE
          )
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  harmonized_matrix <- scaled_matrix[
    mapping$input_model_gene,
    ,
    drop = FALSE
  ]

  rownames(
    harmonized_matrix
  ) <- mapping$W_model_gene

  if (
    !identical(
      rownames(
        harmonized_matrix
      ),
      W_model_genes
    ) ||
      anyDuplicated(
        normalize_symbol_key(
          rownames(
            harmonized_matrix
          )
        )
      )
  ) {
    stop(
      branch_name,
      ": failed to construct a unique frozen-W-ordered expression matrix.",
      call. = FALSE
    )
  }

  mapping[
    ,
    branch := branch_name
  ]

  mapping[
    ,
    input_row_index := match(
      input_model_gene,
      input_genes
    )
  ]

  mapping[
    ,
    W_gene_index := seq_len(.N)
  ]

  setcolorder(
    mapping,
    c(
      "branch",
      "W_gene_index",
      "W_model_gene",
      "input_row_index",
      "input_model_gene",
      "harmonization"
    )
  )

  list(
    matrix = harmonized_matrix,
    mapping = mapping,
    W_model_gene_count = length(
      W_model_genes
    ),
    exact_count = sum(
      mapping$harmonization ==
        "exact_symbol"
    ),
    case_harmonized_count = sum(
      mapping$harmonization ==
        "case_harmonized_to_W"
    ),
    alias_harmonized_count = sum(
      mapping$harmonization ==
        "alias_harmonized_to_W"
    )
  )
}

make_directional_bundle <- function(
    scaled_matrix,
    W
) {
  model_genes <- unique(
    sub(
      "__.*$",
      "",
      rownames(W)
    )
  )

  missing_model_genes <- setdiff(
    model_genes,
    rownames(
      scaled_matrix
    )
  )

  if (length(missing_model_genes) > 0L) {
    return(
      list(
        V = NULL,
        W = NULL,
        retained_features = character(),
        retained_feature_count = 0L,
        supported_gene_count = length(
          intersect(
            model_genes,
            rownames(
              scaled_matrix
            )
          )
        ),
        expected_gene_count = length(
          model_genes
        ),
        missing_model_genes = missing_model_genes,
        missing_W_features = rownames(W),
        zero_variance_W_features = character()
      )
    )
  }

  testdat <- scaled_matrix[
    model_genes,
    ,
    drop = FALSE
  ]

  to_predict_all <- as.matrix(
    NMF::posneg(testdat)
  )

  duplicate_index <- duplicated(
    rownames(
      to_predict_all
    )
  )

  rownames(
    to_predict_all
  )[!duplicate_index] <- paste0(
    rownames(
      to_predict_all
    )[!duplicate_index],
    "__pos"
  )

  rownames(
    to_predict_all
  )[duplicate_index] <- paste0(
    rownames(
      to_predict_all
    )[duplicate_index],
    "__neg"
  )

  row_variance <- apply(
    to_predict_all,
    1L,
    stats::var
  )

  variable_rows <- is.finite(
    row_variance
  ) &
    row_variance > 0

  zero_variance_features <- rownames(
    to_predict_all
  )[
    !variable_rows
  ]

  to_predict <- to_predict_all[
    variable_rows,
    ,
    drop = FALSE
  ]

  features <- intersect(
    rownames(
      to_predict
    ),
    rownames(W)
  )

  missing_W_features <- setdiff(
    rownames(W),
    features
  )

  zero_variance_W_features <- intersect(
    rownames(W),
    zero_variance_features
  )

  W_subset <- W[
    match(
      features,
      rownames(W)
    ),
    ,
    drop = FALSE
  ]

  V <- to_predict[
    match(
      features,
      rownames(
        to_predict
      )
    ),
    ,
    drop = FALSE
  ]

  W_subset[
    is.na(
      W_subset
    )
  ] <- 0

  V[
    is.na(V)
  ] <- 0

  list(
    V = V,
    W = W_subset,
    retained_features = features,
    retained_feature_count = length(
      features
    ),
    supported_gene_count = length(
      model_genes
    ),
    expected_gene_count = length(
      model_genes
    ),
    missing_model_genes = missing_model_genes,
    missing_W_features = missing_W_features,
    zero_variance_W_features = zero_variance_W_features
  )
}

project_raw_h <- function(
    scaled_matrix,
    W,
    branch_name
) {
  if (ncol(scaled_matrix) > NSAMPLE_PER_RUN) {
    stop(
      branch_name,
      ": canonical sample count exceeds NSAMPLE_PER_RUN; this would create multiple internal prediction batches.",
      call. = FALSE
    )
  }

  raw_h <- SpatialEcoTyper::DeconvoluteSE(
    dat = scaled_matrix,
    scale = FALSE,
    W = W,
    nsample.per.run = NSAMPLE_PER_RUN,
    sum2one = FALSE,
    ncores = NCORES
  )

  raw_h <- as.matrix(raw_h)
  storage.mode(raw_h) <- "double"

  if (
    !identical(
      colnames(raw_h),
      COMPONENTS
    )
  ) {
    stop(
      branch_name,
      ": projected components differ from the locked component order.",
      call. = FALSE
    )
  }

  if (!setequal(rownames(raw_h), colnames(scaled_matrix))) {
    stop(
      branch_name,
      ": projected sample IDs differ from input IDs.",
      call. = FALSE
    )
  }

  raw_h <- raw_h[
    colnames(scaled_matrix),
    COMPONENTS,
    drop = FALSE
  ]

  if (
    any(!is.finite(raw_h)) ||
      any(raw_h < 0)
  ) {
    stop(
      branch_name,
      ": raw H contains negative or non-finite values.",
      call. = FALSE
    )
  }

  raw_h
}

calculate_reconstruction_adequacy <- function(
    V,
    V_hat,
    fractions
) {
  if (!identical(dim(V), dim(V_hat))) {
    stop("Observed and reconstructed directional matrices differ in shape.", call. = FALSE)
  }

  V_hat_safe <- pmax(
    V_hat,
    RECONSTRUCTION_EPS
  )

  kl_terms <- V_hat - V
  positive <- V > 0

  kl_terms[positive] <- (
    V[positive] *
      log(
        V[positive] /
          V_hat_safe[positive]
      ) -
      V[positive] +
      V_hat[positive]
  )

  kl_divergence <- colSums(kl_terms)
  input_mass <- colSums(V)
  reconstructed_mass <- colSums(V_hat)

  residual <- V - V_hat
  squared_error <- colSums(residual^2)
  rmse <- sqrt(colMeans(residual^2))
  mae <- colMeans(abs(residual))

  V_centered <- sweep(
    V,
    2L,
    colMeans(V),
    "-"
  )

  V_hat_centered <- sweep(
    V_hat,
    2L,
    colMeans(V_hat),
    "-"
  )

  correlation_denominator <- sqrt(
    colSums(V_centered^2) *
      colSums(V_hat_centered^2)
  )

  pearson <- colSums(
    V_centered *
      V_hat_centered
  ) / correlation_denominator

  total_sum_squares <- colSums(
    V_centered^2
  )

  reconstruction_r2 <- 1 -
    squared_error /
      total_sum_squares

  cosine_denominator <- sqrt(
    colSums(V^2) *
      colSums(V_hat^2)
  )

  cosine_similarity <- colSums(
    V *
      V_hat
  ) / cosine_denominator

  score_entropy <- -rowSums(
    fractions *
      log(
        pmax(
          fractions,
          RECONSTRUCTION_EPS
        )
      )
  ) / log(ncol(fractions))

  top_order <- t(
    apply(
      fractions,
      1L,
      order,
      decreasing = TRUE
    )
  )

  top1_index <- top_order[, 1L]
  top2_index <- top_order[, 2L]

  top1_fraction <- fractions[
    cbind(
      seq_len(nrow(fractions)),
      top1_index
    )
  ]

  top2_fraction <- fractions[
    cbind(
      seq_len(nrow(fractions)),
      top2_index
    )
  ]

  data.table(
    sample_id = colnames(V),
    directional_input_mass = input_mass,
    directional_reconstructed_mass = reconstructed_mass,
    fixed_W_KL_divergence = kl_divergence,
    normalized_KL_per_input_mass = kl_divergence /
      pmax(input_mass, RECONSTRUCTION_EPS),
    reconstruction_RMSE = rmse,
    reconstruction_MAE = mae,
    reconstruction_Pearson = pearson,
    reconstruction_cosine = cosine_similarity,
    reconstruction_R2 = reconstruction_r2,
    normalized_score_entropy = score_entropy,
    dominant_component = COMPONENTS[top1_index],
    dominant_fraction = top1_fraction,
    second_component = COMPONENTS[top2_index],
    second_fraction = top2_fraction,
    top1_minus_top2 = top1_fraction - top2_fraction,
    dominance_class = ifelse(
      top1_fraction >= 0.5,
      "dominant_ge_0.5",
      "mixed_lt_0.5"
    )
  )
}

normalization_sensitivity <- function(
    primary_scores,
    alternative_scores,
    comparison_name,
    evaluable_mask = NULL,
    non_evaluable_reason = NULL
) {
  if (
    !identical(rownames(primary_scores), rownames(alternative_scores)) ||
      !identical(colnames(primary_scores), colnames(alternative_scores))
  ) {
    stop(
      comparison_name,
      ": score matrices are not aligned.",
      call. = FALSE
    )
  }

  n_samples <- nrow(primary_scores)

  if (is.null(evaluable_mask)) {
    evaluable_mask <- rep(TRUE, n_samples)
  }

  evaluable_mask <- as.logical(evaluable_mask)

  if (
    length(evaluable_mask) != n_samples ||
      anyNA(evaluable_mask) ||
      sum(evaluable_mask) < 2L
  ) {
    stop(
      comparison_name,
      ": invalid evaluability mask.",
      call. = FALSE
    )
  }

  if (is.null(non_evaluable_reason)) {
    non_evaluable_reason <- rep(NA_character_, n_samples)
  }

  if (length(non_evaluable_reason) == 1L) {
    non_evaluable_reason <- rep(non_evaluable_reason, n_samples)
  }

  if (length(non_evaluable_reason) != n_samples) {
    stop(
      comparison_name,
      ": non-evaluable reason vector has the wrong length.",
      call. = FALSE
    )
  }

  non_evaluable_reason <- as.character(non_evaluable_reason)
  non_evaluable_reason[evaluable_mask] <- NA_character_

  primary_clr <- clr_matrix(primary_scores)
  alternative_clr <- clr_matrix(alternative_scores)

  absolute_difference <- abs(primary_scores - alternative_scores)

  diagnostic_tvd <- 0.5 * rowSums(absolute_difference)
  diagnostic_max_difference <- apply(absolute_difference, 1L, max)
  diagnostic_aitchison <- sqrt(
    rowSums((primary_clr - alternative_clr)^2)
  )
  diagnostic_dominant_concordant <-
    dominant_component(primary_scores) ==
      dominant_component(alternative_scores)

  sample_table <- data.table(
    sample_id = rownames(primary_scores),
    comparison = comparison_name,
    comparison_evaluable = evaluable_mask,
    non_evaluable_reason = non_evaluable_reason,
    diagnostic_total_variation_distance = diagnostic_tvd,
    diagnostic_maximum_component_absolute_difference =
      diagnostic_max_difference,
    diagnostic_aitchison_distance = diagnostic_aitchison,
    diagnostic_dominant_component_concordant =
      diagnostic_dominant_concordant,
    total_variation_distance = ifelse(
      evaluable_mask,
      diagnostic_tvd,
      NA_real_
    ),
    maximum_component_absolute_difference = ifelse(
      evaluable_mask,
      diagnostic_max_difference,
      NA_real_
    ),
    aitchison_distance = ifelse(
      evaluable_mask,
      diagnostic_aitchison,
      NA_real_
    ),
    dominant_primary = dominant_component(primary_scores),
    dominant_alternative = dominant_component(alternative_scores),
    dominant_component_concordant = ifelse(
      evaluable_mask,
      diagnostic_dominant_concordant,
      NA
    )
  )

  eval_indices <- which(evaluable_mask)

  component_table <- rbindlist(
    lapply(
      COMPONENTS,
      function(component_name) {
        x <- primary_scores[eval_indices, component_name]
        y <- alternative_scores[eval_indices, component_name]
        x_clr <- primary_clr[eval_indices, component_name]
        y_clr <- alternative_clr[eval_indices, component_name]

        data.table(
          comparison = comparison_name,
          component = component_name,
          n_samples_total = n_samples,
          n_samples_evaluable = length(eval_indices),
          fraction_Pearson = stats::cor(x, y, method = "pearson"),
          fraction_Spearman = stats::cor(x, y, method = "spearman"),
          CLR_Pearson = stats::cor(x_clr, y_clr, method = "pearson"),
          median_absolute_fraction_difference = stats::median(abs(x - y)),
          maximum_absolute_fraction_difference = max(abs(x - y))
        )
      }
    )
  )

  list(
    sample = sample_table,
    component = component_table
  )
}
# End of byte-preserved original Stage10 helper function bodies.


# V7P is an independently specified, fixed-W, nonnegative KL coefficient
# estimator. It retains all 500 original directional W rows (including zero
# variance rows). It DOES NOT replace the frozen official reference.
full500_bundle <- function(Z,W) {
  assert(identical(rownames(Z),W_genes),"W base-gene Z order changed")
  v <- as.matrix(NMF::posneg(Z))
  d <- duplicated(rownames(v))
  rownames(v)[!d] <- paste0(rownames(v)[!d],"__pos")
  rownames(v)[d] <- paste0(rownames(v)[d],"__neg")
  ix <- match(rownames(W),rownames(v))
  assert(!anyNA(ix),"One of the original 500 directional W rows is missing")
  v <- v[ix,,drop=FALSE]
  assert(identical(rownames(v),rownames(W))&&
         all(is.finite(v))&&all(v>=0),"Invalid full W directional values")
  vv <- apply(v,1L,var)
  list(V=v,W=W,supported_gene_count=499L,retained_feature_count=500L,
       zero_variance_W_features=rownames(v)[!is.finite(vv)|vv==0])
}
fixed500_project_raw_h <- function(Z,W,branch_name,h_init=NULL) {
  b <- full500_bundle(Z,W)
  V <- b$V
  assert(identical(dim(W),c(500L,10L))&&all(is.finite(W))&&all(W>=0),
         "Frozen W invalid for independent KL optimization")
  ns <- ncol(V)
  if (!is.null(h_init))
    assert(is.matrix(h_init)&&identical(dim(h_init),c(ns,10L))&&
           identical(rownames(h_init),colnames(V))&&
           identical(colnames(h_init),COMPONENTS),
           paste("Frozen archived H start is not sample-aligned",branch_name))
  fits <- parallel::mclapply(seq_len(ns),function(i) {
    v <- as.numeric(V[,i])
    nz <- v>0
    eps <- 1e-12
    objective <- function(h) {
      q <- pmax(as.numeric(W %*% h),eps)
      sum(q-v)+sum(v[nz]*log(v[nz]/q[nz]))
    }
    gradient <- function(h) {
      q <- pmax(as.numeric(W %*% h),eps)
      as.numeric(crossprod(W,1-v/q))
    }
    initial <- if(is.null(h_init)) rep(0.1,10L) else
      pmax(as.numeric(h_init[i,]),1e-7)
    fit <- stats::optim(initial,objective,gradient,method="L-BFGS-B",
      lower=rep(0,10L),
      control=list(maxit=2000L,factr=1e7,pgtol=1e-8))
    c(as.numeric(pmax(fit$par,0)),fit$convergence,fit$value)
  },mc.cores=4L,mc.preschedule=TRUE)
  assert(all(vapply(fits,function(x)is.numeric(x)&&length(x)==12L&&
         all(is.finite(x)),logical(1))),
         paste("KL optimizer produced invalid fit",branch_name))
  fm <- do.call(rbind,fits)
  assert(all(fm[,11]==0),paste("KL optimizer convergence HOLD",branch_name))
  H <- fm[,1:10,drop=FALSE]
  storage.mode(H) <- "double"
  dimnames(H) <- list(colnames(V),COMPONENTS)
  assert(all(H>=0)&&all(is.finite(H))&&all(rowSums(H)>0),
         paste("KL coefficients invalid",branch_name))
  cat("V7P_FIXED500_KL_OPTIMIZATION_PASS",branch_name,
      "SAMPLES",ns,"W_ROWS",nrow(V),
      "CONSTANT_W_ROWS_RETAINED",length(b$zero_variance_W_features),
      "NONZERO_OPTIM_CODES",sum(fm[,11]!=0),"\n")
  H
}

per_pair <- list()
per_gene <- list()
aggregates <- list()
per_branch <- list()
for (branch in BRANCHES) {
  cat("V7P_BEGIN_BRANCH",branch,"\n")
  a <- branch_assets[[branch]]
  params <- read_params(branch)
  arch_X <- read_499(find_frozen(S9,a$expr),gsm,full_atlas=TRUE)
  rec_X <- read_499(file.path(OUT,paste0("499_recount3_",branch,"_LOCAL_ONLY.tsv.gz")),
                    gsm,full_atlas=FALSE)
  assert(identical(rownames(arch_X),rownames(rec_X)),
         paste("Gene row mismatch between original ARCHS4 and recount3",branch))
  az_recalc <- scaled_model(arch_X,params)
  original_arch <- harmonize_expression_to_W(az_recalc,W,gene_audit,branch)
  az <- archiveZ(branch)
  zerror <- max(abs(original_arch$matrix-az))
  cat("V7P_ORIGINAL_ARCHS4_FROZEN_Z_REPLAY",branch,
      "MAX_ABSOLUTE_ERROR",format(zerror,scientific=TRUE,digits=12),"\n")
  assert(zerror<=1e-8,
         paste("Original 143-ARCHS4 normalized baseline not identical to archived Stage10 Z",branch))
  rec_Zinput <- scaled_model(rec_X,params)
  harmonized <- harmonize_expression_to_W(rec_Zinput,W,gene_audit,branch)
  rz <- harmonized$matrix
  assert(identical(rownames(rz),W_genes) && identical(colnames(rz),gsm),
         "Recount3 W gene/paired-sample order mismatch")
  ledger_ids <- which(as.character(frozen_mapping[["branch"]])==branch)
  ledger <- copy(frozen_mapping[ledger_ids])
  setorder(ledger,W_gene_index)
  assert(nrow(ledger)==499L &&
         identical(as.integer(ledger$W_gene_index),seq_len(499L)) &&
         identical(as.character(ledger$W_model_gene),W_genes) &&
         identical(as.character(harmonized$mapping$input_model_gene),
                   as.character(ledger$input_model_gene)) &&
         identical(as.character(original_arch$mapping$input_model_gene),
                   as.character(ledger$input_model_gene)),
         paste("Recount3 or ARCHS4 harmonization differs from original Stage10 ledger",branch))
  cat("V7P_ORIGINAL_STAGE10_GENE_HARMONIZATION_PASS",branch,"499\n")

  ah <- read_scores(find_frozen(S10,a$h))
  af <- read_scores(find_frozen(S10,a$f))
  ferror <- max(abs(score_fractions(ah)-af))
  assert(ferror<=1e-8,
         paste("Original Stage10 H-fraction mismatch",branch))
  cat("V7P_ARCHIVED_ARCHS4_H_FRACTION_PASS",branch,
      "MAX_ABSOLUTE_ERROR",format(ferror,scientific=TRUE,digits=12),"\n")

  # Original archived 9465-context H is the immutable reference.
  # Independent full-500 KL replay is a controlled fidelity check only.
  arch_h500 <- fixed500_project_raw_h(az,W,
    paste0("ARCHS4_ARCHIVE_FULL500_KL_CHECK_",branch),h_init=ah)
  replay_tvd <- rowSums(abs(score_fractions(arch_h500)-af))/2
  cat("V7P_ORIGINAL_ARCHIVE_500_KL_REPLAY",branch,
      "MAX_FRACTION_TVD",format(max(replay_tvd),digits=12,scientific=TRUE),"\n")
  assert(max(replay_tvd)<=5e-4,
         paste("Independent KL full500 does not sufficiently replay archived original H; HOLD",branch))

  rec_dir <- full500_bundle(rz,W)
  cat("V7P_RECOUNT3_DIRECTIONAL_FEATURES",branch,
      "SUPPORTED_MODEL_GENES",rec_dir$supported_gene_count,
      "RETAINED_W_FEATURES",rec_dir$retained_feature_count,
      "ZERO_VARIANCE_W_FEATURES",length(rec_dir$zero_variance_W_features),"\n")
  assert(rec_dir$supported_gene_count==499L &&
         rec_dir$retained_feature_count==nrow(W),
         paste("Recount3 missing/nonvariable original W directional features; HOLD",branch))
  rh <- fixed500_project_raw_h(rz,W,paste0("RECOUNT3_",branch))
  rf <- score_fractions(rh)
  assert(max(abs(rowSums(rf)-1))<SUM_TO_ONE_TOLERANCE,
         paste("Recount3 fraction totals incorrect",branch))
  reconstructed <- rec_dir$W %*% t(rh)
  dimnames(reconstructed) <- dimnames(rec_dir$V)
  adequacy <- calculate_reconstruction_adequacy(rec_dir$V,reconstructed,rf)
  assert(all(is.finite(as.matrix(adequacy[, .(
     normalized_KL_per_input_mass,reconstruction_RMSE,
     reconstruction_MAE,reconstruction_R2,reconstruction_cosine
  )]))),"Nonfinite recount3 original adequacy metrics")
  cat("V7P_RECOUNT3_INDEPENDENT_FIXED500_KL_PASS",branch,"\n")

  # Original archived ARCHS4 H was obtained from the full 9465 atlas. For the
  # 143-subset we align all 500 features without REFITTING or SELECTING W.
  arch_v <- aligned_directional(az)
  arch_vhat <- W %*% t(ah)
  dimnames(arch_vhat) <- dimnames(arch_v)
  arch_adequacy <- calculate_reconstruction_adequacy(arch_v,arch_vhat,af)

  if (branch=="CPM") {
    reference_adequacy <- fread(
      find_frozen(S10,"19_primary_CPM_fixed_W_model_adequacy.tsv.gz"),
      sep="\t",quote="",showProgress=FALSE)
    ia <- match(gsm,as.character(reference_adequacy$sample_id))
    assert(!anyNA(ia),"Original Stage10 CPM adequacy 143 GSM subset absent")
    compare_cols <- c("normalized_KL_per_input_mass","reconstruction_RMSE",
                      "reconstruction_MAE","reconstruction_R2","reconstruction_cosine")
    ref_values <- as.matrix(reference_adequacy[ia,..compare_cols])
    calc_values <- as.matrix(arch_adequacy[,..compare_cols])
    adeq_err <- max(abs(ref_values-calc_values))
    cat("V7P_ARCHS4_FROZEN_PRIMARY_ADEQUACY_REPLAY_MAXABS",
        format(adeq_err,digits=12,scientific=TRUE),"\n")
    assert(is.finite(adeq_err) && adeq_err<=1e-8,
           "Original ARCHS4 143 frozen primary adequacy does not replay")
  }

  tvd <- rowSums(abs(rf-af))/2
  dom_arch <- dominant_component(af)
  dom_rec <- dominant_component(rf)
  joint <- if (branch=="UQ_CPM") uq_arch & uq_rec else rep(TRUE,143L)
  assert(sum(joint)>=2L,"Jointly evaluable UQ matched specimens insufficient")
  used <- summarize(tvd[joint])
  diffZ <- rz-az
  pairing <- data.table(
    sample_id_LOCAL_ONLY=gsm,
    branch=branch,
    n_technical_runs=as.integer(manifest$n_technical_runs),
    jointly_UQ_evaluable=uq_arch & uq_rec,
    assay_fraction_TVD=tvd,
    same_dominant_component=dom_arch==dom_rec,
    original_ARCHS4_dominant=dom_arch,
    projected_recount3_dominant=dom_rec,
    frozen_Z_MAE=colMeans(abs(diffZ)),
    original_ARCHS4_H_total=rowSums(ah),
    recount3_H_total=rowSums(rh),
    original_ARCHS4_normalized_KL=arch_adequacy$normalized_KL_per_input_mass,
    recount3_normalized_KL=adequacy$normalized_KL_per_input_mass,
    original_ARCHS4_reconstruction_RMSE=arch_adequacy$reconstruction_RMSE,
    recount3_reconstruction_RMSE=adequacy$reconstruction_RMSE,
    original_ARCHS4_reconstruction_R2=arch_adequacy$reconstruction_R2,
    recount3_reconstruction_R2=adequacy$reconstruction_R2,
    original_ARCHS4_score_entropy=arch_adequacy$normalized_score_entropy,
    recount3_score_entropy=adequacy$normalized_score_entropy
  )
  for (cmp in COMPONENTS) {
    pairing[[paste0("fraction_delta_",cmp)]] <- rf[,cmp]-af[,cmp]
    pairing[[paste0("raw_H_delta_",cmp)]] <- rh[,cmp]-ah[,cmp]
  }
  per_pair[[branch]] <- pairing
  anno <- c("FAM20C","FCAR","ISLR","MAL2","MATR3","PODNL1",
            "PSME2","SLC39A4","TNNT3")
  per_gene[[branch]] <- data.table(
    branch=branch,W_model_gene=W_genes,
    mean_recount3_minus_ARCHS4_frozen_Z=rowMeans(diffZ),
    median_abs_Z_error=apply(abs(diffZ),1,median),
    known_source_annotation_difference=W_genes %in% anno
  )
  aggregates[[branch]] <- data.table(
    branch=branch,
    n_paired_specimens=143L,
    n_jointly_evaluable=sum(joint),
    median_fraction_TVD=used[["median"]],
    p90_fraction_TVD=used[["p90"]],
    p95_fraction_TVD=used[["p95"]],
    max_fraction_TVD=used[["max"]],
    dominant_component_concordance=mean(dom_arch[joint]==dom_rec[joint]),
    median_frozen_Z_MAE=median(colMeans(abs(diffZ))[joint]),
    median_original_ARCHS4_norm_KL=median(arch_adequacy$normalized_KL_per_input_mass[joint]),
    median_recount3_norm_KL=median(adequacy$normalized_KL_per_input_mass[joint]),
    median_original_ARCHS4_RMSE=median(arch_adequacy$reconstruction_RMSE[joint]),
    median_recount3_RMSE=median(adequacy$reconstruction_RMSE[joint]),
    original_ARCHS4_Z_replay_max_abs=zerror,
    original_archive_full500_KL_replay_max_TVD=max(replay_tvd),
    original_ARCHS4_H_fraction_check_max_abs=ferror,
    recount3_directional_features=rec_dir$retained_feature_count
  )
  cat("V7P_BRANCH_SUMMARY",branch,
      "N",sum(joint),
      "MEDIAN_FRACTION_TVD",format(used[["median"]],digits=7),
      "P95_FRACTION_TVD",format(used[["p95"]],digits=7),
      "DOMINANT_CONCORDANCE",format(mean(dom_arch[joint]==dom_rec[joint]),digits=7),"\n")
  per_branch[[branch]] <- list(
    arch=list(Z=az,H=ah,f=af,adequacy=arch_adequacy),
    recount3=list(Z=rz,H=rh,f=rf,adequacy=adequacy)
  )
}
cat("V7P_ARCHS4_FROZEN_REFERENCE_ALL_THREE_BRANCHES_PASS\n")

# Within-assay normalization sensitivity using the ORIGINAL Stage10 helper,
# which computes TVD, dominant concordance and CLR-related effects.
sens_summary <- list()
for (assay in c("ARCHS4","recount3")) {
  f0 <- per_branch[["CPM"]][[assay]]$f
  for (alternative in c("TMM_CPM","UQ_CPM")) {
    mask <- if (alternative=="UQ_CPM") {
      if (assay=="ARCHS4") uq_arch else uq_rec
    } else rep(TRUE,143L)
    if (sum(mask)<2L) next
    tag <- paste0("CPM_vs_",alternative)
    response <- normalization_sensitivity(
      f0,per_branch[[alternative]][[assay]]$f,
      tag,evaluable_mask=mask
    )
    # Use original function's sample table as local-only evidence.
    fwrite(response$sample,file.path(FINAL,
      paste0("normalization_sensitivity_",assay,"_",alternative,"_LOCAL_ONLY.tsv")),sep="\t")
    s <- rowSums(abs(f0-per_branch[[alternative]][[assay]]$f))/2
    sens_summary[[paste(assay,alternative)]] <- data.table(
      assay=assay,comparison=tag,n_evaluable=sum(mask),
      median_fraction_TVD=median(s[mask]),
      p95_fraction_TVD=unname(quantile(s[mask],.95)),
      dominant_component_concordance=mean(
        dominant_component(f0)[mask]==dominant_component(per_branch[[alternative]][[assay]]$f)[mask]
      )
    )
  }
}
fwrite(rbindlist(sens_summary),
       file.path(FINAL,"aggregate_normalization_robustness.tsv"),sep="\t")
cat("V7P_ORIGINAL_NORMALIZATION_SENSITIVITY_COMPLETE\n")

# Explicit nine-gene counterfactual in frozen Z space (NOT a different W).
anno <- c("FAM20C","FCAR","ISLR","MAL2","MATR3","PODNL1",
          "PSME2","SLC39A4","TNNT3")
assert(sum(W_genes %in% anno)==9L,"Nine annotation-difference genes not found in frozen W")
base <- per_branch[["CPM"]]
hybrid <- base$recount3$Z
affected <- which(W_genes %in% anno)
hybrid[affected,] <- base$arch$Z[affected,,drop=FALSE]
hybrid_directional <- full500_bundle(hybrid,W)
assert(hybrid_directional$retained_feature_count==nrow(W),
       "Nine-gene counterfactual loses original 500 W features; HOLD sensitivity")
Hhybrid <- fixed500_project_raw_h(hybrid,W,"RECOUNT3_NINE_GENE_SWAP")
Fhybrid <- score_fractions(Hhybrid)
full_tvd <- rowSums(abs(base$arch$f-base$recount3$f))/2
swapped_tvd <- rowSums(abs(base$arch$f-Fhybrid))/2
fwrite(data.table(
  sample_id_LOCAL_ONLY=gsm,
  original_CPM_TVD=full_tvd,
  nine_recount3_Z_substituted_with_ARCHS4_TVD=swapped_tvd,
  delta_swapped_minus_original_TVD=swapped_tvd-full_tvd
),file.path(FINAL,"nine_gene_recount3_to_ARCHS4_Z_swap_LOCAL_ONLY.tsv"),sep="\t")
fwrite(data.table(
  sensitivity="recount3_nine_frozen_Z_gene_values_replaced_by_same_sample_ARCHS4",
  n_swapped_genes=9L,n_pairs=143L,
  original_median_TVD=median(full_tvd),
  swapped_median_TVD=median(swapped_tvd),
  median_per_pair_TVD_difference=median(swapped_tvd-full_tvd),
  p95_absolute_TVD_change=unname(quantile(abs(swapped_tvd-full_tvd),.95))
),file.path(FINAL,"aggregate_nine_gene_sensitivity.tsv"),sep="\t")
cat("V7P_NINE_ANNOTATION_DIFFERENCE_SWAP_DIAGNOSTIC_COMPLETE\n")

# Repeat/order sentinel: independent fixed500 KL estimator on 12 paired samples.
# Compare against numerical order/repeat limits originally adopted by Stage10.
r0 <- base$recount3$Z[,gsm[seq_len(12L)],drop=FALSE]
rr1 <- fixed500_project_raw_h(r0,W,"RECOUNT3_CPM_SENTINEL_REPEAT")
rr2 <- fixed500_project_raw_h(r0,W,"RECOUNT3_CPM_SENTINEL_REPEAT2")
repeat_max <- max(abs(score_fractions(rr1)-score_fractions(rr2)))
reversed <- r0[,rev(colnames(r0)),drop=FALSE]
revscore <- score_fractions(fixed500_project_raw_h(reversed,W,"RECOUNT3_CPM_SENTINEL_REVERSED"))
revscore <- revscore[colnames(r0),COMPONENTS,drop=FALSE]
ordered <- score_fractions(rr1)
order_tvd <- rowSums(abs(ordered-revscore))/2
order_median <- median(order_tvd)
order_max <- max(order_tvd)
fwrite(data.table(
   repeat_max_abs_fraction=repeat_max,
   reverse_order_median_TVD=order_median,
   reverse_order_max_TVD=order_max,
   reverse_order_dominant_concordance=mean(
     dominant_component(ordered)==dominant_component(revscore)
   )
),file.path(FINAL,"aggregate_determinism_sentinel.tsv"),sep="\t")
cat("V7P_RECOUNT3_ORDER_SENTINEL",
    "REPEAT_MAX",format(repeat_max,digits=8),
    "REVERSE_MEDIAN_TVD",format(order_median,digits=8),
    "REVERSE_MAX_TVD",format(order_max,digits=8),"\n")
assert(repeat_max<=1e-12 &&
       order_median<=1e-5 &&
       order_max<=1e-3,
       "Recount3 projection repeat/order sentinel outside original Stage10 tolerances")

fwrite(rbindlist(per_pair),file.path(FINAL,"143_paired_transfer_LOCAL_ONLY.tsv"),sep="\t")
fwrite(rbindlist(per_gene),file.path(FINAL,"499_gene_transfer_LOCAL_ONLY.tsv"),sep="\t")
fwrite(rbindlist(aggregates),file.path(FINAL,"aggregate_transfer_metrics.tsv"),sep="\t")
saveRDS(per_branch,file.path(FINAL,"paired_143_frozen_ARCHS4_and_recount3_projections_LOCAL_ONLY.rds"),
        compress="gzip")
cat("V7P_FROZEN_ARCHS4_VS_RECOUNT3_143_PAIRS_499_GENES_COMPLETED_PENDING_REVIEW\n")
cat("ARCHS4_SIDE_ORIGINAL_FROZEN_9465_ATLAS_NOT_RENORMALIZED\n")
cat("RECOUNT3_ALL_454_SOURCE_RUNS_AGGREGATED_WITHIN_143_SRX\n")
cat("COMPUTATIONAL_RESULTS_ARE_TECHNICAL_CALIBRATION_NOT_BIOLOGICAL_VALIDATION\n")
cat("NO_FROZEN_FILES_OR_W_MODIFIED\nCLINICAL_ADMISSIONS 0\n")
