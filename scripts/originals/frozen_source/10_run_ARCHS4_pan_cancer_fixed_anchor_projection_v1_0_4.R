#!/usr/bin/env Rscript

# ==============================================================================
# 10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R
# ==============================================================================
#
# Purpose
# -------
# Project the final locked ARCHS4 solid-tumour atlas onto the frozen
# SpatialEcoTyper v1.0.3 bulk SE basis in one common pan-cancer coordinate
# system.
#
# This stage:
#   1. uses the frozen Stage 09 v1.0.1 alias-resolved 499-gene matrices;
#   2. estimates fixed weighted gene means/SDs from the independent
#      8,962-library cross-sectional anchor; the UQ branch uses the 8,799
#      UQ-evaluable anchors only;
#   3. applies those frozen parameters to all 9,465 canonical libraries;
#   4. calls DeconvoluteSE with scale = FALSE and one prediction run;
#   5. generates CPM-primary, TMM-CPM and UQ-CPM score spaces separately;
#   6. retains raw H and derives sum-to-one compositions from it;
#   7. calculates primary fixed-W reconstruction adequacy;
#   8. quantifies normalization and sample-order sensitivity;
#   9. writes a strict QC contract.
#
# It does NOT:
#   - estimate cancer-specific scaling;
#   - test clinical metadata;
#   - discover residual factors;
#   - call NonSE an unexplained fraction.
#
# ==============================================================================
#
# This v1.0.4 hotfix restores the validated v1.0.2 determinism sentinel and
# QC-row helper that were accidentally omitted from v1.0.3. Projection
# mathematics, inputs, UQ evaluability handling, and output schema are unchanged.
#
options(
  stringsAsFactors = FALSE,
  width = 260,
  warn = 1
)

suppressPackageStartupMessages({
  library(data.table)
  library(SpatialEcoTyper)
  library(NMF)
})

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.", call. = FALSE)
}

# ==============================================================================
# 1. Constants and paths
# ==============================================================================

SCRIPT_VERSION <- "1.0.4"
EXPECTED_PACKAGE_VERSION <- "1.0.3"

PROJECT_HOME <- Sys.getenv(
  "ARCHS4_PROJECT_HOME",
  unset = "/home/tpn890/orestis/ST_ECOTYPES"
)

DATA_ROOT <- Sys.getenv(
  "ARCHS4_V26_ROOT",
  unset = "/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6"
)

INPUT_DIR <- Sys.getenv(
  "ARCHS4_STAGE10_INPUT_DIR",
  unset = file.path(
    DATA_ROOT,
    "10_preprojection_normalization_v1_0_1"
  )
)

OUTPUT_DIR <- Sys.getenv(
  "ARCHS4_STAGE10_OUTPUT_DIR",
  unset = file.path(
    DATA_ROOT,
    "11_pan_cancer_fixed_anchor_projection_v1_0_4"
  )
)

FROZEN_W_FILE <- Sys.getenv(
  "ARCHS4_FROZEN_W_FILE",
  unset = file.path(
    PROJECT_HOME,
    "02_model_assets",
    "official_v1.0.3",
    "Bulk_SE_Recovery_W.rds"
  )
)

INSTALLED_W_FILE <- file.path(
  system.file(
    "extdata",
    package = "SpatialEcoTyper"
  ),
  "Bulk_SE_Recovery_W.rds"
)

CANONICAL_MANIFEST_FILE <- file.path(
  INPUT_DIR,
  "01_canonical_library_manifest.tsv.gz"
)

ANCHOR_MANIFEST_FILE <- file.path(
  INPUT_DIR,
  "02_cross_sectional_anchor_manifest.tsv.gz"
)

ANCHOR_WEIGHT_FILE <- file.path(
  INPUT_DIR,
  "03_cross_sectional_balanced_anchor_weights.tsv"
)

MODEL_GENE_AUDIT_FILE <- file.path(
  INPUT_DIR,
  "08_model_gene_availability.tsv"
)

INPUT_STATUS_FILE <- file.path(
  INPUT_DIR,
  "18_preprojection_normalization_status.json"
)

SAMPLE_NORMALIZATION_QC_FILE <- file.path(
  INPUT_DIR,
  "04_sample_normalization_QC.tsv.gz"
)

BRANCHES <- list(
  CPM = list(
    canonical = file.path(INPUT_DIR, "10_model_gene_CPM_canonical.tsv.gz"),
    anchor = file.path(INPUT_DIR, "13_model_gene_CPM_cross_sectional_anchor.tsv.gz")
  ),
  TMM_CPM = list(
    canonical = file.path(INPUT_DIR, "11_model_gene_TMM_CPM_canonical.tsv.gz"),
    anchor = file.path(INPUT_DIR, "14_model_gene_TMM_CPM_cross_sectional_anchor.tsv.gz")
  ),
  UQ_CPM = list(
    canonical = file.path(INPUT_DIR, "12_model_gene_UQ_CPM_canonical.tsv.gz"),
    anchor = file.path(INPUT_DIR, "15_model_gene_UQ_CPM_cross_sectional_anchor.tsv.gz")
  )
)

COMPONENTS <- c(
  "NonSE",
  paste0("SE", 1:9)
)

EXPECTED_CANONICAL <- 9465L
EXPECTED_ANCHOR <- 8962L
EXPECTED_MODEL_GENES <- 499L
EXPECTED_DIRECTIONAL_FEATURES <- 500L
EXPECTED_CANCERS <- 18L
EXPECTED_UQ_EVALUABLE <- 9302L
EXPECTED_UQ_NOT_EVALUABLE <- 163L
EXPECTED_UQ_ANCHOR <- 8799L

PRIMARY_BRANCH <- "CPM"
NSAMPLE_PER_RUN <- 10000L
NCORES <- 1L
SENTINEL_N <- 512L

SD_FLOOR <- 1e-12
MATRIX_TOLERANCE <- 1e-8
SUM_TO_ONE_TOLERANCE <- 1e-8
REPEAT_TOLERANCE <- 1e-12
ORDER_MEDIAN_TVD_TOLERANCE <- 1e-5
ORDER_MAX_TVD_TOLERANCE <- 1e-3
CLR_FLOOR <- 1e-8
RECONSTRUCTION_EPS <- 1e-12

EARLY_COMMAND_ARGS <- commandArgs(trailingOnly = TRUE)
SELF_TEST_MODE <- "--run-self-tests" %in% EARLY_COMMAND_ARGS

if (!SELF_TEST_MODE) {
  if (
    dir.exists(OUTPUT_DIR) &&
      length(list.files(OUTPUT_DIR, all.files = FALSE, no.. = TRUE)) > 0L
  ) {
    stop(
      "Refusing to overwrite non-empty immutable output directory: ",
      OUTPUT_DIR,
      call. = FALSE
    )
  }

  dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 2. General helpers
# ==============================================================================

write_tsv <- function(x, path, compress = FALSE) {
  dt <- as.data.table(x)
  fwrite(
    dt,
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = if (compress) "gzip" else "none"
  )
  message("Wrote: ", path)
}

write_matrix_rows <- function(matrix_object, row_id_name, path, compress = TRUE) {
  dt <- as.data.table(
    matrix_object,
    keep.rownames = row_id_name
  )
  write_tsv(
    dt,
    path,
    compress = compress
  )
}

read_expression_matrix <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing or empty expression matrix: ", path, call. = FALSE)
  }

  dt <- fread(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    check.names = FALSE,
    showProgress = FALSE
  )

  if (ncol(dt) < 2L) {
    stop("Expression matrix has fewer than two columns: ", path, call. = FALSE)
  }

  genes <- trimws(as.character(dt[[1L]]))

  if (anyNA(genes) || any(genes == "") || anyDuplicated(genes)) {
    stop("Invalid or duplicated gene identifiers in: ", path, call. = FALSE)
  }

  matrix_object <- as.matrix(
    dt[, -1L, with = FALSE]
  )
  storage.mode(matrix_object) <- "double"
  rownames(matrix_object) <- genes

  if (
    any(!is.finite(matrix_object)) ||
      any(matrix_object < 0)
  ) {
    stop(
      "Expression matrix contains negative or non-finite values: ",
      path,
      call. = FALSE
    )
  }

  matrix_object
}

validate_id_order <- function(observed, expected, label) {
  if (!identical(observed, expected)) {
    if (!setequal(observed, expected)) {
      missing_ids <- setdiff(expected, observed)
      extra_ids <- setdiff(observed, expected)
      stop(
        label,
        " ID set mismatch. Missing: ",
        paste(head(missing_ids, 20L), collapse = ", "),
        "; extra: ",
        paste(head(extra_ids, 20L), collapse = ", "),
        call. = FALSE
      )
    }
    stop(label, " IDs are not in the locked manifest order.", call. = FALSE)
  }
}

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

build_pan_cancer_balanced_weights <- function(
    anchor_frame,
    weight_column = "pan_cancer_study_balanced_weight"
) {
  required <- c(
    "sample_id",
    "primary_gse_id",
    "analysis_cancer_class"
  )

  missing <- setdiff(required, names(anchor_frame))
  if (length(missing) > 0L) {
    stop(
      "Cannot build pan-cancer anchor weights; missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  weights <- copy(
    as.data.table(anchor_frame)[
      ,
      .(sample_id, primary_gse_id, analysis_cancer_class)
    ]
  )

  if (
    nrow(weights) < 2L ||
      anyDuplicated(weights$sample_id) ||
      anyNA(weights$sample_id) ||
      any(weights$sample_id == "")
  ) {
    stop("Invalid anchor frame supplied for weight construction.", call. = FALSE)
  }

  n_cancers <- uniqueN(weights$analysis_cancer_class)

  weights[
    ,
    n_studies_in_cancer := uniqueN(primary_gse_id),
    by = analysis_cancer_class
  ]

  weights[
    ,
    n_samples_in_cancer_study := .N,
    by = .(analysis_cancer_class, primary_gse_id)
  ]

  weights[
    ,
    (weight_column) :=
      1 /
        n_cancers /
        n_studies_in_cancer /
        n_samples_in_cancer_study
  ]

  total <- sum(weights[[weight_column]])

  if (!is.finite(total) || abs(total - 1) > 1e-12) {
    stop(
      "Constructed pan-cancer anchor weights sum to ",
      format(total, digits = 16),
      ", not one.",
      call. = FALSE
    )
  }

  setorder(
    weights,
    analysis_cancer_class,
    primary_gse_id,
    sample_id
  )

  weights
}

weighted_anchor_scaling <- function(
    canonical_matrix,
    anchor_matrix,
    anchor_ids,
    weights,
    branch_name
) {
  if (!identical(rownames(canonical_matrix), rownames(anchor_matrix))) {
    stop(
      branch_name,
      ": canonical and anchor gene orders differ.",
      call. = FALSE
    )
  }

  validate_id_order(
    colnames(anchor_matrix),
    anchor_ids,
    paste0(branch_name, " anchor matrix")
  )

  canonical_anchor <- canonical_matrix[
    ,
    anchor_ids,
    drop = FALSE
  ]

  anchor_copy_difference <- max(
    abs(
      canonical_anchor -
        anchor_matrix
    )
  )

  if (
    !is.finite(anchor_copy_difference) ||
      anchor_copy_difference > MATRIX_TOLERANCE
  ) {
    stop(
      branch_name,
      ": anchor values differ from the corresponding canonical subset; max |difference| = ",
      format(anchor_copy_difference, scientific = TRUE),
      call. = FALSE
    )
  }

  if (
    length(weights) != ncol(anchor_matrix) ||
      any(!is.finite(weights)) ||
      any(weights <= 0)
  ) {
    stop(branch_name, ": invalid anchor weights.", call. = FALSE)
  }

  weights <- weights / sum(weights)
  effective_n <- 1 / sum(weights^2)
  unbiased_denominator <- 1 - sum(weights^2)

  if (
    !is.finite(unbiased_denominator) ||
      unbiased_denominator <= 0
  ) {
    stop(
      branch_name,
      ": weighted-SD denominator is non-positive.",
      call. = FALSE
    )
  }

  log_canonical <- log2(canonical_matrix + 1)
  log_anchor <- log2(anchor_matrix + 1)

  weighted_mean <- as.vector(
    log_anchor %*% weights
  )

  centered_anchor <- sweep(
    log_anchor,
    1L,
    weighted_mean,
    "-"
  )

  weighted_variance <- rowSums(
    sweep(
      centered_anchor^2,
      2L,
      weights,
      "*"
    )
  ) / unbiased_denominator

  weighted_sd <- sqrt(
    pmax(
      weighted_variance,
      0
    )
  )

  zero_sd <- !is.finite(weighted_sd) |
    weighted_sd <= SD_FLOOR

  if (any(zero_sd)) {
    stop(
      branch_name,
      ": ",
      sum(zero_sd),
      " model genes have zero or invalid weighted anchor SD: ",
      paste(
        head(
          rownames(anchor_matrix)[zero_sd],
          30L
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  scaled_canonical <- sweep(
    sweep(
      log_canonical,
      1L,
      weighted_mean,
      "-"
    ),
    1L,
    weighted_sd,
    "/"
  )

  scaled_anchor <- scaled_canonical[
    ,
    anchor_ids,
    drop = FALSE
  ]

  scaled_weighted_mean <- as.vector(
    scaled_anchor %*% weights
  )

  scaled_centered <- sweep(
    scaled_anchor,
    1L,
    scaled_weighted_mean,
    "-"
  )

  scaled_weighted_sd <- sqrt(
    rowSums(
      sweep(
        scaled_centered^2,
        2L,
        weights,
        "*"
      )
    ) / unbiased_denominator
  )

  if (
    any(!is.finite(scaled_canonical)) ||
      all(scaled_canonical >= 0)
  ) {
    stop(
      branch_name,
      ": fixed-anchor scaled matrix is invalid or lacks negative values.",
      call. = FALSE
    )
  }

  parameter_table <- data.table(
    model_gene = rownames(canonical_matrix),
    branch = branch_name,
    weighted_anchor_mean_log2 = weighted_mean,
    weighted_anchor_sd_log2 = weighted_sd,
    unweighted_anchor_mean_log2 = rowMeans(log_anchor),
    unweighted_anchor_sd_log2 = apply(log_anchor, 1L, stats::sd),
    weighted_scaled_anchor_mean = scaled_weighted_mean,
    weighted_scaled_anchor_sd = scaled_weighted_sd,
    canonical_scaled_min = apply(scaled_canonical, 1L, min),
    canonical_scaled_max = apply(scaled_canonical, 1L, max),
    effective_anchor_n = effective_n,
    anchor_copy_max_absolute_difference = anchor_copy_difference
  )

  list(
    scaled = scaled_canonical,
    parameters = parameter_table,
    anchor_copy_difference = anchor_copy_difference,
    effective_n = effective_n,
    max_abs_weighted_scaled_mean = max(abs(scaled_weighted_mean)),
    max_abs_weighted_scaled_sd_error = max(abs(scaled_weighted_sd - 1))
  )
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


component_feature_coverage <- function(
    W,
    retained_features,
    branch_name
) {
  retained <- W[
    rownames(W) %in% retained_features,
    ,
    drop = FALSE
  ]

  rbindlist(
    lapply(
      COMPONENTS,
      function(component_name) {
        total_weight <- sum(
          W[, component_name],
          na.rm = TRUE
        )

        retained_weight <- sum(
          retained[, component_name],
          na.rm = TRUE
        )

        data.table(
          branch = branch_name,
          component = component_name,
          total_directional_features = nrow(W),
          retained_directional_features = nrow(retained),
          component_nonzero_features_total = sum(
            W[, component_name] > 0,
            na.rm = TRUE
          ),
          component_nonzero_features_retained = sum(
            retained[, component_name] > 0,
            na.rm = TRUE
          ),
          total_component_weight = total_weight,
          retained_component_weight = retained_weight,
          retained_weight_fraction = if (
            total_weight > 0
          ) {
            retained_weight / total_weight
          } else {
            NA_real_
          }
        )
      }
    )
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

top_three_table <- function(
    fractions
) {
  score_dt <- data.table(
    sample_id = rownames(fractions),
    as.data.table(fractions)
  )

  long <- melt(
    score_dt,
    id.vars = "sample_id",
    variable.name = "component",
    value.name = "projected_fraction"
  )

  setorder(
    long,
    sample_id,
    -projected_fraction,
    component
  )

  long[
    ,
    rank := seq_len(.N),
    by = sample_id
  ][
    rank <= 3L
  ]
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

run_sentinel_validation <- function(
    scaled_matrix,
    W
) {
  sentinel_index <- unique(
    as.integer(
      round(
        seq(
          1,
          ncol(scaled_matrix),
          length.out = min(
            SENTINEL_N,
            ncol(scaled_matrix)
          )
        )
      )
    )
  )

  sentinel_ids <- colnames(scaled_matrix)[
    sentinel_index
  ]

  sentinel <- scaled_matrix[
    ,
    sentinel_ids,
    drop = FALSE
  ]

  directional_normal <- make_directional_bundle(
    sentinel,
    W
  )

  if (
    directional_normal$retained_feature_count !=
      nrow(W)
  ) {
    stop(
      "Sentinel subset retained ",
      directional_normal$retained_feature_count,
      " rather than ",
      nrow(W),
      " directional features.",
      call. = FALSE
    )
  }

  raw_1 <- project_raw_h(
    sentinel,
    W,
    "sentinel_repeat_1"
  )

  raw_2 <- project_raw_h(
    sentinel,
    W,
    "sentinel_repeat_2"
  )

  fraction_1 <- score_fractions(raw_1)
  fraction_2 <- score_fractions(raw_2)

  repeat_max_absolute_difference <- max(
    abs(
      fraction_1 -
        fraction_2
    )
  )

  reversed_sentinel <- sentinel[
    ,
    rev(sentinel_ids),
    drop = FALSE
  ]

  directional_reversed <- make_directional_bundle(
    reversed_sentinel,
    W
  )

  raw_reversed <- project_raw_h(
    reversed_sentinel,
    W,
    "sentinel_reversed"
  )

  fraction_reversed <- score_fractions(
    raw_reversed
  )

  fraction_reversed <- fraction_reversed[
    sentinel_ids,
    COMPONENTS,
    drop = FALSE
  ]

  absolute_difference <- abs(
    fraction_1 -
      fraction_reversed
  )

  order_tvd <- 0.5 *
    rowSums(
      absolute_difference
    )

  dominant_concordance <- mean(
    dominant_component(fraction_1) ==
      dominant_component(fraction_reversed)
  )

  summary_table <- data.table(
    check = c(
      "same_order_repeat",
      "reversed_sample_order"
    ),
    n_sentinel_samples = length(sentinel_ids),
    retained_features_normal = directional_normal$retained_feature_count,
    retained_features_reversed = directional_reversed$retained_feature_count,
    maximum_absolute_fraction_difference = c(
      repeat_max_absolute_difference,
      max(absolute_difference)
    ),
    median_total_variation_distance = c(
      0,
      stats::median(order_tvd)
    ),
    maximum_total_variation_distance = c(
      0,
      max(order_tvd)
    ),
    dominant_component_concordance = c(
      1,
      dominant_concordance
    )
  )

  list(
    summary = summary_table,
    sentinel_ids = sentinel_ids,
    repeat_max_absolute_difference = repeat_max_absolute_difference,
    order_median_tvd = stats::median(order_tvd),
    order_max_tvd = max(order_tvd),
    order_dominant_concordance = dominant_concordance
  )
}

qc_row <- function(
    check,
    passed,
    detail
) {
  data.table(
    check = check,
    status = if (isTRUE(passed)) {
      "PASS"
    } else {
      "FAIL"
    },
    detail = as.character(detail)
  )
}

# ==============================================================================
# 3. Lightweight self-tests
# ==============================================================================

run_self_tests <- function() {
  stopifnot(
    exists("run_sentinel_validation", mode = "function"),
    exists("qc_row", mode = "function")
  )

  qc_helper_test <- qc_row(
    check = "selftest",
    passed = TRUE,
    detail = "helper_present"
  )

  stopifnot(
    nrow(qc_helper_test) == 1L,
    identical(qc_helper_test$status, "PASS")
  )

  harmonization_test_matrix <- matrix(
    c(
      -1,
      1,
      -2,
      2,
      -3,
      3
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(
      c(
        "CXORF36",
        "TP53",
        "C17ORF107"
      ),
      c(
        "S1",
        "S2"
      )
    )
  )

  harmonization_test_W <- matrix(
    1,
    nrow = 3L,
    ncol = 1L,
    dimnames = list(
      c(
        "DIPK2B__pos",
        "TP53__pos",
        "C17orf107__pos"
      ),
      "NonSE"
    )
  )

  harmonization_test_audit <- data.table(
    model_gene = c(
      "CXORF36",
      "TP53",
      "C17ORF107"
    ),
    ARCHS4_source_gene = c(
      "DIPK2B",
      "TP53",
      "C17ORF107"
    ),
    resolution = c(
      "alias_resolved",
      "exact",
      "exact"
    )
  )

  harmonization_test <- harmonize_expression_to_W(
    scaled_matrix = harmonization_test_matrix,
    W = harmonization_test_W,
    model_gene_audit = harmonization_test_audit,
    branch_name = "SELFTEST"
  )

  stopifnot(
    identical(
      rownames(
        harmonization_test$matrix
      ),
      c(
        "DIPK2B",
        "TP53",
        "C17orf107"
      )
    ),
    harmonization_test$alias_harmonized_count ==
      1L,
    harmonization_test$case_harmonized_count ==
      1L,
    harmonization_test$matrix[
      "DIPK2B",
      "S1"
    ] ==
      -1,
    harmonization_test$matrix[
      "C17orf107",
      "S1"
    ] ==
      -3
  )

  canonical <- matrix(
    c(
      1, 2, 4, 8,
      2, 5, 7, 11
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      c("G1", "G2"),
      c("S1", "S2", "S3", "S4")
    )
  )

  anchor <- canonical[
    ,
    c("S1", "S2", "S3"),
    drop = FALSE
  ]

  scaling <- weighted_anchor_scaling(
    canonical_matrix = canonical,
    anchor_matrix = anchor,
    anchor_ids = colnames(anchor),
    weights = c(0.2, 0.3, 0.5),
    branch_name = "SELFTEST"
  )

  stopifnot(
    scaling$max_abs_weighted_scaled_mean < 1e-12,
    scaling$max_abs_weighted_scaled_sd_error < 1e-12,
    identical(
      dim(scaling$scaled),
      dim(canonical)
    )
  )

  raw_h <- matrix(
    c(
      2, 1,
      1, 3
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      c("S1", "S2"),
      c("NonSE", "SE1")
    )
  )

  fractions <- score_fractions(raw_h)

  stopifnot(
    max(
      abs(
        rowSums(fractions) -
          1
      )
    ) <
      1e-12
  )

  identical_sensitivity <- suppressWarnings(
    normalization_sensitivity(
      primary_scores = cbind(
      fractions,
      SE2 = 0,
      SE3 = 0,
      SE4 = 0,
      SE5 = 0,
      SE6 = 0,
      SE7 = 0,
      SE8 = 0,
      SE9 = 0
    )[
      ,
      COMPONENTS
    ],
    alternative_scores = cbind(
      fractions,
      SE2 = 0,
      SE3 = 0,
      SE4 = 0,
      SE5 = 0,
      SE6 = 0,
      SE7 = 0,
      SE8 = 0,
      SE9 = 0
    )[
      ,
      COMPONENTS
    ],
      comparison_name = "SELFTEST"
    )
  )

  stopifnot(
    max(
      identical_sensitivity$sample$total_variation_distance
    ) ==
      0,
    max(
      identical_sensitivity$sample$aitchison_distance
    ) ==
      0
  )

  cat(
    "ALL PAN-CANCER FIXED-ANCHOR PROJECTION V1.0.4 SELF-TESTS PASSED\n"
  )
}

command_args <- commandArgs(
  trailingOnly = TRUE
)

if ("--run-self-tests" %in% command_args) {
  run_self_tests()
  quit(
    save = "no",
    status = 0
  )
}

# ==============================================================================
# 4. Validate software, frozen model and Stage 08 lock
# ==============================================================================

required_files <- c(
  CANONICAL_MANIFEST_FILE,
  ANCHOR_MANIFEST_FILE,
  ANCHOR_WEIGHT_FILE,
  MODEL_GENE_AUDIT_FILE,
  INPUT_STATUS_FILE,
  SAMPLE_NORMALIZATION_QC_FILE,
  FROZEN_W_FILE,
  INSTALLED_W_FILE,
  unlist(
    lapply(
      BRANCHES,
      unlist
    )
  )
)

missing_files <- required_files[
  !file.exists(required_files) |
    file.info(required_files)$size <= 0
]

if (length(missing_files) > 0L) {
  stop(
    "Missing or empty required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

package_version <- as.character(
  packageVersion("SpatialEcoTyper")
)

if (!identical(package_version, EXPECTED_PACKAGE_VERSION)) {
  stop(
    "Expected SpatialEcoTyper ",
    EXPECTED_PACKAGE_VERSION,
    " but found ",
    package_version,
    ".",
    call. = FALSE
  )
}

W <- readRDS(
  FROZEN_W_FILE
)

installed_W <- readRDS(
  INSTALLED_W_FILE
)

if (
  !is.matrix(W) ||
    !identical(
      dim(W),
      c(
        EXPECTED_DIRECTIONAL_FEATURES,
        length(COMPONENTS)
      )
    ) ||
    !identical(
      colnames(W),
      COMPONENTS
    )
) {
  stop("Frozen Bulk_SE_Recovery_W.rds failed structural validation.", call. = FALSE)
}

if (!identical(W, installed_W)) {
  stop(
    "Frozen W is not identical to the model installed in SpatialEcoTyper 1.0.3.",
    call. = FALSE
  )
}

frozen_w_md5 <- unname(
  tools::md5sum(
    FROZEN_W_FILE
  )
)

installed_w_md5 <- unname(
  tools::md5sum(
    INSTALLED_W_FILE
  )
)

input_status <- jsonlite::read_json(
  INPUT_STATUS_FILE,
  simplifyVector = TRUE
)

if (
  !identical(input_status$status, "PASS") ||
    input_status$qc_failures != 0 ||
    !isTRUE(input_status$expression_matrix_exported) ||
    !identical(input_status$projection_run, FALSE) ||
    input_status$canonical_samples != EXPECTED_CANONICAL ||
    input_status$cross_sectional_samples != EXPECTED_ANCHOR ||
    input_status$cancer_classes != EXPECTED_CANCERS ||
    input_status$model_genes_available != EXPECTED_MODEL_GENES ||
    input_status$model_genes_unresolved != 0
) {
  stop(
    "Stage 09 v1.0.1 status is not the required frozen 9,465/8,962, 499/499 PASS state.",
    call. = FALSE
  )
}

canonical_manifest <- fread(
  CANONICAL_MANIFEST_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

anchor_manifest <- fread(
  ANCHOR_MANIFEST_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

anchor_weights <- fread(
  ANCHOR_WEIGHT_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

model_gene_audit <- fread(
  MODEL_GENE_AUDIT_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

sample_normalization_qc <- fread(
  SAMPLE_NORMALIZATION_QC_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

for (
  required_column in c(
    "sample_id",
    "primary_gse_id",
    "provisional_cohort",
    "analysis_cancer_class"
  )
) {
  if (!required_column %in% names(canonical_manifest)) {
    stop(
      "Canonical manifest lacks ",
      required_column,
      ".",
      call. = FALSE
    )
  }
}

if (
  nrow(canonical_manifest) != EXPECTED_CANONICAL ||
    anyDuplicated(canonical_manifest$sample_id)
) {
  stop(
    "Expected ",
    EXPECTED_CANONICAL,
    " unique canonical samples.",
    call. = FALSE
  )
}

if (
  nrow(anchor_manifest) != EXPECTED_ANCHOR ||
    anyDuplicated(anchor_manifest$sample_id)
) {
  stop(
    "Expected ",
    EXPECTED_ANCHOR,
    " unique anchor samples.",
    call. = FALSE
  )
}

if (
  !setequal(
    anchor_manifest$sample_id,
    anchor_weights$sample_id
  ) ||
    anyDuplicated(anchor_weights$sample_id)
) {
  stop(
    "Anchor weights do not match the locked anchor manifest.",
    call. = FALSE
  )
}

anchor_weights <- anchor_weights[
  match(
    anchor_manifest$sample_id,
    sample_id
  )
]

weight_vector <- anchor_weights$pan_cancer_study_balanced_weight

if (
  any(!is.finite(weight_vector)) ||
    any(weight_vector <= 0) ||
    abs(sum(weight_vector) - 1) > 1e-12
) {
  stop("Pan-cancer anchor weights are invalid.", call. = FALSE)
}

if (
  nrow(model_gene_audit) != EXPECTED_MODEL_GENES ||
    sum(model_gene_audit$resolution == "alias_resolved") != 7L ||
    any(model_gene_audit$resolution == "unresolved")
) {
  stop(
    "Model-gene audit is not the expected 492 exact + 7 alias-resolved state.",
    call. = FALSE
  )
}

canonical_ids <- as.character(canonical_manifest$sample_id)
anchor_ids <- as.character(anchor_manifest$sample_id)

if (
  uniqueN(canonical_manifest$analysis_cancer_class) != EXPECTED_CANCERS ||
    !all(anchor_ids %in% canonical_ids)
) {
  stop(
    "Canonical/anchor cancer-class or subset contract failed.",
    call. = FALSE
  )
}

required_sample_qc_columns <- c(
  "sample_id",
  "upper_quartile_fraction",
  "tmm_factor",
  "uq_factor",
  "model_genes_detected_recomputed",
  "stage08_tier",
  "stage07_expression_quality_tier"
)

missing_sample_qc_columns <- setdiff(
  required_sample_qc_columns,
  names(sample_normalization_qc)
)

if (length(missing_sample_qc_columns) > 0L) {
  stop(
    "Stage 09 sample-normalization QC lacks columns: ",
    paste(missing_sample_qc_columns, collapse = ", "),
    call. = FALSE
  )
}

if (
  nrow(sample_normalization_qc) != EXPECTED_CANONICAL ||
    anyDuplicated(sample_normalization_qc$sample_id) ||
    !identical(
      as.character(sample_normalization_qc$sample_id),
      canonical_ids
    )
) {
  stop(
    "Stage 09 sample-normalization QC is not in exact canonical order.",
    call. = FALSE
  )
}

uq_evaluable_mask <-
  is.finite(sample_normalization_qc$upper_quartile_fraction) &
    sample_normalization_qc$upper_quartile_fraction > 0

if (
  sum(uq_evaluable_mask) != EXPECTED_UQ_EVALUABLE ||
    sum(!uq_evaluable_mask) != EXPECTED_UQ_NOT_EVALUABLE
) {
  stop(
    "Unexpected UQ evaluability partition: evaluable=",
    sum(uq_evaluable_mask),
    "; non-evaluable=",
    sum(!uq_evaluable_mask),
    call. = FALSE
  )
}

uq_evaluable_ids <- canonical_ids[uq_evaluable_mask]
uq_anchor_manifest <- anchor_manifest[sample_id %in% uq_evaluable_ids]

if (nrow(uq_anchor_manifest) != EXPECTED_UQ_ANCHOR) {
  stop(
    "Expected ",
    EXPECTED_UQ_ANCHOR,
    " UQ-evaluable cross-sectional anchors; observed ",
    nrow(uq_anchor_manifest),
    ".",
    call. = FALSE
  )
}

uq_anchor_weights <- build_pan_cancer_balanced_weights(
  uq_anchor_manifest
)

setkey(uq_anchor_weights, sample_id)
uq_anchor_weights <- uq_anchor_weights[
  data.table(sample_id = as.character(uq_anchor_manifest$sample_id)),
  on = "sample_id"
]

if (anyNA(uq_anchor_weights$pan_cancer_study_balanced_weight)) {
  stop("Failed to align UQ-evaluable anchor weights.", call. = FALSE)
}

branch_anchor_ids <- list(
  CPM = anchor_ids,
  TMM_CPM = anchor_ids,
  UQ_CPM = as.character(uq_anchor_manifest$sample_id)
)

branch_weight_vectors <- list(
  CPM = weight_vector,
  TMM_CPM = weight_vector,
  UQ_CPM = uq_anchor_weights$pan_cancer_study_balanced_weight
)

projection_evaluability <- copy(
  sample_normalization_qc[
    ,
    .(
      sample_id,
      primary_gse_id,
      provisional_cohort,
      analysis_cancer_class,
      stage08_tier,
      stage07_expression_quality_tier,
      model_genes_detected_recomputed,
      upper_quartile_fraction,
      tmm_factor,
      uq_factor
    )
  ]
)

projection_evaluability[
  ,
  `:=`(
    CPM_projection_evaluable = TRUE,
    TMM_projection_evaluable = TRUE,
    UQ_projection_evaluable = uq_evaluable_mask,
    UQ_non_evaluable_reason = fifelse(
      uq_evaluable_mask,
      NA_character_,
      "upper_quartile_fraction_nonpositive_stage09_fallback"
    )
  )
]

input_provenance <- data.table(
  role = c(
    "stage09_status",
    "canonical_manifest",
    "anchor_manifest",
    "anchor_weights",
    "sample_normalization_qc",
    "model_gene_audit",
    "frozen_W",
    "installed_W"
  ),
  path = c(
    INPUT_STATUS_FILE,
    CANONICAL_MANIFEST_FILE,
    ANCHOR_MANIFEST_FILE,
    ANCHOR_WEIGHT_FILE,
    SAMPLE_NORMALIZATION_QC_FILE,
    MODEL_GENE_AUDIT_FILE,
    FROZEN_W_FILE,
    INSTALLED_W_FILE
  )
)

input_provenance[
  ,
  `:=`(
    size_bytes = file.info(path)$size,
    md5 = unname(tools::md5sum(path))
  )
]

write_tsv(
  input_provenance,
  file.path(OUTPUT_DIR, "00_input_provenance.tsv")
)

write_tsv(
  projection_evaluability,
  file.path(
    OUTPUT_DIR,
    "02A_projection_normalization_evaluability.tsv.gz"
  ),
  compress = TRUE
)

write_tsv(
  uq_anchor_weights,
  file.path(
    OUTPUT_DIR,
    "02B_UQ_evaluable_cross_sectional_balanced_anchor_weights.tsv.gz"
  ),
  compress = TRUE
)

# ==============================================================================
# 5. Load, scale and project each normalization branch
# ==============================================================================

message(
  "Loading, fixed-anchor scaling and projecting CPM/TMM/UQ branches."
)

branch_results <- list()
scaling_summary_rows <- list()
coverage_rows <- list()
W_harmonization_rows <- list()

for (branch_name in names(BRANCHES)) {
  message(
    "--------------------------------------------------------------------------"
  )
  message(
    "BRANCH: ",
    branch_name
  )

  canonical_matrix <- read_expression_matrix(
    BRANCHES[[branch_name]]$canonical
  )

  anchor_matrix_all <- read_expression_matrix(
    BRANCHES[[branch_name]]$anchor
  )

  current_anchor_ids <- branch_anchor_ids[[branch_name]]
  current_weight_vector <- branch_weight_vectors[[branch_name]]

  anchor_matrix <- anchor_matrix_all[
    ,
    current_anchor_ids,
    drop = FALSE
  ]

  if (
    nrow(canonical_matrix) != EXPECTED_MODEL_GENES ||
      nrow(anchor_matrix) != EXPECTED_MODEL_GENES
  ) {
    stop(
      branch_name,
      ": expected ",
      EXPECTED_MODEL_GENES,
      " model genes.",
      call. = FALSE
    )
  }

  validate_id_order(
    colnames(canonical_matrix),
    canonical_ids,
    paste0(branch_name, " canonical matrix")
  )

  scaling <- weighted_anchor_scaling(
    canonical_matrix = canonical_matrix,
    anchor_matrix = anchor_matrix,
    anchor_ids = current_anchor_ids,
    weights = current_weight_vector,
    branch_name = branch_name
  )

  harmonized <- harmonize_expression_to_W(
    scaled_matrix = scaling$scaled,
    W = W,
    model_gene_audit = model_gene_audit,
    branch_name = branch_name
  )

  projection_scaled <- harmonized$matrix

  W_harmonization_rows[[branch_name]] <- harmonized$mapping

  directional <- make_directional_bundle(
    projection_scaled,
    W
  )

  if (
    length(
      directional$missing_model_genes
    ) >
      0L ||
      directional$supported_gene_count !=
        directional$expected_gene_count ||
      directional$retained_feature_count !=
        nrow(W)
  ) {
    stop(
      branch_name,
      ": frozen-W feature recovery failed after symbol harmonization.\n",
      "  W base genes expected: ",
      directional$expected_gene_count,
      "\n",
      "  W base genes supported: ",
      directional$supported_gene_count,
      "\n",
      "  W directional features expected: ",
      nrow(W),
      "\n",
      "  W directional features retained: ",
      directional$retained_feature_count,
      "\n",
      "  Missing W base genes: ",
      if (
        length(
          directional$missing_model_genes
        ) ==
          0L
      ) {
        "none"
      } else {
        paste(
          directional$missing_model_genes,
          collapse = ", "
        )
      },
      "\n",
      "  Missing W directional features: ",
      if (
        length(
          directional$missing_W_features
        ) ==
          0L
      ) {
        "none"
      } else {
        paste(
          directional$missing_W_features,
          collapse = ", "
        )
      },
      "\n",
      "  Zero-variance W directional features: ",
      if (
        length(
          directional$zero_variance_W_features
        ) ==
          0L
      ) {
        "none"
      } else {
        paste(
          directional$zero_variance_W_features,
          collapse = ", "
        )
      },
      call. = FALSE
    )
  }

  raw_h <- project_raw_h(
    projection_scaled,
    W,
    branch_name
  )

  fractions <- score_fractions(
    raw_h
  )

  if (
    max(
      abs(
        rowSums(fractions) -
          1
      )
    ) >
      SUM_TO_ONE_TOLERANCE
  ) {
    stop(
      branch_name,
      ": projected fractions do not sum to one.",
      call. = FALSE
    )
  }

  scaling_file <- file.path(
    OUTPUT_DIR,
    paste0(
      switch(
        branch_name,
        CPM = "03",
        TMM_CPM = "04",
        UQ_CPM = "05"
      ),
      "_pan_cancer_anchor_scaling_",
      branch_name,
      ".tsv.gz"
    )
  )

  scaled_rds_file <- file.path(
    OUTPUT_DIR,
    paste0(
      switch(
        branch_name,
        CPM = "06",
        TMM_CPM = "07",
        UQ_CPM = "08"
      ),
      "_fixed_anchor_scaled_input_",
      branch_name,
      ".rds"
    )
  )

  raw_h_file <- file.path(
    OUTPUT_DIR,
    paste0(
      switch(
        branch_name,
        CPM = "09",
        TMM_CPM = "11",
        UQ_CPM = "13"
      ),
      "_raw_H_",
      branch_name,
      ".tsv.gz"
    )
  )

  fraction_file <- file.path(
    OUTPUT_DIR,
    paste0(
      switch(
        branch_name,
        CPM = "10",
        TMM_CPM = "12",
        UQ_CPM = "14"
      ),
      "_fraction_scores_",
      branch_name,
      ".tsv.gz"
    )
  )

  write_tsv(
    scaling$parameters,
    scaling_file,
    compress = TRUE
  )

  saveRDS(
    projection_scaled,
    scaled_rds_file,
    compress = "gzip"
  )
  message(
    "Wrote: ",
    scaled_rds_file
  )

  write_matrix_rows(
    raw_h,
    "sample_id",
    raw_h_file,
    compress = TRUE
  )

  write_matrix_rows(
    fractions,
    "sample_id",
    fraction_file,
    compress = TRUE
  )

  coverage <- component_feature_coverage(
    W,
    directional$retained_features,
    branch_name
  )

  coverage_rows[[branch_name]] <- coverage

  scaling_summary_rows[[branch_name]] <- data.table(
    branch = branch_name,
    canonical_samples = ncol(canonical_matrix),
    anchor_samples = ncol(anchor_matrix),
    anchor_weight_scheme = if (branch_name == "UQ_CPM") {
      "UQ_evaluable_cancer_study_balanced"
    } else {
      "all_anchor_cancer_study_balanced"
    },
    input_model_genes = nrow(canonical_matrix),
    frozen_W_model_genes = directional$expected_gene_count,
    exact_symbol_rows_to_W = harmonized$exact_count,
    case_harmonized_rows_to_W = harmonized$case_harmonized_count,
    alias_harmonized_rows_to_W = harmonized$alias_harmonized_count,
    retained_directional_features = directional$retained_feature_count,
    effective_anchor_n = scaling$effective_n,
    anchor_copy_max_absolute_difference = scaling$anchor_copy_difference,
    maximum_absolute_weighted_scaled_anchor_mean = scaling$max_abs_weighted_scaled_mean,
    maximum_absolute_weighted_scaled_anchor_SD_error = scaling$max_abs_weighted_scaled_sd_error,
    scaled_global_minimum = min(scaling$scaled),
    scaled_global_maximum = max(scaling$scaled),
    maximum_fraction_row_sum_error = max(
      abs(
        rowSums(fractions) -
          1
      )
    )
  )

  branch_results[[branch_name]] <- list(
    scaled = projection_scaled,
    raw_h = raw_h,
    fractions = fractions,
    directional = directional
  )

  rm(
    canonical_matrix,
    anchor_matrix,
    anchor_matrix_all
  )
  invisible(
    gc()
  )
}

scaling_summary <- rbindlist(
  scaling_summary_rows,
  use.names = TRUE
)

feature_coverage <- rbindlist(
  coverage_rows,
  use.names = TRUE
)

write_tsv(
  canonical_manifest,
  file.path(
    OUTPUT_DIR,
    "01_canonical_library_manifest.tsv.gz"
  ),
  compress = TRUE
)

write_tsv(
  anchor_weights,
  file.path(
    OUTPUT_DIR,
    "02_cross_sectional_balanced_anchor_weights.tsv.gz"
  ),
  compress = TRUE
)

write_tsv(
  scaling_summary,
  file.path(
    OUTPUT_DIR,
    "15_scaling_and_projection_summary.tsv"
  )
)

write_tsv(
  feature_coverage,
  file.path(
    OUTPUT_DIR,
    "16_component_directional_feature_coverage.tsv"
  )
)


W_harmonization_audit <- rbindlist(
  W_harmonization_rows,
  use.names = TRUE
)

write_tsv(
  W_harmonization_audit,
  file.path(
    OUTPUT_DIR,
    "16A_frozen_W_gene_symbol_harmonization.tsv"
  )
)

# ==============================================================================
# 6. Primary CPM metadata, composition and model adequacy outputs
# ==============================================================================

primary <- branch_results[[PRIMARY_BRANCH]]

primary_scores <- data.table(
  sample_id = rownames(primary$fractions),
  as.data.table(primary$fractions)
)

canonical_ordered <- canonical_manifest[
  match(
    primary_scores$sample_id,
    sample_id
  )
]

primary_scores_with_metadata <- cbind(
  canonical_ordered,
  primary_scores[
    ,
    -1L,
    with = FALSE
  ]
)

write_tsv(
  primary_scores_with_metadata,
  file.path(
    OUTPUT_DIR,
    "17_primary_CPM_fraction_scores_with_metadata.tsv.gz"
  ),
  compress = TRUE
)

top_three <- top_three_table(
  primary$fractions
)

write_tsv(
  top_three,
  file.path(
    OUTPUT_DIR,
    "18_primary_CPM_top3_components_per_sample.tsv.gz"
  ),
  compress = TRUE
)

message(
  "Calculating fixed-W reconstruction adequacy for the primary CPM branch."
)

primary_V <- primary$directional$V
primary_W <- primary$directional$W
primary_V_hat <- primary_W %*%
  t(
    primary$raw_h
  )

rownames(primary_V_hat) <- rownames(primary_W)
colnames(primary_V_hat) <- rownames(primary$raw_h)

primary_adequacy <- calculate_reconstruction_adequacy(
  V = primary_V,
  V_hat = primary_V_hat,
  fractions = primary$fractions
)

primary_adequacy <- merge(
  primary_adequacy,
  canonical_manifest[
    ,
    .(
      sample_id,
      primary_gse_id,
      provisional_cohort,
      analysis_cancer_class
    )
  ],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

primary_adequacy <- primary_adequacy[
  match(
    canonical_ids,
    sample_id
  )
]

write_tsv(
  primary_adequacy,
  file.path(
    OUTPUT_DIR,
    "19_primary_CPM_fixed_W_model_adequacy.tsv.gz"
  ),
  compress = TRUE
)

feature_metadata <- data.table(
  directional_feature = rownames(primary_W),
  model_gene = sub(
    "__.*$",
    "",
    rownames(primary_W)
  ),
  direction = sub(
    "^.*__",
    "",
    rownames(primary_W)
  ),
  row_index = seq_len(nrow(primary_W))
)

write_tsv(
  feature_metadata,
  file.path(
    OUTPUT_DIR,
    "20_primary_CPM_directional_feature_metadata.tsv"
  )
)

primary_reconstruction <- list(
  description = paste(
    "Primary CPM fixed-anchor directional feature matrices.",
    "Observed and reconstructed matrices are nonnegative SpatialEcoTyper",
    "pos/neg feature-space values; residual = observed - reconstructed."
  ),
  directional_feature_metadata = feature_metadata,
  sample_ids = colnames(primary_V),
  observed_directional_features = primary_V,
  reconstructed_directional_features = primary_V_hat,
  signed_directional_residual = primary_V -
    primary_V_hat
)

saveRDS(
  primary_reconstruction,
  file.path(
    OUTPUT_DIR,
    "21_primary_CPM_directional_reconstruction_and_residual.rds"
  ),
  compress = "gzip"
)
message(
  "Wrote: ",
  file.path(
    OUTPUT_DIR,
    "21_primary_CPM_directional_reconstruction_and_residual.rds"
  )
)

# ==============================================================================
# 7. Normalization sensitivity
# ==============================================================================

message(
  "Calculating CPM-versus-TMM/UQ compositional sensitivity."
)

uq_non_evaluable_reason <- ifelse(
  uq_evaluable_mask,
  NA_character_,
  "upper_quartile_fraction_nonpositive_stage09_fallback"
)

sensitivity_objects <- list(
  CPM_vs_TMM_CPM = normalization_sensitivity(
    primary$fractions,
    branch_results$TMM_CPM$fractions,
    "CPM_vs_TMM_CPM"
  ),
  CPM_vs_UQ_CPM = normalization_sensitivity(
    primary$fractions,
    branch_results$UQ_CPM$fractions,
    "CPM_vs_UQ_CPM",
    evaluable_mask = uq_evaluable_mask,
    non_evaluable_reason = uq_non_evaluable_reason
  )
)

sensitivity_by_sample <- rbindlist(
  lapply(
    sensitivity_objects,
    `[[`,
    "sample"
  ),
  use.names = TRUE
)

sensitivity_by_sample <- merge(
  sensitivity_by_sample,
  projection_evaluability,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

sensitivity_by_sample[
  ,
  `:=`(
    canonical_order_index = match(sample_id, canonical_ids),
    comparison_order_index = match(
      comparison,
      c("CPM_vs_TMM_CPM", "CPM_vs_UQ_CPM")
    )
  )
]

setorder(
  sensitivity_by_sample,
  canonical_order_index,
  comparison_order_index
)

sensitivity_by_sample[
  ,
  c("canonical_order_index", "comparison_order_index") := NULL
]

sensitivity_by_component <- rbindlist(
  lapply(
    sensitivity_objects,
    `[[`,
    "component"
  ),
  use.names = TRUE
)

sensitivity_summary <- sensitivity_by_sample[
  ,
  .(
    n_samples_total = .N,
    n_samples_evaluable = sum(comparison_evaluable),
    n_samples_not_evaluable = sum(!comparison_evaluable),
    median_total_variation_distance = stats::median(
      total_variation_distance,
      na.rm = TRUE
    ),
    q95_total_variation_distance = stats::quantile(
      total_variation_distance,
      0.95,
      names = FALSE,
      na.rm = TRUE
    ),
    maximum_total_variation_distance = max(
      total_variation_distance,
      na.rm = TRUE
    ),
    median_aitchison_distance = stats::median(
      aitchison_distance,
      na.rm = TRUE
    ),
    q95_aitchison_distance = stats::quantile(
      aitchison_distance,
      0.95,
      names = FALSE,
      na.rm = TRUE
    ),
    maximum_aitchison_distance = max(
      aitchison_distance,
      na.rm = TRUE
    ),
    dominant_component_concordance_fraction = mean(
      dominant_component_concordant,
      na.rm = TRUE
    ),
    diagnostic_median_TVD_all_samples = stats::median(
      diagnostic_total_variation_distance
    ),
    diagnostic_maximum_TVD_all_samples = max(
      diagnostic_total_variation_distance
    )
  ),
  by = comparison
]

write_tsv(
  sensitivity_by_sample,
  file.path(
    OUTPUT_DIR,
    "22_normalization_sensitivity_by_sample.tsv.gz"
  ),
  compress = TRUE
)

write_tsv(
  sensitivity_by_component,
  file.path(
    OUTPUT_DIR,
    "23_normalization_sensitivity_by_component.tsv"
  )
)

write_tsv(
  sensitivity_summary,
  file.path(
    OUTPUT_DIR,
    "24_normalization_sensitivity_summary.tsv"
  )
)

# ==============================================================================
# 8. Determinism and sample-order sentinel
# ==============================================================================

message(
  "Running primary CPM determinism and sample-order sentinel."
)

sentinel <- run_sentinel_validation(
  primary$scaled,
  W
)

write_tsv(
  sentinel$summary,
  file.path(
    OUTPUT_DIR,
    "25_primary_CPM_determinism_and_order_sentinel.tsv"
  )
)

write_tsv(
  data.table(
    sentinel_order = seq_along(
      sentinel$sentinel_ids
    ),
    sample_id = sentinel$sentinel_ids
  ),
  file.path(
    OUTPUT_DIR,
    "26_primary_CPM_sentinel_sample_ledger.tsv"
  )
)

# ==============================================================================
# 9. QC contract, status and report
# ==============================================================================

adequacy_core_columns <- c(
  "fixed_W_KL_divergence",
  "normalized_KL_per_input_mass",
  "reconstruction_RMSE",
  "reconstruction_MAE",
  "reconstruction_Pearson",
  "reconstruction_cosine",
  "reconstruction_R2",
  "normalized_score_entropy"
)

qc <- rbindlist(
  list(
    qc_row(
      "SpatialEcoTyper_version",
      identical(
        package_version,
        EXPECTED_PACKAGE_VERSION
      ),
      paste0(
        "observed=",
        package_version,
        ";expected=",
        EXPECTED_PACKAGE_VERSION
      )
    ),
    qc_row(
      "frozen_W_identical_to_installed_W",
      identical(
        W,
        installed_W
      ),
      paste0(
        "frozen_md5=",
        frozen_w_md5,
        ";installed_md5=",
        installed_w_md5
      )
    ),
    qc_row(
      "frozen_W_dimensions",
      identical(
        dim(W),
        c(
          EXPECTED_DIRECTIONAL_FEATURES,
          length(COMPONENTS)
        )
      ),
      paste(
        dim(W),
        collapse = "x"
      )
    ),
    qc_row(
      "stage09_frozen_preprojection_status",
      identical(
        input_status$status,
        "PASS"
      ) &&
        input_status$qc_failures == 0 &&
        isTRUE(input_status$expression_matrix_exported) &&
        input_status$canonical_samples == EXPECTED_CANONICAL &&
        input_status$cross_sectional_samples == EXPECTED_ANCHOR &&
        input_status$cancer_classes == EXPECTED_CANCERS &&
        input_status$model_genes_available == EXPECTED_MODEL_GENES &&
        input_status$model_genes_unresolved == 0,
      paste0(
        "status=",
        input_status$status,
        ";available=",
        input_status$model_genes_available,
        ";unresolved=",
        input_status$model_genes_unresolved
      )
    ),
    qc_row(
      "canonical_sample_count",
      nrow(canonical_manifest) ==
        EXPECTED_CANONICAL,
      paste0(
        "observed=",
        nrow(canonical_manifest),
        ";expected=",
        EXPECTED_CANONICAL
      )
    ),
    qc_row(
      "anchor_sample_count",
      nrow(anchor_manifest) ==
        EXPECTED_ANCHOR,
      paste0(
        "observed=",
        nrow(anchor_manifest),
        ";expected=",
        EXPECTED_ANCHOR
      )
    ),
    qc_row(
      "anchor_weights_sum_one",
      abs(
        sum(weight_vector) -
          1
      ) <=
        1e-12,
      paste0(
        "sum=",
        format(
          sum(weight_vector),
          digits = 16
        )
      )
    ),
    qc_row(
      "single_internal_prediction_run",
      EXPECTED_CANONICAL <=
        NSAMPLE_PER_RUN,
      paste0(
        "samples=",
        EXPECTED_CANONICAL,
        ";nsample_per_run=",
        NSAMPLE_PER_RUN
      )
    ),
    qc_row(
      "all_branches_have_499_stage08_input_rows",
      all(
        scaling_summary$input_model_genes ==
          EXPECTED_MODEL_GENES
      ),
      paste(
        scaling_summary$input_model_genes,
        collapse = ","
      )
    ),
    qc_row(
      "all_W_symbol_harmonizations_accounted_for",
      all(
        scaling_summary$exact_symbol_rows_to_W +
          scaling_summary$case_harmonized_rows_to_W +
          scaling_summary$alias_harmonized_rows_to_W ==
          scaling_summary$frozen_W_model_genes
      ),
      paste0(
        "exact=",
        paste(
          scaling_summary$exact_symbol_rows_to_W,
          collapse = ","
        ),
        ";case=",
        paste(
          scaling_summary$case_harmonized_rows_to_W,
          collapse = ","
        ),
        ";alias=",
        paste(
          scaling_summary$alias_harmonized_rows_to_W,
          collapse = ","
        )
      )
    ),
    qc_row(
      "expected_case_harmonized_W_gene_count",
      all(
        scaling_summary$case_harmonized_rows_to_W ==
          6L
      ),
      paste(
        scaling_summary$case_harmonized_rows_to_W,
        collapse = ","
      )
    ),
    qc_row(
      "all_frozen_W_genes_resolved",
      all(
        scaling_summary$frozen_W_model_genes ==
          length(
            unique(
              sub(
                "__.*$",
                "",
                rownames(W)
              )
            )
          )
      ),
      paste(
        scaling_summary$frozen_W_model_genes,
        collapse = ","
      )
    ),
    qc_row(
      "all_branches_retain_all_frozen_W_directional_features",
      all(
        scaling_summary$retained_directional_features ==
          nrow(W)
      ),
      paste(
        scaling_summary$retained_directional_features,
        collapse = ","
      )
    ),
    qc_row(
      "anchor_copy_matches_canonical_subset",
      all(
        scaling_summary$anchor_copy_max_absolute_difference <=
          MATRIX_TOLERANCE
      ),
      paste0(
        "max=",
        format(
          max(
            scaling_summary$anchor_copy_max_absolute_difference
          ),
          scientific = TRUE
        )
      )
    ),
    qc_row(
      "weighted_scaled_anchor_means_zero",
      all(
        scaling_summary$maximum_absolute_weighted_scaled_anchor_mean <=
          1e-10
      ),
      paste0(
        "max=",
        format(
          max(
            scaling_summary$maximum_absolute_weighted_scaled_anchor_mean
          ),
          scientific = TRUE
        )
      )
    ),
    qc_row(
      "weighted_scaled_anchor_SDs_one",
      all(
        scaling_summary$maximum_absolute_weighted_scaled_anchor_SD_error <=
          1e-10
      ),
      paste0(
        "max_error=",
        format(
          max(
            scaling_summary$maximum_absolute_weighted_scaled_anchor_SD_error
          ),
          scientific = TRUE
        )
      )
    ),
    qc_row(
      "all_fraction_scores_sum_to_one",
      all(
        scaling_summary$maximum_fraction_row_sum_error <=
          SUM_TO_ONE_TOLERANCE
      ),
      paste0(
        "max_error=",
        format(
          max(
            scaling_summary$maximum_fraction_row_sum_error
          ),
          scientific = TRUE
        )
      )
    ),
    qc_row(
      "primary_adequacy_metrics_finite",
      all(
        is.finite(
          as.matrix(
            primary_adequacy[
              ,
              ..adequacy_core_columns
            ]
          )
        )
      ),
      paste0(
        "rows=",
        nrow(primary_adequacy)
      )
    ),
    qc_row(
      "branch_anchor_sample_counts",
      identical(
        as.integer(scaling_summary$anchor_samples),
        c(EXPECTED_ANCHOR, EXPECTED_ANCHOR, EXPECTED_UQ_ANCHOR)
      ),
      paste0(
        paste(scaling_summary$branch, scaling_summary$anchor_samples, sep = "="),
        collapse = ";"
      )
    ),
    qc_row(
      "UQ_evaluability_partition",
      sum(uq_evaluable_mask) == EXPECTED_UQ_EVALUABLE &&
        sum(!uq_evaluable_mask) == EXPECTED_UQ_NOT_EVALUABLE,
      paste0(
        "evaluable=",
        sum(uq_evaluable_mask),
        ";non_evaluable=",
        sum(!uq_evaluable_mask)
      )
    ),
    qc_row(
      "normalization_sensitivity_complete",
      nrow(sensitivity_by_sample) == 2L * EXPECTED_CANONICAL &&
        nrow(sensitivity_by_sample[
          comparison == "CPM_vs_TMM_CPM" & comparison_evaluable == TRUE
        ]) == EXPECTED_CANONICAL &&
        nrow(sensitivity_by_sample[
          comparison == "CPM_vs_UQ_CPM" & comparison_evaluable == TRUE
        ]) == EXPECTED_UQ_EVALUABLE &&
        nrow(sensitivity_by_sample[
          comparison == "CPM_vs_UQ_CPM" & comparison_evaluable == FALSE
        ]) == EXPECTED_UQ_NOT_EVALUABLE &&
        all(is.finite(
          sensitivity_by_sample$diagnostic_total_variation_distance
        )) &&
        all(is.finite(
          sensitivity_by_sample$diagnostic_aitchison_distance
        )) &&
        all(is.finite(
          sensitivity_by_sample[comparison_evaluable == TRUE]$total_variation_distance
        )) &&
        all(is.finite(
          sensitivity_by_sample[comparison_evaluable == TRUE]$aitchison_distance
        )) &&
        all(is.na(
          sensitivity_by_sample[comparison_evaluable == FALSE]$total_variation_distance
        )) &&
        all(is.na(
          sensitivity_by_sample[comparison_evaluable == FALSE]$aitchison_distance
        )) &&
        all(is.finite(
          as.matrix(
            sensitivity_by_component[
              ,
              .(
                fraction_Pearson,
                fraction_Spearman,
                CLR_Pearson,
                median_absolute_fraction_difference,
                maximum_absolute_fraction_difference
              )
            ]
          )
        )),
      paste0(
        "rows=",
        nrow(sensitivity_by_sample),
        ";TMM_evaluable=",
        sensitivity_by_sample[comparison == "CPM_vs_TMM_CPM", sum(comparison_evaluable)],
        ";UQ_evaluable=",
        sensitivity_by_sample[comparison == "CPM_vs_UQ_CPM", sum(comparison_evaluable)]
      )
    ),
    qc_row(
      "same_order_repeat_deterministic",
      sentinel$repeat_max_absolute_difference <=
        REPEAT_TOLERANCE,
      paste0(
        "max_abs_difference=",
        format(
          sentinel$repeat_max_absolute_difference,
          scientific = TRUE
        ),
        ";tolerance=",
        REPEAT_TOLERANCE
      )
    ),
    qc_row(
      "reversed_order_median_TVD_stable",
      sentinel$order_median_tvd <=
        ORDER_MEDIAN_TVD_TOLERANCE,
      paste0(
        "median_TVD=",
        format(
          sentinel$order_median_tvd,
          scientific = TRUE
        ),
        ";tolerance=",
        ORDER_MEDIAN_TVD_TOLERANCE
      )
    ),
    qc_row(
      "reversed_order_maximum_TVD_stable",
      sentinel$order_max_tvd <=
        ORDER_MAX_TVD_TOLERANCE,
      paste0(
        "maximum_TVD=",
        format(
          sentinel$order_max_tvd,
          scientific = TRUE
        ),
        ";tolerance=",
        ORDER_MAX_TVD_TOLERANCE
      )
    ),
    qc_row(
      "reversed_order_dominant_component_concordance",
      sentinel$order_dominant_concordance ==
        1,
      paste0(
        "concordance=",
        sentinel$order_dominant_concordance
      )
    )
  ),
  use.names = TRUE
)

write_tsv(
  qc,
  file.path(
    OUTPUT_DIR,
    "27_pan_cancer_projection_QC_contract.tsv"
  )
)

failed_qc <- qc[
  status ==
    "FAIL"
]

status <- list(
  script_version = SCRIPT_VERSION,
  generated_at_utc = format(
    Sys.time(),
    tz = "UTC",
    usetz = TRUE
  ),
  status = if (
    nrow(failed_qc) ==
      0L
  ) {
    "PASS"
  } else {
    "FAIL"
  },
  SpatialEcoTyper_version = package_version,
  frozen_W_md5 = frozen_w_md5,
  canonical_samples = nrow(canonical_manifest),
  cross_sectional_anchor_samples = nrow(anchor_manifest),
  stage09_input_model_genes = EXPECTED_MODEL_GENES,
  source_stage09_dir = INPUT_DIR,
  source_stage09_script_version = input_status$script_version,
  cancer_classes = EXPECTED_CANCERS,
  UQ_evaluable_samples = EXPECTED_UQ_EVALUABLE,
  UQ_non_evaluable_samples = EXPECTED_UQ_NOT_EVALUABLE,
  UQ_scaling_anchor_samples = EXPECTED_UQ_ANCHOR,
  frozen_W_model_genes = length(
    unique(
      sub(
        "__.*$",
        "",
        rownames(W)
      )
    )
  ),
  directional_features = nrow(W),
  projection_branches = names(BRANCHES),
  primary_branch = PRIMARY_BRANCH,
  nsample_per_run = NSAMPLE_PER_RUN,
  effective_anchor_n = scaling_summary$effective_anchor_n[
    scaling_summary$branch ==
      PRIMARY_BRANCH
  ],
  sentinel_repeat_max_absolute_difference = sentinel$repeat_max_absolute_difference,
  sentinel_reversed_order_median_TVD = sentinel$order_median_tvd,
  sentinel_reversed_order_maximum_TVD = sentinel$order_max_tvd,
  sentinel_reversed_order_dominant_concordance = sentinel$order_dominant_concordance,
  qc_failures = nrow(failed_qc)
)

jsonlite::write_json(
  status,
  path = file.path(
    OUTPUT_DIR,
    "28_pan_cancer_projection_status.json"
  ),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
message(
  "Wrote: ",
  file.path(
    OUTPUT_DIR,
    "28_pan_cancer_projection_status.json"
  )
)

report_lines <- c(
  "ARCHS4 PAN-CANCER FIXED-ANCHOR SPATIAL ECOTYPER PROJECTION V1.0.4",
  paste(
    rep(
      "=",
      88
    ),
    collapse = ""
  ),
  "",
  paste0(
    "Status: ",
    status$status
  ),
  paste0(
    "SpatialEcoTyper version: ",
    package_version
  ),
  paste0(
    "Frozen W: ",
    nrow(W),
    " directional features x ",
    ncol(W),
    " components"
  ),
  paste0(
    "Canonical libraries: ",
    nrow(canonical_manifest)
  ),
  paste0(
    "Cross-sectional anchor libraries: ",
    nrow(anchor_manifest)
  ),
  paste0(
    "Effective weighted anchor N: ",
    format(
      status$effective_anchor_n,
      digits = 6
    )
  ),
  paste0(
    "Prediction batching: one run; nsample.per.run=",
    NSAMPLE_PER_RUN
  ),
  "",
  "COORDINATE SYSTEM",
  "-----------------",
  "For each normalization branch, expression was transformed as log2(abundance + 1).",
  "Gene means and weighted unbiased SDs were estimated only from the locked",
  "cross-sectional anchor using the cancer/study-balanced weights. The same",
  "frozen parameters were then applied to all canonical samples.",
  "",
  "SpatialEcoTyper was called with:",
  "  scale = FALSE",
  "  sum2one = FALSE",
  paste0(
    "  nsample.per.run = ",
    NSAMPLE_PER_RUN
  ),
  paste0(
    "  ncores = ",
    NCORES
  ),
  "  frozen official Bulk_SE_Recovery_W.rds",
  "",
  "The raw fixed-W coefficient matrix H was retained. Sum-to-one NonSE/SE1-SE9",
  "compositions were derived directly from raw H.",
  "",
  "PRIMARY AND SENSITIVITY BRANCHES",
  "--------------------------------",
  "Primary: CPM",
  "Sensitivity: TMM-CPM and upper-quartile CPM",
  paste0("UQ robustness evaluable: ", EXPECTED_UQ_EVALUABLE, " samples"),
  paste0("UQ fallback/non-evaluable: ", EXPECTED_UQ_NOT_EVALUABLE, " samples"),
  paste0("UQ scaling anchor excludes fallback profiles: ", EXPECTED_UQ_ANCHOR, " anchors"),
  "The three branches use separate anchor scaling parameters and must remain",
  "separate score spaces. CPM and TMM use all 8,962 anchors. UQ uses only",
  "the 8,799 anchors with a positive Stage 09 upper-quartile fraction.",
  "UQ projections are retained for all samples, but the 163 fallback profiles",
  "are explicitly non-evaluable for CPM-versus-UQ robustness claims.",
  "",
  "SCALING AND PROJECTION SUMMARY",
  "------------------------------",
  paste(
    capture.output(
      print(
        scaling_summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "NORMALIZATION SENSITIVITY SUMMARY",
  "---------------------------------",
  paste(
    capture.output(
      print(
        sensitivity_summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "PRIMARY MODEL ADEQUACY",
  "----------------------",
  "Adequacy is measured in the nonnegative directional feature space used by",
  "the frozen model. NonSE is a published model component and is not used as",
  "an unexplained fraction or goodness-of-fit measure.",
  "",
  paste(
    capture.output(
      print(
        primary_adequacy[
          ,
          lapply(
            .SD,
            function(x) {
              c(
                minimum = min(x),
                median = stats::median(x),
                maximum = max(x)
              )
            }
          ),
          .SDcols = c(
            "normalized_KL_per_input_mass",
            "reconstruction_RMSE",
            "reconstruction_Pearson",
            "reconstruction_cosine",
            "reconstruction_R2",
            "normalized_score_entropy",
            "dominant_fraction"
          )
        ],
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "DETERMINISM AND ORDER SENTINEL",
  "------------------------------",
  paste(
    capture.output(
      print(
        sentinel$summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "BOUNDARY OF THIS STAGE",
  "----------------------",
  "This stage establishes the shared pan-cancer score space and its technical",
  "adequacy. Cancer-specific score spaces, residual-factor discovery, clinical",
  "covariate eligibility and Layer 2/Layer 3 inference are separate later stages.",
  "",
  "QC CONTRACT",
  "-----------",
  paste(
    capture.output(
      print(
        qc,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  )
)

report_file <- file.path(
  OUTPUT_DIR,
  "ARCHS4_pan_cancer_fixed_anchor_projection_v1_report.txt"
)

writeLines(
  report_lines,
  con = report_file
)
message(
  "Wrote: ",
  report_file
)

if (nrow(failed_qc) > 0L) {
  print(
    failed_qc
  )
  stop(
    "Pan-cancer projection QC contract failed.",
    call. = FALSE
  )
}

cat(
  "\n",
  paste(
    rep(
      "=",
      88
    ),
    collapse = ""
  ),
  "\n",
  "ARCHS4 PAN-CANCER FIXED-ANCHOR PROJECTION V1.0.4: PASS\n",
  paste(
    rep(
      "=",
      88
    ),
    collapse = ""
  ),
  "\n",
  "Canonical libraries:       ",
  nrow(canonical_manifest),
  "\n",
  "Cross-sectional anchor:    ",
  nrow(anchor_manifest),
  "\n",
  "Stage 09 input rows:       ",
  EXPECTED_MODEL_GENES,
  "\n",
  "Frozen W model genes:      ",
  length(
    unique(
      sub(
        "__.*$",
        "",
        rownames(W)
      )
    )
  ),
  "\n",
  "Directional features:      ",
  nrow(W),
  "\n",
  "Primary score space:       CPM fixed-anchor\n",
  "Sensitivity score spaces: TMM-CPM, UQ-CPM fixed-anchor\n",
  "Output:                    ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
