#!/usr/bin/env Rscript

# ==============================================================================
# 12_run_ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3.R
# ==============================================================================
#
# Purpose
# -------
# Project the locked ARCHS4 solid-tumour atlas onto the frozen
# SpatialEcoTyper v1.0.3 bulk SE basis using a separate, study-balanced
# fixed-anchor coordinate system for each cancer class.
#
# Why this stage is required
# --------------------------
# The published SEs and their bulk-recovery basis are pan-cancer constructs,
# but the published bulk workflow standardized expression within cancer type.
# The preceding shared pan-cancer stage intentionally used one common ruler
# across cancers. This stage provides the complementary within-cancer ruler.
#
# For each cancer and normalization branch:
#   1. derive study-balanced weights from the locked cross-sectional samples;
#   2. estimate weighted log2(abundance + 1) gene means and unbiased SDs;
#   3. apply those frozen parameters to every canonical sample in that cancer;
#   4. harmonize Stage 09 gene names to the exact frozen-W namespace;
#   5. project with fixed W, scale=FALSE and sum2one=FALSE;
#   6. retain raw H, sum-to-one compositions, fixed-W adequacy and residuals.
#
# The script also compares:
#   - CPM versus TMM-CPM and UQ-CPM within cancer-specific score spaces;
#   - shared pan-cancer versus cancer-specific CPM coordinates;
#   - pan-cancer versus cancer-specific fixed-W adequacy and residual structure.
#
# Important
# ---------
# Cancer-specific fractions are comparable among samples within the same cancer.
# They are not placed on one common cross-cancer coordinate system.
#
# NonSE is a published component and is never treated as the residual.
#
# ==============================================================================

options(
  stringsAsFactors = FALSE,
  width = 260,
  warn = 1
)

suppressPackageStartupMessages({
  library(data.table)
  library(SpatialEcoTyper)
  library(NMF)
  library(ggplot2)
  library(scales)
})

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.", call. = FALSE)
}

# ==============================================================================
# 1. Constants and paths
# ==============================================================================

SCRIPT_VERSION <- "1.0.3"
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
  "ARCHS4_STAGE12_INPUT_DIR",
  unset = file.path(
    DATA_ROOT,
    "10_preprojection_normalization_v1_0_1"
  )
)

PAN_INPUT_DIR <- Sys.getenv(
  "ARCHS4_STAGE12_PAN_INPUT_DIR",
  unset = file.path(
    DATA_ROOT,
    "11_pan_cancer_fixed_anchor_projection_v1_0_4"
  )
)

LAYER1_INPUT_DIR <- Sys.getenv(
  "ARCHS4_STAGE12_LAYER1_INPUT_DIR",
  unset = file.path(
    DATA_ROOT,
    "12_layer1_adequacy_robustness_atlas_v1_0_2_1"
  )
)

OUTPUT_DIR <- Sys.getenv(
  "ARCHS4_STAGE12_OUTPUT_DIR",
  unset = file.path(
    DATA_ROOT,
    "13_cancer_specific_fixed_anchor_projection_v1_0_3"
  )
)

FIGURE_DIR <- file.path(
  OUTPUT_DIR,
  "figures"
)

SOURCE_DATA_DIR <- file.path(
  FIGURE_DIR,
  "source_data"
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

MODEL_GENE_AUDIT_FILE <- file.path(
  INPUT_DIR,
  "08_model_gene_availability.tsv"
)

INPUT_STATUS_FILE <- file.path(
  INPUT_DIR,
  "18_preprojection_normalization_status.json"
)

PAN_STATUS_FILE <- file.path(
  PAN_INPUT_DIR,
  "28_pan_cancer_projection_status.json"
)

LAYER1_STATUS_FILE <- file.path(
  LAYER1_INPUT_DIR,
  "13_layer1_adequacy_robustness_status.json"
)

SAMPLE_NORMALIZATION_QC_FILE <- file.path(
  INPUT_DIR,
  "04_sample_normalization_QC.tsv.gz"
)

PAN_QC_FILE <- file.path(
  PAN_INPUT_DIR,
  "27_pan_cancer_projection_QC_contract.tsv"
)

LAYER1_SAMPLE_ATLAS_FILE <- file.path(
  LAYER1_INPUT_DIR,
  "01_sample_level_adequacy_robustness_atlas.tsv.gz"
)

LAYER1_QC_FILE <- file.path(
  LAYER1_INPUT_DIR,
  "12_layer1_adequacy_robustness_QC_contract.tsv"
)

PAN_SCORE_FILE <- file.path(
  PAN_INPUT_DIR,
  "17_primary_CPM_fraction_scores_with_metadata.tsv.gz"
)

PAN_ADEQUACY_FILE <- file.path(
  PAN_INPUT_DIR,
  "19_primary_CPM_fixed_W_model_adequacy.tsv.gz"
)

PAN_RECONSTRUCTION_FILE <- file.path(
  PAN_INPUT_DIR,
  "21_primary_CPM_directional_reconstruction_and_residual.rds"
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
  paste0(
    "SE",
    1:9
  )
)

EXPECTED_CANONICAL <- 9465L
EXPECTED_ANCHOR <- 8962L
EXPECTED_MODEL_GENES <- 499L
EXPECTED_DIRECTIONAL_FEATURES <- 500L
EXPECTED_CANCERS <- 18L
EXPECTED_UQ_EVALUABLE <- 9302L
EXPECTED_UQ_NOT_EVALUABLE <- 163L
EXPECTED_UQ_ANCHOR <- 8799L
EXPECTED_STAGE09_VERSION <- "1.0.1"
EXPECTED_PAN_VERSION <- "1.0.4"
EXPECTED_LAYER1_VERSION <- "1.0.2.1"

PRIMARY_BRANCH <- "CPM"
NSAMPLE_PER_RUN <- 10000L
NCORES <- 1L

SD_FLOOR <- 1e-12
MATRIX_TOLERANCE <- 1e-8
SUM_TO_ONE_TOLERANCE <- 1e-8
REPEAT_TOLERANCE <- 1e-12
ORDER_MEDIAN_TVD_TOLERANCE <- 1e-5
ORDER_MAX_TVD_TOLERANCE <- 1e-3
CLR_FLOOR <- 1e-8
RECONSTRUCTION_EPS <- 1e-12

MIN_ANCHOR_SAMPLES_PER_CANCER <- 20L
MIN_RETAINED_DIRECTIONAL_FEATURES <- 450L
MIN_COMPONENT_RETAINED_WEIGHT_FRACTION <- 0.90

PERSISTENT_RESIDUAL_ABS_MEAN_THRESHOLD <- 0.25
PERSISTENT_RESIDUAL_SIGN_CONCORDANCE <- 0.75
PAN_RESIDUAL_MIN_ABS_FOR_ATTENUATION <- 0.10

DPI <- 600
NATURE_FONT <- "Helvetica"

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

  dir.create(
    FIGURE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    SOURCE_DATA_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ==============================================================================
# 2. General helpers
# ==============================================================================

write_tsv <- function(
    x,
    path,
    compress = FALSE
) {
  fwrite(
    as.data.table(x),
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = if (compress) {
      "gzip"
    } else {
      "none"
    }
  )

  message(
    "Wrote: ",
    path
  )
}

write_matrix_rows <- function(
    matrix_object,
    row_id_name,
    path,
    compress = TRUE
) {
  write_tsv(
    as.data.table(
      matrix_object,
      keep.rownames = row_id_name
    ),
    path,
    compress = compress
  )
}

check_columns <- function(
    x,
    required_columns,
    label
) {
  missing_columns <- setdiff(
    required_columns,
    names(x)
  )

  if (length(missing_columns) > 0L) {
    stop(
      label,
      " lacks required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
}

qc_row <- function(check, passed, detail) {
  data.table(
    check = check,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    detail = as.character(detail)
  )
}

theme_nature <- function(
    base_size = 7,
    base_family = NATURE_FONT
) {
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      text = element_text(
        colour = "black"
      ),
      panel.grid = element_blank(),
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      ),
      axis.line = element_line(
        linewidth = 0.3,
        colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.3,
        colour = "black"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "plain",
        colour = "black"
      ),
      plot.title = element_text(
        face = "plain",
        colour = "black",
        hjust = 0
      ),
      plot.subtitle = element_text(
        face = "plain",
        colour = "grey20",
        hjust = 0
      ),
      plot.caption = element_text(
        face = "plain",
        colour = "grey20",
        hjust = 0
      ),
      legend.title = element_text(
        face = "plain"
      )
    )
}

save_plot <- function(
    plot_object,
    stem,
    width_mm,
    height_mm
) {
  pdf_path <- file.path(
    FIGURE_DIR,
    paste0(
      stem,
      ".pdf"
    )
  )

  png_path <- file.path(
    FIGURE_DIR,
    paste0(
      stem,
      ".png"
    )
  )

  ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = png_path,
    plot = plot_object,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = DPI,
    bg = "white"
  )

  message(
    "Wrote: ",
    pdf_path
  )

  message(
    "Wrote: ",
    png_path
  )
}

build_cancer_anchor_weights <- function(
    anchor_manifest
) {
  check_columns(
    anchor_manifest,
    c(
      "sample_id",
      "primary_gse_id",
      "analysis_cancer_class"
    ),
    "Anchor manifest"
  )

  weights <- copy(
    anchor_manifest[
      ,
      .(
        sample_id,
        primary_gse_id,
        provisional_cohort,
        analysis_cancer_class
      )
    ]
  )

  weights[
    ,
    samples_in_study_within_cancer :=
      .N,
    by = .(
      analysis_cancer_class,
      primary_gse_id
    )
  ]

  weights[
    ,
    studies_in_cancer :=
      uniqueN(
        primary_gse_id
      ),
    by = analysis_cancer_class
  ]

  weights[
    ,
    cancer_specific_study_balanced_weight :=
      1 /
        studies_in_cancer /
        samples_in_study_within_cancer
  ]

  weights[
    ,
    cancer_specific_study_balanced_weight :=
      cancer_specific_study_balanced_weight /
        sum(
          cancer_specific_study_balanced_weight
        ),
    by = analysis_cancer_class
  ]

  weights[
    ,
    effective_anchor_n :=
      1 /
        sum(
          cancer_specific_study_balanced_weight^2
        ),
    by = analysis_cancer_class
  ]

  setorder(
    weights,
    analysis_cancer_class,
    primary_gse_id,
    sample_id
  )

  weights
}

weighted_cancer_scaling <- function(
    canonical_matrix,
    anchor_matrix,
    canonical_ids,
    anchor_ids,
    weights,
    cancer_name,
    branch_name
) {
  if (!identical(
    rownames(
      canonical_matrix
    ),
    rownames(
      anchor_matrix
    )
  )) {
    stop(
      branch_name,
      "/",
      cancer_name,
      ": canonical and anchor gene orders differ.",
      call. = FALSE
    )
  }

  validate_id_order(
    colnames(
      canonical_matrix
    ),
    canonical_ids,
    paste0(
      branch_name,
      "/",
      cancer_name,
      " canonical matrix"
    )
  )

  validate_id_order(
    colnames(
      anchor_matrix
    ),
    anchor_ids,
    paste0(
      branch_name,
      "/",
      cancer_name,
      " anchor matrix"
    )
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
    !is.finite(
      anchor_copy_difference
    ) ||
      anchor_copy_difference >
        MATRIX_TOLERANCE
  ) {
    stop(
      branch_name,
      "/",
      cancer_name,
      ": anchor values differ from canonical subset; max |difference|=",
      format(
        anchor_copy_difference,
        scientific = TRUE
      ),
      call. = FALSE
    )
  }

  if (
    length(weights) !=
      ncol(
        anchor_matrix
      ) ||
      any(
        !is.finite(
          weights
        )
      ) ||
      any(
        weights <=
          0
      )
  ) {
    stop(
      branch_name,
      "/",
      cancer_name,
      ": invalid anchor weights.",
      call. = FALSE
    )
  }

  weights <- weights /
    sum(weights)

  unbiased_denominator <- 1 -
    sum(
      weights^2
    )

  if (
    !is.finite(
      unbiased_denominator
    ) ||
      unbiased_denominator <=
        0
  ) {
    stop(
      branch_name,
      "/",
      cancer_name,
      ": weighted-SD denominator is non-positive.",
      call. = FALSE
    )
  }

  log_canonical <- log2(
    canonical_matrix +
      1
  )

  log_anchor <- log2(
    anchor_matrix +
      1
  )

  weighted_mean <- as.vector(
    log_anchor %*%
      weights
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
  ) /
    unbiased_denominator

  weighted_sd <- sqrt(
    pmax(
      weighted_variance,
      0
    )
  )

  zero_sd <- !is.finite(
    weighted_sd
  ) |
    weighted_sd <=
      SD_FLOOR

  scaled <- matrix(
    0,
    nrow = nrow(
      canonical_matrix
    ),
    ncol = ncol(
      canonical_matrix
    ),
    dimnames = dimnames(
      canonical_matrix
    )
  )

  variable <- !zero_sd

  if (any(variable)) {
    scaled[
      variable,
      ] <- sweep(
        sweep(
          log_canonical[
            variable,
            ,
            drop = FALSE
          ],
          1L,
          weighted_mean[
            variable
          ],
          "-"
        ),
        1L,
        weighted_sd[
          variable
        ],
        "/"
      )
  }

  if (
    any(
      !is.finite(
        scaled
      )
    ) ||
      all(
        scaled ==
          0
      )
  ) {
    stop(
      branch_name,
      "/",
      cancer_name,
      ": scaled matrix is non-finite or entirely zero.",
      call. = FALSE
    )
  }

  scaled_anchor <- scaled[
    ,
    anchor_ids,
    drop = FALSE
  ]

  scaled_weighted_mean <- rep(
    NA_real_,
    nrow(
      scaled
    )
  )

  scaled_weighted_sd <- rep(
    NA_real_,
    nrow(
      scaled
    )
  )

  if (any(variable)) {
    scaled_weighted_mean[
      variable
    ] <- as.vector(
      scaled_anchor[
        variable,
        ,
        drop = FALSE
      ] %*%
        weights
    )

    scaled_centered <- sweep(
      scaled_anchor[
        variable,
        ,
        drop = FALSE
      ],
      1L,
      scaled_weighted_mean[
        variable
      ],
      "-"
    )

    scaled_weighted_sd[
      variable
    ] <- sqrt(
      rowSums(
        sweep(
          scaled_centered^2,
          2L,
          weights,
          "*"
        )
      ) /
        unbiased_denominator
    )
  }

  parameter_table <- data.table(
    analysis_cancer_class = cancer_name,
    branch = branch_name,
    model_gene = rownames(
      canonical_matrix
    ),
    weighted_anchor_mean_log2 = weighted_mean,
    weighted_anchor_sd_log2 = weighted_sd,
    zero_or_invalid_anchor_SD = zero_sd,
    weighted_scaled_anchor_mean = scaled_weighted_mean,
    weighted_scaled_anchor_sd = scaled_weighted_sd,
    canonical_scaled_min = apply(
      scaled,
      1L,
      min
    ),
    canonical_scaled_max = apply(
      scaled,
      1L,
      max
    ),
    effective_anchor_n = 1 /
      sum(
        weights^2
      ),
    anchor_copy_max_absolute_difference = anchor_copy_difference
  )

  list(
    scaled = scaled,
    parameters = parameter_table,
    zero_sd_genes = rownames(
      canonical_matrix
    )[
      zero_sd
    ],
    effective_n = 1 /
      sum(
        weights^2
      ),
    anchor_copy_difference = anchor_copy_difference,
    max_abs_weighted_scaled_mean = if (
      any(variable)
    ) {
      max(
        abs(
          scaled_weighted_mean[
            variable
          ]
        )
      )
    } else {
      NA_real_
    },
    max_abs_weighted_scaled_sd_error = if (
      any(variable)
    ) {
      max(
        abs(
          scaled_weighted_sd[
            variable
          ] -
            1
        )
      )
    } else {
      NA_real_
    }
  )
}

coordinate_shift <- function(
    pan_scores,
    cancer_specific_scores
) {
  if (
    !identical(
      rownames(
        pan_scores
      ),
      rownames(
        cancer_specific_scores
      )
    ) ||
      !identical(
        colnames(
          pan_scores
        ),
        colnames(
          cancer_specific_scores
        )
      )
  ) {
    stop(
      "Pan-cancer and cancer-specific score matrices are not aligned.",
      call. = FALSE
    )
  }

  absolute_difference <- abs(
    pan_scores -
      cancer_specific_scores
  )

  data.table(
    sample_id = rownames(
      pan_scores
    ),
    pan_to_cancer_specific_TVD =
      0.5 *
        rowSums(
          absolute_difference
        ),
    pan_to_cancer_specific_max_component_difference =
      apply(
        absolute_difference,
        1L,
        max
      ),
    pan_dominant_component =
      dominant_component(
        pan_scores
      ),
    cancer_specific_dominant_component =
      dominant_component(
        cancer_specific_scores
      ),
    dominant_component_concordant =
      dominant_component(
        pan_scores
      ) ==
        dominant_component(
          cancer_specific_scores
        )
  )
}

run_cancer_order_sentinel <- function(
    scaled_matrix,
    W,
    cancer_name
) {
  raw_1 <- project_raw_h(
    scaled_matrix,
    W,
    paste0(
      cancer_name,
      "_repeat_1"
    )
  )

  raw_2 <- project_raw_h(
    scaled_matrix,
    W,
    paste0(
      cancer_name,
      "_repeat_2"
    )
  )

  fraction_1 <- score_fractions(
    raw_1
  )

  fraction_2 <- score_fractions(
    raw_2
  )

  repeat_max <- max(
    abs(
      fraction_1 -
        fraction_2
    )
  )

  reversed_matrix <- scaled_matrix[
    ,
    rev(
      colnames(
        scaled_matrix
      )
    ),
    drop = FALSE
  ]

  raw_reversed <- project_raw_h(
    reversed_matrix,
    W,
    paste0(
      cancer_name,
      "_reversed"
    )
  )

  fraction_reversed <- score_fractions(
    raw_reversed
  )

  fraction_reversed <- fraction_reversed[
    rownames(
      fraction_1
    ),
    ,
    drop = FALSE
  ]

  absolute_difference <- abs(
    fraction_1 -
      fraction_reversed
  )

  data.table(
    analysis_cancer_class = cancer_name,
    n_samples = nrow(
      fraction_1
    ),
    same_order_repeat_max_absolute_difference = repeat_max,
    reversed_order_median_TVD =
      median(
        0.5 *
          rowSums(
            absolute_difference
          )
      ),
    reversed_order_maximum_TVD =
      max(
        0.5 *
          rowSums(
            absolute_difference
          )
      ),
    reversed_order_dominant_concordance =
      mean(
        dominant_component(
          fraction_1
        ) ==
          dominant_component(
            fraction_reversed
          )
      )
  )
}

# ==============================================================================
# 3. Reused validated projection helpers
# ==============================================================================

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
      ": Stage 09 matrix contains case-insensitive duplicate row symbols: ",
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
          " maps ambiguously to Stage 09 rows: ",
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
      " frozen W genes cannot be resolved from the Stage 09 matrix: ",
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
      ": at least one Stage 09 row maps to more than one frozen W gene:\n",
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

# ==============================================================================
# 4. Self-tests
# ==============================================================================

run_self_tests <- function() {
  test_manifest <- data.table(
    sample_id = paste0(
      "S",
      1:6
    ),
    primary_gse_id = c(
      "G1",
      "G1",
      "G2",
      "G3",
      "G3",
      "G3"
    ),
    provisional_cohort = "test",
    analysis_cancer_class = c(
      "A",
      "A",
      "A",
      "B",
      "B",
      "B"
    )
  )

  test_weights <- build_cancer_anchor_weights(
    test_manifest
  )

  stopifnot(
    max(
      abs(
        test_weights[
          ,
          sum(
            cancer_specific_study_balanced_weight
          ),
          by = analysis_cancer_class
        ]$V1 -
          1
      )
    ) <=
      1e-15,
    max(
      abs(
        test_weights[
          analysis_cancer_class ==
            "A" &
            primary_gse_id ==
            "G1",
          sum(
            cancer_specific_study_balanced_weight
          )
        ] -
          0.5
      )
    ) <=
      1e-15
  )

  canonical <- matrix(
    c(
      0, 0, 0, 0,
      1, 2, 3, 4
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      c(
        "G_ZERO",
        "G_VAR"
      ),
      paste0(
        "S",
        1:4
      )
    )
  )

  anchor <- canonical[
    ,
    1:3,
    drop = FALSE
  ]

  scaling <- weighted_cancer_scaling(
    canonical_matrix = canonical,
    anchor_matrix = anchor,
    canonical_ids = colnames(
      canonical
    ),
    anchor_ids = colnames(
      anchor
    ),
    weights = rep(
      1 / 3,
      3
    ),
    cancer_name = "TEST",
    branch_name = "CPM"
  )

  stopifnot(
    identical(
      scaling$zero_sd_genes,
      "G_ZERO"
    ),
    all(
      scaling$scaled[
        "G_ZERO",
        ] ==
          0
    ),
    max(
      abs(
        mean(
          scaling$scaled[
            "G_VAR",
            1:3
          ]
        )
      )
    ) <=
      1e-15
  )

  test_scores <- matrix(
    c(
      0.4, 0.6, rep(0, 8),
      0.3, 0.7, rep(0, 8),
      0.2, 0.8, rep(0, 8)
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(
      c("S1", "S2", "S3"),
      COMPONENTS
    )
  )

  sensitivity_test <- suppressWarnings(
    normalization_sensitivity(
      primary_scores = test_scores,
      alternative_scores = test_scores,
      comparison_name = "SELFTEST_UQ",
      evaluable_mask = c(TRUE, TRUE, FALSE),
      non_evaluable_reason = c(NA, NA, "fallback")
    )
  )

  stopifnot(
    nrow(sensitivity_test$sample) == 3L,
    sum(sensitivity_test$sample$comparison_evaluable) == 2L,
    is.na(sensitivity_test$sample[sample_id == "S3", total_variation_distance]),
    sensitivity_test$sample[sample_id == "S3", diagnostic_total_variation_distance] == 0,
    all(sensitivity_test$sample[comparison_evaluable == TRUE, total_variation_distance] == 0)
  )

  cat(
    "ALL CANCER-SPECIFIC FIXED-ANCHOR PROJECTION V1.0.3 SELF-TESTS PASSED\n"
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
# 5. Validate inputs and locked assets
# ==============================================================================

required_files <- c(
  FROZEN_W_FILE,
  INSTALLED_W_FILE,
  CANONICAL_MANIFEST_FILE,
  ANCHOR_MANIFEST_FILE,
  MODEL_GENE_AUDIT_FILE,
  SAMPLE_NORMALIZATION_QC_FILE,
  INPUT_STATUS_FILE,
  PAN_STATUS_FILE,
  PAN_QC_FILE,
  LAYER1_STATUS_FILE,
  LAYER1_SAMPLE_ATLAS_FILE,
  LAYER1_QC_FILE,
  PAN_SCORE_FILE,
  PAN_ADEQUACY_FILE,
  PAN_RECONSTRUCTION_FILE,
  unlist(
    lapply(
      BRANCHES,
      unlist
    ),
    use.names = FALSE
  )
)

missing_files <- required_files[
  !file.exists(
    required_files
  ) |
    file.info(
      required_files
    )$size <=
      0
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

if (
  as.character(
    packageVersion(
      "SpatialEcoTyper"
    )
  ) !=
    EXPECTED_PACKAGE_VERSION
) {
  stop(
    "SpatialEcoTyper version mismatch.",
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
        length(
          COMPONENTS
        )
      )
    ) ||
    !identical(
      colnames(W),
      COMPONENTS
    ) ||
    !identical(
      W,
      installed_W
    )
) {
  stop(
    "Frozen W failed structural or installed-model identity validation.",
    call. = FALSE
  )
}

input_status <- jsonlite::read_json(
  INPUT_STATUS_FILE,
  simplifyVector = TRUE
)

pan_status <- jsonlite::read_json(
  PAN_STATUS_FILE,
  simplifyVector = TRUE
)

layer1_status <- jsonlite::read_json(
  LAYER1_STATUS_FILE,
  simplifyVector = TRUE
)

if (
  !identical(input_status$status, "PASS") ||
    !identical(input_status$script_version, EXPECTED_STAGE09_VERSION) ||
    input_status$qc_failures != 0 ||
    input_status$canonical_samples != EXPECTED_CANONICAL ||
    input_status$cross_sectional_samples != EXPECTED_ANCHOR ||
    input_status$cancer_classes != EXPECTED_CANCERS ||
    input_status$model_genes_available != EXPECTED_MODEL_GENES ||
    input_status$model_genes_unresolved != 0 ||
    !identical(pan_status$status, "PASS") ||
    !identical(pan_status$script_version, EXPECTED_PAN_VERSION) ||
    pan_status$qc_failures != 0 ||
    pan_status$canonical_samples != EXPECTED_CANONICAL ||
    pan_status$cross_sectional_anchor_samples != EXPECTED_ANCHOR ||
    pan_status$cancer_classes != EXPECTED_CANCERS ||
    pan_status$UQ_evaluable_samples != EXPECTED_UQ_EVALUABLE ||
    pan_status$UQ_non_evaluable_samples != EXPECTED_UQ_NOT_EVALUABLE ||
    pan_status$UQ_scaling_anchor_samples != EXPECTED_UQ_ANCHOR ||
    pan_status$directional_features != EXPECTED_DIRECTIONAL_FEATURES ||
    !identical(layer1_status$status, "PASS") ||
    !identical(layer1_status$script_version, EXPECTED_LAYER1_VERSION) ||
    layer1_status$qc_failures != 0 ||
    layer1_status$canonical_samples != EXPECTED_CANONICAL ||
    layer1_status$cancer_classes != EXPECTED_CANCERS ||
    layer1_status$UQ_evaluable_samples != EXPECTED_UQ_EVALUABLE ||
    layer1_status$UQ_non_evaluable_samples != EXPECTED_UQ_NOT_EVALUABLE
) {
  stop(
    "Prerequisite stages are not the required frozen Stage 09 v1.0.1, Stage 10 v1.0.4 and Stage 11 v1.0.2.1 zero-failure PASS states.",
    call. = FALSE
  )
}

pan_qc <- fread(
  PAN_QC_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

layer1_qc <- fread(
  LAYER1_QC_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

check_columns(pan_qc, c("check", "status", "detail"), "Stage 10 QC")
check_columns(layer1_qc, c("check", "status", "detail"), "Stage 11 QC")

if (
  nrow(pan_qc) == 0L || any(pan_qc$status != "PASS") ||
    nrow(layer1_qc) == 0L || any(layer1_qc$status != "PASS")
) {
  stop("Stage 10 or Stage 11 QC contract contains non-PASS rows.", call. = FALSE)
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

layer1_sample_atlas <- fread(
  LAYER1_SAMPLE_ATLAS_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

check_columns(
  canonical_manifest,
  c(
    "sample_id",
    "primary_gse_id",
    "provisional_cohort",
    "analysis_cancer_class"
  ),
  "Canonical manifest"
)

check_columns(
  anchor_manifest,
  c(
    "sample_id",
    "primary_gse_id",
    "provisional_cohort",
    "analysis_cancer_class"
  ),
  "Anchor manifest"
)

if (
  nrow(
    canonical_manifest
  ) !=
    EXPECTED_CANONICAL ||
    nrow(
      anchor_manifest
    ) !=
      EXPECTED_ANCHOR ||
    anyDuplicated(
      canonical_manifest$sample_id
    ) ||
    anyDuplicated(
      anchor_manifest$sample_id
    ) ||
    !all(
      anchor_manifest$sample_id %in%
        canonical_manifest$sample_id
    )
) {
  stop(
    "Canonical or anchor manifest failed locked-count/subset validation.",
    call. = FALSE
  )
}

cancer_levels <- unique(
  canonical_manifest$analysis_cancer_class
)

if (
  length(
    cancer_levels
  ) !=
    EXPECTED_CANCERS ||
    !setequal(
      cancer_levels,
      unique(
        anchor_manifest$analysis_cancer_class
      )
    )
) {
  stop(
    "Expected the same ",
    EXPECTED_CANCERS,
    " cancer classes in canonical and anchor manifests.",
    call. = FALSE
  )
}

check_columns(
  sample_normalization_qc,
  c(
    "sample_id",
    "primary_gse_id",
    "analysis_cancer_class",
    "upper_quartile_fraction",
    "tmm_factor",
    "uq_factor",
    "model_genes_detected_recomputed",
    "stage08_tier",
    "stage07_expression_quality_tier"
  ),
  "Stage 09 sample-normalization QC"
)

check_columns(
  layer1_sample_atlas,
  c(
    "sample_id",
    "UQ_evaluable",
    "projection_evidence_tier",
    "adequacy_composite_percentile",
    "TMM_stable",
    "UQ_robustness_state",
    "recommended_use"
  ),
  "Stage 11 sample-level atlas"
)

if (
  nrow(sample_normalization_qc) != EXPECTED_CANONICAL ||
    anyDuplicated(sample_normalization_qc$sample_id) ||
    !identical(as.character(sample_normalization_qc$sample_id), as.character(canonical_manifest$sample_id)) ||
    nrow(layer1_sample_atlas) != EXPECTED_CANONICAL ||
    anyDuplicated(layer1_sample_atlas$sample_id) ||
    !identical(as.character(layer1_sample_atlas$sample_id), as.character(canonical_manifest$sample_id))
) {
  stop(
    "Stage 09 normalization QC or Stage 11 sample atlas is not in exact canonical order.",
    call. = FALSE
  )
}

if (
  nrow(model_gene_audit) != EXPECTED_MODEL_GENES ||
    sum(model_gene_audit$resolution == "alias_resolved") != 7L ||
    any(model_gene_audit$resolution == "unresolved")
) {
  stop(
    "Model-gene audit is not the required 492 exact + 7 alias-resolved, zero-unresolved state.",
    call. = FALSE
  )
}

canonical_ids <- as.character(canonical_manifest$sample_id)
anchor_ids <- as.character(anchor_manifest$sample_id)

uq_evaluable_mask <-
  is.finite(sample_normalization_qc$upper_quartile_fraction) &
    sample_normalization_qc$upper_quartile_fraction > 0

layer1_sample_atlas[, UQ_evaluable := as.logical(UQ_evaluable)]

if (
  sum(uq_evaluable_mask) != EXPECTED_UQ_EVALUABLE ||
    sum(!uq_evaluable_mask) != EXPECTED_UQ_NOT_EVALUABLE ||
    sum(layer1_sample_atlas$UQ_evaluable) != EXPECTED_UQ_EVALUABLE ||
    !identical(layer1_sample_atlas$UQ_evaluable, uq_evaluable_mask)
) {
  stop(
    "Stage 09 and Stage 11 UQ evaluability partitions are not the required identical 9,302/163 split.",
    call. = FALSE
  )
}

uq_evaluable_ids <- canonical_ids[uq_evaluable_mask]
uq_anchor_manifest <- anchor_manifest[sample_id %in% uq_evaluable_ids]

if (
  nrow(uq_anchor_manifest) != EXPECTED_UQ_ANCHOR ||
    uniqueN(uq_anchor_manifest$analysis_cancer_class) != EXPECTED_CANCERS
) {
  stop(
    "UQ-evaluable cancer-specific anchor is not the required 8,799-sample, 18-cancer set.",
    call. = FALSE
  )
}

cancer_weights <- build_cancer_anchor_weights(anchor_manifest)
uq_cancer_weights <- build_cancer_anchor_weights(uq_anchor_manifest)

weight_sums <- rbindlist(
  list(
    cancer_weights[
      ,
      .(weight_sum = sum(cancer_specific_study_balanced_weight)),
      by = analysis_cancer_class
    ][, branch := "CPM_TMM"],
    uq_cancer_weights[
      ,
      .(weight_sum = sum(cancer_specific_study_balanced_weight)),
      by = analysis_cancer_class
    ][, branch := "UQ_CPM"]
  ),
  use.names = TRUE
)

if (max(abs(weight_sums$weight_sum - 1)) > 1e-12) {
  stop("Within-cancer anchor weights do not sum to one.", call. = FALSE)
}

anchor_counts_by_branch <- rbindlist(
  list(
    anchor_manifest[
      ,
      .(
        anchor_samples = .N,
        anchor_studies = uniqueN(primary_gse_id)
      ),
      by = analysis_cancer_class
    ][, branch := "CPM"],
    anchor_manifest[
      ,
      .(
        anchor_samples = .N,
        anchor_studies = uniqueN(primary_gse_id)
      ),
      by = analysis_cancer_class
    ][, branch := "TMM_CPM"],
    uq_anchor_manifest[
      ,
      .(
        anchor_samples = .N,
        anchor_studies = uniqueN(primary_gse_id)
      ),
      by = analysis_cancer_class
    ][, branch := "UQ_CPM"]
  ),
  use.names = TRUE
)

setcolorder(
  anchor_counts_by_branch,
  c("analysis_cancer_class", "branch", "anchor_samples", "anchor_studies")
)

if (
  any(anchor_counts_by_branch$anchor_samples < MIN_ANCHOR_SAMPLES_PER_CANCER) ||
    sum(anchor_counts_by_branch[branch == "CPM", anchor_samples]) != EXPECTED_ANCHOR ||
    sum(anchor_counts_by_branch[branch == "TMM_CPM", anchor_samples]) != EXPECTED_ANCHOR ||
    sum(anchor_counts_by_branch[branch == "UQ_CPM", anchor_samples]) != EXPECTED_UQ_ANCHOR
) {
  stop(
    "At least one branch/cancer anchor is too small or branch anchor totals are incorrect.",
    call. = FALSE
  )
}

branch_anchor_manifests <- list(
  CPM = anchor_manifest,
  TMM_CPM = anchor_manifest,
  UQ_CPM = uq_anchor_manifest
)

branch_weight_tables <- list(
  CPM = cancer_weights,
  TMM_CPM = cancer_weights,
  UQ_CPM = uq_cancer_weights
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

layer1_metadata <- layer1_sample_atlas[
  ,
  .(
    sample_id,
    layer1_projection_evidence_tier = projection_evidence_tier,
    layer1_adequacy_composite_percentile = adequacy_composite_percentile,
    layer1_TMM_stable = TMM_stable,
    layer1_UQ_evaluable = UQ_evaluable,
    layer1_UQ_robustness_state = UQ_robustness_state,
    layer1_recommended_use = recommended_use
  )
]

input_provenance <- data.table(
  role = c(
    "stage09_status",
    "canonical_manifest",
    "anchor_manifest",
    "sample_normalization_qc",
    "model_gene_audit",
    "stage10_status",
    "stage10_QC",
    "stage10_scores",
    "stage10_adequacy",
    "stage10_reconstruction",
    "stage11_status",
    "stage11_QC",
    "stage11_sample_atlas",
    "frozen_W",
    "installed_W"
  ),
  path = c(
    INPUT_STATUS_FILE,
    CANONICAL_MANIFEST_FILE,
    ANCHOR_MANIFEST_FILE,
    SAMPLE_NORMALIZATION_QC_FILE,
    MODEL_GENE_AUDIT_FILE,
    PAN_STATUS_FILE,
    PAN_QC_FILE,
    PAN_SCORE_FILE,
    PAN_ADEQUACY_FILE,
    PAN_RECONSTRUCTION_FILE,
    LAYER1_STATUS_FILE,
    LAYER1_QC_FILE,
    LAYER1_SAMPLE_ATLAS_FILE,
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
  canonical_manifest,
  file.path(OUTPUT_DIR, "01_canonical_library_manifest.tsv.gz"),
  compress = TRUE
)

write_tsv(
  cancer_weights,
  file.path(OUTPUT_DIR, "02_cancer_specific_cross_sectional_anchor_weights.tsv.gz"),
  compress = TRUE
)

write_tsv(
  projection_evaluability,
  file.path(OUTPUT_DIR, "02A_projection_normalization_evaluability.tsv.gz"),
  compress = TRUE
)

write_tsv(
  uq_cancer_weights,
  file.path(OUTPUT_DIR, "02B_UQ_evaluable_cancer_specific_anchor_weights.tsv.gz"),
  compress = TRUE
)

write_tsv(
  anchor_counts_by_branch,
  file.path(OUTPUT_DIR, "02C_anchor_count_summary_by_cancer_branch.tsv")
)

# Load pan-cancer references for direct comparison.
pan_scores_dt <- fread(
  PAN_SCORE_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

pan_adequacy <- fread(
  PAN_ADEQUACY_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

pan_reconstruction <- readRDS(
  PAN_RECONSTRUCTION_FILE
)

check_columns(
  pan_scores_dt,
  c(
    "sample_id",
    COMPONENTS
  ),
  "Pan-cancer score table"
)

check_columns(
  pan_adequacy,
  c(
    "sample_id",
    "normalized_KL_per_input_mass",
    "reconstruction_Pearson",
    "reconstruction_cosine",
    "reconstruction_R2"
  ),
  "Pan-cancer adequacy table"
)

pan_scores_dt <- pan_scores_dt[
  match(
    canonical_ids,
    sample_id
  )
]

pan_adequacy <- pan_adequacy[
  match(
    canonical_ids,
    sample_id
  )
]

if (
  anyNA(
    pan_scores_dt$sample_id
  ) ||
    anyNA(
      pan_adequacy$sample_id
    )
) {
  stop(
    "Pan-cancer comparison tables do not cover all canonical samples.",
    call. = FALSE
  )
}

pan_score_matrix <- as.matrix(
  pan_scores_dt[
    ,
    ..COMPONENTS
  ]
)

rownames(
  pan_score_matrix
) <- pan_scores_dt$sample_id

pan_residual_matrix <- pan_reconstruction$signed_directional_residual

if (
  !is.matrix(
    pan_residual_matrix
  ) ||
    !setequal(
      colnames(
        pan_residual_matrix
      ),
      canonical_ids
    )
) {
  stop(
    "Pan-cancer residual matrix failed validation.",
    call. = FALSE
  )
}

pan_residual_matrix <- pan_residual_matrix[
  ,
  canonical_ids,
  drop = FALSE
]

# ==============================================================================
# 6. Load branch matrices once
# ==============================================================================

branch_input <- list()

for (branch_name in names(BRANCHES)) {
  canonical_matrix <- read_expression_matrix(
    BRANCHES[[branch_name]]$canonical
  )

  anchor_matrix <- read_expression_matrix(
    BRANCHES[[branch_name]]$anchor
  )

  if (
    nrow(
      canonical_matrix
    ) !=
      EXPECTED_MODEL_GENES ||
    nrow(
      anchor_matrix
    ) !=
      EXPECTED_MODEL_GENES
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
    colnames(
      canonical_matrix
    ),
    canonical_ids,
    paste0(
      branch_name,
      " canonical matrix"
    )
  )

  validate_id_order(
    colnames(
      anchor_matrix
    ),
    anchor_ids,
    paste0(
      branch_name,
      " anchor matrix"
    )
  )

  branch_input[[branch_name]] <- list(
    canonical = canonical_matrix,
    anchor = anchor_matrix
  )
}

# ==============================================================================
# 7. Cancer-specific scaling and projection
# ==============================================================================

message(
  "Running cancer-specific fixed-anchor projections."
)

branch_results <- list()
scaling_parameters <- list()
scaling_summary <- list()
coverage_rows <- list()
zero_sd_rows <- list()
harmonization_rows <- list()
sentinel_rows <- list()

W_model_genes <- unique(
  sub(
    "__.*$",
    "",
    rownames(W)
  )
)

for (branch_name in names(BRANCHES)) {
  message(
    "================================================================================"
  )
  message(
    "BRANCH: ",
    branch_name
  )

  full_scaled <- matrix(
    NA_real_,
    nrow = length(
      W_model_genes
    ),
    ncol = EXPECTED_CANONICAL,
    dimnames = list(
      W_model_genes,
      canonical_ids
    )
  )

  raw_h_full <- matrix(
    NA_real_,
    nrow = EXPECTED_CANONICAL,
    ncol = length(
      COMPONENTS
    ),
    dimnames = list(
      canonical_ids,
      COMPONENTS
    )
  )

  fractions_full <- raw_h_full

  branch_adequacy_rows <- list()
  branch_observed <- matrix(
    NA_real_,
    nrow = nrow(W),
    ncol = EXPECTED_CANONICAL,
    dimnames = list(
      rownames(W),
      canonical_ids
    )
  )

  branch_reconstructed <- branch_observed
  branch_residual <- branch_observed
  branch_feature_available <- matrix(
    FALSE,
    nrow = nrow(W),
    ncol = length(
      cancer_levels
    ),
    dimnames = list(
      rownames(W),
      cancer_levels
    )
  )

  for (cancer_name in cancer_levels) {
    message(
      "  ",
      cancer_name
    )

    cancer_canonical_ids <- canonical_manifest[
      analysis_cancer_class ==
        cancer_name,
      sample_id
    ]

    current_anchor_manifest <- branch_anchor_manifests[[branch_name]]
    current_weight_table <- branch_weight_tables[[branch_name]]

    cancer_anchor_ids <- current_anchor_manifest[
      analysis_cancer_class ==
        cancer_name,
      sample_id
    ]

    canonical_sub <- branch_input[[branch_name]]$canonical[
      ,
      cancer_canonical_ids,
      drop = FALSE
    ]

    anchor_sub <- branch_input[[branch_name]]$anchor[
      ,
      cancer_anchor_ids,
      drop = FALSE
    ]

    weight_sub <- current_weight_table[
      analysis_cancer_class ==
        cancer_name
    ][
      match(
        cancer_anchor_ids,
        sample_id
      ),
      cancer_specific_study_balanced_weight
    ]

    scaling <- weighted_cancer_scaling(
      canonical_matrix = canonical_sub,
      anchor_matrix = anchor_sub,
      canonical_ids = cancer_canonical_ids,
      anchor_ids = cancer_anchor_ids,
      weights = weight_sub,
      cancer_name = cancer_name,
      branch_name = branch_name
    )

    harmonized <- harmonize_expression_to_W(
      scaled_matrix = scaling$scaled,
      W = W,
      model_gene_audit = model_gene_audit,
      branch_name = paste0(
        branch_name,
        "/",
        cancer_name
      )
    )

    scaled_sub <- harmonized$matrix

    directional <- make_directional_bundle(
      scaled_sub,
      W
    )

    coverage <- component_feature_coverage(
      W = W,
      retained_features = directional$retained_features,
      branch_name = branch_name
    )

    coverage[
      ,
      analysis_cancer_class :=
        cancer_name
    ]

    setcolorder(
      coverage,
      c(
        "analysis_cancer_class",
        "branch",
        setdiff(
          names(
            coverage
          ),
          c(
            "analysis_cancer_class",
            "branch"
          )
        )
      )
    )

    if (
      directional$retained_feature_count <
        MIN_RETAINED_DIRECTIONAL_FEATURES ||
        min(
          coverage$retained_weight_fraction,
          na.rm = TRUE
        ) <
          MIN_COMPONENT_RETAINED_WEIGHT_FRACTION
    ) {
      stop(
        branch_name,
        "/",
        cancer_name,
        ": insufficient frozen-W coverage. Retained features=",
        directional$retained_feature_count,
        "; minimum component weight fraction=",
        format(
          min(
            coverage$retained_weight_fraction,
            na.rm = TRUE
          ),
          digits = 5
        ),
        call. = FALSE
      )
    }

    raw_h <- project_raw_h(
      scaled_sub,
      W,
      paste0(
        branch_name,
        "/",
        cancer_name
      )
    )

    fractions <- score_fractions(
      raw_h
    )

    directional_W <- directional$W
    V <- directional$V

    V_hat <- directional_W %*%
      t(
        raw_h
      )

    adequacy <- calculate_reconstruction_adequacy(
      V = V,
      V_hat = V_hat,
      fractions = fractions
    )

    adequacy[
      ,
      `:=`(
        analysis_cancer_class = cancer_name,
        branch = branch_name,
        retained_directional_features =
          directional$retained_feature_count
      )
    ]

    full_scaled[
      ,
      cancer_canonical_ids
    ] <- scaled_sub

    raw_h_full[
      cancer_canonical_ids,
      ] <- raw_h

    fractions_full[
      cancer_canonical_ids,
      ] <- fractions

    branch_adequacy_rows[[cancer_name]] <- adequacy

    retained_features <- directional$retained_features

    branch_observed[
      retained_features,
      cancer_canonical_ids
    ] <- V

    branch_reconstructed[
      retained_features,
      cancer_canonical_ids
    ] <- V_hat

    branch_residual[
      retained_features,
      cancer_canonical_ids
    ] <- V -
      V_hat

    branch_feature_available[
      retained_features,
      cancer_name
    ] <- TRUE

    scaling_parameters[[paste(branch_name, cancer_name, sep = "::")]] <- scaling$parameters

    scaling_summary[[paste(branch_name, cancer_name, sep = "::")]] <- data.table(
      analysis_cancer_class = cancer_name,
      branch = branch_name,
      canonical_samples =
        length(
          cancer_canonical_ids
        ),
      anchor_samples =
        length(
          cancer_anchor_ids
        ),
      anchor_studies =
        uniqueN(
          current_anchor_manifest[
            analysis_cancer_class ==
              cancer_name,
            primary_gse_id
          ]
        ),
      anchor_weight_scheme = if (branch_name == "UQ_CPM") {
        "UQ_evaluable_within_cancer_study_balanced"
      } else {
        "all_anchor_within_cancer_study_balanced"
      },
      effective_anchor_n =
        scaling$effective_n,
      zero_SD_model_genes =
        length(
          scaling$zero_sd_genes
        ),
      retained_directional_features =
        directional$retained_feature_count,
      minimum_component_retained_weight_fraction =
        min(
          coverage$retained_weight_fraction,
          na.rm = TRUE
        ),
      anchor_copy_max_absolute_difference =
        scaling$anchor_copy_difference,
      maximum_absolute_weighted_scaled_anchor_mean =
        scaling$max_abs_weighted_scaled_mean,
      maximum_absolute_weighted_scaled_anchor_SD_error =
        scaling$max_abs_weighted_scaled_sd_error,
      maximum_fraction_row_sum_error =
        max(
          abs(
            rowSums(
              fractions
            ) -
              1
          )
        ),
      exact_symbol_rows_to_W =
        harmonized$exact_count,
      case_harmonized_rows_to_W =
        harmonized$case_harmonized_count,
      alias_harmonized_rows_to_W =
        harmonized$alias_harmonized_count
    )

    coverage_rows[[paste(branch_name, cancer_name, sep = "::")]] <- coverage

    zero_sd_rows[[paste(branch_name, cancer_name, sep = "::")]] <- data.table(
      analysis_cancer_class = cancer_name,
      branch = branch_name,
      model_gene =
        scaling$zero_sd_genes
    )

    harmonization_table <- copy(
      harmonized$mapping
    )

    harmonization_table[
      ,
      analysis_cancer_class :=
        cancer_name
    ]

    harmonization_rows[[paste(branch_name, cancer_name, sep = "::")]] <- harmonization_table

    if (
      branch_name ==
        PRIMARY_BRANCH
    ) {
      sentinel_rows[[cancer_name]] <- run_cancer_order_sentinel(
        scaled_matrix = scaled_sub,
        W = W,
        cancer_name = cancer_name
      )
    }
  }

  if (
    any(
      !is.finite(
        full_scaled
      )
    ) ||
      any(
        !is.finite(
          raw_h_full
        )
      ) ||
      any(
        !is.finite(
          fractions_full
        )
      )
  ) {
    stop(
      branch_name,
      ": combined matrices contain non-finite values.",
      call. = FALSE
    )
  }

  branch_adequacy <- rbindlist(
    branch_adequacy_rows,
    use.names = TRUE
  )

  branch_adequacy <- branch_adequacy[
    match(
      canonical_ids,
      sample_id
    )
  ]

  branch_results[[branch_name]] <- list(
    scaled = full_scaled,
    raw_h = raw_h_full,
    fractions = fractions_full,
    adequacy = branch_adequacy,
    observed = branch_observed,
    reconstructed = branch_reconstructed,
    residual = branch_residual,
    feature_available = branch_feature_available
  )
}

# ==============================================================================
# 8. Write primary projection products
# ==============================================================================

scaling_parameter_table <- rbindlist(
  scaling_parameters,
  use.names = TRUE,
  fill = TRUE
)

scaling_summary_table <- rbindlist(
  scaling_summary,
  use.names = TRUE,
  fill = TRUE
)

coverage_table <- rbindlist(
  coverage_rows,
  use.names = TRUE,
  fill = TRUE
)

zero_sd_table <- rbindlist(
  zero_sd_rows,
  use.names = TRUE,
  fill = TRUE
)

if (
  is.null(
    zero_sd_table
  ) ||
    ncol(
      zero_sd_table
    ) ==
      0L
) {
  zero_sd_table <- data.table(
    analysis_cancer_class = character(),
    branch = character(),
    model_gene = character()
  )
}

harmonization_table <- rbindlist(
  harmonization_rows,
  use.names = TRUE,
  fill = TRUE
)

sentinel_table <- rbindlist(
  sentinel_rows,
  use.names = TRUE,
  fill = TRUE
)

write_tsv(
  scaling_parameter_table,
  file.path(
    OUTPUT_DIR,
    "03_cancer_specific_scaling_parameters_all_branches.tsv.gz"
  ),
  compress = TRUE
)

write_tsv(
  scaling_summary_table,
  file.path(
    OUTPUT_DIR,
    "04_cancer_specific_scaling_projection_summary.tsv"
  )
)

write_tsv(
  coverage_table,
  file.path(
    OUTPUT_DIR,
    "05_cancer_specific_component_feature_coverage.tsv"
  )
)

write_tsv(
  zero_sd_table,
  file.path(
    OUTPUT_DIR,
    "06_zero_SD_model_genes_by_cancer_branch.tsv"
  )
)

write_tsv(
  harmonization_table,
  file.path(
    OUTPUT_DIR,
    "07_frozen_W_gene_symbol_harmonization.tsv.gz"
  ),
  compress = TRUE
)

branch_output_numbers <- list(
  CPM = c(
    raw = "08",
    fraction = "09",
    scaled = "14"
  ),
  TMM_CPM = c(
    raw = "10",
    fraction = "11",
    scaled = "15"
  ),
  UQ_CPM = c(
    raw = "12",
    fraction = "13",
    scaled = "16"
  )
)

for (branch_name in names(BRANCHES)) {
  numbers <- branch_output_numbers[[branch_name]]

  raw_dt <- data.table(
    sample_id = canonical_ids,
    branch_results[[branch_name]]$raw_h
  )

  fraction_dt <- data.table(
    sample_id = canonical_ids,
    branch_results[[branch_name]]$fractions
  )

  write_tsv(
    raw_dt,
    file.path(
      OUTPUT_DIR,
      paste0(
        numbers[["raw"]],
        "_raw_H_",
        branch_name,
        ".tsv.gz"
      )
    ),
    compress = TRUE
  )

  write_tsv(
    fraction_dt,
    file.path(
      OUTPUT_DIR,
      paste0(
        numbers[["fraction"]],
        "_fraction_scores_",
        branch_name,
        ".tsv.gz"
      )
    ),
    compress = TRUE
  )

  saveRDS(
    branch_results[[branch_name]]$scaled,
    file.path(
      OUTPUT_DIR,
      paste0(
        numbers[["scaled"]],
        "_cancer_specific_scaled_input_",
        branch_name,
        ".rds"
      )
    ),
    compress = "gzip"
  )
}

primary_scores <- branch_results[[PRIMARY_BRANCH]]$fractions

primary_score_dt <- data.table(
  sample_id = canonical_ids,
  primary_scores
)

primary_score_metadata <- merge(
  canonical_manifest,
  primary_score_dt,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

primary_score_metadata <- merge(
  primary_score_metadata,
  layer1_metadata,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

primary_score_metadata <- primary_score_metadata[
  match(
    canonical_ids,
    sample_id
  )
]

primary_score_metadata[
  ,
  score_space_id :=
    paste0(
      "cancer_specific_",
      analysis_cancer_class
    )
]

write_tsv(
  primary_score_metadata,
  file.path(
    OUTPUT_DIR,
    "17_primary_CPM_cancer_specific_fraction_scores_with_metadata.tsv.gz"
  ),
  compress = TRUE
)

primary_adequacy <- branch_results[[PRIMARY_BRANCH]]$adequacy

primary_adequacy <- merge(
  canonical_manifest,
  primary_adequacy,
  by = c(
    "sample_id",
    "analysis_cancer_class"
  ),
  all.y = TRUE,
  sort = FALSE
)

primary_adequacy <- merge(
  primary_adequacy,
  layer1_metadata,
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
    "18_primary_CPM_cancer_specific_fixed_W_model_adequacy.tsv.gz"
  ),
  compress = TRUE
)

saveRDS(
  list(
    directional_feature_names = rownames(W),
    sample_ids = canonical_ids,
    observed_directional_features =
      branch_results[[PRIMARY_BRANCH]]$observed,
    reconstructed_directional_features =
      branch_results[[PRIMARY_BRANCH]]$reconstructed,
    signed_directional_residual =
      branch_results[[PRIMARY_BRANCH]]$residual,
    feature_available_by_cancer =
      branch_results[[PRIMARY_BRANCH]]$feature_available
  ),
  file.path(
    OUTPUT_DIR,
    "19_primary_CPM_cancer_specific_reconstruction_and_residual.rds"
  ),
  compress = "gzip"
)

# ==============================================================================
# 9. Normalization sensitivity within cancer-specific spaces
# ==============================================================================

uq_non_evaluable_reason <- ifelse(
  uq_evaluable_mask,
  NA_character_,
  "upper_quartile_fraction_nonpositive_stage09_fallback"
)

tmm_sensitivity <- normalization_sensitivity(
  primary_scores,
  branch_results$TMM_CPM$fractions,
  "CPM_vs_TMM_CPM"
)

uq_sensitivity <- normalization_sensitivity(
  primary_scores,
  branch_results$UQ_CPM$fractions,
  "CPM_vs_UQ_CPM",
  evaluable_mask = uq_evaluable_mask,
  non_evaluable_reason = uq_non_evaluable_reason
)

normalization_sample <- rbindlist(
  list(
    tmm_sensitivity$sample,
    uq_sensitivity$sample
  ),
  use.names = TRUE
)

normalization_sample <- merge(
  normalization_sample,
  projection_evaluability,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

normalization_sample <- merge(
  normalization_sample,
  layer1_metadata,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

normalization_sample[
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
  normalization_sample,
  canonical_order_index,
  comparison_order_index
)

normalization_sample[
  ,
  c("canonical_order_index", "comparison_order_index") := NULL
]

normalization_component <- rbindlist(
  lapply(
    cancer_levels,
    function(cancer_name) {
      cancer_ids <- canonical_manifest[
        analysis_cancer_class == cancer_name,
        sample_id
      ]

      cancer_uq_mask <- uq_evaluable_mask[
        match(cancer_ids, canonical_ids)
      ]

      cancer_uq_reason <- uq_non_evaluable_reason[
        match(cancer_ids, canonical_ids)
      ]

      tmm_component <- normalization_sensitivity(
        primary_scores[cancer_ids, , drop = FALSE],
        branch_results$TMM_CPM$fractions[cancer_ids, , drop = FALSE],
        "CPM_vs_TMM_CPM"
      )$component

      uq_component <- normalization_sensitivity(
        primary_scores[cancer_ids, , drop = FALSE],
        branch_results$UQ_CPM$fractions[cancer_ids, , drop = FALSE],
        "CPM_vs_UQ_CPM",
        evaluable_mask = cancer_uq_mask,
        non_evaluable_reason = cancer_uq_reason
      )$component

      out <- rbindlist(
        list(tmm_component, uq_component),
        use.names = TRUE
      )

      out[, analysis_cancer_class := cancer_name]

      setcolorder(
        out,
        c(
          "analysis_cancer_class",
          setdiff(names(out), "analysis_cancer_class")
        )
      )

      out
    }
  ),
  use.names = TRUE
)

normalization_summary <- normalization_sample[
  ,
  .(
    n_samples_total = .N,
    n_samples_evaluable = sum(comparison_evaluable),
    n_samples_not_evaluable = sum(!comparison_evaluable),
    median_TVD = stats::median(total_variation_distance, na.rm = TRUE),
    q95_TVD = stats::quantile(
      total_variation_distance,
      0.95,
      names = FALSE,
      na.rm = TRUE
    ),
    maximum_TVD = max(total_variation_distance, na.rm = TRUE),
    median_Aitchison = stats::median(aitchison_distance, na.rm = TRUE),
    q95_Aitchison = stats::quantile(
      aitchison_distance,
      0.95,
      names = FALSE,
      na.rm = TRUE
    ),
    maximum_Aitchison = max(aitchison_distance, na.rm = TRUE),
    dominant_component_concordance = mean(
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
  by = .(
    comparison,
    analysis_cancer_class
  )
]

write_tsv(
  normalization_sample,
  file.path(
    OUTPUT_DIR,
    "20_normalization_sensitivity_by_sample.tsv.gz"
  ),
  compress = TRUE
)

write_tsv(
  normalization_component,
  file.path(
    OUTPUT_DIR,
    "21_normalization_sensitivity_by_component.tsv"
  )
)

write_tsv(
  normalization_summary,
  file.path(
    OUTPUT_DIR,
    "22_normalization_sensitivity_by_cancer.tsv"
  )
)

# ==============================================================================
# 10. Shared pan-cancer versus cancer-specific comparison
# ==============================================================================

coordinate_shift_table <- coordinate_shift(
  pan_scores = pan_score_matrix,
  cancer_specific_scores = primary_scores
)

coordinate_shift_table <- merge(
  coordinate_shift_table,
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

coordinate_shift_table <- merge(
  coordinate_shift_table,
  layer1_metadata,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

coordinate_shift_table <- coordinate_shift_table[
  match(canonical_ids, sample_id)
]

coordinate_shift_summary <- coordinate_shift_table[
  ,
  .(
    n_samples = .N,
    median_TVD =
      median(
        pan_to_cancer_specific_TVD
      ),
    q95_TVD =
      quantile(
        pan_to_cancer_specific_TVD,
        0.95,
        names = FALSE
      ),
    maximum_TVD =
      max(
        pan_to_cancer_specific_TVD
      ),
    dominant_component_concordance =
      mean(
        dominant_component_concordant
      )
  ),
  by = analysis_cancer_class
][
  order(
    -median_TVD
  )
]

write_tsv(
  coordinate_shift_table,
  file.path(
    OUTPUT_DIR,
    "23_pan_vs_cancer_specific_coordinate_shift_by_sample.tsv.gz"
  ),
  compress = TRUE
)

write_tsv(
  coordinate_shift_summary,
  file.path(
    OUTPUT_DIR,
    "24_pan_vs_cancer_specific_coordinate_shift_by_cancer.tsv"
  )
)

adequacy_comparison <- merge(
  pan_adequacy[
    ,
    .(
      sample_id,
      pan_normalized_KL =
        normalized_KL_per_input_mass,
      pan_reconstruction_Pearson =
        reconstruction_Pearson,
      pan_reconstruction_cosine =
        reconstruction_cosine,
      pan_reconstruction_R2 =
        reconstruction_R2
    )
  ],
  primary_adequacy[
    ,
    .(
      sample_id,
      analysis_cancer_class,
      cancer_specific_normalized_KL =
        normalized_KL_per_input_mass,
      cancer_specific_reconstruction_Pearson =
        reconstruction_Pearson,
      cancer_specific_reconstruction_cosine =
        reconstruction_cosine,
      cancer_specific_reconstruction_R2 =
        reconstruction_R2
    )
  ],
  by = "sample_id",
  all = FALSE,
  sort = FALSE
)

adequacy_comparison[
  ,
  `:=`(
    delta_normalized_KL =
      cancer_specific_normalized_KL -
        pan_normalized_KL,
    delta_reconstruction_Pearson =
      cancer_specific_reconstruction_Pearson -
        pan_reconstruction_Pearson,
    delta_reconstruction_cosine =
      cancer_specific_reconstruction_cosine -
        pan_reconstruction_cosine,
    delta_reconstruction_R2 =
      cancer_specific_reconstruction_R2 -
        pan_reconstruction_R2
  )
]

adequacy_comparison[
  ,
  all_four_metrics_improved :=
    delta_normalized_KL <
      0 &
      delta_reconstruction_Pearson >
        0 &
      delta_reconstruction_cosine >
        0 &
      delta_reconstruction_R2 >
        0
]

adequacy_comparison <- merge(
  adequacy_comparison,
  layer1_metadata,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

adequacy_comparison <- adequacy_comparison[
  match(canonical_ids, sample_id)
]

adequacy_comparison_summary <- adequacy_comparison[
  ,
  .(
    n_samples = .N,
    median_delta_normalized_KL =
      median(
        delta_normalized_KL
      ),
    median_delta_reconstruction_Pearson =
      median(
        delta_reconstruction_Pearson
      ),
    median_delta_reconstruction_cosine =
      median(
        delta_reconstruction_cosine
      ),
    median_delta_reconstruction_R2 =
      median(
        delta_reconstruction_R2
      ),
    fraction_all_four_metrics_improved =
      mean(
        all_four_metrics_improved
      )
  ),
  by = analysis_cancer_class
][
  order(
    median_delta_normalized_KL
  )
]

write_tsv(
  adequacy_comparison,
  file.path(
    OUTPUT_DIR,
    "25_pan_vs_cancer_specific_adequacy_delta_by_sample.tsv.gz"
  ),
  compress = TRUE
)

write_tsv(
  adequacy_comparison_summary,
  file.path(
    OUTPUT_DIR,
    "26_pan_vs_cancer_specific_adequacy_delta_by_cancer.tsv"
  )
)

# ==============================================================================
# 11. Residual attenuation and persistence
# ==============================================================================

cs_residual <- branch_results[[PRIMARY_BRANCH]]$residual

if (
  !setequal(
    rownames(
      cs_residual
    ),
    rownames(
      pan_residual_matrix
    )
  )
) {
  stop(
    "Pan-cancer and cancer-specific residual feature sets differ.",
    call. = FALSE
  )
}

cs_residual <- cs_residual[
  rownames(
    pan_residual_matrix
  ),
  ,
  drop = FALSE
]

residual_comparison_rows <- list()

for (cancer_name in cancer_levels) {
  sample_ids <- canonical_manifest[
    analysis_cancer_class ==
      cancer_name,
    sample_id
  ]

  pan_sub <- pan_residual_matrix[
    ,
    sample_ids,
    drop = FALSE
  ]

  cs_sub <- cs_residual[
    ,
    sample_ids,
    drop = FALSE
  ]

  available <- rowSums(
    !is.na(
      cs_sub
    )
  ) >
    0

  pan_mean <- rowMeans(
    pan_sub,
    na.rm = TRUE
  )

  cs_mean <- rep(
    NA_real_,
    nrow(
      cs_sub
    )
  )

  cs_mean_abs <- rep(
    NA_real_,
    nrow(
      cs_sub
    )
  )

  cs_positive_fraction <- rep(
    NA_real_,
    nrow(
      cs_sub
    )
  )

  cs_negative_fraction <- rep(
    NA_real_,
    nrow(
      cs_sub
    )
  )

  cs_mean[
    available
  ] <- rowMeans(
    cs_sub[
      available,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  cs_mean_abs[
    available
  ] <- rowMeans(
    abs(
      cs_sub[
        available,
        ,
        drop = FALSE
      ]
    ),
    na.rm = TRUE
  )

  cs_positive_fraction[
    available
  ] <- rowMeans(
    cs_sub[
      available,
      ,
      drop = FALSE
    ] >
      0,
    na.rm = TRUE
  )

  cs_negative_fraction[
    available
  ] <- rowMeans(
    cs_sub[
      available,
      ,
      drop = FALSE
    ] <
      0,
    na.rm = TRUE
  )

  residual_comparison_rows[[cancer_name]] <- data.table(
    analysis_cancer_class = cancer_name,
    directional_feature = rownames(
      pan_residual_matrix
    ),
    pan_mean_signed_residual = pan_mean,
    cancer_specific_mean_signed_residual = cs_mean,
    cancer_specific_mean_absolute_residual = cs_mean_abs,
    cancer_specific_positive_fraction = cs_positive_fraction,
    cancer_specific_negative_fraction = cs_negative_fraction,
    feature_available = available
  )
}

residual_comparison <- rbindlist(
  residual_comparison_rows,
  use.names = TRUE
)

residual_comparison[
  ,
  absolute_pan_mean :=
    abs(
      pan_mean_signed_residual
    )
]

residual_comparison[
  ,
  absolute_cancer_specific_mean :=
    abs(
      cancer_specific_mean_signed_residual
    )
]

residual_comparison[
  ,
  absolute_mean_residual_reduction :=
    absolute_pan_mean -
      absolute_cancer_specific_mean
]

residual_comparison[
  ,
  attenuation_fraction := fifelse(
    absolute_pan_mean >=
      PAN_RESIDUAL_MIN_ABS_FOR_ATTENUATION,
    1 -
      absolute_cancer_specific_mean /
        absolute_pan_mean,
    NA_real_
  )
]

residual_comparison[
  ,
  cancer_specific_sign_concordance := fifelse(
    feature_available,
    pmax(
      cancer_specific_positive_fraction,
      cancer_specific_negative_fraction,
      na.rm = TRUE
    ),
    NA_real_
  )
]

residual_comparison[
  ,
  persistent_cancer_specific_residual :=
    feature_available &
      absolute_cancer_specific_mean >=
        PERSISTENT_RESIDUAL_ABS_MEAN_THRESHOLD &
      cancer_specific_sign_concordance >=
        PERSISTENT_RESIDUAL_SIGN_CONCORDANCE
]

write_tsv(
  residual_comparison,
  file.path(
    OUTPUT_DIR,
    "27_pan_vs_cancer_specific_residual_attenuation.tsv.gz"
  ),
  compress = TRUE
)

persistent_residual_ranking <- residual_comparison[
  feature_available ==
    TRUE
][
  ,
  .(
    cancers_persistent =
      sum(
        persistent_cancer_specific_residual
      ),
    maximum_absolute_cancer_specific_mean =
      max(
        absolute_cancer_specific_mean,
        na.rm = TRUE
      ),
    median_absolute_cancer_specific_mean =
      median(
        absolute_cancer_specific_mean,
        na.rm = TRUE
      ),
    median_attenuation_fraction =
      median(
        attenuation_fraction,
        na.rm = TRUE
      )
  ),
  by = directional_feature
][
  order(
    -cancers_persistent,
    -maximum_absolute_cancer_specific_mean
  )
]

write_tsv(
  persistent_residual_ranking,
  file.path(
    OUTPUT_DIR,
    "28_cancer_specific_persistent_residual_feature_ranking.tsv"
  )
)

residual_attenuation_summary <- residual_comparison[
  absolute_pan_mean >=
    PAN_RESIDUAL_MIN_ABS_FOR_ATTENUATION
][
  ,
  .(
    n_evaluable_features = .N,
    median_attenuation_fraction =
      median(
        attenuation_fraction,
        na.rm = TRUE
      ),
    q25_attenuation_fraction =
      quantile(
        attenuation_fraction,
        0.25,
        na.rm = TRUE,
        names = FALSE
      ),
    q75_attenuation_fraction =
      quantile(
        attenuation_fraction,
        0.75,
        na.rm = TRUE,
        names = FALSE
      ),
    fraction_attenuated_at_least_half =
      mean(
        attenuation_fraction >=
          0.5,
        na.rm = TRUE
      ),
    persistent_features =
      sum(
        persistent_cancer_specific_residual,
        na.rm = TRUE
      )
  ),
  by = analysis_cancer_class
][
  order(
    -median_attenuation_fraction
  )
]

write_tsv(
  residual_attenuation_summary,
  file.path(
    OUTPUT_DIR,
    "29_residual_attenuation_summary_by_cancer.tsv"
  )
)

# ==============================================================================
# 12. Determinism, figures and source data
# ==============================================================================

write_tsv(
  sentinel_table,
  file.path(
    OUTPUT_DIR,
    "30_primary_CPM_determinism_and_order_sentinel_by_cancer.tsv"
  )
)

# Figure 1: adequacy changes.
adequacy_plot_source <- merge(
  adequacy_comparison,
  canonical_manifest[
    ,
    .(
      sample_id,
      analysis_cancer_class_manifest =
        analysis_cancer_class
    )
  ],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

adequacy_plot_source[
  ,
  analysis_cancer_class :=
    fifelse(
      is.na(
        analysis_cancer_class
      ),
      analysis_cancer_class_manifest,
      analysis_cancer_class
    )
]

cancer_order <- adequacy_comparison_summary[
  order(
    median_delta_reconstruction_R2
  ),
  analysis_cancer_class
]

adequacy_plot_source[
  ,
  cancer_label := factor(
    analysis_cancer_class,
    levels = cancer_order
  )
]

write_tsv(
  adequacy_plot_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_1_adequacy_delta_source.tsv.gz"
  ),
  compress = TRUE
)

figure_1 <- ggplot(
  adequacy_plot_source,
  aes(
    x = delta_reconstruction_R2,
    y = cancer_label
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.3,
    linetype = "dashed",
    colour = "grey55"
  ) +
  geom_boxplot(
    width = 0.58,
    linewidth = 0.35,
    outlier.shape = NA,
    fill = "white",
    colour = "black"
  ) +
  labs(
    title = "Cancer-specific scaling changes fixed-W reconstruction adequacy",
    subtitle = "Positive values indicate higher reconstruction R² than in the shared pan-cancer coordinate system.",
    x = expression(Delta*" reconstruction R"^2*" (cancer-specific - pan-cancer)"),
    y = NULL
  ) +
  theme_nature(
    base_size = 7
  )

save_plot(
  figure_1,
  "FIGURE_1_pan_vs_cancer_specific_reconstruction_R2",
  width_mm = 183,
  height_mm = 112
)

# Figure 2: coordinate shift.
coordinate_plot_source <- copy(
  coordinate_shift_table
)

coordinate_plot_source[
  ,
  cancer_label := factor(
    analysis_cancer_class,
    levels = rev(
      coordinate_shift_summary$analysis_cancer_class
    )
  )
]

write_tsv(
  coordinate_plot_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_2_coordinate_shift_source.tsv.gz"
  ),
  compress = TRUE
)

figure_2 <- ggplot(
  coordinate_plot_source,
  aes(
    x = pan_to_cancer_specific_TVD,
    y = cancer_label
  )
) +
  geom_boxplot(
    width = 0.58,
    linewidth = 0.35,
    outlier.shape = NA,
    fill = "white",
    colour = "black"
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    )
  ) +
  labs(
    title = "Shared and cancer-specific coordinate systems yield distinct SE compositions",
    subtitle = "TVD is a descriptive coordinate-system shift, not a normalization-error metric.",
    x = "pan-cancer versus cancer-specific total variation distance",
    y = NULL
  ) +
  theme_nature(
    base_size = 7
  )

save_plot(
  figure_2,
  "FIGURE_2_pan_vs_cancer_specific_coordinate_shift",
  width_mm = 183,
  height_mm = 112
)

# Figure 3: residual attenuation.
residual_plot_source <- residual_comparison[
  feature_available ==
    TRUE &
    absolute_pan_mean >=
      PAN_RESIDUAL_MIN_ABS_FOR_ATTENUATION
]

write_tsv(
  residual_plot_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_3_residual_attenuation_source.tsv.gz"
  ),
  compress = TRUE
)

figure_3 <- ggplot(
  residual_plot_source,
  aes(
    x = absolute_pan_mean,
    y = absolute_cancer_specific_mean
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.3,
    linetype = "dashed",
    colour = "grey55"
  ) +
  geom_point(
    size = 0.55,
    alpha = 0.35,
    colour = "black"
  ) +
  facet_wrap(
    ~ analysis_cancer_class,
    ncol = 4,
    scales = "free"
  ) +
  labs(
    title = "Cancer-specific scaling separates tissue-baseline residuals from persistent model mismatch",
    subtitle = "Points below the diagonal have smaller cancer-mean residual magnitude after within-cancer anchoring.",
    x = "absolute mean residual in shared pan-cancer space",
    y = "absolute mean residual in cancer-specific space"
  ) +
  theme_nature(
    base_size = 6
  )

save_plot(
  figure_3,
  "FIGURE_3_residual_attenuation_by_cancer",
  width_mm = 183,
  height_mm = 150
)

# Figure 4: top persistent residual heatmap.
top_feature_table <- persistent_residual_ranking[
  cancers_persistent >
    0 |
    maximum_absolute_cancer_specific_mean >=
      PERSISTENT_RESIDUAL_ABS_MEAN_THRESHOLD
]

if (nrow(top_feature_table) == 0L) {
  top_feature_table <- persistent_residual_ranking
}

top_features <- head(
  top_feature_table,
  30L
)$directional_feature

heatmap_source <- residual_comparison[
  directional_feature %in%
    top_features &
    feature_available ==
      TRUE
]

feature_order <- persistent_residual_ranking[
  directional_feature %in%
    top_features,
  directional_feature
]

heatmap_source[
  ,
  directional_feature :=
    factor(
      directional_feature,
      levels = rev(
        feature_order
      )
    )
]

heatmap_source[
  ,
  analysis_cancer_class :=
    factor(
      analysis_cancer_class,
      levels = cancer_levels
    )
]

write_tsv(
  heatmap_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_4_persistent_residual_heatmap_source.tsv"
  )
)

heat_limit <- quantile(
  abs(
    heatmap_source$cancer_specific_mean_signed_residual
  ),
  0.98,
  na.rm = TRUE,
  names = FALSE
)

if (
  !is.finite(
    heat_limit
  ) ||
    heat_limit <=
      0
) {
  heat_limit <- 1e-6
}

figure_4 <- ggplot(
  heatmap_source,
  aes(
    x = analysis_cancer_class,
    y = directional_feature,
    fill = cancer_specific_mean_signed_residual
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -heat_limit,
      heat_limit
    ),
    oob = squish
  ) +
  labs(
    title = "Residual programs persisting after cancer-specific anchoring",
    subtitle = "Positive values are under-reconstructed directional features; negative values are over-reconstructed.",
    x = NULL,
    y = "directional feature",
    fill = "mean\nresidual"
  ) +
  theme_nature(
    base_size = 6
  ) +
  theme(
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 5.5
    )
  )

save_plot(
  figure_4,
  "FIGURE_4_persistent_cancer_specific_residuals",
  width_mm = 183,
  height_mm = 132
)

# ==============================================================================
# 13. QC contract, report and status
# ==============================================================================

n_tmm_evaluable <- normalization_sample[
  comparison == "CPM_vs_TMM_CPM" & comparison_evaluable == TRUE,
  .N
]

n_uq_evaluable <- normalization_sample[
  comparison == "CPM_vs_UQ_CPM" & comparison_evaluable == TRUE,
  .N
]

n_uq_not_evaluable <- normalization_sample[
  comparison == "CPM_vs_UQ_CPM" & comparison_evaluable == FALSE,
  .N
]

formal_sensitivity_columns <- c(
  "total_variation_distance",
  "maximum_component_absolute_difference",
  "aitchison_distance"
)

diagnostic_sensitivity_columns <- c(
  "diagnostic_total_variation_distance",
  "diagnostic_maximum_component_absolute_difference",
  "diagnostic_aitchison_distance"
)

qc <- rbindlist(
  list(
    qc_row(
      "prerequisite_frozen_stages_PASS",
      identical(input_status$status, "PASS") &&
        identical(input_status$script_version, EXPECTED_STAGE09_VERSION) &&
        identical(pan_status$status, "PASS") &&
        identical(pan_status$script_version, EXPECTED_PAN_VERSION) &&
        identical(layer1_status$status, "PASS") &&
        identical(layer1_status$script_version, EXPECTED_LAYER1_VERSION) &&
        all(pan_qc$status == "PASS") &&
        all(layer1_qc$status == "PASS"),
      paste0(
        "stage09=", input_status$status, "/", input_status$script_version,
        ";stage10=", pan_status$status, "/", pan_status$script_version,
        ";stage11=", layer1_status$status, "/", layer1_status$script_version
      )
    ),
    qc_row(
      "frozen_W_identical_to_installed_W",
      identical(W, installed_W),
      paste0(
        "frozen_md5=", unname(tools::md5sum(FROZEN_W_FILE)),
        ";installed_md5=", unname(tools::md5sum(INSTALLED_W_FILE))
      )
    ),
    qc_row(
      "within_cancer_anchor_weights_sum_one",
      max(abs(weight_sums$weight_sum - 1)) <= 1e-12,
      paste0(
        "max_error=",
        format(max(abs(weight_sums$weight_sum - 1)), scientific = TRUE)
      )
    ),
    qc_row(
      "branch_anchor_counts_exact",
      sum(anchor_counts_by_branch[branch == "CPM", anchor_samples]) == EXPECTED_ANCHOR &&
        sum(anchor_counts_by_branch[branch == "TMM_CPM", anchor_samples]) == EXPECTED_ANCHOR &&
        sum(anchor_counts_by_branch[branch == "UQ_CPM", anchor_samples]) == EXPECTED_UQ_ANCHOR &&
        min(anchor_counts_by_branch$anchor_samples) >= MIN_ANCHOR_SAMPLES_PER_CANCER,
      paste0(
        "CPM=", sum(anchor_counts_by_branch[branch == "CPM", anchor_samples]),
        ";TMM=", sum(anchor_counts_by_branch[branch == "TMM_CPM", anchor_samples]),
        ";UQ=", sum(anchor_counts_by_branch[branch == "UQ_CPM", anchor_samples]),
        ";minimum_cancer_branch_anchor=", min(anchor_counts_by_branch$anchor_samples)
      )
    ),
    qc_row(
      "UQ_evaluability_partition_exact",
      sum(uq_evaluable_mask) == EXPECTED_UQ_EVALUABLE &&
        sum(!uq_evaluable_mask) == EXPECTED_UQ_NOT_EVALUABLE &&
        nrow(uq_anchor_manifest) == EXPECTED_UQ_ANCHOR &&
        identical(layer1_sample_atlas$UQ_evaluable, uq_evaluable_mask),
      paste0(
        "evaluable=", sum(uq_evaluable_mask),
        ";non_evaluable=", sum(!uq_evaluable_mask),
        ";UQ_anchor=", nrow(uq_anchor_manifest)
      )
    ),
    qc_row(
      "single_internal_prediction_run_per_cancer",
      max(scaling_summary_table$canonical_samples) <= NSAMPLE_PER_RUN,
      paste0(
        "maximum_cancer_samples=", max(scaling_summary_table$canonical_samples),
        ";nsample_per_run=", NSAMPLE_PER_RUN
      )
    ),
    qc_row(
      "scaling_summary_complete",
      nrow(scaling_summary_table) == EXPECTED_CANCERS * length(BRANCHES) &&
        all(scaling_summary_table[, sum(canonical_samples), by = branch]$V1 == EXPECTED_CANONICAL) &&
        sum(scaling_summary_table[branch == "CPM", anchor_samples]) == EXPECTED_ANCHOR &&
        sum(scaling_summary_table[branch == "TMM_CPM", anchor_samples]) == EXPECTED_ANCHOR &&
        sum(scaling_summary_table[branch == "UQ_CPM", anchor_samples]) == EXPECTED_UQ_ANCHOR,
      paste0(
        "rows=", nrow(scaling_summary_table),
        ";expected=", EXPECTED_CANCERS * length(BRANCHES)
      )
    ),
    qc_row(
      "frozen_W_symbol_harmonization_expected",
      all(scaling_summary_table$exact_symbol_rows_to_W == 493L) &&
        all(scaling_summary_table$case_harmonized_rows_to_W == 6L) &&
        all(scaling_summary_table$alias_harmonized_rows_to_W == 0L),
      paste0(
        "exact_range=", min(scaling_summary_table$exact_symbol_rows_to_W), "-", max(scaling_summary_table$exact_symbol_rows_to_W),
        ";case_range=", min(scaling_summary_table$case_harmonized_rows_to_W), "-", max(scaling_summary_table$case_harmonized_rows_to_W),
        ";alias_range=", min(scaling_summary_table$alias_harmonized_rows_to_W), "-", max(scaling_summary_table$alias_harmonized_rows_to_W)
      )
    ),
    qc_row(
      "all_score_rows_complete",
      all(
        vapply(
          branch_results,
          function(x) {
            nrow(x$fractions) == EXPECTED_CANONICAL &&
              all(is.finite(x$fractions)) &&
              all(x$fractions >= 0)
          },
          logical(1)
        )
      ),
      paste0("samples=", EXPECTED_CANONICAL, ";branches=", length(branch_results))
    ),
    qc_row(
      "all_fraction_scores_sum_to_one",
      max(scaling_summary_table$maximum_fraction_row_sum_error) <= SUM_TO_ONE_TOLERANCE,
      paste0(
        "max_error=",
        format(max(scaling_summary_table$maximum_fraction_row_sum_error), scientific = TRUE)
      )
    ),
    qc_row(
      "minimum_directional_feature_coverage",
      min(scaling_summary_table$retained_directional_features) >= MIN_RETAINED_DIRECTIONAL_FEATURES,
      paste0(
        "minimum=", min(scaling_summary_table$retained_directional_features),
        ";threshold=", MIN_RETAINED_DIRECTIONAL_FEATURES
      )
    ),
    qc_row(
      "minimum_component_weight_coverage",
      min(scaling_summary_table$minimum_component_retained_weight_fraction) >=
        MIN_COMPONENT_RETAINED_WEIGHT_FRACTION,
      paste0(
        "minimum=",
        format(min(scaling_summary_table$minimum_component_retained_weight_fraction), digits = 6),
        ";threshold=", MIN_COMPONENT_RETAINED_WEIGHT_FRACTION
      )
    ),
    qc_row(
      "weighted_anchor_scaling_valid",
      max(scaling_summary_table$anchor_copy_max_absolute_difference) <= MATRIX_TOLERANCE &&
        max(
          scaling_summary_table$maximum_absolute_weighted_scaled_anchor_mean,
          na.rm = TRUE
        ) <= 1e-8 &&
        max(
          scaling_summary_table$maximum_absolute_weighted_scaled_anchor_SD_error,
          na.rm = TRUE
        ) <= 1e-8,
      paste0(
        "copy_max=",
        format(max(scaling_summary_table$anchor_copy_max_absolute_difference), scientific = TRUE),
        ";mean_max=",
        format(
          max(
            scaling_summary_table$maximum_absolute_weighted_scaled_anchor_mean,
            na.rm = TRUE
          ),
          scientific = TRUE
        ),
        ";sd_error_max=",
        format(
          max(
            scaling_summary_table$maximum_absolute_weighted_scaled_anchor_SD_error,
            na.rm = TRUE
          ),
          scientific = TRUE
        )
      )
    ),
    qc_row(
      "layer1_metadata_complete",
      nrow(layer1_metadata) == EXPECTED_CANONICAL &&
        !anyNA(layer1_metadata$layer1_projection_evidence_tier) &&
        !anyNA(primary_score_metadata$layer1_projection_evidence_tier) &&
        !anyNA(primary_adequacy$layer1_projection_evidence_tier),
      paste0("rows=", nrow(layer1_metadata), ";expected=", EXPECTED_CANONICAL)
    ),
    qc_row(
      "primary_adequacy_complete",
      nrow(primary_adequacy) == EXPECTED_CANONICAL &&
        all(is.finite(primary_adequacy$normalized_KL_per_input_mass)) &&
        all(is.finite(primary_adequacy$reconstruction_R2)),
      paste0("rows=", nrow(primary_adequacy))
    ),
    qc_row(
      "normalization_sensitivity_partition_complete",
      nrow(normalization_sample) == 2L * EXPECTED_CANONICAL &&
        n_tmm_evaluable == EXPECTED_CANONICAL &&
        n_uq_evaluable == EXPECTED_UQ_EVALUABLE &&
        n_uq_not_evaluable == EXPECTED_UQ_NOT_EVALUABLE,
      paste0(
        "rows=", nrow(normalization_sample),
        ";TMM_evaluable=", n_tmm_evaluable,
        ";UQ_evaluable=", n_uq_evaluable,
        ";UQ_non_evaluable=", n_uq_not_evaluable
      )
    ),
    qc_row(
      "normalization_formal_diagnostic_NA_contract",
      all(
        is.finite(
          as.matrix(normalization_sample[, ..diagnostic_sensitivity_columns])
        )
      ) &&
        all(
          is.finite(
            as.matrix(
              normalization_sample[
                comparison_evaluable == TRUE,
                ..formal_sensitivity_columns
              ]
            )
          )
        ) &&
        all(
          is.na(
            as.matrix(
              normalization_sample[
                comparison_evaluable == FALSE,
                ..formal_sensitivity_columns
              ]
            )
          )
        ) &&
        all(
          is.na(
            normalization_sample[
              comparison_evaluable == FALSE,
              dominant_component_concordant
            ]
          )
        ),
      paste0("formal_UQ_rows=", n_uq_evaluable, ";NA_UQ_rows=", n_uq_not_evaluable)
    ),
    qc_row(
      "coordinate_comparison_complete",
      nrow(coordinate_shift_table) == EXPECTED_CANONICAL &&
        nrow(adequacy_comparison) == EXPECTED_CANONICAL,
      paste0(
        "score_shift_rows=", nrow(coordinate_shift_table),
        ";adequacy_rows=", nrow(adequacy_comparison)
      )
    ),
    qc_row(
      "residual_comparison_complete",
      nrow(residual_comparison) == EXPECTED_DIRECTIONAL_FEATURES * EXPECTED_CANCERS,
      paste0(
        "rows=", nrow(residual_comparison),
        ";expected=", EXPECTED_DIRECTIONAL_FEATURES * EXPECTED_CANCERS
      )
    ),
    qc_row(
      "same_order_repeat_deterministic_all_cancers",
      nrow(sentinel_table) == EXPECTED_CANCERS &&
        max(sentinel_table$same_order_repeat_max_absolute_difference) <= REPEAT_TOLERANCE,
      paste0(
        "cancers=", nrow(sentinel_table),
        ";max=",
        format(
          max(sentinel_table$same_order_repeat_max_absolute_difference),
          scientific = TRUE
        )
      )
    ),
    qc_row(
      "reversed_order_stable_all_cancers",
      max(sentinel_table$reversed_order_median_TVD) <= ORDER_MEDIAN_TVD_TOLERANCE &&
        max(sentinel_table$reversed_order_maximum_TVD) <= ORDER_MAX_TVD_TOLERANCE &&
        min(sentinel_table$reversed_order_dominant_concordance) == 1,
      paste0(
        "max_median_TVD=",
        format(max(sentinel_table$reversed_order_median_TVD), scientific = TRUE),
        ";max_TVD=",
        format(max(sentinel_table$reversed_order_maximum_TVD), scientific = TRUE),
        ";min_concordance=",
        min(sentinel_table$reversed_order_dominant_concordance)
      )
    )
  ),
  use.names = TRUE
)

write_tsv(
  qc,
  file.path(
    OUTPUT_DIR,
    "31_cancer_specific_projection_QC_contract.tsv"
  )
)

failed_qc <- qc[
  status ==
    "FAIL"
]

report_lines <- c(
  "ARCHS4 CANCER-SPECIFIC FIXED-ANCHOR SPATIAL ECOTYPER PROJECTION V1.0.3",
  paste(
    rep(
      "=",
      96
    ),
    collapse = ""
  ),
  "",
  paste0(
    "Status: ",
    if (
      nrow(
        failed_qc
      ) ==
        0L
    ) {
      "PASS"
    } else {
      "FAIL"
    }
  ),
  paste0(
    "SpatialEcoTyper version: ",
    as.character(
      packageVersion(
        "SpatialEcoTyper"
      )
    )
  ),
  paste0(
    "Frozen W: ",
    nrow(W),
    " directional features x ",
    ncol(W),
    " components"
  ),
  paste0(
    "Stage 09 input version: ",
    input_status$script_version
  ),
  paste0(
    "Pan-cancer projection version: ",
    pan_status$script_version
  ),
  paste0(
    "Layer 1 atlas version: ",
    layer1_status$script_version
  ),
  paste0(
    "Canonical samples: ",
    EXPECTED_CANONICAL
  ),
  paste0(
    "Cross-sectional anchors for CPM/TMM: ",
    EXPECTED_ANCHOR
  ),
  paste0(
    "UQ-evaluable anchors: ",
    EXPECTED_UQ_ANCHOR
  ),
  paste0(
    "UQ robustness evaluable/non-evaluable: ",
    EXPECTED_UQ_EVALUABLE,
    "/",
    EXPECTED_UQ_NOT_EVALUABLE
  ),
  paste0(
    "Cancer-specific score spaces: ",
    EXPECTED_CANCERS
  ),
  paste0(
    "Prediction batching: one internal run per cancer and branch; nsample.per.run=",
    NSAMPLE_PER_RUN
  ),
  "",
  "COORDINATE SYSTEM",
  "-----------------",
  "Each cancer class has its own study-balanced fixed anchor. Means and weighted",
  "unbiased SDs are estimated from locked cross-sectional samples in that cancer",
  "and applied to all canonical samples from the same cancer. These scores support",
  "within-cancer analysis and are not one shared cross-cancer coordinate system.",
  "CPM and TMM use all 8,962 cross-sectional anchors. UQ scaling is re-estimated",
  "within each cancer using only the 8,799 anchors with a positive Stage 09 upper-",
  "quartile fraction. All 9,465 samples retain diagnostic UQ projections, but the",
  "163 fallback profiles are non-evaluable for formal CPM-versus-UQ robustness.",
  "Layer 1 evidence fields are attached as metadata only; this stage projects all",
  "canonical samples and does not filter on Tier A/B/C.",
  "",
  "ANCHOR COUNTS BY CANCER AND BRANCH",
  "----------------------------------",
  paste(
    capture.output(
      print(
        anchor_counts_by_branch,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "ZERO-VARIANCE MODEL GENES",
  "-------------------------",
  "A model gene with zero anchor SD within a cancer is set to zero in that cancer",
  "and removed by the official directional-feature filtering. It is never divided",
  "by zero or silently assigned an artificial variance.",
  "",
  "SCALING AND FEATURE-COVERAGE SUMMARY",
  "------------------------------------",
  paste(
    capture.output(
      print(
        scaling_summary_table,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "NORMALIZATION SENSITIVITY BY CANCER",
  "-----------------------------------",
  paste(
    capture.output(
      print(
        normalization_summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "PAN-CANCER VERSUS CANCER-SPECIFIC SCORE SHIFT",
  "---------------------------------------------",
  "This TVD describes the change of coordinate system. It must not be interpreted",
  "as a technical normalization error.",
  paste(
    capture.output(
      print(
        coordinate_shift_summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "PAN-CANCER VERSUS CANCER-SPECIFIC ADEQUACY",
  "------------------------------------------",
  paste(
    capture.output(
      print(
        adequacy_comparison_summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "RESIDUAL ATTENUATION",
  "--------------------",
  "Positive attenuation indicates that the cancer-mean residual becomes smaller",
  "after within-cancer anchoring. Persistent residuals remain candidates for",
  "within-cancer or recurrent TME structure, but are not new spatial ecotypes",
  "without independent spatial validation.",
  paste(
    capture.output(
      print(
        residual_attenuation_summary,
        row.names = FALSE
      )
    ),
    collapse = "\n"
  ),
  "",
  "INTERPRETATION BOUNDARY",
  "-----------------------",
  "Cancer-specific scores are intended for within-cancer analysis. Their absolute",
  "values must not be compared directly across cancers because each cancer uses its",
  "own fixed-anchor mean and SD. Pan-versus-cancer-specific TVD quantifies coordinate",
  "shift, not technical error. Persistent residuals remain hypotheses rather than",
  "new ecotypes without independent spatial validation.",
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
  "ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3_report.txt"
)

writeLines(
  report_lines,
  con = report_file
)

message(
  "Wrote: ",
  report_file
)

status <- list(
  script_version = SCRIPT_VERSION,
  generated_at_utc = format(
    Sys.time(),
    tz = "UTC",
    usetz = TRUE
  ),
  status = if (nrow(failed_qc) == 0L) "PASS" else "FAIL",
  SpatialEcoTyper_version = as.character(packageVersion("SpatialEcoTyper")),
  frozen_W_md5 = unname(tools::md5sum(FROZEN_W_FILE)),
  source_stage09_dir = INPUT_DIR,
  source_stage09_script_version = input_status$script_version,
  source_pan_projection_dir = PAN_INPUT_DIR,
  source_pan_projection_script_version = pan_status$script_version,
  source_layer1_dir = LAYER1_INPUT_DIR,
  source_layer1_script_version = layer1_status$script_version,
  canonical_samples = EXPECTED_CANONICAL,
  cross_sectional_anchor_samples_CPM_TMM = EXPECTED_ANCHOR,
  UQ_evaluable_samples = EXPECTED_UQ_EVALUABLE,
  UQ_non_evaluable_samples = EXPECTED_UQ_NOT_EVALUABLE,
  UQ_scaling_anchor_samples = EXPECTED_UQ_ANCHOR,
  cancer_specific_score_spaces = EXPECTED_CANCERS,
  projection_branches = names(BRANCHES),
  primary_branch = PRIMARY_BRANCH,
  nsample_per_run = NSAMPLE_PER_RUN,
  maximum_cancer_samples = max(scaling_summary_table$canonical_samples),
  minimum_cancer_branch_anchor_samples = min(anchor_counts_by_branch$anchor_samples),
  minimum_retained_directional_features =
    min(scaling_summary_table$retained_directional_features),
  minimum_component_retained_weight_fraction =
    min(scaling_summary_table$minimum_component_retained_weight_fraction),
  zero_SD_cancer_branch_gene_rows = nrow(zero_sd_table),
  persistent_residual_features =
    sum(persistent_residual_ranking$cancers_persistent > 0),
  normalization_sensitivity_rows = nrow(normalization_sample),
  TMM_evaluable_samples = n_tmm_evaluable,
  UQ_formally_evaluable_samples = n_uq_evaluable,
  UQ_formally_non_evaluable_samples = n_uq_not_evaluable,
  sentinel_same_order_max_absolute_difference =
    max(sentinel_table$same_order_repeat_max_absolute_difference),
  sentinel_reversed_order_maximum_median_TVD =
    max(sentinel_table$reversed_order_median_TVD),
  sentinel_reversed_order_maximum_TVD =
    max(sentinel_table$reversed_order_maximum_TVD),
  sentinel_reversed_order_minimum_dominant_concordance =
    min(sentinel_table$reversed_order_dominant_concordance),
  qc_failures = nrow(failed_qc)
)

status_file <- file.path(
  OUTPUT_DIR,
  "32_cancer_specific_projection_status.json"
)

jsonlite::write_json(
  status,
  path = status_file,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)

message(
  "Wrote: ",
  status_file
)

if (nrow(failed_qc) > 0L) {
  print(
    failed_qc
  )

  stop(
    "Cancer-specific projection QC contract failed.",
    call. = FALSE
  )
}

cat(
  "\n",
  paste(
    rep(
      "=",
      96
    ),
    collapse = ""
  ),
  "\n",
  "ARCHS4 CANCER-SPECIFIC FIXED-ANCHOR PROJECTION V1.0.3: PASS\n",
  paste(
    rep(
      "=",
      96
    ),
    collapse = ""
  ),
  "\n",
  "Canonical samples:                    ",
  EXPECTED_CANONICAL,
  "\n",
  "Cancer-specific score spaces:         ",
  EXPECTED_CANCERS,
  "\n",
  "Minimum retained features:            ",
  status$minimum_retained_directional_features,
  "\n",
  "Minimum component weight coverage:    ",
  format(
    status$minimum_component_retained_weight_fraction,
    digits = 6
  ),
  "\n",
  "Zero-SD cancer/branch/gene rows:       ",
  status$zero_SD_cancer_branch_gene_rows,
  "\n",
  "Persistent residual features:         ",
  status$persistent_residual_features,
  "\n",
  "Output:                               ",
  OUTPUT_DIR,
  "\n",
  sep = ""
)
