#!/usr/bin/env Rscript

# ==============================================================================
# 11_build_ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1.R
# ==============================================================================
#
# Purpose
# -------
# Convert the validated Stage 10 pan-cancer fixed-anchor SpatialEcoTyper
# projection into a Layer 1 atlas of:
#   - relative fixed-W reconstruction adequacy;
#   - normalization robustness;
#   - explicit UQ evaluability;
#   - extreme fixed-anchor z-score burden;
#   - study- and cancer-level transferability;
#   - directional residual structure;
#   - conservative evidence tiers for later analyses.
#
# This is a descriptive/QC stage. It performs no cancer-comparison hypothesis
# tests and no clinical-covariate analyses.
#
# Important interpretation
# ------------------------
# Adequacy tiers are empirical ranks within this ARCHS4 atlas. They are not
# externally calibrated validity thresholds. Continuous metrics remain primary.
#
# NonSE is a published component and is never interpreted as unexplained
# expression. The true signed residual is V - W H.
#
# UQ handling
# -----------
# Stage 10 projects all 9,465 samples in the UQ branch, but 163 profiles used a
# Stage 09 fallback because their raw upper-quartile fraction was non-positive.
# Those profiles retain diagnostic UQ projections but are formally
# UQ-non-evaluable. They cannot enter Tier A and are not counted in formal UQ
# robustness denominators.
#
# ==============================================================================

options(
  stringsAsFactors = FALSE,
  width = 260,
  warn = 1
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.", call. = FALSE)
}

# ==============================================================================
# 1. Paths and constants
# ==============================================================================

SCRIPT_VERSION <- "1.0.2.1"
EXPECTED_PROJECTION_VERSION <- "1.0.4"

DATA_ROOT <- Sys.getenv(
  "ARCHS4_V26_ROOT",
  unset = "/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6"
)

INPUT_DIR <- Sys.getenv(
  "ARCHS4_STAGE11_INPUT_DIR",
  unset = file.path(
    DATA_ROOT,
    "11_pan_cancer_fixed_anchor_projection_v1_0_4"
  )
)

OUTPUT_DIR <- Sys.getenv(
  "ARCHS4_STAGE11_OUTPUT_DIR",
  unset = file.path(
    DATA_ROOT,
    "12_layer1_adequacy_robustness_atlas_v1_0_2_1"
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

  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(SOURCE_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
}

SCORE_FILE <- file.path(
  INPUT_DIR,
  "17_primary_CPM_fraction_scores_with_metadata.tsv.gz"
)

ADEQUACY_FILE <- file.path(
  INPUT_DIR,
  "19_primary_CPM_fixed_W_model_adequacy.tsv.gz"
)

SENSITIVITY_FILE <- file.path(
  INPUT_DIR,
  "22_normalization_sensitivity_by_sample.tsv.gz"
)

SCALED_INPUT_FILE <- file.path(
  INPUT_DIR,
  "06_fixed_anchor_scaled_input_CPM.rds"
)

RECONSTRUCTION_FILE <- file.path(
  INPUT_DIR,
  "21_primary_CPM_directional_reconstruction_and_residual.rds"
)

PROJECTION_QC_FILE <- file.path(
  INPUT_DIR,
  "27_pan_cancer_projection_QC_contract.tsv"
)

PROJECTION_STATUS_FILE <- file.path(
  INPUT_DIR,
  "28_pan_cancer_projection_status.json"
)

COMPONENTS <- c(
  "NonSE",
  paste0("SE", 1:9)
)

EXPECTED_SAMPLES <- 9465L
EXPECTED_CANCERS <- 18L
EXPECTED_GENES <- 499L
EXPECTED_DIRECTIONAL_FEATURES <- 500L
EXPECTED_UQ_EVALUABLE <- 9302L
EXPECTED_UQ_NOT_EVALUABLE <- 163L

# Transparent normalization-robustness thresholds inherited from the validated
# previous Layer 1 analysis. They are descriptive policy thresholds, not model
# fitting parameters.
TMM_STABLE_TVD <- 0.10
UQ_STABLE_TVD <- 0.15
TMM_HIGH_SENSITIVITY_TVD <- 0.25
UQ_HIGH_SENSITIVITY_TVD <- 0.35

# Extreme fixed-anchor z-score diagnostics.
EXTREME_Z_THRESHOLD <- 8
SEVERE_Z_THRESHOLD <- 10
SEVERE_Z_FEATURE_COUNT <- 5L

# Study summaries below this sample count remain descriptive but are not used
# for study-outlier designation.
MIN_STUDY_N_FOR_OUTLIER <- 5L
MIN_STUDIES_FOR_WITHIN_CANCER_OUTLIER <- 3L
STUDY_ROBUST_Z_THRESHOLD <- -2.5

# Residual-feature summaries.
TOP_RESIDUAL_FEATURES <- 30L
TOP_STUDIES_FOR_TABLE <- 50L

DPI <- 600
NATURE_W_SINGLE_MM <- 89
NATURE_W_DOUBLE_MM <- 183
NATURE_FONT <- "Helvetica"

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

empirical_percentile <- function(
    x,
    higher_is_better = TRUE
) {
  n_valid <- sum(
    is.finite(x)
  )

  if (n_valid < 2L) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  ranks <- frank(
    x,
    ties.method = "average",
    na.last = "keep"
  )

  percentile <- (
    ranks -
      1
  ) / (
    n_valid -
      1
  )

  if (!higher_is_better) {
    percentile <- 1 -
      percentile
  }

  percentile
}

robust_z <- function(x) {
  x <- as.numeric(x)

  centre <- stats::median(
    x,
    na.rm = TRUE
  )

  scale_value <- stats::mad(
    x,
    center = centre,
    constant = 1.4826,
    na.rm = TRUE
  )

  if (
    !is.finite(scale_value) ||
      scale_value <= 0
  ) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  (
    x -
      centre
  ) /
    scale_value
}

format_number <- function(
    x,
    digits = 3
) {
  formatC(
    x,
    format = "f",
    digits = digits
  )
}

finite_median <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  stats::median(x[is.finite(x)])
}

finite_quantile <- function(x, probability) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) {
    return(NA_real_)
  }
  as.numeric(
    stats::quantile(
      x[is.finite(x)],
      probs = probability,
      names = FALSE
    )
  )
}

logical_mean_na_safe <- function(x) {
  x <- as.logical(x)
  if (!any(!is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
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
        colour = "black",
        size = base_size
      ),
      axis.title = element_text(
        colour = "black",
        size = base_size
      ),
      axis.line = element_line(
        linewidth = 0.3,
        colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.3,
        colour = "black"
      ),
      axis.ticks.length = grid::unit(
        2,
        "pt"
      ),
      plot.title = element_text(
        face = "plain",
        colour = "black",
        size = base_size,
        hjust = 0
      ),
      plot.subtitle = element_text(
        face = "plain",
        colour = "grey20",
        size = base_size - 1,
        hjust = 0
      ),
      plot.caption = element_text(
        face = "plain",
        colour = "grey20",
        size = base_size - 1,
        hjust = 0
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "plain",
        colour = "black",
        size = base_size
      ),
      legend.title = element_text(
        size = base_size
      ),
      legend.text = element_text(
        size = base_size - 1
      ),
      legend.key.size = grid::unit(
        0.30,
        "cm"
      ),
      plot.margin = margin(
        3,
        3,
        3,
        3
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

matrix_sample_metrics <- function(
    matrix_object
) {
  abs_matrix <- abs(
    matrix_object
  )

  data.table(
    sample_id = colnames(
      matrix_object
    ),
    maximum_absolute_fixed_anchor_z = apply(
      abs_matrix,
      2L,
      max
    ),
    q99_absolute_fixed_anchor_z = apply(
      abs_matrix,
      2L,
      stats::quantile,
      probs = 0.99,
      names = FALSE
    ),
    n_features_abs_z_ge_6 = colSums(
      abs_matrix >=
        6
    ),
    n_features_abs_z_ge_8 = colSums(
      abs_matrix >=
        EXTREME_Z_THRESHOLD
    ),
    n_features_abs_z_ge_10 = colSums(
      abs_matrix >=
        SEVERE_Z_THRESHOLD
    )
  )
}

matrix_gene_metrics <- function(
    matrix_object
) {
  abs_matrix <- abs(
    matrix_object
  )

  data.table(
    model_gene = rownames(
      matrix_object
    ),
    maximum_absolute_fixed_anchor_z = apply(
      abs_matrix,
      1L,
      max
    ),
    q99_absolute_fixed_anchor_z = apply(
      abs_matrix,
      1L,
      stats::quantile,
      probs = 0.99,
      names = FALSE
    ),
    fraction_samples_abs_z_ge_6 = rowMeans(
      abs_matrix >=
        6
    ),
    fraction_samples_abs_z_ge_8 = rowMeans(
      abs_matrix >=
        EXTREME_Z_THRESHOLD
    ),
    fraction_samples_abs_z_ge_10 = rowMeans(
      abs_matrix >=
        SEVERE_Z_THRESHOLD
    )
  )
}

# ==============================================================================
# 3. Self-tests
# ==============================================================================

run_self_tests <- function() {
  x <- c(1, 2, 3, 4)

  p_good <- empirical_percentile(x, higher_is_better = TRUE)
  p_reverse <- empirical_percentile(x, higher_is_better = FALSE)

  expected_good <- c(0, 1 / 3, 2 / 3, 1)
  expected_reverse <- c(1, 2 / 3, 1 / 3, 0)

  stopifnot(
    max(abs(p_good - expected_good)) <= 1e-15,
    max(abs(p_reverse - expected_reverse)) <= 1e-15,
    max(abs(p_good + p_reverse - 1)) <= 1e-15
  )

  z <- robust_z(c(1, 2, 3, 4, 100))
  stopifnot(length(z) == 5L, is.finite(z[[5L]]), z[[5L]] > 10)

  test_matrix <- matrix(
    c(-1, 2, -8, 9, 0, 10),
    nrow = 2L,
    dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
  )

  metrics <- matrix_sample_metrics(test_matrix)
  stopifnot(
    nrow(metrics) == 3L,
    metrics[sample_id == "S3", maximum_absolute_fixed_anchor_z] == 10
  )

  # UQ-non-evaluable values must not enter formal UQ summaries.
  stopifnot(
    is.na(finite_median(c(NA_real_, NA_real_))),
    finite_median(c(NA_real_, 0.1, 0.3)) == 0.2,
    logical_mean_na_safe(c(NA, TRUE, FALSE)) == 0.5
  )

  cat("ALL LAYER 1 ADEQUACY/ROBUSTNESS ATLAS V1.0.2.1 SELF-TESTS PASSED\n")
}

command_args <- commandArgs(trailingOnly = TRUE)

if ("--run-self-tests" %in% command_args) {
  run_self_tests()
  quit(save = "no", status = 0)
}

# ==============================================================================
# 4. Load and validate the locked projection
# ==============================================================================

required_files <- c(
  SCORE_FILE,
  ADEQUACY_FILE,
  SENSITIVITY_FILE,
  SCALED_INPUT_FILE,
  RECONSTRUCTION_FILE,
  PROJECTION_QC_FILE,
  PROJECTION_STATUS_FILE
)

missing_files <- required_files[
  !file.exists(required_files) |
    file.info(required_files)$size <= 0
]

if (length(missing_files) > 0L) {
  stop(
    "Missing or empty required files:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

projection_status <- jsonlite::read_json(
  PROJECTION_STATUS_FILE,
  simplifyVector = TRUE
)

projection_qc <- fread(
  PROJECTION_QC_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

if (
  !identical(projection_status$status, "PASS") ||
    projection_status$qc_failures != 0 ||
    !identical(projection_status$script_version, EXPECTED_PROJECTION_VERSION) ||
    projection_status$canonical_samples != EXPECTED_SAMPLES ||
    projection_status$cancer_classes != EXPECTED_CANCERS ||
    projection_status$stage09_input_model_genes != EXPECTED_GENES ||
    projection_status$directional_features != EXPECTED_DIRECTIONAL_FEATURES ||
    projection_status$UQ_evaluable_samples != EXPECTED_UQ_EVALUABLE ||
    projection_status$UQ_non_evaluable_samples != EXPECTED_UQ_NOT_EVALUABLE
) {
  stop(
    "Stage 10 projection is not the required frozen v1.0.4 PASS state.",
    call. = FALSE
  )
}

check_columns(projection_qc, c("check", "status", "detail"), "Projection QC")
if (nrow(projection_qc) == 0L || any(projection_qc$status != "PASS")) {
  stop("Stage 10 projection QC contract contains non-PASS rows.", call. = FALSE)
}

scores <- fread(
  SCORE_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

adequacy <- fread(
  ADEQUACY_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

sensitivity_long <- fread(
  SENSITIVITY_FILE,
  sep = "\t",
  header = TRUE,
  quote = "",
  showProgress = FALSE
)

scaled_input <- readRDS(SCALED_INPUT_FILE)
reconstruction <- readRDS(RECONSTRUCTION_FILE)

check_columns(
  scores,
  c(
    "sample_id",
    "primary_gse_id",
    "provisional_cohort",
    "analysis_cancer_class",
    "stage08_tier",
    "stage07_expression_quality_tier",
    COMPONENTS
  ),
  "Primary score table"
)

check_columns(
  adequacy,
  c(
    "sample_id",
    "primary_gse_id",
    "provisional_cohort",
    "analysis_cancer_class",
    "normalized_KL_per_input_mass",
    "reconstruction_RMSE",
    "reconstruction_MAE",
    "reconstruction_Pearson",
    "reconstruction_cosine",
    "reconstruction_R2",
    "normalized_score_entropy",
    "dominant_component",
    "dominant_fraction",
    "top1_minus_top2",
    "dominance_class"
  ),
  "Primary adequacy table"
)

check_columns(
  sensitivity_long,
  c(
    "sample_id",
    "comparison",
    "comparison_evaluable",
    "non_evaluable_reason",
    "diagnostic_total_variation_distance",
    "diagnostic_maximum_component_absolute_difference",
    "diagnostic_aitchison_distance",
    "diagnostic_dominant_component_concordant",
    "total_variation_distance",
    "maximum_component_absolute_difference",
    "aitchison_distance",
    "dominant_primary",
    "dominant_alternative",
    "dominant_component_concordant",
    "UQ_projection_evaluable"
  ),
  "Normalization-sensitivity table"
)

if (
  nrow(scores) != EXPECTED_SAMPLES ||
    nrow(adequacy) != EXPECTED_SAMPLES ||
    uniqueN(scores$analysis_cancer_class) != EXPECTED_CANCERS ||
    anyDuplicated(scores$sample_id) ||
    anyDuplicated(adequacy$sample_id)
) {
  stop(
    "Expected ", EXPECTED_SAMPLES,
    " unique score and adequacy rows across ", EXPECTED_CANCERS,
    " cancer classes.",
    call. = FALSE
  )
}

if (!identical(scores$sample_id, adequacy$sample_id)) {
  stop("Score and adequacy sample orders differ.", call. = FALSE)
}

if (
  nrow(sensitivity_long) != 2L * EXPECTED_SAMPLES ||
    any(sensitivity_long[, .N, by = sample_id]$N != 2L) ||
    !setequal(
      unique(sensitivity_long$comparison),
      c("CPM_vs_TMM_CPM", "CPM_vs_UQ_CPM")
    )
) {
  stop("Expected exactly two normalization-sensitivity rows per sample.", call. = FALSE)
}

sensitivity_long[, comparison_evaluable := as.logical(comparison_evaluable)]
sensitivity_long[, dominant_component_concordant := as.logical(dominant_component_concordant)]
sensitivity_long[, diagnostic_dominant_component_concordant := as.logical(diagnostic_dominant_component_concordant)]

n_tmm_evaluable <- sensitivity_long[
  comparison == "CPM_vs_TMM_CPM" & comparison_evaluable == TRUE,
  .N
]
n_uq_evaluable <- sensitivity_long[
  comparison == "CPM_vs_UQ_CPM" & comparison_evaluable == TRUE,
  .N
]
n_uq_not_evaluable <- sensitivity_long[
  comparison == "CPM_vs_UQ_CPM" & comparison_evaluable == FALSE,
  .N
]

if (
  n_tmm_evaluable != EXPECTED_SAMPLES ||
    n_uq_evaluable != EXPECTED_UQ_EVALUABLE ||
    n_uq_not_evaluable != EXPECTED_UQ_NOT_EVALUABLE
) {
  stop(
    "Unexpected normalization-evaluability partition: TMM=", n_tmm_evaluable,
    ";UQ_evaluable=", n_uq_evaluable,
    ";UQ_non_evaluable=", n_uq_not_evaluable,
    call. = FALSE
  )
}

formal_columns <- c(
  "total_variation_distance",
  "maximum_component_absolute_difference",
  "aitchison_distance"
)

diagnostic_columns <- c(
  "diagnostic_total_variation_distance",
  "diagnostic_maximum_component_absolute_difference",
  "diagnostic_aitchison_distance"
)

if (
  any(!is.finite(as.matrix(sensitivity_long[, ..diagnostic_columns]))) ||
    any(!is.finite(as.matrix(sensitivity_long[comparison_evaluable == TRUE, ..formal_columns]))) ||
    any(!is.na(as.matrix(sensitivity_long[comparison_evaluable == FALSE, ..formal_columns]))) ||
    any(!is.na(sensitivity_long[comparison_evaluable == FALSE, dominant_component_concordant]))
) {
  stop("Formal/diagnostic normalization-sensitivity NA contract failed.", call. = FALSE)
}

if (
  !is.matrix(scaled_input) ||
    nrow(scaled_input) != EXPECTED_GENES ||
    ncol(scaled_input) != EXPECTED_SAMPLES ||
    !identical(colnames(scaled_input), scores$sample_id) ||
    any(!is.finite(scaled_input))
) {
  stop("Primary fixed-anchor scaled matrix failed validation.", call. = FALSE)
}

required_reconstruction_items <- c(
  "directional_feature_metadata",
  "sample_ids",
  "observed_directional_features",
  "reconstructed_directional_features",
  "signed_directional_residual"
)

missing_reconstruction_items <- setdiff(
  required_reconstruction_items,
  names(reconstruction)
)

if (length(missing_reconstruction_items) > 0L) {
  stop(
    "Reconstruction RDS lacks: ",
    paste(missing_reconstruction_items, collapse = ", "),
    call. = FALSE
  )
}

residual_matrix <- reconstruction$signed_directional_residual

if (
  !is.matrix(residual_matrix) ||
    nrow(residual_matrix) != EXPECTED_DIRECTIONAL_FEATURES ||
    ncol(residual_matrix) != EXPECTED_SAMPLES ||
    !identical(colnames(residual_matrix), scores$sample_id) ||
    any(!is.finite(residual_matrix))
) {
  stop("Signed directional residual matrix failed validation.", call. = FALSE)
}

input_provenance <- data.table(
  role = c(
    "projection_scores",
    "projection_adequacy",
    "projection_sensitivity",
    "fixed_anchor_scaled_input",
    "directional_reconstruction_residual",
    "projection_QC_contract",
    "projection_status"
  ),
  path = c(
    SCORE_FILE,
    ADEQUACY_FILE,
    SENSITIVITY_FILE,
    SCALED_INPUT_FILE,
    RECONSTRUCTION_FILE,
    PROJECTION_QC_FILE,
    PROJECTION_STATUS_FILE
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

policy_thresholds <- data.table(
  policy = c(
    "TMM_stable_TVD",
    "UQ_stable_TVD",
    "TMM_high_sensitivity_TVD",
    "UQ_high_sensitivity_TVD",
    "extreme_fixed_anchor_z",
    "severe_fixed_anchor_z",
    "severe_z_feature_count",
    "Tier_A_min_adequacy_percentile",
    "Tier_B_min_adequacy_percentile",
    "minimum_study_n_for_outlier",
    "minimum_studies_for_within_cancer_outlier",
    "study_robust_z_threshold"
  ),
  value = c(
    TMM_STABLE_TVD,
    UQ_STABLE_TVD,
    TMM_HIGH_SENSITIVITY_TVD,
    UQ_HIGH_SENSITIVITY_TVD,
    EXTREME_Z_THRESHOLD,
    SEVERE_Z_THRESHOLD,
    SEVERE_Z_FEATURE_COUNT,
    0.25,
    0.10,
    MIN_STUDY_N_FOR_OUTLIER,
    MIN_STUDIES_FOR_WITHIN_CANCER_OUTLIER,
    STUDY_ROBUST_Z_THRESHOLD
  ),
  interpretation = c(
    "TMM stable only when TVD is at or below this threshold and dominant component is unchanged",
    "UQ stable only among formally UQ-evaluable profiles when TVD is at or below this threshold and dominant component is unchanged",
    "TMM TVD at or above this threshold contributes to the highly-sensitive flag",
    "Formally evaluable UQ TVD at or above this threshold contributes to the highly-sensitive flag",
    "Absolute fixed-anchor z threshold used for burden counts",
    "Absolute fixed-anchor z threshold used for severe single-feature burden",
    "Number of features at abs(z)>=8 sufficient for the severe-z flag",
    "Minimum empirical adequacy percentile for Tier A",
    "Minimum empirical adequacy percentile for Tier B",
    "Studies below this sample count remain descriptive only",
    "Minimum number of studies required to compute within-cancer robust-z outliers",
    "Within-cancer study median-adequacy robust-z threshold for poor-fit outlier designation"
  )
)

write_tsv(
  policy_thresholds,
  file.path(OUTPUT_DIR, "00B_layer1_policy_thresholds.tsv")
)

# ==============================================================================
# 5. Build the sample-level adequacy atlas
# ==============================================================================

message(
  "Building empirical reconstruction-adequacy percentiles."
)

sample_atlas <- copy(
  adequacy
)

layer1_metadata_candidates <- c(
  "stage08_tier",
  "stage07_expression_quality_tier",
  "expression_qc_warning_codes",
  "model_genes_detected",
  "overall_model_detected_weight_fraction",
  "min_component_detected_weight_fraction",
  "library_sum_archs4_pseudocounts"
)

layer1_metadata_columns <- intersect(
  layer1_metadata_candidates,
  names(scores)
)

if (length(layer1_metadata_columns) > 0L) {
  layer1_metadata <- scores[
    ,
    c("sample_id", layer1_metadata_columns),
    with = FALSE
  ]

  sample_atlas <- merge(
    sample_atlas,
    layer1_metadata,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  sample_atlas <- sample_atlas[
    match(scores$sample_id, sample_id)
  ]
}

sample_atlas[
  ,
  fit_percentile_low_KL :=
    empirical_percentile(
      normalized_KL_per_input_mass,
      higher_is_better = FALSE
    )
]

sample_atlas[
  ,
  fit_percentile_high_Pearson :=
    empirical_percentile(
      reconstruction_Pearson,
      higher_is_better = TRUE
    )
]

sample_atlas[
  ,
  fit_percentile_high_cosine :=
    empirical_percentile(
      reconstruction_cosine,
      higher_is_better = TRUE
    )
]

sample_atlas[
  ,
  fit_percentile_high_R2 :=
    empirical_percentile(
      reconstruction_R2,
      higher_is_better = TRUE
    )
]

sample_atlas[
  ,
  adequacy_composite_percentile :=
    rowMeans(
      cbind(
        fit_percentile_low_KL,
        fit_percentile_high_Pearson,
        fit_percentile_high_cosine,
        fit_percentile_high_R2
      )
    )
]

sample_atlas[
  ,
  adequacy_percentile_rank :=
    empirical_percentile(
      adequacy_composite_percentile,
      higher_is_better = TRUE
    )
]

sample_atlas[
  ,
  adequacy_percentile_within_cancer :=
    empirical_percentile(
      adequacy_composite_percentile,
      higher_is_better = TRUE
    ),
  by = analysis_cancer_class
]

sample_atlas[
  ,
  relative_fit_tier :=
    fifelse(
      adequacy_percentile_rank >=
        0.75,
      "Q4_relatively_higher_fit",
      fifelse(
        adequacy_percentile_rank <
          0.25,
        "Q1_relatively_lower_fit",
        "Q2_Q3_intermediate_fit"
      )
    )
]

sample_atlas[
  ,
  global_lowest_10pct_fit :=
    adequacy_percentile_rank <=
      0.10
]

sample_atlas[
  ,
  within_cancer_lowest_10pct_fit :=
    adequacy_percentile_within_cancer <=
      0.10
]

# ==============================================================================
# 6. Add normalization sensitivity with explicit UQ evaluability
# ==============================================================================

message(
  "Reshaping normalization sensitivity and assigning explicit evaluability states."
)

sensitivity_wide <- dcast(
  sensitivity_long[
    ,
    .(
      sample_id,
      comparison,
      comparison_evaluable,
      non_evaluable_reason,
      diagnostic_total_variation_distance,
      diagnostic_maximum_component_absolute_difference,
      diagnostic_aitchison_distance,
      diagnostic_dominant_component_concordant,
      total_variation_distance,
      maximum_component_absolute_difference,
      aitchison_distance,
      dominant_component_concordant
    )
  ],
  sample_id ~ comparison,
  value.var = c(
    "comparison_evaluable",
    "non_evaluable_reason",
    "diagnostic_total_variation_distance",
    "diagnostic_maximum_component_absolute_difference",
    "diagnostic_aitchison_distance",
    "diagnostic_dominant_component_concordant",
    "total_variation_distance",
    "maximum_component_absolute_difference",
    "aitchison_distance",
    "dominant_component_concordant"
  )
)

rename_map <- c(
  "comparison_evaluable_CPM_vs_TMM_CPM" = "TMM_evaluable",
  "comparison_evaluable_CPM_vs_UQ_CPM" = "UQ_evaluable",
  "non_evaluable_reason_CPM_vs_TMM_CPM" = "TMM_non_evaluable_reason",
  "non_evaluable_reason_CPM_vs_UQ_CPM" = "UQ_non_evaluable_reason",
  "diagnostic_total_variation_distance_CPM_vs_TMM_CPM" = "TMM_diagnostic_total_variation_distance",
  "diagnostic_total_variation_distance_CPM_vs_UQ_CPM" = "UQ_diagnostic_total_variation_distance",
  "diagnostic_maximum_component_absolute_difference_CPM_vs_TMM_CPM" = "TMM_diagnostic_maximum_component_absolute_difference",
  "diagnostic_maximum_component_absolute_difference_CPM_vs_UQ_CPM" = "UQ_diagnostic_maximum_component_absolute_difference",
  "diagnostic_aitchison_distance_CPM_vs_TMM_CPM" = "TMM_diagnostic_aitchison_distance",
  "diagnostic_aitchison_distance_CPM_vs_UQ_CPM" = "UQ_diagnostic_aitchison_distance",
  "diagnostic_dominant_component_concordant_CPM_vs_TMM_CPM" = "TMM_diagnostic_dominant_component_concordant",
  "diagnostic_dominant_component_concordant_CPM_vs_UQ_CPM" = "UQ_diagnostic_dominant_component_concordant",
  "total_variation_distance_CPM_vs_TMM_CPM" = "TMM_total_variation_distance",
  "total_variation_distance_CPM_vs_UQ_CPM" = "UQ_total_variation_distance",
  "maximum_component_absolute_difference_CPM_vs_TMM_CPM" = "TMM_maximum_component_absolute_difference",
  "maximum_component_absolute_difference_CPM_vs_UQ_CPM" = "UQ_maximum_component_absolute_difference",
  "aitchison_distance_CPM_vs_TMM_CPM" = "TMM_aitchison_distance",
  "aitchison_distance_CPM_vs_UQ_CPM" = "UQ_aitchison_distance",
  "dominant_component_concordant_CPM_vs_TMM_CPM" = "TMM_dominant_component_concordant",
  "dominant_component_concordant_CPM_vs_UQ_CPM" = "UQ_dominant_component_concordant"
)

for (old_name in names(rename_map)) {
  if (!old_name %in% names(sensitivity_wide)) {
    stop("Missing reshaped sensitivity column: ", old_name, call. = FALSE)
  }
  setnames(sensitivity_wide, old_name, rename_map[[old_name]])
}

sensitivity_wide[, TMM_evaluable := as.logical(TMM_evaluable)]
sensitivity_wide[, UQ_evaluable := as.logical(UQ_evaluable)]
sensitivity_wide[, TMM_dominant_component_concordant := as.logical(TMM_dominant_component_concordant)]
sensitivity_wide[, UQ_dominant_component_concordant := as.logical(UQ_dominant_component_concordant)]

if (
  anyNA(sensitivity_wide$TMM_evaluable) ||
    any(!sensitivity_wide$TMM_evaluable) ||
    sum(sensitivity_wide$UQ_evaluable) != EXPECTED_UQ_EVALUABLE ||
    sum(!sensitivity_wide$UQ_evaluable) != EXPECTED_UQ_NOT_EVALUABLE
) {
  stop("Wide normalization-evaluability partition failed.", call. = FALSE)
}

sensitivity_wide[
  ,
  TMM_stable :=
    TMM_total_variation_distance <= TMM_STABLE_TVD &
      TMM_dominant_component_concordant
]

sensitivity_wide[
  ,
  UQ_stable :=
    UQ_evaluable &
      !is.na(UQ_total_variation_distance) &
      UQ_total_variation_distance <= UQ_STABLE_TVD &
      UQ_dominant_component_concordant == TRUE
]

sensitivity_wide[
  ,
  UQ_robustness_state := fifelse(
    !UQ_evaluable,
    "not_evaluable",
    fifelse(
      UQ_stable,
      "stable",
      fifelse(
        UQ_total_variation_distance >= UQ_HIGH_SENSITIVITY_TVD |
          UQ_dominant_component_concordant == FALSE,
        "unstable",
        "intermediate"
      )
    )
  )
]

sensitivity_wide[
  ,
  maximum_normalization_TVD := fifelse(
    UQ_evaluable,
    pmax(TMM_total_variation_distance, UQ_total_variation_distance),
    TMM_total_variation_distance
  )
]

sensitivity_wide[
  ,
  highly_normalization_sensitive :=
    TMM_total_variation_distance >= TMM_HIGH_SENSITIVITY_TVD |
      (
        UQ_evaluable &
          UQ_total_variation_distance >= UQ_HIGH_SENSITIVITY_TVD
      ) |
      (
        UQ_evaluable &
          TMM_dominant_component_concordant == FALSE &
          UQ_dominant_component_concordant == FALSE
      )
]

sensitivity_wide[
  ,
  normalization_robustness_tier := fifelse(
    TMM_stable & UQ_evaluable & UQ_stable,
    "robust_TMM_and_UQ",
    fifelse(
      highly_normalization_sensitive,
      "highly_sensitive",
      fifelse(
        TMM_stable & !UQ_evaluable,
        "robust_TMM_UQ_not_evaluable",
        fifelse(
          TMM_stable,
          "robust_TMM_only_UQ_not_stable",
          fifelse(
            !UQ_evaluable,
            "intermediate_TMM_UQ_not_evaluable",
            "intermediate_sensitivity"
          )
        )
      )
    )
  )
]

sample_atlas <- merge(
  sample_atlas,
  sensitivity_wide,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

sample_atlas <- sample_atlas[match(scores$sample_id, sample_id)]

if (
  anyNA(sample_atlas$normalization_robustness_tier) ||
    any(!is.finite(sample_atlas$TMM_total_variation_distance)) ||
    any(!is.finite(sample_atlas$maximum_normalization_TVD)) ||
    any(!is.na(sample_atlas[UQ_evaluable == FALSE, UQ_total_variation_distance])) ||
    any(!is.finite(sample_atlas[UQ_evaluable == TRUE, UQ_total_variation_distance]))
) {
  stop("At least one sample violates the normalization-sensitivity contract.", call. = FALSE)
}

# ==============================================================================
# 7. Add extreme-z diagnostics
# ==============================================================================

message(
  "Calculating primary fixed-anchor z-score burden."
)

z_sample_metrics <- matrix_sample_metrics(
  scaled_input
)

z_gene_metrics <- matrix_gene_metrics(
  scaled_input
)

z_sample_metrics[
  ,
  severe_extreme_z_flag :=
    maximum_absolute_fixed_anchor_z >=
      SEVERE_Z_THRESHOLD |
      n_features_abs_z_ge_8 >=
        SEVERE_Z_FEATURE_COUNT
]

z_sample_metrics[
  ,
  extreme_z_burden_percentile :=
    empirical_percentile(
      maximum_absolute_fixed_anchor_z,
      higher_is_better = TRUE
    )
]

sample_atlas <- merge(
  sample_atlas,
  z_sample_metrics,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

sample_atlas <- sample_atlas[
  match(
    scores$sample_id,
    sample_id
  )
]

# ==============================================================================
# 8. Conservative evidence tiers for later use
# ==============================================================================

sample_atlas[
  ,
  projection_evidence_tier := fifelse(
    adequacy_percentile_rank >= 0.25 &
      TMM_stable &
      UQ_evaluable &
      UQ_stable &
      !severe_extreme_z_flag,
    "Tier_A_fit_and_both_normalizations",
    fifelse(
      adequacy_percentile_rank >= 0.10 &
        TMM_stable &
        !highly_normalization_sensitive &
        (
          !UQ_evaluable |
            UQ_robustness_state %in% c("stable", "intermediate")
        ) &
        !severe_extreme_z_flag,
      "Tier_B_fit_and_TMM_supported",
      "Tier_C_exploratory_or_sensitive"
    )
  )
]

sample_atlas[
  ,
  projection_evidence_tier_detail := fifelse(
    projection_evidence_tier == "Tier_A_fit_and_both_normalizations",
    "Tier_A_UQ_evaluable_and_stable",
    fifelse(
      projection_evidence_tier == "Tier_B_fit_and_TMM_supported" & !UQ_evaluable,
      "Tier_B_TMM_supported_UQ_not_evaluable",
      fifelse(
        projection_evidence_tier == "Tier_B_fit_and_TMM_supported",
        paste0(
          "Tier_B_TMM_supported_UQ_",
          UQ_robustness_state
        ),
        "Tier_C_exploratory_or_sensitive"
      )
    )
  )
]

sample_atlas[
  ,
  recommended_use := fifelse(
    projection_evidence_tier == "Tier_A_fit_and_both_normalizations",
    "primary_layer2_evidence",
    fifelse(
      projection_evidence_tier == "Tier_B_fit_and_TMM_supported",
      "sensitivity_supported_evidence",
      "exploratory_only"
    )
  )
]

sample_atlas[
  ,
  technical_caution_reasons := vapply(
    seq_len(.N),
    function(i) {
      reasons <- character()

      if (global_lowest_10pct_fit[[i]]) {
        reasons <- c(reasons, "global_low_fit_tail")
      }
      if (within_cancer_lowest_10pct_fit[[i]]) {
        reasons <- c(reasons, "within_cancer_low_fit_tail")
      }
      if (highly_normalization_sensitive[[i]]) {
        reasons <- c(reasons, "normalization_sensitive")
      }
      if (!UQ_evaluable[[i]]) {
        reasons <- c(reasons, "UQ_robustness_not_evaluable")
      }
      if (severe_extreme_z_flag[[i]]) {
        reasons <- c(reasons, "extreme_fixed_anchor_z")
      }

      if (length(reasons) == 0L) {
        "none"
      } else {
        paste(reasons, collapse = ";")
      }
    },
    FUN.VALUE = character(1)
  )
]

if (any(sample_atlas[UQ_evaluable == FALSE, projection_evidence_tier] == "Tier_A_fit_and_both_normalizations")) {
  stop("UQ-non-evaluable samples were incorrectly assigned to Tier A.", call. = FALSE)
}

write_tsv(
  sample_atlas,
  file.path(OUTPUT_DIR, "01_sample_level_adequacy_robustness_atlas.tsv.gz"),
  compress = TRUE
)

write_tsv(
  z_gene_metrics[
    order(-maximum_absolute_fixed_anchor_z)
  ],
  file.path(OUTPUT_DIR, "02_model_gene_extreme_z_diagnostics.tsv")
)

# ==============================================================================
# 9. Study-level and cancer-level summaries
# ==============================================================================

message("Building study- and cancer-level transferability summaries.")

study_by_cohort_summary <- sample_atlas[
  ,
  .(
    n_samples = .N,
    n_UQ_evaluable = sum(UQ_evaluable),
    n_UQ_not_evaluable = sum(!UQ_evaluable),
    UQ_evaluable_fraction = mean(UQ_evaluable),
    median_adequacy_composite = median(adequacy_composite_percentile),
    q10_adequacy_composite = quantile(adequacy_composite_percentile, 0.10, names = FALSE),
    median_normalized_KL = median(normalized_KL_per_input_mass),
    median_reconstruction_Pearson = median(reconstruction_Pearson),
    median_reconstruction_cosine = median(reconstruction_cosine),
    median_reconstruction_R2 = median(reconstruction_R2),
    median_TMM_TVD = finite_median(TMM_total_variation_distance),
    median_UQ_TVD_evaluable = finite_median(UQ_total_variation_distance),
    TMM_dominant_concordance = logical_mean_na_safe(TMM_dominant_component_concordant),
    UQ_dominant_concordance_evaluable = logical_mean_na_safe(UQ_dominant_component_concordant),
    TMM_stable_fraction = mean(TMM_stable),
    UQ_stable_fraction_evaluable = if (sum(UQ_evaluable) > 0L) mean(UQ_stable[UQ_evaluable]) else NA_real_,
    Tier_A_fraction = mean(projection_evidence_tier == "Tier_A_fit_and_both_normalizations"),
    Tier_B_fraction = mean(projection_evidence_tier == "Tier_B_fit_and_TMM_supported"),
    Tier_C_fraction = mean(projection_evidence_tier == "Tier_C_exploratory_or_sensitive"),
    global_low_fit_tail_fraction = mean(global_lowest_10pct_fit),
    highly_sensitive_fraction = mean(highly_normalization_sensitive),
    severe_extreme_z_fraction = mean(severe_extreme_z_flag)
  ),
  by = .(
    analysis_cancer_class,
    primary_gse_id,
    provisional_cohort
  )
]

setorder(
  study_by_cohort_summary,
  analysis_cancer_class,
  primary_gse_id,
  provisional_cohort
)

write_tsv(
  study_by_cohort_summary,
  file.path(
    OUTPUT_DIR,
    "03A_study_by_provisional_cohort_transferability_summary.tsv"
  )
)

study_summary <- sample_atlas[
  ,
  .(
    n_samples = .N,
    n_Gold = sum(provisional_cohort == "GOLD"),
    n_Silver = sum(provisional_cohort == "SILVER"),
    Gold_fraction = mean(provisional_cohort == "GOLD"),
    Silver_fraction = mean(provisional_cohort == "SILVER"),
    provisional_cohorts_present = paste(
      sort(unique(as.character(provisional_cohort))),
      collapse = ";"
    ),
    n_UQ_evaluable = sum(UQ_evaluable),
    n_UQ_not_evaluable = sum(!UQ_evaluable),
    UQ_evaluable_fraction = mean(UQ_evaluable),
    median_adequacy_composite = median(adequacy_composite_percentile),
    q10_adequacy_composite = quantile(adequacy_composite_percentile, 0.10, names = FALSE),
    median_normalized_KL = median(normalized_KL_per_input_mass),
    median_reconstruction_Pearson = median(reconstruction_Pearson),
    median_reconstruction_cosine = median(reconstruction_cosine),
    median_reconstruction_R2 = median(reconstruction_R2),
    median_TMM_TVD = finite_median(TMM_total_variation_distance),
    median_UQ_TVD_evaluable = finite_median(UQ_total_variation_distance),
    TMM_dominant_concordance = logical_mean_na_safe(TMM_dominant_component_concordant),
    UQ_dominant_concordance_evaluable = logical_mean_na_safe(UQ_dominant_component_concordant),
    TMM_stable_fraction = mean(TMM_stable),
    UQ_stable_fraction_evaluable = if (sum(UQ_evaluable) > 0L) mean(UQ_stable[UQ_evaluable]) else NA_real_,
    Tier_A_fraction = mean(projection_evidence_tier == "Tier_A_fit_and_both_normalizations"),
    Tier_B_fraction = mean(projection_evidence_tier == "Tier_B_fit_and_TMM_supported"),
    Tier_C_fraction = mean(projection_evidence_tier == "Tier_C_exploratory_or_sensitive"),
    global_low_fit_tail_fraction = mean(global_lowest_10pct_fit),
    highly_sensitive_fraction = mean(highly_normalization_sensitive),
    severe_extreme_z_fraction = mean(severe_extreme_z_flag)
  ),
  by = .(
    analysis_cancer_class,
    primary_gse_id
  )
]

study_summary[
  ,
  n_studies_in_cancer := uniqueN(primary_gse_id),
  by = analysis_cancer_class
]

study_summary[
  ,
  adequacy_robust_z_within_cancer := if (.N >= MIN_STUDIES_FOR_WITHIN_CANCER_OUTLIER) {
    robust_z(median_adequacy_composite)
  } else {
    rep(NA_real_, .N)
  },
  by = analysis_cancer_class
]

study_summary[
  ,
  poor_fit_study_outlier :=
    n_samples >= MIN_STUDY_N_FOR_OUTLIER &
      !is.na(adequacy_robust_z_within_cancer) &
      adequacy_robust_z_within_cancer <= STUDY_ROBUST_Z_THRESHOLD
]

study_summary[
  ,
  study_projection_status := fifelse(
    n_samples < MIN_STUDY_N_FOR_OUTLIER,
    "small_study_descriptive_only",
    fifelse(
      poor_fit_study_outlier,
      "within_cancer_poor_fit_outlier",
      "no_within_cancer_fit_outlier"
    )
  )
]

setorder(study_summary, analysis_cancer_class, median_adequacy_composite)

write_tsv(
  study_summary,
  file.path(
    OUTPUT_DIR,
    "03_study_level_transferability_summary.tsv"
  )
)

cancer_summary <- sample_atlas[
  ,
  .(
    n_samples = .N,
    n_studies = uniqueN(primary_gse_id),
    n_UQ_evaluable = sum(UQ_evaluable),
    n_UQ_not_evaluable = sum(!UQ_evaluable),
    UQ_evaluable_fraction = mean(UQ_evaluable),
    median_adequacy_composite = median(adequacy_composite_percentile),
    q10_adequacy_composite = quantile(adequacy_composite_percentile, 0.10, names = FALSE),
    q90_adequacy_composite = quantile(adequacy_composite_percentile, 0.90, names = FALSE),
    median_normalized_KL = median(normalized_KL_per_input_mass),
    median_reconstruction_Pearson = median(reconstruction_Pearson),
    median_reconstruction_cosine = median(reconstruction_cosine),
    median_reconstruction_R2 = median(reconstruction_R2),
    median_score_entropy = median(normalized_score_entropy),
    median_dominant_fraction = median(dominant_fraction),
    median_TMM_TVD = finite_median(TMM_total_variation_distance),
    q95_TMM_TVD = finite_quantile(TMM_total_variation_distance, 0.95),
    median_UQ_TVD_evaluable = finite_median(UQ_total_variation_distance),
    q95_UQ_TVD_evaluable = finite_quantile(UQ_total_variation_distance, 0.95),
    TMM_dominant_concordance = logical_mean_na_safe(TMM_dominant_component_concordant),
    UQ_dominant_concordance_evaluable = logical_mean_na_safe(UQ_dominant_component_concordant),
    TMM_stable_fraction = mean(TMM_stable),
    UQ_stable_fraction_evaluable = if (sum(UQ_evaluable) > 0L) mean(UQ_stable[UQ_evaluable]) else NA_real_,
    robust_both_fraction_all_samples = mean(normalization_robustness_tier == "robust_TMM_and_UQ"),
    robust_both_fraction_UQ_evaluable = if (sum(UQ_evaluable) > 0L) mean((TMM_stable & UQ_stable)[UQ_evaluable]) else NA_real_,
    highly_sensitive_fraction = mean(highly_normalization_sensitive),
    Tier_A_fraction = mean(projection_evidence_tier == "Tier_A_fit_and_both_normalizations"),
    Tier_B_fraction = mean(projection_evidence_tier == "Tier_B_fit_and_TMM_supported"),
    Tier_C_fraction = mean(projection_evidence_tier == "Tier_C_exploratory_or_sensitive"),
    severe_extreme_z_fraction = mean(severe_extreme_z_flag)
  ),
  by = analysis_cancer_class
]

setorder(cancer_summary, -median_adequacy_composite)
write_tsv(cancer_summary, file.path(OUTPUT_DIR, "04_cancer_level_transferability_summary.tsv"))

# ==============================================================================
# 10. Residual-feature atlas
# ==============================================================================

message(
  "Summarizing signed directional residual structure."
)

feature_metadata <- copy(
  as.data.table(
    reconstruction$directional_feature_metadata
  )
)

check_columns(
  feature_metadata,
  c(
    "directional_feature",
    "model_gene",
    "direction",
    "row_index"
  ),
  "Directional feature metadata"
)

if (
  !identical(
    feature_metadata$directional_feature,
    rownames(
      residual_matrix
    )
  )
) {
  stop(
    "Directional feature metadata and residual-matrix rows differ.",
    call. = FALSE
  )
}

overall_residual_summary <- data.table(
  directional_feature = rownames(
    residual_matrix
  ),
  model_gene = feature_metadata$model_gene,
  direction = feature_metadata$direction,
  mean_signed_residual = rowMeans(
    residual_matrix
  ),
  median_signed_residual = apply(
    residual_matrix,
    1L,
    stats::median
  ),
  mean_absolute_residual = rowMeans(
    abs(
      residual_matrix
    )
  ),
  q90_absolute_residual = apply(
    abs(
      residual_matrix
    ),
    1L,
    stats::quantile,
    probs = 0.90,
    names = FALSE
  ),
  fraction_samples_under_reconstructed = rowMeans(
    residual_matrix >
      0
  ),
  fraction_samples_over_reconstructed = rowMeans(
    residual_matrix <
      0
  )
)

setorder(
  overall_residual_summary,
  -mean_absolute_residual
)

write_tsv(
  overall_residual_summary,
  file.path(
    OUTPUT_DIR,
    "05_overall_directional_residual_feature_summary.tsv"
  )
)

cancer_levels <- unique(
  scores$analysis_cancer_class
)

cancer_residual_summary <- rbindlist(
  lapply(
    cancer_levels,
    function(cancer_name) {
      sample_index <- which(
        scores$analysis_cancer_class ==
          cancer_name
      )

      residual_sub <- residual_matrix[
        ,
        sample_index,
        drop = FALSE
      ]

      data.table(
        analysis_cancer_class = cancer_name,
        directional_feature = rownames(
          residual_sub
        ),
        model_gene = feature_metadata$model_gene,
        direction = feature_metadata$direction,
        n_samples = length(
          sample_index
        ),
        mean_signed_residual = rowMeans(
          residual_sub
        ),
        mean_absolute_residual = rowMeans(
          abs(
            residual_sub
          )
        ),
        fraction_samples_under_reconstructed = rowMeans(
          residual_sub >
            0
        )
      )
    }
  ),
  use.names = TRUE
)

write_tsv(
  cancer_residual_summary,
  file.path(
    OUTPUT_DIR,
    "06_cancer_directional_residual_feature_summary.tsv.gz"
  ),
  compress = TRUE
)

feature_selection <- cancer_residual_summary[
  ,
  .(
    maximum_absolute_cancer_mean_signed_residual =
      max(
        abs(
          mean_signed_residual
        )
      ),
    between_cancer_SD_signed_residual =
      stats::sd(
        mean_signed_residual
      ),
    global_mean_absolute_residual =
      mean(
        mean_absolute_residual
      )
  ),
  by = .(
    directional_feature,
    model_gene,
    direction
  )
][
  order(
    -maximum_absolute_cancer_mean_signed_residual,
    -between_cancer_SD_signed_residual
  )
]

top_residual_features <- head(
  feature_selection,
  TOP_RESIDUAL_FEATURES
)

write_tsv(
  feature_selection,
  file.path(
    OUTPUT_DIR,
    "07_residual_feature_cancer_specificity_ranking.tsv"
  )
)

residual_heatmap_source <- merge(
  cancer_residual_summary,
  top_residual_features[
    ,
    .(
      directional_feature,
      feature_rank = seq_len(.N)
    )
  ],
  by = "directional_feature",
  all = FALSE,
  sort = FALSE
)

write_tsv(
  residual_heatmap_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_5_residual_heatmap_source.tsv"
  )
)

# ==============================================================================
# 11. Essential figures
# ==============================================================================

message(
  "Generating Nature-style Layer 1 figures."
)

cancer_order <- cancer_summary$analysis_cancer_class

sample_atlas[
  ,
  cancer_label := factor(
    analysis_cancer_class,
    levels = rev(
      cancer_order
    )
  )
]

# Figure 1: empirical model adequacy distributions.
figure_1_source <- sample_atlas[
  ,
  .(
    sample_id,
    analysis_cancer_class,
    cancer_label,
    adequacy_composite_percentile,
    relative_fit_tier
  )
]

write_tsv(
  figure_1_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_1_model_adequacy_by_cancer_source.tsv.gz"
  ),
  compress = TRUE
)

figure_1 <- ggplot(
  figure_1_source,
  aes(
    x = adequacy_composite_percentile,
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
  geom_point(
    data = cancer_summary,
    aes(
      x = median_adequacy_composite,
      y = factor(
        analysis_cancer_class,
        levels = rev(
          cancer_order
        )
      )
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.8,
    stroke = 0.25,
    fill = "black",
    colour = "black"
  ) +
  geom_vline(
    xintercept = c(
      0.25,
      0.75
    ),
    linetype = "dashed",
    linewidth = 0.25,
    colour = "grey65"
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.5,
      0.75,
      1
    )
  ) +
  labs(
    title = "Fixed-W reconstruction adequacy differs across cancer classes",
    subtitle = "The composite is the mean of empirical good-direction percentiles for KL, Pearson, cosine similarity and R².",
    x = "relative adequacy composite (0 = lower, 1 = higher)",
    y = NULL,
    caption = "Black points show cancer medians. Dashed lines mark atlas-wide quartiles; these are relative ranks, not biological validity cutoffs."
  ) +
  theme_nature(
    base_size = 7
  )

save_plot(
  figure_1,
  "FIGURE_1_model_adequacy_by_cancer",
  width_mm = NATURE_W_DOUBLE_MM,
  height_mm = 118
)

# Figure 2: normalization TVD by cancer.
normalization_plot_source <- rbindlist(
  list(
    sample_atlas[
      ,
      .(
        sample_id,
        analysis_cancer_class,
        comparison = "CPM vs TMM-CPM",
        comparison_evaluable = TRUE,
        total_variation_distance = TMM_total_variation_distance
      )
    ],
    sample_atlas[
      UQ_evaluable == TRUE,
      .(
        sample_id,
        analysis_cancer_class,
        comparison = "CPM vs UQ-CPM",
        comparison_evaluable = TRUE,
        total_variation_distance = UQ_total_variation_distance
      )
    ]
  )
)

normalization_plot_source[
  ,
  cancer_label := factor(
    analysis_cancer_class,
    levels = rev(cancer_order)
  )
]

write_tsv(
  normalization_plot_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_2_normalization_sensitivity_by_cancer_source.tsv.gz"
  ),
  compress = TRUE
)

figure_2 <- ggplot(
  normalization_plot_source,
  aes(
    x = total_variation_distance,
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
  facet_wrap(~ comparison, nrow = 1, scales = "fixed") +
  geom_vline(
    data = data.table(
      comparison = c("CPM vs TMM-CPM", "CPM vs UQ-CPM"),
      threshold = c(TMM_STABLE_TVD, UQ_STABLE_TVD)
    ),
    aes(xintercept = threshold),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.30,
    colour = "grey45"
  ) +
  scale_x_continuous(
    breaks = seq(0, 0.8, 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(xlim = c(0, 0.85)) +
  labs(
    title = "Normalization sensitivity is heterogeneous across cancers",
    subtitle = "TMM uses all 9,465 samples; formal upper-quartile comparisons use the 9,302 evaluable profiles only.",
    x = "total variation distance",
    y = NULL,
    caption = "Dashed lines show the prespecified stable thresholds: 0.10 for TMM and 0.15 for upper-quartile normalization."
  ) +
  theme_nature(base_size = 7)

save_plot(
  figure_2,
  "FIGURE_2_normalization_sensitivity_by_cancer",
  width_mm = NATURE_W_DOUBLE_MM,
  height_mm = 118
)

# Figure 3: evidence-tier composition by cancer.# Figure 3: evidence-tier composition by cancer.
tier_plot_source <- sample_atlas[
  ,
  .N,
  by = .(
    analysis_cancer_class,
    projection_evidence_tier
  )
][
  ,
  fraction := N /
    sum(N),
  by = analysis_cancer_class
]

tier_plot_source[
  ,
  analysis_cancer_class := factor(
    analysis_cancer_class,
    levels = rev(
      cancer_order
    )
  )
]

tier_plot_source[
  ,
  projection_evidence_tier := factor(
    projection_evidence_tier,
    levels = c(
      "Tier_A_fit_and_both_normalizations",
      "Tier_B_fit_and_TMM_supported",
      "Tier_C_exploratory_or_sensitive"
    ),
    labels = c(
      "Tier A",
      "Tier B",
      "Tier C"
    )
  )
]

write_tsv(
  tier_plot_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_3_projection_evidence_tiers_by_cancer_source.tsv"
  )
)

figure_3 <- ggplot(
  tier_plot_source,
  aes(
    x = fraction,
    y = analysis_cancer_class,
    fill = projection_evidence_tier
  )
) +
  geom_col(
    width = 0.72,
    colour = "black",
    linewidth = 0.15
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    labels = label_percent(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  scale_fill_manual(
    values = c(
      "Tier A" = "#FFFFFF",
      "Tier B" = "#BDBDBD",
      "Tier C" = "#4D4D4D"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Projection evidence tiers vary across cancer classes",
    subtitle = "Tier A requires evaluable agreement with both alternatives; Tier B permits TMM-supported profiles with intermediate or non-evaluable UQ robustness.",
    x = "fraction of samples",
    y = NULL,
    fill = NULL,
    caption = "UQ-non-evaluable profiles cannot enter Tier A. Tier C is an evidence-quality label, not a biological exclusion label."
  ) +
  theme_nature(
    base_size = 7
  ) +
  theme(
    legend.position = "top"
  )

save_plot(
  figure_3,
  "FIGURE_3_projection_evidence_tiers_by_cancer",
  width_mm = NATURE_W_DOUBLE_MM,
  height_mm = 118
)

# Figure 4: adequacy versus normalization sensitivity.
figure_4_source <- sample_atlas[
  ,
  .(
    sample_id,
    analysis_cancer_class,
    primary_gse_id,
    adequacy_composite_percentile,
    maximum_normalization_TVD,
    UQ_evaluable,
    projection_evidence_tier,
    severe_extreme_z_flag
  )
]

write_tsv(
  figure_4_source,
  file.path(
    SOURCE_DATA_DIR,
    "FIGURE_4_adequacy_vs_normalization_sensitivity_source.tsv.gz"
  ),
  compress = TRUE
)

figure_4 <- ggplot(
  figure_4_source,
  aes(
    x = adequacy_composite_percentile,
    y = maximum_normalization_TVD
  )
) +
  geom_vline(
    xintercept = 0.25,
    linewidth = 0.25,
    linetype = "dashed",
    colour = "grey55"
  ) +
  geom_point(
    aes(shape = UQ_evaluable),
    size = 0.65,
    alpha = 0.30,
    colour = "black"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "UQ evaluable", `FALSE` = "UQ not evaluable")
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  coord_cartesian(ylim = c(0, 0.85)) +
  labs(
    title = "Model adequacy and normalization robustness are distinct properties",
    subtitle = "For UQ-non-evaluable profiles, maximum TVD is defined by the complete CPM-versus-TMM comparison.",
    x = "relative adequacy composite",
    y = "maximum formally evaluable normalization TVD",
    shape = NULL,
    caption = "The vertical guide marks the minimum adequacy percentile used for Tier A; continuous metrics remain primary."
  ) +
  theme_nature(base_size = 7) +
  theme(legend.position = "top")

save_plot(
  figure_4,
  "FIGURE_4_adequacy_vs_normalization_sensitivity",
  width_mm = NATURE_W_SINGLE_MM,
  height_mm = 90
)

# Figure 5: top cancer-specific directional residuals.# Figure 5: top cancer-specific directional residuals.
residual_heatmap_source[
  ,
  cancer_label := factor(
    analysis_cancer_class,
    levels = cancer_order
  )
]

feature_order <- top_residual_features$directional_feature

residual_heatmap_source[
  ,
  feature_label := factor(
    directional_feature,
    levels = rev(
      feature_order
    )
  )
]

residual_limit <- stats::quantile(
  abs(
    residual_heatmap_source$mean_signed_residual
  ),
  probs = 0.98,
  names = FALSE
)

figure_5 <- ggplot(
  residual_heatmap_source,
  aes(
    x = cancer_label,
    y = feature_label,
    fill = mean_signed_residual
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -residual_limit,
      residual_limit
    ),
    oob = squish
  ) +
  labs(
    title = "Directional residual features reveal cancer-specific under- and over-reconstruction",
    subtitle = "Positive residuals indicate observed directional signal stronger than W×H; negative residuals indicate model over-reconstruction.",
    x = NULL,
    y = "directional model feature",
    fill = "mean\nresidual",
    caption = "Shown are the 30 features with the largest cancer-specific mean signed residuals."
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
    ),
    legend.position = "right"
  )

save_plot(
  figure_5,
  "FIGURE_5_cancer_specific_directional_residuals",
  width_mm = NATURE_W_DOUBLE_MM,
  height_mm = 132
)

# ==============================================================================
# 12. Focused outlier and diagnostic tables
# ==============================================================================

poor_fit_studies <- study_summary[
  poor_fit_study_outlier ==
    TRUE
][
  order(
    adequacy_robust_z_within_cancer
  )
]

write_tsv(
  poor_fit_studies,
  file.path(
    OUTPUT_DIR,
    "08_within_cancer_poor_fit_study_outliers.tsv"
  )
)

worst_samples <- sample_atlas[
  order(
    adequacy_composite_percentile,
    -maximum_normalization_TVD
  )
][
  1:min(
    100L,
    .N
  )
]

write_tsv(
  worst_samples,
  file.path(
    OUTPUT_DIR,
    "09_lowest_adequacy_samples_with_robustness.tsv"
  )
)

normalization_sensitive_samples <- sample_atlas[
  highly_normalization_sensitive ==
    TRUE
][
  order(
    -maximum_normalization_TVD
  )
]

write_tsv(
  normalization_sensitive_samples,
  file.path(
    OUTPUT_DIR,
    "10_highly_normalization_sensitive_samples.tsv.gz"
  ),
  compress = TRUE
)

extreme_z_samples <- sample_atlas[
  severe_extreme_z_flag ==
    TRUE
][
  order(
    -maximum_absolute_fixed_anchor_z,
    -n_features_abs_z_ge_8
  )
]

write_tsv(
  extreme_z_samples,
  file.path(
    OUTPUT_DIR,
    "11_extreme_fixed_anchor_z_samples.tsv"
  )
)

# ==============================================================================
# 13. QC contract and report
# ==============================================================================

qc <- rbindlist(
  list(
    data.table(
      check = "projection_input_status_PASS_v1_0_4",
      status = if (
        identical(projection_status$status, "PASS") &&
          projection_status$qc_failures == 0 &&
          identical(projection_status$script_version, EXPECTED_PROJECTION_VERSION)
      ) "PASS" else "FAIL",
      detail = paste0(
        "status=", projection_status$status,
        ";version=", projection_status$script_version,
        ";qc_failures=", projection_status$qc_failures
      )
    ),
    data.table(
      check = "projection_input_QC_all_PASS",
      status = if (nrow(projection_qc) > 0L && all(projection_qc$status == "PASS")) "PASS" else "FAIL",
      detail = paste0("rows=", nrow(projection_qc), ";failures=", sum(projection_qc$status != "PASS"))
    ),
    data.table(
      check = "sample_atlas_complete",
      status = if (
        nrow(sample_atlas) == EXPECTED_SAMPLES &&
          !anyDuplicated(sample_atlas$sample_id) &&
          uniqueN(sample_atlas$analysis_cancer_class) == EXPECTED_CANCERS
      ) "PASS" else "FAIL",
      detail = paste0(
        "rows=", nrow(sample_atlas),
        ";cancers=", uniqueN(sample_atlas$analysis_cancer_class)
      )
    ),
    data.table(
      check = "adequacy_composite_finite_bounded",
      status = if (
        all(is.finite(sample_atlas$adequacy_composite_percentile)) &&
          all(sample_atlas$adequacy_composite_percentile >= 0 & sample_atlas$adequacy_composite_percentile <= 1)
      ) "PASS" else "FAIL",
      detail = paste0(
        "range=", format_number(min(sample_atlas$adequacy_composite_percentile)),
        "-", format_number(max(sample_atlas$adequacy_composite_percentile))
      )
    ),
    data.table(
      check = "TMM_metrics_complete_all_samples",
      status = if (
        all(sample_atlas$TMM_evaluable) &&
          all(is.finite(sample_atlas$TMM_total_variation_distance)) &&
          all(!is.na(sample_atlas$TMM_dominant_component_concordant))
      ) "PASS" else "FAIL",
      detail = paste0("evaluable=", sum(sample_atlas$TMM_evaluable), ";expected=", EXPECTED_SAMPLES)
    ),
    data.table(
      check = "UQ_evaluability_partition_exact",
      status = if (
        sum(sample_atlas$UQ_evaluable) == EXPECTED_UQ_EVALUABLE &&
          sum(!sample_atlas$UQ_evaluable) == EXPECTED_UQ_NOT_EVALUABLE
      ) "PASS" else "FAIL",
      detail = paste0(
        "evaluable=", sum(sample_atlas$UQ_evaluable),
        ";non_evaluable=", sum(!sample_atlas$UQ_evaluable)
      )
    ),
    data.table(
      check = "UQ_formal_metrics_only_when_evaluable",
      status = if (
        all(is.finite(sample_atlas[UQ_evaluable == TRUE, UQ_total_variation_distance])) &&
          all(!is.na(sample_atlas[UQ_evaluable == TRUE, UQ_dominant_component_concordant])) &&
          all(is.na(sample_atlas[UQ_evaluable == FALSE, UQ_total_variation_distance])) &&
          all(is.na(sample_atlas[UQ_evaluable == FALSE, UQ_dominant_component_concordant]))
      ) "PASS" else "FAIL",
      detail = paste0(
        "formal_rows=", sum(sample_atlas$UQ_evaluable),
        ";NA_rows=", sum(!sample_atlas$UQ_evaluable)
      )
    ),
    data.table(
      check = "maximum_normalization_TVD_complete",
      status = if (all(is.finite(sample_atlas$maximum_normalization_TVD))) "PASS" else "FAIL",
      detail = paste0("finite=", sum(is.finite(sample_atlas$maximum_normalization_TVD)))
    ),
    data.table(
      check = "Tier_A_requires_UQ_evaluable_and_stable",
      status = if (
        all(
          sample_atlas[
            projection_evidence_tier == "Tier_A_fit_and_both_normalizations",
            UQ_evaluable & UQ_stable & TMM_stable
          ]
        )
      ) "PASS" else "FAIL",
      detail = paste0(
        "Tier_A=", sum(sample_atlas$projection_evidence_tier == "Tier_A_fit_and_both_normalizations"),
        ";Tier_A_UQ_non_evaluable=", sum(
          sample_atlas$projection_evidence_tier == "Tier_A_fit_and_both_normalizations" & !sample_atlas$UQ_evaluable
        )
      )
    ),
    data.table(
      check = "Tier_B_UQ_state_allowed",
      status = if (
        all(
          sample_atlas[
            projection_evidence_tier == "Tier_B_fit_and_TMM_supported",
            !UQ_evaluable |
              UQ_robustness_state %in% c("stable", "intermediate")
          ]
        )
      ) "PASS" else "FAIL",
      detail = paste0(
        "Tier_B=", sum(sample_atlas$projection_evidence_tier == "Tier_B_fit_and_TMM_supported"),
        ";Tier_B_UQ_unstable=", sum(
          sample_atlas$projection_evidence_tier == "Tier_B_fit_and_TMM_supported" &
            sample_atlas$UQ_robustness_state == "unstable"
        )
      )
    ),
    data.table(
      check = "fixed_anchor_z_metrics_complete",
      status = if (
        all(is.finite(sample_atlas$maximum_absolute_fixed_anchor_z)) &&
          nrow(z_gene_metrics) == EXPECTED_GENES
      ) "PASS" else "FAIL",
      detail = paste0("samples=", nrow(z_sample_metrics), ";genes=", nrow(z_gene_metrics))
    ),
    data.table(
      check = "residual_feature_summary_complete",
      status = if (
        nrow(overall_residual_summary) == EXPECTED_DIRECTIONAL_FEATURES &&
          nrow(cancer_residual_summary) == EXPECTED_DIRECTIONAL_FEATURES * EXPECTED_CANCERS
      ) "PASS" else "FAIL",
      detail = paste0(
        "overall_features=", nrow(overall_residual_summary),
        ";cancer_feature_rows=", nrow(cancer_residual_summary)
      )
    ),
    data.table(
      check = "study_and_cancer_summaries_complete",
      status = if (
        nrow(study_summary) ==
          uniqueN(sample_atlas[, paste(analysis_cancer_class, primary_gse_id, sep = "|")]) &&
          nrow(study_by_cohort_summary) ==
            uniqueN(sample_atlas[, paste(analysis_cancer_class, primary_gse_id, provisional_cohort, sep = "|")]) &&
          nrow(cancer_summary) == EXPECTED_CANCERS &&
          !anyDuplicated(study_summary[, .(analysis_cancer_class, primary_gse_id)])
      ) "PASS" else "FAIL",
      detail = paste0(
        "cancer_studies=", nrow(study_summary),
        ";study_tier_strata=", nrow(study_by_cohort_summary),
        ";cancers=", nrow(cancer_summary)
      )
    ),
    data.table(
      check = "evidence_tiers_exhaustive",
      status = if (
        sum(table(sample_atlas$projection_evidence_tier)) == EXPECTED_SAMPLES &&
          !anyNA(sample_atlas$projection_evidence_tier)
      ) "PASS" else "FAIL",
      detail = paste(
        names(table(sample_atlas$projection_evidence_tier)),
        as.integer(table(sample_atlas$projection_evidence_tier)),
        sep = "=",
        collapse = ";"
      )
    )
  ),
  use.names = TRUE
)

write_tsv(
  qc,
  file.path(OUTPUT_DIR, "12_layer1_adequacy_robustness_QC_contract.tsv")
)

failed_qc <- qc[status == "FAIL"]

tier_counts <- table(sample_atlas$projection_evidence_tier)

get_tier_count <- function(label) {
  if (label %in% names(tier_counts)) as.integer(tier_counts[[label]]) else 0L
}

global_summary <- sample_atlas[
  ,
  .(
    n_samples = .N,
    n_UQ_evaluable = sum(UQ_evaluable),
    n_UQ_not_evaluable = sum(!UQ_evaluable),
    median_adequacy_composite = median(adequacy_composite_percentile),
    global_lowest_10pct_fit_fraction = mean(global_lowest_10pct_fit),
    TMM_stable_fraction = mean(TMM_stable),
    UQ_stable_fraction_evaluable = mean(UQ_stable[UQ_evaluable]),
    robust_both_fraction_all_samples = mean(TMM_stable & UQ_stable),
    robust_both_fraction_UQ_evaluable = mean((TMM_stable & UQ_stable)[UQ_evaluable]),
    highly_normalization_sensitive_fraction = mean(highly_normalization_sensitive),
    severe_extreme_z_fraction = mean(severe_extreme_z_flag),
    Tier_A_fraction = mean(projection_evidence_tier == "Tier_A_fit_and_both_normalizations"),
    Tier_B_fraction = mean(projection_evidence_tier == "Tier_B_fit_and_TMM_supported"),
    Tier_C_fraction = mean(projection_evidence_tier == "Tier_C_exploratory_or_sensitive"),
    n_poor_fit_study_outliers = nrow(poor_fit_studies)
  )
]

report_lines <- c(
  "ARCHS4 LAYER 1 ADEQUACY AND ROBUSTNESS ATLAS V1.0.2.1",
  paste(rep("=", 96), collapse = ""),
  "",
  paste0("Status: ", if (nrow(failed_qc) == 0L) "PASS" else "FAIL"),
  paste0("Input projection version: ", projection_status$script_version),
  paste0("Canonical samples: ", nrow(sample_atlas)),
  paste0("Cancer classes: ", uniqueN(sample_atlas$analysis_cancer_class)),
  paste0("Unique GEO studies: ", uniqueN(sample_atlas$primary_gse_id)),
  paste0("Cancer-study combinations: ", nrow(study_summary)),
  paste0("Study-by-Gold/Silver strata retained separately: ", nrow(study_by_cohort_summary)),
  paste0("UQ robustness evaluable: ", sum(sample_atlas$UQ_evaluable)),
  paste0("UQ robustness non-evaluable: ", sum(!sample_atlas$UQ_evaluable)),
  "",
  "INTERPRETATION BOUNDARY",
  "-----------------------",
  "This stage is descriptive and diagnostic. Adequacy tiers are empirical ranks",
  "within the current ARCHS4 atlas and are not externally calibrated validity",
  "thresholds. Continuous adequacy and sensitivity metrics remain primary.",
  "",
  "ADEQUACY COMPOSITE",
  "------------------",
  "The composite is the unweighted mean of four empirical good-direction",
  "percentiles: lower normalized KL, higher Pearson, higher cosine similarity",
  "and higher reconstruction R2. Q1/Q2-Q3/Q4 labels describe relative atlas",
  "position only.",
  "",
  "NORMALIZATION ROBUSTNESS",
  "------------------------",
  paste0("TMM stable: TVD <= ", TMM_STABLE_TVD, " and dominant component unchanged."),
  paste0("UQ stable: TVD <= ", UQ_STABLE_TVD, " and dominant component unchanged, among UQ-evaluable samples only."),
  paste0("Highly sensitive: TMM TVD >= ", TMM_HIGH_SENSITIVITY_TVD, ", evaluable UQ TVD >= ", UQ_HIGH_SENSITIVITY_TVD, ", or both evaluable alternative dominant components differ."),
  "The 163 Stage 09 UQ-fallback profiles retain diagnostic projections but are",
  "not classified as UQ stable or unstable and cannot enter Tier A.",
  "",
  "PROJECTION EVIDENCE TIERS",
  "-------------------------",
  "Tier A: adequacy percentile >= 0.25, stable under TMM and evaluable UQ, and",
  "no severe fixed-anchor z-score flag.",
  "Tier B: adequacy percentile >= 0.10, TMM stable, not highly normalization-",
  "sensitive, and no severe fixed-anchor z-score flag. When UQ is evaluable it",
  "must be stable or intermediate; UQ may also be formally non-evaluable.",
  "Tier C: exploratory, low-relative-adequacy, normalization-sensitive, or severe",
  "fixed-anchor-z evidence.",
  "",
  "GLOBAL SUMMARY",
  "--------------",
  paste(capture.output(print(global_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "CANCER-LEVEL SUMMARY",
  "--------------------",
  paste(capture.output(print(cancer_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "WITHIN-CANCER POOR-FIT STUDY OUTLIERS",
  "-------------------------------------",
  if (nrow(poor_fit_studies) == 0L) {
    "None under the prespecified robust-z rule."
  } else {
    paste(capture.output(print(head(poor_fit_studies, TOP_STUDIES_FOR_TABLE), row.names = FALSE)), collapse = "\n")
  },
  "",
  "RESIDUAL INTERPRETATION",
  "-----------------------",
  "Positive signed residual means the observed directional feature is stronger",
  "than W x H reconstructs. Negative residual means W x H reconstructs more",
  "signal than observed. NonSE is not used as the unexplained fraction.",
  "",
  "QC CONTRACT",
  "-----------",
  paste(capture.output(print(qc, row.names = FALSE)), collapse = "\n")
)

report_file <- file.path(
  OUTPUT_DIR,
  "ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1_report.txt"
)
writeLines(report_lines, con = report_file)
message("Wrote: ", report_file)

status <- list(
  script_version = SCRIPT_VERSION,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = if (nrow(failed_qc) == 0L) "PASS" else "FAIL",
  input_projection_status = projection_status$status,
  input_projection_script_version = projection_status$script_version,
  source_projection_dir = normalizePath(INPUT_DIR, winslash = "/", mustWork = TRUE),
  canonical_samples = nrow(sample_atlas),
  cancer_classes = uniqueN(sample_atlas$analysis_cancer_class),
  GEO_studies = uniqueN(sample_atlas$primary_gse_id),
  cancer_study_combinations = nrow(study_summary),
  study_tier_strata = nrow(study_by_cohort_summary),
  UQ_evaluable_samples = sum(sample_atlas$UQ_evaluable),
  UQ_non_evaluable_samples = sum(!sample_atlas$UQ_evaluable),
  Tier_A_samples = get_tier_count("Tier_A_fit_and_both_normalizations"),
  Tier_B_samples = get_tier_count("Tier_B_fit_and_TMM_supported"),
  Tier_C_samples = get_tier_count("Tier_C_exploratory_or_sensitive"),
  Tier_A_UQ_non_evaluable_samples = sum(
    sample_atlas$projection_evidence_tier == "Tier_A_fit_and_both_normalizations" &
      !sample_atlas$UQ_evaluable
  ),
  poor_fit_study_outliers = nrow(poor_fit_studies),
  qc_failures = nrow(failed_qc)
)

status_file <- file.path(
  OUTPUT_DIR,
  "13_layer1_adequacy_robustness_status.json"
)
jsonlite::write_json(
  status,
  path = status_file,
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
message("Wrote: ", status_file)

if (nrow(failed_qc) > 0L) {
  print(failed_qc)
  stop("Layer 1 adequacy/robustness QC contract failed.", call. = FALSE)
}

cat(
  "\n", paste(rep("=", 96), collapse = ""), "\n",
  "ARCHS4 LAYER 1 ADEQUACY AND ROBUSTNESS ATLAS V1.0.2.1: PASS\n",
  paste(rep("=", 96), collapse = ""), "\n",
  "Samples:                  ", nrow(sample_atlas), "\n",
  "Cancer classes:           ", uniqueN(sample_atlas$analysis_cancer_class), "\n",
  "GEO studies:              ", uniqueN(sample_atlas$primary_gse_id), "\n",
  "UQ evaluable:             ", status$UQ_evaluable_samples, "\n",
  "UQ non-evaluable:         ", status$UQ_non_evaluable_samples, "\n",
  "Tier A samples:           ", status$Tier_A_samples, "\n",
  "Tier B samples:           ", status$Tier_B_samples, "\n",
  "Tier C samples:           ", status$Tier_C_samples, "\n",
  "Poor-fit study outliers:  ", status$poor_fit_study_outliers, "\n",
  "Output:                   ", OUTPUT_DIR, "\n",
  sep = ""
)

