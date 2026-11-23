#!/usr/bin/env Rscript

# ==============================================================================
# 13_build_ARCHS4_projection_use_analysis_ready_atlas_v1_0_3.R
# ==============================================================================
#
# Purpose
# -------
# Build the final combined, analysis-ready projection atlas after the frozen:
#   Stage 09 expression/normalization release,
#   Stage 10 shared pan-cancer projection,
#   Stage 11 Layer 1 adequacy/robustness atlas, and
#   Stage 12 cancer-specific projection.
#
# This stage performs no new SpatialEcoTyper projection. It:
#   1. joins shared pan-cancer and cancer-specific score spaces;
#   2. preserves Layer 1 as the authoritative evidence system for cross-cancer
#      analyses;
#   3. derives a parallel within-cancer evidence system from cancer-specific
#      reconstruction adequacy, normalization sensitivity and fixed-anchor-z
#      burden;
#   4. carries pan-cancer and cancer-specific poor-fit study cautions;
#   5. writes strict primary and supported-use manifests for each score space;
#   6. exports a long score table and explicit score-space policy contract;
#   7. retains persistent residuals as hypotheses, not newly discovered
#      ecotypes.
#
# The two coordinate systems remain separate:
#   - pan-cancer scores support comparisons across cancer classes;
#   - cancer-specific scores support analyses within one cancer class only.
#
# ==============================================================================

options(stringsAsFactors = FALSE, width = 260, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required.", call. = FALSE)
}

# ==============================================================================
# 1. Constants and paths
# ==============================================================================

SCRIPT_VERSION <- "1.0.3"
EXPECTED_STAGE09_VERSION <- "1.0.1"
EXPECTED_STAGE10_VERSION <- "1.0.4"
EXPECTED_STAGE11_VERSION <- "1.0.2.1"
EXPECTED_STAGE12_VERSION <- "1.0.3"

DATA_ROOT <- Sys.getenv(
  "ARCHS4_V26_ROOT",
  unset = "/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6"
)

STAGE09_DIR <- Sys.getenv(
  "ARCHS4_STAGE13_STAGE09_DIR",
  unset = file.path(DATA_ROOT, "10_preprojection_normalization_v1_0_1")
)

PAN_DIR <- Sys.getenv(
  "ARCHS4_STAGE13_PAN_DIR",
  unset = file.path(DATA_ROOT, "11_pan_cancer_fixed_anchor_projection_v1_0_4")
)

LAYER1_DIR <- Sys.getenv(
  "ARCHS4_STAGE13_LAYER1_DIR",
  unset = file.path(DATA_ROOT, "12_layer1_adequacy_robustness_atlas_v1_0_2_1")
)

CANCER_SPECIFIC_DIR <- Sys.getenv(
  "ARCHS4_STAGE13_CANCER_SPECIFIC_DIR",
  unset = file.path(DATA_ROOT, "13_cancer_specific_fixed_anchor_projection_v1_0_3")
)

OUTPUT_DIR <- Sys.getenv(
  "ARCHS4_STAGE13_OUTPUT_DIR",
  unset = file.path(DATA_ROOT, "14_projection_use_analysis_ready_atlas_v1_0_3")
)

COMPONENTS <- c("NonSE", paste0("SE", 1:9))
EXPECTED_SAMPLES <- 9465L
EXPECTED_CANCERS <- 18L
EXPECTED_MODEL_GENES <- 499L
EXPECTED_UQ_EVALUABLE <- 9302L
EXPECTED_UQ_NOT_EVALUABLE <- 163L
EXPECTED_PAN_TIER_A <- 4200L
EXPECTED_PAN_TIER_B <- 1026L
EXPECTED_PAN_TIER_C <- 4239L
EXPECTED_PAN_POOR_FIT_STUDIES <- 5L
EXPECTED_PERSISTENT_FEATURES <- 49L

TMM_STABLE_TVD <- 0.10
UQ_STABLE_TVD <- 0.15
TMM_HIGH_SENSITIVITY_TVD <- 0.25
UQ_HIGH_SENSITIVITY_TVD <- 0.35
EXTREME_Z_THRESHOLD <- 8
SEVERE_Z_THRESHOLD <- 10
SEVERE_Z_FEATURE_COUNT <- 5L
MIN_STUDY_N_FOR_OUTLIER <- 5L
MIN_STUDIES_FOR_WITHIN_CANCER_OUTLIER <- 3L
STUDY_ROBUST_Z_THRESHOLD <- -2.5

SELF_TEST_MODE <- "--run-self-tests" %in% commandArgs(trailingOnly = TRUE)

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
}

# Input files
STAGE09_STATUS_FILE <- file.path(STAGE09_DIR, "18_preprojection_normalization_status.json")
STAGE09_QC_FILE <- file.path(STAGE09_DIR, "17_preprojection_normalization_QC_contract.tsv")

PAN_STATUS_FILE <- file.path(PAN_DIR, "28_pan_cancer_projection_status.json")
PAN_QC_FILE <- file.path(PAN_DIR, "27_pan_cancer_projection_QC_contract.tsv")
PAN_SCORE_FILE <- file.path(PAN_DIR, "17_primary_CPM_fraction_scores_with_metadata.tsv.gz")

LAYER1_STATUS_FILE <- file.path(LAYER1_DIR, "13_layer1_adequacy_robustness_status.json")
LAYER1_QC_FILE <- file.path(LAYER1_DIR, "12_layer1_adequacy_robustness_QC_contract.tsv")
LAYER1_SAMPLE_FILE <- file.path(LAYER1_DIR, "01_sample_level_adequacy_robustness_atlas.tsv.gz")
LAYER1_STUDY_FILE <- file.path(LAYER1_DIR, "03_study_level_transferability_summary.tsv")
LAYER1_OUTLIER_FILE <- file.path(LAYER1_DIR, "08_within_cancer_poor_fit_study_outliers.tsv")

CS_STATUS_FILE <- file.path(CANCER_SPECIFIC_DIR, "32_cancer_specific_projection_status.json")
CS_QC_FILE <- file.path(CANCER_SPECIFIC_DIR, "31_cancer_specific_projection_QC_contract.tsv")
CS_SCORE_FILE <- file.path(CANCER_SPECIFIC_DIR, "17_primary_CPM_cancer_specific_fraction_scores_with_metadata.tsv.gz")
CS_ADEQUACY_FILE <- file.path(CANCER_SPECIFIC_DIR, "18_primary_CPM_cancer_specific_fixed_W_model_adequacy.tsv.gz")
CS_SENSITIVITY_FILE <- file.path(CANCER_SPECIFIC_DIR, "20_normalization_sensitivity_by_sample.tsv.gz")
CS_SCALED_INPUT_FILE <- file.path(CANCER_SPECIFIC_DIR, "14_cancer_specific_scaled_input_CPM.rds")
CS_SHIFT_FILE <- file.path(CANCER_SPECIFIC_DIR, "23_pan_vs_cancer_specific_coordinate_shift_by_sample.tsv.gz")
CS_ADEQUACY_DELTA_FILE <- file.path(CANCER_SPECIFIC_DIR, "25_pan_vs_cancer_specific_adequacy_delta_by_sample.tsv.gz")
CS_RESIDUAL_ATTENUATION_FILE <- file.path(CANCER_SPECIFIC_DIR, "27_pan_vs_cancer_specific_residual_attenuation.tsv.gz")
CS_PERSISTENT_FEATURE_FILE <- file.path(CANCER_SPECIFIC_DIR, "28_cancer_specific_persistent_residual_feature_ranking.tsv")
CS_ATTENUATION_SUMMARY_FILE <- file.path(CANCER_SPECIFIC_DIR, "29_residual_attenuation_summary_by_cancer.tsv")

# ==============================================================================
# 2. Helpers
# ==============================================================================

write_tsv <- function(x, path, compress = FALSE) {
  fwrite(
    as.data.table(x),
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = if (compress) "gzip" else "none"
  )
  message("Wrote: ", path)
}

check_columns <- function(x, required_columns, label) {
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0L) {
    stop(
      label,
      " lacks required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

empirical_percentile <- function(x, higher_is_better = TRUE) {
  x <- as.numeric(x)
  n_valid <- sum(is.finite(x))
  if (n_valid < 2L) {
    return(rep(NA_real_, length(x)))
  }
  ranks <- frank(x, ties.method = "average", na.last = "keep")
  percentile <- (ranks - 1) / (n_valid - 1)
  if (!higher_is_better) {
    percentile <- 1 - percentile
  }
  percentile
}

robust_z <- function(x) {
  x <- as.numeric(x)
  centre <- stats::median(x, na.rm = TRUE)
  scale_value <- stats::mad(
    x,
    center = centre,
    constant = 1.4826,
    na.rm = TRUE
  )
  if (!is.finite(scale_value) || scale_value <= 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - centre) / scale_value
}

finite_median <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) return(NA_real_)
  stats::median(x[is.finite(x)])
}

logical_mean_na_safe <- function(x) {
  x <- as.logical(x)
  if (!any(!is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

qc_row <- function(check, passed, detail) {
  data.table(
    check = check,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    detail = as.character(detail)
  )
}

normalization_to_wide <- function(x, prefix) {
  required <- c(
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
    "dominant_component_concordant"
  )
  check_columns(x, required, paste0(prefix, " normalization sensitivity"))

  use <- copy(x[, ..required])
  use[, comparison_evaluable := as.logical(comparison_evaluable)]
  use[, dominant_component_concordant := as.logical(dominant_component_concordant)]
  use[, diagnostic_dominant_component_concordant :=
        as.logical(diagnostic_dominant_component_concordant)]

  wide <- dcast(
    use,
    sample_id ~ comparison,
    value.var = setdiff(required, c("sample_id", "comparison"))
  )

  rename_map <- c(
    comparison_evaluable_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_evaluable"),
    comparison_evaluable_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_evaluable"),
    non_evaluable_reason_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_non_evaluable_reason"),
    non_evaluable_reason_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_non_evaluable_reason"),
    diagnostic_total_variation_distance_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_diagnostic_TVD"),
    diagnostic_total_variation_distance_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_diagnostic_TVD"),
    diagnostic_maximum_component_absolute_difference_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_diagnostic_max_component_difference"),
    diagnostic_maximum_component_absolute_difference_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_diagnostic_max_component_difference"),
    diagnostic_aitchison_distance_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_diagnostic_aitchison_distance"),
    diagnostic_aitchison_distance_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_diagnostic_aitchison_distance"),
    diagnostic_dominant_component_concordant_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_diagnostic_dominant_concordant"),
    diagnostic_dominant_component_concordant_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_diagnostic_dominant_concordant"),
    total_variation_distance_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_TVD"),
    total_variation_distance_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_TVD"),
    maximum_component_absolute_difference_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_max_component_difference"),
    maximum_component_absolute_difference_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_max_component_difference"),
    aitchison_distance_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_aitchison_distance"),
    aitchison_distance_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_aitchison_distance"),
    dominant_component_concordant_CPM_vs_TMM_CPM = paste0(prefix, "_TMM_dominant_concordant"),
    dominant_component_concordant_CPM_vs_UQ_CPM = paste0(prefix, "_UQ_dominant_concordant")
  )

  for (old_name in names(rename_map)) {
    if (!old_name %in% names(wide)) {
      stop("Missing reshaped normalization column: ", old_name, call. = FALSE)
    }
    setnames(wide, old_name, rename_map[[old_name]])
  }
  wide
}

assign_normalization_states <- function(dt, prefix) {
  tmm_eval <- paste0(prefix, "_TMM_evaluable")
  uq_eval <- paste0(prefix, "_UQ_evaluable")
  tmm_tvd <- paste0(prefix, "_TMM_TVD")
  uq_tvd <- paste0(prefix, "_UQ_TVD")
  tmm_dom <- paste0(prefix, "_TMM_dominant_concordant")
  uq_dom <- paste0(prefix, "_UQ_dominant_concordant")

  dt[, (tmm_eval) := as.logical(get(tmm_eval))]
  dt[, (uq_eval) := as.logical(get(uq_eval))]
  dt[, (tmm_dom) := as.logical(get(tmm_dom))]
  dt[, (uq_dom) := as.logical(get(uq_dom))]

  dt[, (paste0(prefix, "_TMM_stable")) :=
       get(tmm_eval) & is.finite(get(tmm_tvd)) &
       get(tmm_tvd) <= TMM_STABLE_TVD & get(tmm_dom) == TRUE]

  dt[, (paste0(prefix, "_UQ_stable")) :=
       get(uq_eval) & is.finite(get(uq_tvd)) &
       get(uq_tvd) <= UQ_STABLE_TVD & get(uq_dom) == TRUE]

  dt[, (paste0(prefix, "_UQ_robustness_state")) :=
       fifelse(
         !get(uq_eval),
         "not_evaluable",
         fifelse(
           get(paste0(prefix, "_UQ_stable")),
           "stable",
           fifelse(
             get(uq_tvd) >= UQ_HIGH_SENSITIVITY_TVD | get(uq_dom) == FALSE,
             "unstable",
             "intermediate"
           )
         )
       )]

  dt[, (paste0(prefix, "_maximum_normalization_TVD")) :=
       fifelse(
         get(uq_eval),
         pmax(get(tmm_tvd), get(uq_tvd)),
         get(tmm_tvd)
       )]

  dt[, (paste0(prefix, "_highly_normalization_sensitive")) :=
       get(tmm_tvd) >= TMM_HIGH_SENSITIVITY_TVD |
       (get(uq_eval) & get(uq_tvd) >= UQ_HIGH_SENSITIVITY_TVD) |
       (get(uq_eval) & get(tmm_dom) == FALSE & get(uq_dom) == FALSE)]

  invisible(dt)
}

assign_evidence_tier <- function(
    dt,
    adequacy_rank_column,
    prefix,
    severe_z_column
) {
  tmm_stable <- paste0(prefix, "_TMM_stable")
  uq_eval <- paste0(prefix, "_UQ_evaluable")
  uq_stable <- paste0(prefix, "_UQ_stable")
  uq_state <- paste0(prefix, "_UQ_robustness_state")
  high_sensitive <- paste0(prefix, "_highly_normalization_sensitive")
  tier_column <- paste0(prefix, "_projection_evidence_tier")
  detail_column <- paste0(prefix, "_projection_evidence_tier_detail")
  use_column <- paste0(prefix, "_recommended_use")

  dt[, (tier_column) := fifelse(
    get(adequacy_rank_column) >= 0.25 &
      get(tmm_stable) &
      get(uq_eval) &
      get(uq_stable) &
      !get(severe_z_column),
    "Tier_A_fit_and_both_normalizations",
    fifelse(
      get(adequacy_rank_column) >= 0.10 &
        get(tmm_stable) &
        !get(high_sensitive) &
        (!get(uq_eval) | get(uq_state) %in% c("stable", "intermediate")) &
        !get(severe_z_column),
      "Tier_B_fit_and_TMM_supported",
      "Tier_C_exploratory_or_sensitive"
    )
  )]

  dt[, (detail_column) := fifelse(
    get(tier_column) == "Tier_A_fit_and_both_normalizations",
    "Tier_A_UQ_evaluable_and_stable",
    fifelse(
      get(tier_column) == "Tier_B_fit_and_TMM_supported" & !get(uq_eval),
      "Tier_B_TMM_supported_UQ_not_evaluable",
      fifelse(
        get(tier_column) == "Tier_B_fit_and_TMM_supported",
        paste0("Tier_B_TMM_supported_UQ_", get(uq_state)),
        "Tier_C_exploratory_or_sensitive"
      )
    )
  )]

  dt[, (use_column) := fifelse(
    get(tier_column) == "Tier_A_fit_and_both_normalizations",
    "primary_evidence",
    fifelse(
      get(tier_column) == "Tier_B_fit_and_TMM_supported",
      "sensitivity_supported_evidence",
      "exploratory_only"
    )
  )]

  invisible(dt)
}

# ==============================================================================
# 3. Self-tests
# ==============================================================================

run_self_tests <- function() {
  x <- c(1, 2, 3, 4)
  stopifnot(
    max(abs(empirical_percentile(x, TRUE) - c(0, 1/3, 2/3, 1))) < 1e-15,
    max(abs(empirical_percentile(x, FALSE) - c(1, 2/3, 1/3, 0))) < 1e-15
  )

  z <- robust_z(c(1, 2, 3, 4, 100))
  stopifnot(length(z) == 5L, z[[5L]] > 10)

  test <- data.table(
    sample_id = c("S1", "S2", "S3"),
    cs_TMM_evaluable = TRUE,
    cs_UQ_evaluable = c(TRUE, FALSE, TRUE),
    cs_TMM_TVD = c(0.01, 0.02, 0.40),
    cs_UQ_TVD = c(0.02, NA_real_, 0.50),
    cs_TMM_dominant_concordant = c(TRUE, TRUE, FALSE),
    cs_UQ_dominant_concordant = c(TRUE, NA, FALSE),
    cs_adequacy_rank = c(0.8, 0.8, 0.8),
    cs_severe_z = FALSE
  )
  assign_normalization_states(test, "cs")
  assign_evidence_tier(test, "cs_adequacy_rank", "cs", "cs_severe_z")
  stopifnot(
    test[sample_id == "S1", cs_projection_evidence_tier] ==
      "Tier_A_fit_and_both_normalizations",
    test[sample_id == "S2", cs_projection_evidence_tier] ==
      "Tier_B_fit_and_TMM_supported",
    test[sample_id == "S3", cs_projection_evidence_tier] ==
      "Tier_C_exploratory_or_sensitive",
    !any(test[cs_UQ_evaluable == FALSE, cs_projection_evidence_tier] ==
           "Tier_A_fit_and_both_normalizations")
  )

  cat("ALL PROJECTION-USE ANALYSIS-READY ATLAS V1.0.3 SELF-TESTS PASSED\n")
}

if (SELF_TEST_MODE) {
  run_self_tests()
  quit(save = "no", status = 0)
}

# ==============================================================================
# 4. Validate frozen inputs
# ==============================================================================

required_files <- c(
  STAGE09_STATUS_FILE, STAGE09_QC_FILE,
  PAN_STATUS_FILE, PAN_QC_FILE, PAN_SCORE_FILE,
  LAYER1_STATUS_FILE, LAYER1_QC_FILE, LAYER1_SAMPLE_FILE,
  LAYER1_STUDY_FILE, LAYER1_OUTLIER_FILE,
  CS_STATUS_FILE, CS_QC_FILE, CS_SCORE_FILE, CS_ADEQUACY_FILE,
  CS_SENSITIVITY_FILE, CS_SCALED_INPUT_FILE, CS_SHIFT_FILE,
  CS_ADEQUACY_DELTA_FILE, CS_RESIDUAL_ATTENUATION_FILE,
  CS_PERSISTENT_FEATURE_FILE, CS_ATTENUATION_SUMMARY_FILE
)

missing_files <- required_files[
  !file.exists(required_files) | file.info(required_files)$size <= 0
]
if (length(missing_files) > 0L) {
  stop(
    "Missing or empty required files:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

stage09_status <- jsonlite::read_json(STAGE09_STATUS_FILE, simplifyVector = TRUE)
pan_status <- jsonlite::read_json(PAN_STATUS_FILE, simplifyVector = TRUE)
layer1_status <- jsonlite::read_json(LAYER1_STATUS_FILE, simplifyVector = TRUE)
cs_status <- jsonlite::read_json(CS_STATUS_FILE, simplifyVector = TRUE)

status_contract <- list(
  stage09 = identical(stage09_status$status, "PASS") &&
    stage09_status$qc_failures == 0 &&
    identical(stage09_status$script_version, EXPECTED_STAGE09_VERSION),
  stage10 = identical(pan_status$status, "PASS") &&
    pan_status$qc_failures == 0 &&
    identical(pan_status$script_version, EXPECTED_STAGE10_VERSION),
  stage11 = identical(layer1_status$status, "PASS") &&
    layer1_status$qc_failures == 0 &&
    identical(layer1_status$script_version, EXPECTED_STAGE11_VERSION),
  stage12 = identical(cs_status$status, "PASS") &&
    cs_status$qc_failures == 0 &&
    identical(cs_status$script_version, EXPECTED_STAGE12_VERSION)
)
if (!all(unlist(status_contract))) {
  stop("At least one frozen prerequisite status/version contract failed.", call. = FALSE)
}

if (
  pan_status$canonical_samples != EXPECTED_SAMPLES ||
    layer1_status$canonical_samples != EXPECTED_SAMPLES ||
    cs_status$canonical_samples != EXPECTED_SAMPLES ||
    pan_status$cancer_classes != EXPECTED_CANCERS ||
    layer1_status$cancer_classes != EXPECTED_CANCERS ||
    cs_status$cancer_specific_score_spaces != EXPECTED_CANCERS ||
    pan_status$stage09_input_model_genes != EXPECTED_MODEL_GENES ||
    pan_status$UQ_evaluable_samples != EXPECTED_UQ_EVALUABLE ||
    pan_status$UQ_non_evaluable_samples != EXPECTED_UQ_NOT_EVALUABLE
) {
  stop("Frozen prerequisite dimensions are not the expected release state.", call. = FALSE)
}

qc_files <- c(STAGE09_QC_FILE, PAN_QC_FILE, LAYER1_QC_FILE, CS_QC_FILE)
for (path in qc_files) {
  qc_input <- fread(path, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  check_columns(qc_input, c("check", "status", "detail"), basename(path))
  if (nrow(qc_input) == 0L || any(qc_input$status != "PASS")) {
    stop("Upstream QC contract has non-PASS rows: ", path, call. = FALSE)
  }
}

# ==============================================================================
# 5. Load source tables and align sample axes
# ==============================================================================

pan_scores <- fread(PAN_SCORE_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
layer1 <- fread(LAYER1_SAMPLE_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
pan_study <- fread(LAYER1_STUDY_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
pan_outliers <- fread(LAYER1_OUTLIER_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
cs_scores <- fread(CS_SCORE_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
cs_adequacy <- fread(CS_ADEQUACY_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
cs_sensitivity_long <- fread(CS_SENSITIVITY_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
coordinate_shift <- fread(CS_SHIFT_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
adequacy_delta <- fread(CS_ADEQUACY_DELTA_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
residual_attenuation <- fread(CS_RESIDUAL_ATTENUATION_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
persistent_feature_rank <- fread(CS_PERSISTENT_FEATURE_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
attenuation_summary <- fread(CS_ATTENUATION_SUMMARY_FILE, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
cs_scaled_input <- readRDS(CS_SCALED_INPUT_FILE)

for (object_name in c("pan_scores", "layer1", "cs_scores", "cs_adequacy", "coordinate_shift", "adequacy_delta")) {
  object <- get(object_name)
  if (nrow(object) != EXPECTED_SAMPLES || anyDuplicated(object$sample_id)) {
    stop(object_name, " does not contain exactly 9,465 unique samples.", call. = FALSE)
  }
}

check_columns(pan_scores, c(
  "sample_id", "primary_gse_id", "analysis_cancer_class", "provisional_cohort",
  COMPONENTS
), "Pan-cancer scores")
check_columns(cs_scores, c(
  "sample_id", "analysis_cancer_class", "score_space_id", COMPONENTS
), "Cancer-specific scores")
check_columns(layer1, c(
  "sample_id", "projection_evidence_tier", "adequacy_composite_percentile",
  "adequacy_percentile_rank", "TMM_stable", "UQ_evaluable",
  "UQ_robustness_state", "recommended_use", "severe_extreme_z_flag",
  "technical_caution_reasons"
), "Layer 1 sample atlas")
check_columns(cs_adequacy, c(
  "sample_id", "analysis_cancer_class", "normalized_KL_per_input_mass",
  "reconstruction_Pearson", "reconstruction_cosine", "reconstruction_R2",
  "dominant_component", "dominant_fraction"
), "Cancer-specific adequacy")
check_columns(coordinate_shift, c(
  "sample_id", "pan_to_cancer_specific_TVD",
  "pan_to_cancer_specific_max_component_difference",
  "pan_dominant_component", "cancer_specific_dominant_component",
  "dominant_component_concordant"
), "Coordinate shift")
check_columns(adequacy_delta, c(
  "sample_id", "delta_normalized_KL", "delta_reconstruction_Pearson",
  "delta_reconstruction_cosine", "delta_reconstruction_R2",
  "all_four_metrics_improved"
), "Adequacy delta")

canonical_ids <- as.character(pan_scores$sample_id)
for (x in list(layer1, cs_scores, cs_adequacy, coordinate_shift, adequacy_delta)) {
  if (!setequal(canonical_ids, as.character(x$sample_id))) {
    stop("At least one sample-level input has a different sample-ID set.", call. = FALSE)
  }
}

if (
  !is.matrix(cs_scaled_input) ||
    nrow(cs_scaled_input) != EXPECTED_MODEL_GENES ||
    ncol(cs_scaled_input) != EXPECTED_SAMPLES ||
    !setequal(colnames(cs_scaled_input), canonical_ids) ||
    any(!is.finite(cs_scaled_input))
) {
  stop("Cancer-specific fixed-anchor scaled CPM matrix failed validation.", call. = FALSE)
}
cs_scaled_input <- cs_scaled_input[, canonical_ids, drop = FALSE]

# ==============================================================================
# 6. Input provenance and policy contract
# ==============================================================================

input_provenance <- data.table(
  role = c(
    "stage09_status", "stage09_QC",
    "pan_status", "pan_QC", "pan_scores",
    "layer1_status", "layer1_QC", "layer1_samples", "layer1_studies", "layer1_outliers",
    "cancer_specific_status", "cancer_specific_QC", "cancer_specific_scores",
    "cancer_specific_adequacy", "cancer_specific_sensitivity",
    "cancer_specific_scaled_CPM", "coordinate_shift", "adequacy_delta",
    "residual_attenuation", "persistent_residual_ranking", "attenuation_summary"
  ),
  path = c(
    STAGE09_STATUS_FILE, STAGE09_QC_FILE,
    PAN_STATUS_FILE, PAN_QC_FILE, PAN_SCORE_FILE,
    LAYER1_STATUS_FILE, LAYER1_QC_FILE, LAYER1_SAMPLE_FILE, LAYER1_STUDY_FILE, LAYER1_OUTLIER_FILE,
    CS_STATUS_FILE, CS_QC_FILE, CS_SCORE_FILE, CS_ADEQUACY_FILE, CS_SENSITIVITY_FILE,
    CS_SCALED_INPUT_FILE, CS_SHIFT_FILE, CS_ADEQUACY_DELTA_FILE,
    CS_RESIDUAL_ATTENUATION_FILE, CS_PERSISTENT_FEATURE_FILE, CS_ATTENUATION_SUMMARY_FILE
  )
)
input_provenance[, `:=`(
  size_bytes = file.info(path)$size,
  md5 = unname(tools::md5sum(path))
)]
write_tsv(input_provenance, file.path(OUTPUT_DIR, "00_input_provenance.tsv"))

policy <- data.table(
  analysis_scope = c(
    "cross_cancer_primary",
    "cross_cancer_sensitivity_supported",
    "within_cancer_primary",
    "within_cancer_sensitivity_supported",
    "persistent_residual_hypothesis"
  ),
  authoritative_score_space = c(
    "pan_cancer_shared",
    "pan_cancer_shared",
    "cancer_specific_<analysis_cancer_class>",
    "cancer_specific_<analysis_cancer_class>",
    "cancer_specific_residual_feature_space"
  ),
  cross_cancer_comparable = c(TRUE, TRUE, FALSE, FALSE, FALSE),
  minimum_evidence = c(
    "pan Tier A",
    "pan Tier A or B",
    "cancer-specific Tier A",
    "cancer-specific Tier A or B",
    "persistent residual flag; independent validation required"
  ),
  strict_study_rule = c(
    "exclude Stage 11 pan-cancer poor-fit study outliers",
    "exclude Stage 11 pan-cancer poor-fit study outliers",
    "exclude cancer-specific poor-fit study outliers",
    "exclude cancer-specific poor-fit study outliers",
    "do not interpret as a new ecotype"
  ),
  interpretation = c(
    "Primary comparisons among cancer classes use only shared pan-cancer coordinates.",
    "Tier B may support sensitivity analyses but should not replace Tier A primary results.",
    "Within-cancer analyses use the cancer-specific score space for that cancer only.",
    "Cancer-specific Tier B may support sensitivity analyses within one cancer.",
    "Persistent residuals are model-mismatch hypotheses, not automatically novel spatial ecotypes."
  )
)
write_tsv(policy, file.path(OUTPUT_DIR, "01_projection_use_policy.tsv"))

# ==============================================================================
# 7. Build combined master table
# ==============================================================================

metadata_candidates <- c(
  "expression_column_index", "expression_column_number", "sample_id",
  "primary_gse_id", "all_gse_ids", "analysis_cancer_class", "provisional_cohort",
  "stage08_tier", "stage08_patient_group", "stage08_biological_condition_signature",
  "stage08_condition_group_id", "stage08_group_decision", "stage08_cross_sectional_anchor",
  "stage07_expression_quality_tier", "expression_qc_warning_codes",
  "model_genes_detected", "overall_model_detected_weight_fraction",
  "min_component_detected_weight_fraction", "library_sum_archs4_pseudocounts",
  "longitudinal_or_repeated_context", "stage05_independence_group",
  "title", "source_name_ch1", "characteristics_ch1"
)
metadata_columns <- intersect(metadata_candidates, names(pan_scores))
master <- copy(pan_scores[, c(metadata_columns, COMPONENTS), with = FALSE])
setnames(master, COMPONENTS, paste0("pan_", COMPONENTS))

cs_score_sub <- copy(cs_scores[, c("sample_id", "score_space_id", COMPONENTS), with = FALSE])
setnames(cs_score_sub, COMPONENTS, paste0("cs_", COMPONENTS))
master <- merge(master, cs_score_sub, by = "sample_id", all.x = TRUE, sort = FALSE)
master <- master[match(canonical_ids, sample_id)]

layer1_derived_columns <- c(
  "sample_id",
  "fit_percentile_low_KL", "fit_percentile_high_Pearson",
  "fit_percentile_high_cosine", "fit_percentile_high_R2",
  "adequacy_composite_percentile", "adequacy_percentile_rank",
  "adequacy_percentile_within_cancer", "relative_fit_tier",
  "global_lowest_10pct_fit", "within_cancer_lowest_10pct_fit",
  "TMM_evaluable", "UQ_evaluable", "TMM_total_variation_distance",
  "UQ_total_variation_distance", "TMM_dominant_component_concordant",
  "UQ_dominant_component_concordant", "TMM_stable", "UQ_stable",
  "UQ_robustness_state", "maximum_normalization_TVD",
  "highly_normalization_sensitive", "normalization_robustness_tier",
  "maximum_absolute_fixed_anchor_z", "q99_absolute_fixed_anchor_z",
  "n_features_abs_z_ge_6", "n_features_abs_z_ge_8", "n_features_abs_z_ge_10",
  "severe_extreme_z_flag", "projection_evidence_tier",
  "projection_evidence_tier_detail", "recommended_use", "technical_caution_reasons"
)
check_columns(layer1, layer1_derived_columns, "Layer 1 derived fields")
pan_evidence <- copy(layer1[, ..layer1_derived_columns])
setnames(
  pan_evidence,
  setdiff(names(pan_evidence), "sample_id"),
  paste0("pan_", setdiff(names(pan_evidence), "sample_id"))
)
master <- merge(master, pan_evidence, by = "sample_id", all.x = TRUE, sort = FALSE)
master <- master[match(canonical_ids, sample_id)]

# Cancer-specific adequacy percentiles are calculated within each cancer.
cs_fit <- copy(cs_adequacy[, .(
  sample_id,
  analysis_cancer_class,
  cs_normalized_KL_per_input_mass = normalized_KL_per_input_mass,
  cs_reconstruction_Pearson = reconstruction_Pearson,
  cs_reconstruction_cosine = reconstruction_cosine,
  cs_reconstruction_R2 = reconstruction_R2,
  cs_normalized_score_entropy = normalized_score_entropy,
  cs_dominant_component = dominant_component,
  cs_dominant_fraction = dominant_fraction,
  cs_second_component = second_component,
  cs_second_fraction = second_fraction,
  cs_top1_minus_top2 = top1_minus_top2,
  cs_dominance_class = dominance_class
)])

cs_fit[, cs_fit_percentile_low_KL :=
         empirical_percentile(cs_normalized_KL_per_input_mass, FALSE),
       by = analysis_cancer_class]
cs_fit[, cs_fit_percentile_high_Pearson :=
         empirical_percentile(cs_reconstruction_Pearson, TRUE),
       by = analysis_cancer_class]
cs_fit[, cs_fit_percentile_high_cosine :=
         empirical_percentile(cs_reconstruction_cosine, TRUE),
       by = analysis_cancer_class]
cs_fit[, cs_fit_percentile_high_R2 :=
         empirical_percentile(cs_reconstruction_R2, TRUE),
       by = analysis_cancer_class]
cs_fit[, cs_adequacy_composite_within_cancer := rowMeans(cbind(
  cs_fit_percentile_low_KL,
  cs_fit_percentile_high_Pearson,
  cs_fit_percentile_high_cosine,
  cs_fit_percentile_high_R2
))]
cs_fit[, cs_adequacy_percentile_rank_within_cancer :=
         empirical_percentile(cs_adequacy_composite_within_cancer, TRUE),
       by = analysis_cancer_class]

cs_sensitivity_wide <- normalization_to_wide(cs_sensitivity_long, "cs")
assign_normalization_states(cs_sensitivity_wide, "cs")

cs_abs <- abs(cs_scaled_input)
cs_z <- data.table(
  sample_id = colnames(cs_scaled_input),
  cs_maximum_absolute_fixed_anchor_z = apply(cs_abs, 2L, max),
  cs_q99_absolute_fixed_anchor_z = apply(
    cs_abs, 2L, stats::quantile, probs = 0.99, names = FALSE
  ),
  cs_n_features_abs_z_ge_6 = colSums(cs_abs >= 6),
  cs_n_features_abs_z_ge_8 = colSums(cs_abs >= EXTREME_Z_THRESHOLD),
  cs_n_features_abs_z_ge_10 = colSums(cs_abs >= SEVERE_Z_THRESHOLD)
)
cs_z[, cs_severe_extreme_z_flag :=
       cs_maximum_absolute_fixed_anchor_z >= SEVERE_Z_THRESHOLD |
       cs_n_features_abs_z_ge_8 >= SEVERE_Z_FEATURE_COUNT]

master <- merge(master, cs_fit[, analysis_cancer_class := NULL], by = "sample_id", all.x = TRUE, sort = FALSE)
master <- merge(master, cs_sensitivity_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
master <- merge(master, cs_z, by = "sample_id", all.x = TRUE, sort = FALSE)
master <- master[match(canonical_ids, sample_id)]

assign_evidence_tier(
  master,
  adequacy_rank_column = "cs_adequacy_percentile_rank_within_cancer",
  prefix = "cs",
  severe_z_column = "cs_severe_extreme_z_flag"
)

# Coordinate shift and adequacy response fields.
shift_sub <- copy(coordinate_shift[, .(
  sample_id,
  pan_to_cancer_specific_TVD,
  pan_to_cancer_specific_max_component_difference,
  pan_dominant_component,
  cancer_specific_dominant_component,
  pan_to_cancer_specific_dominant_concordant = dominant_component_concordant
)])
delta_sub <- copy(adequacy_delta[, .(
  sample_id,
  delta_normalized_KL,
  delta_reconstruction_Pearson,
  delta_reconstruction_cosine,
  delta_reconstruction_R2,
  all_four_metrics_improved
)])
master <- merge(master, shift_sub, by = "sample_id", all.x = TRUE, sort = FALSE)
master <- merge(master, delta_sub, by = "sample_id", all.x = TRUE, sort = FALSE)
master <- master[match(canonical_ids, sample_id)]

# ==============================================================================
# 8. Study cautions and analysis-use flags
# ==============================================================================

check_columns(pan_study, c(
  "analysis_cancer_class", "primary_gse_id", "n_samples",
  "poor_fit_study_outlier", "study_projection_status",
  "adequacy_robust_z_within_cancer"
), "Layer 1 study summary")

pan_study_flags <- copy(pan_study[, .(
  analysis_cancer_class,
  primary_gse_id,
  pan_study_n_samples = n_samples,
  pan_poor_fit_study_outlier = as.logical(poor_fit_study_outlier),
  pan_study_projection_status = study_projection_status,
  pan_study_adequacy_robust_z_within_cancer = adequacy_robust_z_within_cancer
)])
master <- merge(
  master,
  pan_study_flags,
  by = c("analysis_cancer_class", "primary_gse_id"),
  all.x = TRUE,
  sort = FALSE
)
master <- master[match(canonical_ids, sample_id)]

cs_study_summary <- master[, .(
  n_samples = .N,
  n_Gold = sum(provisional_cohort == "GOLD", na.rm = TRUE),
  n_Silver = sum(provisional_cohort == "SILVER", na.rm = TRUE),
  median_cs_adequacy_composite = median(cs_adequacy_composite_within_cancer),
  q10_cs_adequacy_composite = as.numeric(stats::quantile(
    cs_adequacy_composite_within_cancer, 0.10, names = FALSE
  )),
  median_cs_normalized_KL = median(cs_normalized_KL_per_input_mass),
  median_cs_reconstruction_Pearson = median(cs_reconstruction_Pearson),
  median_cs_reconstruction_cosine = median(cs_reconstruction_cosine),
  median_cs_reconstruction_R2 = median(cs_reconstruction_R2),
  median_cs_TMM_TVD = finite_median(cs_TMM_TVD),
  median_cs_UQ_TVD_evaluable = finite_median(cs_UQ_TVD),
  cs_TMM_stable_fraction = mean(cs_TMM_stable),
  cs_UQ_stable_fraction_evaluable = if (sum(cs_UQ_evaluable) > 0L) {
    mean(cs_UQ_stable[cs_UQ_evaluable])
  } else {
    NA_real_
  },
  cs_Tier_A_fraction = mean(cs_projection_evidence_tier == "Tier_A_fit_and_both_normalizations"),
  cs_Tier_B_fraction = mean(cs_projection_evidence_tier == "Tier_B_fit_and_TMM_supported"),
  cs_Tier_C_fraction = mean(cs_projection_evidence_tier == "Tier_C_exploratory_or_sensitive"),
  median_pan_to_cs_TVD = median(pan_to_cancer_specific_TVD),
  fraction_all_four_adequacy_metrics_improved = mean(all_four_metrics_improved),
  cs_severe_extreme_z_fraction = mean(cs_severe_extreme_z_flag)
), by = .(analysis_cancer_class, primary_gse_id)]

cs_study_summary[, n_studies_in_cancer := .N, by = analysis_cancer_class]
cs_study_summary[, cs_study_adequacy_robust_z_within_cancer :=
                   if (.N >= MIN_STUDIES_FOR_WITHIN_CANCER_OUTLIER) {
                     robust_z(median_cs_adequacy_composite)
                   } else {
                     rep(NA_real_, .N)
                   },
                 by = analysis_cancer_class]
cs_study_summary[, cs_poor_fit_study_outlier :=
                   n_samples >= MIN_STUDY_N_FOR_OUTLIER &
                   !is.na(cs_study_adequacy_robust_z_within_cancer) &
                   cs_study_adequacy_robust_z_within_cancer <= STUDY_ROBUST_Z_THRESHOLD]
cs_study_summary[, cs_study_projection_status := fifelse(
  n_samples < MIN_STUDY_N_FOR_OUTLIER,
  "small_study_descriptive_only",
  fifelse(
    cs_poor_fit_study_outlier,
    "within_cancer_poor_fit_outlier",
    "no_within_cancer_fit_outlier"
  )
)]

# Hotfix v1.0.2: assert the final cancer-specific study fields before merging.
# Hotfix v1.0.3: use explicit logical comparisons in data.table row filters.
check_columns(cs_study_summary, c(
  "analysis_cancer_class", "primary_gse_id", "n_samples",
  "cs_poor_fit_study_outlier", "cs_study_projection_status",
  "cs_study_adequacy_robust_z_within_cancer"
), "Cancer-specific study summary")

master <- merge(
  master,
  cs_study_summary[, .(
    analysis_cancer_class,
    primary_gse_id,
    cs_study_n_samples = n_samples,
    cs_poor_fit_study_outlier,
    cs_study_projection_status,
    cs_study_adequacy_robust_z_within_cancer
  )],
  by = c("analysis_cancer_class", "primary_gse_id"),
  all.x = TRUE,
  sort = FALSE
)
master <- master[match(canonical_ids, sample_id)]

master[, cross_cancer_use_class := fcase(
  pan_projection_evidence_tier == "Tier_A_fit_and_both_normalizations", "primary",
  pan_projection_evidence_tier == "Tier_B_fit_and_TMM_supported", "sensitivity_supported",
  default = "exploratory"
)]
master[, within_cancer_use_class := fcase(
  cs_projection_evidence_tier == "Tier_A_fit_and_both_normalizations", "primary",
  cs_projection_evidence_tier == "Tier_B_fit_and_TMM_supported", "sensitivity_supported",
  default = "exploratory"
)]

master[, heterogeneous_other_solid_tumour_caution :=
         analysis_cancer_class == "other_solid_tumour"]
master[, cross_cancer_primary_eligible_strict :=
         cross_cancer_use_class == "primary" & !pan_poor_fit_study_outlier]
master[, cross_cancer_supported_eligible_strict :=
         cross_cancer_use_class %in% c("primary", "sensitivity_supported") &
         !pan_poor_fit_study_outlier]
master[, within_cancer_primary_eligible_strict :=
         within_cancer_use_class == "primary" &
         !cs_poor_fit_study_outlier &
         !heterogeneous_other_solid_tumour_caution]
master[, within_cancer_supported_eligible_strict :=
         within_cancer_use_class %in% c("primary", "sensitivity_supported") &
         !cs_poor_fit_study_outlier &
         !heterogeneous_other_solid_tumour_caution]

master[, recommended_cross_cancer_score_space := "pan_cancer_shared"]
master[, recommended_within_cancer_score_space := score_space_id]
master[, cancer_specific_scores_cross_cancer_comparable := FALSE]
master[, pan_scores_cross_cancer_comparable := TRUE]

master[, combined_analysis_caution_reasons := vapply(seq_len(.N), function(i) {
  reasons <- character()
  if (isTRUE(pan_poor_fit_study_outlier[[i]])) reasons <- c(reasons, "pan_poor_fit_study")
  if (isTRUE(cs_poor_fit_study_outlier[[i]])) reasons <- c(reasons, "cancer_specific_poor_fit_study")
  if (isTRUE(pan_highly_normalization_sensitive[[i]])) reasons <- c(reasons, "pan_normalization_sensitive")
  if (isTRUE(cs_highly_normalization_sensitive[[i]])) reasons <- c(reasons, "cancer_specific_normalization_sensitive")
  if (isTRUE(pan_severe_extreme_z_flag[[i]])) reasons <- c(reasons, "pan_extreme_fixed_anchor_z")
  if (isTRUE(cs_severe_extreme_z_flag[[i]])) reasons <- c(reasons, "cancer_specific_extreme_fixed_anchor_z")
  if (identical(pan_UQ_evaluable[[i]], FALSE)) reasons <- c(reasons, "UQ_not_evaluable")
  if (isTRUE(heterogeneous_other_solid_tumour_caution[[i]])) reasons <- c(reasons, "heterogeneous_other_solid_tumour")
  if (length(reasons) == 0L) "none" else paste(reasons, collapse = ";")
}, FUN.VALUE = character(1))]

# ==============================================================================
# 9. Write master, core and long score tables
# ==============================================================================

write_tsv(
  master,
  file.path(OUTPUT_DIR, "02_projection_analysis_ready_master.tsv.gz"),
  compress = TRUE
)

core_columns <- c(
  "sample_id", "primary_gse_id", "all_gse_ids", "analysis_cancer_class",
  "provisional_cohort", "stage08_tier", "stage08_patient_group",
  "stage08_condition_group_id", "stage08_cross_sectional_anchor",
  "stage07_expression_quality_tier", "expression_qc_warning_codes",
  "model_genes_detected", "overall_model_detected_weight_fraction",
  "min_component_detected_weight_fraction", "library_sum_archs4_pseudocounts",
  paste0("pan_", COMPONENTS), paste0("cs_", COMPONENTS), "score_space_id",
  "pan_adequacy_composite_percentile", "pan_adequacy_percentile_rank",
  "pan_TMM_stable", "pan_UQ_evaluable", "pan_UQ_robustness_state",
  "pan_maximum_normalization_TVD", "pan_severe_extreme_z_flag",
  "pan_projection_evidence_tier", "pan_recommended_use",
  "cs_adequacy_composite_within_cancer", "cs_adequacy_percentile_rank_within_cancer",
  "cs_TMM_stable", "cs_UQ_evaluable", "cs_UQ_robustness_state",
  "cs_maximum_normalization_TVD", "cs_severe_extreme_z_flag",
  "cs_projection_evidence_tier", "cs_recommended_use",
  "pan_to_cancer_specific_TVD", "pan_to_cancer_specific_dominant_concordant",
  "delta_normalized_KL", "delta_reconstruction_Pearson",
  "delta_reconstruction_cosine", "delta_reconstruction_R2",
  "all_four_metrics_improved",
  "pan_poor_fit_study_outlier", "cs_poor_fit_study_outlier",
  "cross_cancer_use_class", "within_cancer_use_class",
  "cross_cancer_primary_eligible_strict", "cross_cancer_supported_eligible_strict",
  "within_cancer_primary_eligible_strict", "within_cancer_supported_eligible_strict",
  "heterogeneous_other_solid_tumour_caution",
  "recommended_cross_cancer_score_space", "recommended_within_cancer_score_space",
  "combined_analysis_caution_reasons"
)
core_columns <- intersect(core_columns, names(master))
core <- master[, ..core_columns]
write_tsv(
  core,
  file.path(OUTPUT_DIR, "03_projection_analysis_ready_core.tsv.gz"),
  compress = TRUE
)

long_id_columns <- c(
  "sample_id", "primary_gse_id", "analysis_cancer_class", "provisional_cohort",
  "cross_cancer_use_class", "within_cancer_use_class",
  "cross_cancer_primary_eligible_strict", "cross_cancer_supported_eligible_strict",
  "within_cancer_primary_eligible_strict", "within_cancer_supported_eligible_strict"
)

pan_long <- melt(
  master[, c(long_id_columns, paste0("pan_", COMPONENTS)), with = FALSE],
  id.vars = long_id_columns,
  measure.vars = paste0("pan_", COMPONENTS),
  variable.name = "component",
  value.name = "projected_fraction"
)
pan_long[, `:=`(
  component = sub("^pan_", "", as.character(component)),
  score_space_scope = "pan_cancer_shared",
  score_space_id = "pan_cancer_shared",
  cross_cancer_comparable = TRUE
)]

cs_long <- melt(
  master[, c(long_id_columns, "score_space_id", paste0("cs_", COMPONENTS)), with = FALSE],
  id.vars = c(long_id_columns, "score_space_id"),
  measure.vars = paste0("cs_", COMPONENTS),
  variable.name = "component",
  value.name = "projected_fraction"
)
cs_long[, `:=`(
  component = sub("^cs_", "", as.character(component)),
  score_space_scope = "cancer_specific",
  cross_cancer_comparable = FALSE
)]

score_long <- rbindlist(list(pan_long, cs_long), use.names = TRUE, fill = TRUE)
setcolorder(score_long, c(
  "sample_id", "primary_gse_id", "analysis_cancer_class", "provisional_cohort",
  "score_space_scope", "score_space_id", "cross_cancer_comparable",
  "component", "projected_fraction",
  setdiff(names(score_long), c(
    "sample_id", "primary_gse_id", "analysis_cancer_class", "provisional_cohort",
    "score_space_scope", "score_space_id", "cross_cancer_comparable",
    "component", "projected_fraction"
  ))
))
write_tsv(
  score_long,
  file.path(OUTPUT_DIR, "04_projection_scores_long.tsv.gz"),
  compress = TRUE
)

# Strict manifests retain the relevant score columns and audit fields.
pan_manifest_columns <- intersect(c(
  "sample_id", "primary_gse_id", "analysis_cancer_class", "provisional_cohort",
  "stage08_patient_group", "stage08_condition_group_id",
  paste0("pan_", COMPONENTS),
  "pan_projection_evidence_tier", "cross_cancer_use_class",
  "pan_poor_fit_study_outlier", "combined_analysis_caution_reasons"
), names(master))
cs_manifest_columns <- intersect(c(
  "sample_id", "primary_gse_id", "analysis_cancer_class", "provisional_cohort",
  "stage08_patient_group", "stage08_condition_group_id", "score_space_id",
  paste0("cs_", COMPONENTS),
  "cs_projection_evidence_tier", "within_cancer_use_class",
  "cs_poor_fit_study_outlier", "combined_analysis_caution_reasons"
), names(master))

manifest_filter_columns <- c(
  "cross_cancer_primary_eligible_strict",
  "cross_cancer_supported_eligible_strict",
  "within_cancer_primary_eligible_strict",
  "within_cancer_supported_eligible_strict"
)
check_columns(master, manifest_filter_columns, "Master manifest eligibility flags")
if (
  any(!vapply(master[, ..manifest_filter_columns], is.logical, logical(1))) ||
  anyNA(master[, ..manifest_filter_columns])
) {
  stop(
    "Manifest eligibility flags must be complete logical columns before filtering.",
    call. = FALSE
  )
}

pan_primary_manifest <- master[cross_cancer_primary_eligible_strict == TRUE, ..pan_manifest_columns]
pan_supported_manifest <- master[cross_cancer_supported_eligible_strict == TRUE, ..pan_manifest_columns]
cs_primary_manifest <- master[within_cancer_primary_eligible_strict == TRUE, ..cs_manifest_columns]
cs_supported_manifest <- master[within_cancer_supported_eligible_strict == TRUE, ..cs_manifest_columns]

write_tsv(pan_primary_manifest, file.path(OUTPUT_DIR, "05_pan_cancer_primary_manifest.tsv.gz"), TRUE)
write_tsv(pan_supported_manifest, file.path(OUTPUT_DIR, "06_pan_cancer_supported_manifest.tsv.gz"), TRUE)
write_tsv(cs_primary_manifest, file.path(OUTPUT_DIR, "07_within_cancer_primary_manifest.tsv.gz"), TRUE)
write_tsv(cs_supported_manifest, file.path(OUTPUT_DIR, "08_within_cancer_supported_manifest.tsv.gz"), TRUE)

# ==============================================================================
# 10. Study/cancer summaries and persistent residual hypotheses
# ==============================================================================

study_use_summary <- merge(
  pan_study,
  cs_study_summary,
  by = c("analysis_cancer_class", "primary_gse_id"),
  all = TRUE,
  suffixes = c("_pan", "_cs"),
  sort = FALSE
)

strict_study_counts <- master[, .(
  cross_cancer_primary_strict_samples = sum(cross_cancer_primary_eligible_strict),
  cross_cancer_supported_strict_samples = sum(cross_cancer_supported_eligible_strict),
  within_cancer_primary_strict_samples = sum(within_cancer_primary_eligible_strict),
  within_cancer_supported_strict_samples = sum(within_cancer_supported_eligible_strict)
), by = .(analysis_cancer_class, primary_gse_id)]
study_use_summary <- merge(
  study_use_summary,
  strict_study_counts,
  by = c("analysis_cancer_class", "primary_gse_id"),
  all.x = TRUE,
  sort = FALSE
)
setorder(study_use_summary, analysis_cancer_class, primary_gse_id)
write_tsv(study_use_summary, file.path(OUTPUT_DIR, "09_study_projection_use_summary.tsv"))

cancer_use_summary <- master[, .(
  n_samples = .N,
  n_studies = uniqueN(primary_gse_id),
  n_UQ_evaluable = sum(pan_UQ_evaluable),
  pan_Tier_A = sum(pan_projection_evidence_tier == "Tier_A_fit_and_both_normalizations"),
  pan_Tier_B = sum(pan_projection_evidence_tier == "Tier_B_fit_and_TMM_supported"),
  pan_Tier_C = sum(pan_projection_evidence_tier == "Tier_C_exploratory_or_sensitive"),
  cs_Tier_A = sum(cs_projection_evidence_tier == "Tier_A_fit_and_both_normalizations"),
  cs_Tier_B = sum(cs_projection_evidence_tier == "Tier_B_fit_and_TMM_supported"),
  cs_Tier_C = sum(cs_projection_evidence_tier == "Tier_C_exploratory_or_sensitive"),
  cross_cancer_primary_strict = sum(cross_cancer_primary_eligible_strict),
  cross_cancer_supported_strict = sum(cross_cancer_supported_eligible_strict),
  within_cancer_primary_strict = sum(within_cancer_primary_eligible_strict),
  within_cancer_supported_strict = sum(within_cancer_supported_eligible_strict),
  pan_poor_fit_study_samples = sum(pan_poor_fit_study_outlier),
  cs_poor_fit_study_samples = sum(cs_poor_fit_study_outlier),
  median_pan_adequacy = median(pan_adequacy_composite_percentile),
  median_cs_adequacy_within_cancer = median(cs_adequacy_composite_within_cancer),
  median_pan_to_cs_TVD = median(pan_to_cancer_specific_TVD),
  dominant_coordinate_concordance = mean(pan_to_cancer_specific_dominant_concordant),
  median_delta_normalized_KL = median(delta_normalized_KL),
  median_delta_reconstruction_Pearson = median(delta_reconstruction_Pearson),
  median_delta_reconstruction_cosine = median(delta_reconstruction_cosine),
  median_delta_reconstruction_R2 = median(delta_reconstruction_R2),
  fraction_all_four_metrics_improved = mean(all_four_metrics_improved),
  pan_normalization_sensitive_fraction = mean(pan_highly_normalization_sensitive),
  cs_normalization_sensitive_fraction = mean(cs_highly_normalization_sensitive),
  pan_severe_z_fraction = mean(pan_severe_extreme_z_flag),
  cs_severe_z_fraction = mean(cs_severe_extreme_z_flag)
), by = analysis_cancer_class]

cancer_use_summary <- merge(
  cancer_use_summary,
  attenuation_summary,
  by = "analysis_cancer_class",
  all.x = TRUE,
  sort = FALSE
)
setorder(cancer_use_summary, -within_cancer_primary_strict)
write_tsv(cancer_use_summary, file.path(OUTPUT_DIR, "10_cancer_projection_use_summary.tsv"))

check_columns(residual_attenuation, c(
  "analysis_cancer_class", "directional_feature",
  "pan_mean_signed_residual", "cancer_specific_mean_signed_residual",
  "attenuation_fraction", "persistent_cancer_specific_residual"
), "Residual attenuation")

persistent_pairs <- residual_attenuation[
  as.logical(persistent_cancer_specific_residual) == TRUE
]
persistent_pairs <- merge(
  persistent_pairs,
  persistent_feature_rank,
  by = "directional_feature",
  all.x = TRUE,
  sort = FALSE
)
persistent_pairs[, attenuation_evaluable := is.finite(attenuation_fraction)]
persistent_pairs[, residual_interpretation := fifelse(
  cancer_specific_mean_signed_residual > 0,
  "positive_under_reconstructed",
  fifelse(
    cancer_specific_mean_signed_residual < 0,
    "negative_over_reconstructed",
    "zero_mean_residual"
  )
)]
persistent_pairs[, absolute_cancer_specific_mean_residual := abs(cancer_specific_mean_signed_residual)]
setorder(
  persistent_pairs,
  -absolute_cancer_specific_mean_residual,
  analysis_cancer_class,
  directional_feature
)
write_tsv(
  persistent_pairs,
  file.path(OUTPUT_DIR, "11_persistent_residual_cancer_feature_candidates.tsv")
)

# ==============================================================================
# 11. Derived-field dictionary
# ==============================================================================

derived_dictionary <- data.table(
  field = c(
    "pan_*", "cs_*", "score_space_id",
    "pan_projection_evidence_tier", "cs_projection_evidence_tier",
    "cross_cancer_primary_eligible_strict", "cross_cancer_supported_eligible_strict",
    "within_cancer_primary_eligible_strict", "within_cancer_supported_eligible_strict",
    "pan_to_cancer_specific_TVD", "all_four_metrics_improved",
    "pan_poor_fit_study_outlier", "cs_poor_fit_study_outlier",
    "combined_analysis_caution_reasons"
  ),
  source = c(
    "Stage 10 shared projection", "Stage 12 cancer-specific projection", "Stage 12",
    "Stage 11 Layer 1", "Stage 13 derived from Stage 12 metrics",
    "Stage 13", "Stage 13", "Stage 13", "Stage 13",
    "Stage 12", "Stage 12", "Stage 11", "Stage 13", "Stage 13"
  ),
  interpretation = c(
    "Shared pan-cancer NonSE/SE1-SE9 scores; comparable across cancer classes.",
    "Cancer-specific NonSE/SE1-SE9 scores; comparable only within the same cancer class.",
    "Explicit identifier of the cancer-specific coordinate system.",
    "Authoritative evidence tier for pan-cancer use.",
    "Parallel evidence tier for within-cancer use, based on within-cancer adequacy, normalization robustness and cancer-specific z burden.",
    "Tier A pan sample outside a Stage 11 poor-fit study.",
    "Tier A/B pan sample outside a Stage 11 poor-fit study.",
    "Cancer-specific Tier A sample outside a cancer-specific poor-fit study and outside the heterogeneous other-solid-tumour class.",
    "Cancer-specific Tier A/B sample outside a cancer-specific poor-fit study and outside the heterogeneous other-solid-tumour class.",
    "Coordinate change between shared and cancer-specific scores; not a technical error metric.",
    "TRUE only when cancer-specific scaling improves KL, Pearson, cosine and R2 simultaneously.",
    "Stage 11 within-cancer poor-fit study flag in the shared coordinate system.",
    "Stage 13 analogous poor-fit study flag in the cancer-specific coordinate system.",
    "Semicolon-delimited cautions; 'none' means no listed technical/use caution."
  )
)
write_tsv(derived_dictionary, file.path(OUTPUT_DIR, "01A_derived_field_dictionary.tsv"))

# ==============================================================================
# 12. QC contract, status and report
# ==============================================================================

pan_score_matrix <- as.matrix(master[, paste0("pan_", COMPONENTS), with = FALSE])
cs_score_matrix <- as.matrix(master[, paste0("cs_", COMPONENTS), with = FALSE])

pan_tier_counts <- table(master$pan_projection_evidence_tier)
cs_tier_counts <- table(master$cs_projection_evidence_tier)

qc <- rbindlist(list(
  qc_row(
    "frozen_prerequisite_stages_PASS",
    all(unlist(status_contract)),
    paste0(
      "stage09=", stage09_status$status, "/", stage09_status$script_version,
      ";stage10=", pan_status$status, "/", pan_status$script_version,
      ";stage11=", layer1_status$status, "/", layer1_status$script_version,
      ";stage12=", cs_status$status, "/", cs_status$script_version
    )
  ),
  qc_row(
    "master_sample_axis_complete",
    nrow(master) == EXPECTED_SAMPLES && !anyDuplicated(master$sample_id) &&
      uniqueN(master$analysis_cancer_class) == EXPECTED_CANCERS,
    paste0("rows=", nrow(master), ";cancers=", uniqueN(master$analysis_cancer_class))
  ),
  qc_row(
    "pan_scores_finite_nonnegative_sum_one",
    all(is.finite(pan_score_matrix)) && all(pan_score_matrix >= 0) &&
      max(abs(rowSums(pan_score_matrix) - 1)) <= 1e-8,
    paste0("max_sum_error=", format(max(abs(rowSums(pan_score_matrix) - 1)), scientific = TRUE))
  ),
  qc_row(
    "cancer_specific_scores_finite_nonnegative_sum_one",
    all(is.finite(cs_score_matrix)) && all(cs_score_matrix >= 0) &&
      max(abs(rowSums(cs_score_matrix) - 1)) <= 1e-8,
    paste0("max_sum_error=", format(max(abs(rowSums(cs_score_matrix) - 1)), scientific = TRUE))
  ),
  qc_row(
    "pan_evidence_partition_exact",
    unname(pan_tier_counts["Tier_A_fit_and_both_normalizations"]) == EXPECTED_PAN_TIER_A &&
      unname(pan_tier_counts["Tier_B_fit_and_TMM_supported"]) == EXPECTED_PAN_TIER_B &&
      unname(pan_tier_counts["Tier_C_exploratory_or_sensitive"]) == EXPECTED_PAN_TIER_C,
    paste0(
      "A=", pan_tier_counts["Tier_A_fit_and_both_normalizations"],
      ";B=", pan_tier_counts["Tier_B_fit_and_TMM_supported"],
      ";C=", pan_tier_counts["Tier_C_exploratory_or_sensitive"]
    )
  ),
  qc_row(
    "cancer_specific_evidence_partition_exhaustive",
    sum(cs_tier_counts) == EXPECTED_SAMPLES && length(cs_tier_counts) == 3L,
    paste(capture.output(print(cs_tier_counts)), collapse = " ")
  ),
  qc_row(
    "UQ_evaluability_exact_in_both_spaces",
    sum(master$pan_UQ_evaluable) == EXPECTED_UQ_EVALUABLE &&
      sum(!master$pan_UQ_evaluable) == EXPECTED_UQ_NOT_EVALUABLE &&
      sum(master$cs_UQ_evaluable) == EXPECTED_UQ_EVALUABLE &&
      sum(!master$cs_UQ_evaluable) == EXPECTED_UQ_NOT_EVALUABLE,
    paste0(
      "pan=", sum(master$pan_UQ_evaluable), "/", sum(!master$pan_UQ_evaluable),
      ";cs=", sum(master$cs_UQ_evaluable), "/", sum(!master$cs_UQ_evaluable)
    )
  ),
  qc_row(
    "UQ_non_evaluable_never_Tier_A",
    !any(master[pan_UQ_evaluable == FALSE, pan_projection_evidence_tier] ==
           "Tier_A_fit_and_both_normalizations") &&
      !any(master[cs_UQ_evaluable == FALSE, cs_projection_evidence_tier] ==
             "Tier_A_fit_and_both_normalizations"),
    "checked pan and cancer-specific evidence tiers"
  ),
  qc_row(
    "cancer_specific_z_metrics_complete",
    all(is.finite(master$cs_maximum_absolute_fixed_anchor_z)) &&
      all(master$cs_n_features_abs_z_ge_8 >= 0),
    paste0("samples=", sum(is.finite(master$cs_maximum_absolute_fixed_anchor_z)))
  ),
  qc_row(
    "pan_poor_fit_study_count_exact",
    uniqueN(pan_outliers[, paste(analysis_cancer_class, primary_gse_id)]) ==
      EXPECTED_PAN_POOR_FIT_STUDIES,
    paste0("unique_outliers=", uniqueN(pan_outliers[, paste(analysis_cancer_class, primary_gse_id)]))
  ),
  qc_row(
    "study_summaries_complete",
    nrow(cs_study_summary) == uniqueN(master[, paste(analysis_cancer_class, primary_gse_id)]) &&
      sum(cs_study_summary$n_samples) == EXPECTED_SAMPLES,
    paste0("cancer_studies=", nrow(cs_study_summary), ";samples=", sum(cs_study_summary$n_samples))
  ),
  qc_row(
    "strict_manifests_are_valid_subsets",
    !anyDuplicated(pan_primary_manifest$sample_id) &&
      !anyDuplicated(pan_supported_manifest$sample_id) &&
      !anyDuplicated(cs_primary_manifest$sample_id) &&
      !anyDuplicated(cs_supported_manifest$sample_id) &&
      all(pan_primary_manifest$sample_id %in% canonical_ids) &&
      all(cs_primary_manifest$sample_id %in% canonical_ids),
    paste0(
      "pan_primary=", nrow(pan_primary_manifest),
      ";pan_supported=", nrow(pan_supported_manifest),
      ";cs_primary=", nrow(cs_primary_manifest),
      ";cs_supported=", nrow(cs_supported_manifest)
    )
  ),
  qc_row(
    "strict_manifests_exclude_relevant_outlier_studies",
    all(master[sample_id %in% pan_primary_manifest$sample_id, !pan_poor_fit_study_outlier]) &&
      all(master[sample_id %in% cs_primary_manifest$sample_id, !cs_poor_fit_study_outlier]) &&
      !any(master[sample_id %in% cs_primary_manifest$sample_id, heterogeneous_other_solid_tumour_caution]),
    "pan primary excludes pan outliers; within-cancer primary excludes cs outliers and heterogeneous other-solid-tumour"
  ),
  qc_row(
    "long_score_table_complete",
    nrow(score_long) == EXPECTED_SAMPLES * length(COMPONENTS) * 2L &&
      all(is.finite(score_long$projected_fraction)) &&
      all(score_long$projected_fraction >= 0),
    paste0("rows=", nrow(score_long), ";expected=", EXPECTED_SAMPLES * 20L)
  ),
  qc_row(
    "persistent_residual_contract",
    uniqueN(persistent_pairs$directional_feature) == EXPECTED_PERSISTENT_FEATURES &&
      nrow(persistent_pairs) > 0L,
    paste0(
      "unique_features=", uniqueN(persistent_pairs$directional_feature),
      ";cancer_feature_pairs=", nrow(persistent_pairs),
      ";attenuation_evaluable=", sum(persistent_pairs$attenuation_evaluable)
    )
  ),
  qc_row(
    "score_space_use_policy_explicit",
    all(master$recommended_cross_cancer_score_space == "pan_cancer_shared") &&
      all(master$cancer_specific_scores_cross_cancer_comparable == FALSE) &&
      all(master$pan_scores_cross_cancer_comparable == TRUE),
    "pan scores are cross-cancer comparable; cancer-specific scores are within-cancer only"
  )
), use.names = TRUE)

write_tsv(qc, file.path(OUTPUT_DIR, "12_projection_use_QC_contract.tsv"))

qc_failures <- sum(qc$status != "PASS")
status <- list(
  script_version = SCRIPT_VERSION,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = if (qc_failures == 0L) "PASS" else "FAIL",
  source_stage09_version = stage09_status$script_version,
  source_stage10_version = pan_status$script_version,
  source_stage11_version = layer1_status$script_version,
  source_stage12_version = cs_status$script_version,
  canonical_samples = nrow(master),
  cancer_classes = uniqueN(master$analysis_cancer_class),
  pan_Tier_A_samples = unname(pan_tier_counts["Tier_A_fit_and_both_normalizations"]),
  pan_Tier_B_samples = unname(pan_tier_counts["Tier_B_fit_and_TMM_supported"]),
  pan_Tier_C_samples = unname(pan_tier_counts["Tier_C_exploratory_or_sensitive"]),
  cancer_specific_Tier_A_samples = unname(cs_tier_counts["Tier_A_fit_and_both_normalizations"]),
  cancer_specific_Tier_B_samples = unname(cs_tier_counts["Tier_B_fit_and_TMM_supported"]),
  cancer_specific_Tier_C_samples = unname(cs_tier_counts["Tier_C_exploratory_or_sensitive"]),
  pan_primary_strict_samples = nrow(pan_primary_manifest),
  pan_supported_strict_samples = nrow(pan_supported_manifest),
  within_cancer_primary_strict_samples = nrow(cs_primary_manifest),
  within_cancer_supported_strict_samples = nrow(cs_supported_manifest),
  pan_poor_fit_studies = uniqueN(pan_outliers[, paste(analysis_cancer_class, primary_gse_id)]),
  cancer_specific_poor_fit_studies = sum(cs_study_summary$cs_poor_fit_study_outlier),
  UQ_evaluable_samples = sum(master$pan_UQ_evaluable),
  UQ_non_evaluable_samples = sum(!master$pan_UQ_evaluable),
  persistent_residual_unique_features = uniqueN(persistent_pairs$directional_feature),
  persistent_residual_cancer_feature_pairs = nrow(persistent_pairs),
  persistent_residual_attenuation_evaluable_pairs = sum(persistent_pairs$attenuation_evaluable),
  qc_failures = qc_failures
)
jsonlite::write_json(
  status,
  path = file.path(OUTPUT_DIR, "13_projection_use_status.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)
message("Wrote: ", file.path(OUTPUT_DIR, "13_projection_use_status.json"))

report_path <- file.path(
  OUTPUT_DIR,
  "ARCHS4_projection_use_analysis_ready_atlas_v1_0_3_report.txt"
)
report_lines <- c(
  "ARCHS4 PROJECTION-USE AND ANALYSIS-READY ATLAS V1.0.3",
  paste(rep("=", 96), collapse = ""),
  "",
  paste0("Status: ", status$status),
  paste0("Canonical samples: ", status$canonical_samples),
  paste0("Cancer classes: ", status$cancer_classes),
  paste0("Frozen inputs: Stage 09 v", status$source_stage09_version,
         "; Stage 10 v", status$source_stage10_version,
         "; Stage 11 v", status$source_stage11_version,
         "; Stage 12 v", status$source_stage12_version),
  "",
  "SCORE-SPACE POLICY",
  "------------------",
  "Shared pan-cancer scores are authoritative for comparisons across cancer classes.",
  "Cancer-specific scores are authoritative for analyses within one cancer class only.",
  "Cancer-specific score values must not be compared directly across cancers.",
  "Persistent residuals remain hypotheses and are not newly discovered ecotypes without independent validation.",
  "",
  "EVIDENCE SYSTEMS",
  "----------------",
  paste0("Pan-cancer Tier A/B/C: ", status$pan_Tier_A_samples, "/",
         status$pan_Tier_B_samples, "/", status$pan_Tier_C_samples),
  paste0("Cancer-specific Tier A/B/C: ", status$cancer_specific_Tier_A_samples, "/",
         status$cancer_specific_Tier_B_samples, "/", status$cancer_specific_Tier_C_samples),
  "The pan-cancer tier is the frozen Stage 11 evidence tier.",
  "The cancer-specific tier is calculated within each cancer using the same adequacy and normalization policy plus cancer-specific fixed-anchor-z burden.",
  "",
  "STRICT ANALYSIS MANIFESTS",
  "-------------------------",
  paste0("Cross-cancer primary: ", status$pan_primary_strict_samples),
  paste0("Cross-cancer primary + sensitivity-supported: ", status$pan_supported_strict_samples),
  paste0("Within-cancer primary: ", status$within_cancer_primary_strict_samples),
  paste0("Within-cancer primary + sensitivity-supported: ", status$within_cancer_supported_strict_samples),
  "Strict manifests exclude the relevant poor-fit study outliers. The within-cancer strict manifests also exclude the heterogeneous other_solid_tumour category.",
  "",
  "STUDY CAUTIONS",
  "--------------",
  paste0("Pan-cancer poor-fit studies: ", status$pan_poor_fit_studies),
  paste0("Cancer-specific poor-fit studies: ", status$cancer_specific_poor_fit_studies),
  "Study outlier flags are technical/use cautions, not automatic biological exclusions from exploratory work.",
  "",
  "PERSISTENT RESIDUAL HYPOTHESES",
  "------------------------------",
  paste0("Unique features: ", status$persistent_residual_unique_features),
  paste0("Cancer-feature pairs: ", status$persistent_residual_cancer_feature_pairs),
  paste0("Attenuation-evaluable pairs: ", status$persistent_residual_attenuation_evaluable_pairs),
  "",
  "QC CONTRACT",
  "-----------",
  paste(capture.output(print(qc)), collapse = "\n")
)
writeLines(report_lines, con = report_path)
message("Wrote: ", report_path)

if (qc_failures > 0L) {
  stop("Projection-use atlas completed with ", qc_failures, " QC failure(s).", call. = FALSE)
}

cat("\n")
cat(paste(rep("=", 96), collapse = ""), "\n", sep = "")
cat("ARCHS4 PROJECTION-USE ANALYSIS-READY ATLAS V1.0.3: PASS\n")
cat("Output: ", OUTPUT_DIR, "\n", sep = "")
cat("Samples: ", nrow(master), "\n", sep = "")
cat("Pan strict primary/supported: ", nrow(pan_primary_manifest), "/", nrow(pan_supported_manifest), "\n", sep = "")
cat("Within-cancer strict primary/supported: ", nrow(cs_primary_manifest), "/", nrow(cs_supported_manifest), "\n", sep = "")
cat(paste(rep("=", 96), collapse = ""), "\n", sep = "")
