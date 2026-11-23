# FUSCC frozen-W methods: direct source inspection and transfer contract

**Date:** 2026-11-20. **Status:** Source-methods review; not a successful recount3 fixed-W transfer and not a release authorization.

## Sources inspected

The user-supplied `FUSCC_ORIGINAL_PROJECTION_SOURCE.tar.gz` contains all six exact original scripts: Stage 06 compatibility; Stage 09 expression preparation v1.0.1; Stage 10 pan-cancer projection v1.0.4; Stage 11 adequacy v1.0.2.1; Stage 12 cancer-specific projection v1.0.3; Stage 13 atlas construction v1.0.3, plus `Bulk_SE_Recovery_W.rds`, harmonization and README. The RDS SHA-256 is `39a96b598ec74ccbe67bea49254ac94ae568ca22ec162241ab7849976b52cdcf`, matching the user's earlier copy. The original code is complete; no rerun was performed.

## The exact projection contract

1. **Full-gene expression preprocessing.** Stage 09 uppercases/trims gene symbols, sums all source rows with a shared symbol and resolves seven historical symbols: `CXORF36→DIPK2B`, `FAM71A→GARIN4`, `GPR1→CMKLR2`, `H2AFY→MACROH2A1`, `MRVI1→IRAG1`, `RARRES3→PLAAT4`, `WARS→WARS1`. It requires 499 distinct model genes, 492 direct plus seven alias-resolved; raw Ensembl-row selection is not interchangeable with symbol-group aggregation. Source: Stage 09, `MODEL_GENE_ALIASES`, `build_gene_maps`, `collapse_batch`.
2. **Normalizations.** CPM uses the abundance sum across **all valid-symbol gene groups**, not the 499 model genes. Stage 09 also estimates TMM factors on full-symbol-group abundances with log-ratio trim 0.30 and sum trim 0.05, and upper-quartile normalization with explicit evaluability/fallback. These are specified for ARCHS4 pseudocounts. Recount3 G026 input is gene-body **base-pair coverage**: using coverage as counts is an assay-transfer assumption to be tested rather than assumed equivalent.
3. **Frozen scaling.** Stage 10 pan-cancer uses gene-wise weighted anchor means and SDs of `log2(CPM+1)` from the 8,962 frozen cross-sectional anchors. The original SD uses denominator `1 - sum(w^2)`, not the population-variance denominator. Three normalizations have separate frozen anchor scaling tables. Stage 12 uses **different breast-cancer-specific** means/SDs derived from the original cancer anchors (and zero-SD handling). These parameters must be read from existing frozen files; never recomputed using FUSCC or the 143 bridge pairs.
4. **Feature alignment and estimator.** W is 500×10 (`NonSE`, `SE1`–`SE9`), with 499 genes and two direction-labelled features for PPIB. Original Stage 10 checks matching W symbol order, uses `NMF::posneg` for directional features and checks for zero-variance features, then invokes `SpatialEcoTyper::DeconvoluteSE(dat=scaled_matrix, scale=FALSE, W=W, nsample.per.run=10000L, sum2one=FALSE, ncores=1L)`. Input to `DeconvoluteSE` is the **499×samples signed standardized matrix**, not a pre-subset 465-gene model. Raw H coefficients are nonnegative and normalized afterward within sample to 10 fractions. `NonSE` is a model component, not residual.
5. **Reconstruction and robustness.** Build `V` using the original directional-feature transformation; reconstruction is `W %*% H` and residual `V-WH`. Stage 11 evaluates continuous KL, Pearson/cosine/R² and sample-level TMM/UQ sensitivity. Original Tier A/B/C are **within-ARCHS4 empirical grades**, not calibrated cross-assay probabilities.
6. **Score-space distinction.** Pan-cancer fractions use Stage 10 anchored parameters; cancer-specific FUSCC fractions must use the frozen **breast-cancer Stage 12** parameters. Lung-study pair calibration can validate the **pan-cancer** transfer (and, separately, lung-specific transfer) but cannot on its own empirically prove the breast-specific transfer is bias-free.

## The earlier 34-gene count was NOT a full-gene exclusion list

The prior Stage 5 one-to-one Ensembl/symbol join could compare only 465 genes; 23 genes had multiple recount3 symbol candidates and 11 had ARCHS4-side alias/annotation difficulties. The original Stage 09 contract **requires aggregating all source rows for a valid symbol** and resolving seven aliases, instead of selecting one Ensembl ID. The user's prior complete 499-gene membership audit, documented in `README_MANUAL_CURATION.md`, reports **490/499 identical ARCHS4/recount3 Ensembl membership multisets**, with **nine known annotation-row-membership differences**: `FAM20C`, `FCAR`, `ISLR`, `MAL2`, `MATR3`, `PODNL1`, `PSME2`, `SLC39A4`, `TNNT3`. Those nine discrepancies are not nine missing W genes. The original full-membership audit is documented in the README but was not independently regenerated from both complete annotation axes in this archive.

The nine genes account for the following shares of total W loading (sum of nonnegative weights) for each component:

| Component | Share of W loading from nine genes |
|---|---:|
| NonSE | 4.12% |
| SE1 | 1.64% |
| SE2 | 2.42% |
| SE3 | 2.02% |
| SE4 | 0.08% |
| SE5 | 3.85% |
| SE6 | 1.85% |
| SE7 | 1.63% |
| SE8 | 0.89% |
| SE9 | 0.00% |

These are **exposure indicators**, not expected score errors. All nine genes should remain in the intended 499-gene input, with a dedicated source-membership sensitivity study; do not selectively drop them or refit W.

The earlier exploratory 143-specimen matched-lung bridge found 465 same-gene expression profiles with high median correlations, but did **not** demonstrate a complete 499-gene matched fixed-W H/fraction transfer. Several of the nine genes appear in the earlier 465 subset with measurable log2CPM shifts; see `02_nine_gene_bridge_observations.tsv`.

## Files required from frozen results (not part of source-code archive)

- Stage 09: `10_preprojection_normalization_v1_0_1/08_model_gene_availability.tsv` and `04_sample_normalization_QC.tsv.gz` (reference contracts).
- Stage 10: `11_pan_cancer_fixed_anchor_projection_v1_0_4/03_pan_cancer_anchor_scaling_CPM.tsv.gz`, plus `04_...TMM_CPM` and `05_...UQ_CPM`; `09_raw_H_CPM.tsv.gz` and `17_primary_CPM_fraction_scores_with_metadata.tsv.gz` to verify replay/score identity.
- Stage 12: `13_cancer_specific_fixed_anchor_projection_v1_0_3/03_cancer_specific_scaling_parameters_all_branches.tsv.gz`, and `17_primary_CPM_cancer_specific_fraction_scores_with_metadata.tsv.gz`.
- Recount3: `26_recount3_external_primary_v1/stage2_output/studies/SRP074349/selected_run_coverage_G026.h5` for same-specimen calibration and `.../SRP157974/selected_run_coverage_G026.h5` for FUSCC; bridge pairs in `stage5_output/03_bridge_165_availability.tsv`; source G026 row annotations in `stage5_output/00_gene_row_mapping.tsv`.

**Safeguards:** No original Stage 09–13 scripts should be rerun against historical output directories. Prespecify calibration success criteria before evaluating FUSCC ecotype outcomes. Validate source-to-coverage conversion/gene-length effects, complete 499-gene transfer, raw-H agreement, fraction TVD/bias and dominance, directional reconstruction, three normalization branches, and the nine annotation differences. The 334 FUSCC libraries remain provisional; the 26 expression-outlier libraries remain on technical hold pending separate sensitivity.

The accompanying `03_exact_cluster_transfer_contract_READONLY.py` probes existing frozen scaling/output tables and opens both G026 H5 files **only for metadata**, so the actual dataset layout is observed rather than guessed. It writes nothing and performs no deconvolution.
