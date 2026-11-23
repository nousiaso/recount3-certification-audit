# ARCHS4 SpatialEcoTyper — original execution contract audit

**Audit version:** v1.1, 2026-11-19; corrected and source-adjudicated against the actual frozen Stage 11 script  
**Scope:** Frozen Stage 09 v1.0.1 → Stage 10 v1.0.4 → Stage 11 v1.0.2.1 → Stage 12 v1.0.3 → Stage 13 v1.0.3.  
**Scientific rule:** The completed ARCHS4 cohort and analyses are frozen. This audit establishes their reference contract for future recount3 adaptation; it is not a proposal to rerun or change ARCHS4 analyses.

## Verification state and provenance

- Stage 09 and 10 source files were obtained from the user's **uploaded frozen-source archive** and inspected directly.
- Stage 11, 12, and 13 source files were retrieved from the user's **file Library** and reviewed. Stage 09, 10, 12 and 13 reviewed source hashes matched the original frozen copies in the user's 2026-11-19 cluster audit. Stage 11 initially differed by one executable argument name; the actual frozen Stage 11 script was then uploaded and compared byte-for-byte. This v1.1 version uses the actual frozen Stage 11 SHA as authoritative.
- The 2026-11-19 cluster run of the v1 output validator reported PASS status and zero QC failures in Stage 09–13, all expected core rows/scaling counts/tier counts, and exact Stage 13 pan/cancer-specific score replay (max absolute difference 0). Its only HOLD was the mismatched Stage 11 reference hash, which is now explained. **The v1.1 validator has not yet been run on the cluster**; it adds the corrected SHA and an exact check of the five study-level outliers.
- Exact historical **job-log execution provenance** is not established solely from current frozen source hashes and output consistency. This certification establishes inspected source and output contract, not cryptographic proof of which bytes the historical scheduler executed.
- The W gene harmonization/aggregation rule is additionally supported by the earlier 499-gene source-row membership audit: 490 identical G026/ARCHS4 row-ID multisets; 9 known differences. This is feature identity, **not assay equivalence**.


### Stage 11 source-provenance adjudication (added v1.1)

- Reviewed Library copy: SHA256 `cec6fd5143af5ceaecedbcd6b567a8687c13bc9ffca7fea9049b61ed17b77153`.
- Uploaded **actual frozen** copy: SHA256 `529551d7ecce253b1ab8626956afd076690631bf7d8b51cfcf0c613440a603f5`.
- Both 71,616 bytes and 2,728 lines; **one changed line / two changed bytes**: line 289 `stats::mad(..., centre = centre, ...)` in the reviewed copy versus `stats::mad(..., center = centre, ...)` in the actual frozen source. R `stats::mad()` expects `center`; the frozen source is the valid call.
- Only `robust_z()` is affected. This function is invoked by the Stage 11 self-test and by the within-cancer study-median adequacy outlier calculation (`adequacy_robust_z_within_cancer`, threshold -2.5, minimum five samples per study and at least three studies per cancer). The sample-level 4,200/1,026/4,239 tier rule and Stage 10/12 deconvolution are unchanged.
- A file with the typo would fail if the function were invoked; do not describe the prior reviewed copy as equivalent runnable code. The frozen source has the correct call.
- The observed frozen output status was PASS with five reported poor-fit study outliers; v1.1 validator adds a specific cross-check of the exact five names from the frozen Stage 11 README. No rerun of Stage 11 is indicated.

## 1. Input preparation — Stage 09

**Script:** `09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py`.

Inputs: authoritative Stage 08 patient/canonical lock, HDF5 `human_gene_v2.latest.h5`, and official frozen W feature list. Strict 9,465 canonical samples, 8,962 cross-sectional anchors, 18 cancer classes and 499 model genes. Original HDF5 contains 67,186 gene rows and 1,093,742 sample columns. The script verifies exact HDF5 sample-accession order against locked zero-based column indices.

Gene handling:
1. Clean non-finite/negative expression values; fail on any negative/nonfinite in the locked input.
2. Uppercase and trim HDF5 `meta/genes/symbol`; omit blank/invalid names.
3. **Sum all source rows within each normalized gene symbol** with `np.add.at`, not select a single canonical Ensembl row.
4. Seven explicit historical aliases when an old W name is absent: `CXORF36→DIPK2B`, `FAM71A→GARIN4`, `GPR1→CMKLR2`, `H2AFY→MACROH2A1`, `MRVI1→IRAG1`, `RARRES3→PLAAT4`, `WARS→WARS1`; refuse ambiguous concurrent old/new symbols.

For sample `s`, the denominator is the sum of abundance over **all valid-symbol collapsed genes**, not just 499 W features:

`L_s = Σ_{g in all valid collapsed symbols} X_gs`.

Primary `CPM_gs = 10^6 X_gs / L_s`.

The TMM branch computes an in-cohort reference from the mean upper-quartile fraction, then performs trimmed M-value/abundance-weighted effective library size estimation (`--tmm-logratio-trim=0.30`, `--tmm-sum-trim=0.05`); factors are geometrically centered to one, and `TMM_CPM_gs = 10^6 X_gs/(L_s f_s^TMM)`.

The UQ branch defines `q_s = 75th_percentile(all valid-symbol collapsed genes) / L_s`, imputes nonpositive q by the positive-sample median for *diagnostic projections*, geometrically centers q factors and calculates `UQ_CPM_gs = 10^6 X_gs/(L_s f_s^UQ)`. **The 163 UQ fallback profiles are non-evaluable for formal UQ robustness, despite having UQ diagnostic scores.**

Stage 09 expression numbers are ARCHS4/Kallisto **pseudocount abundances**, not ordinary reads or TPM. This matters when recount3 contributes gene-body base-pair coverage.

Main frozen outputs: `10/11/12_model_gene_{CPM,TMM_CPM,UQ_CPM}_canonical.tsv.gz` (actual files use explicit numbered names); `13/14/15` cross-sectional anchor matrices; `04_sample_normalization_QC.tsv.gz`; `18_preprojection_normalization_status.json`.

## 2. Shared pan-cancer fixed-anchor projection — Stage 10

**Script:** `10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R`.

For pan-cancer anchor sample `s` belonging to cancer `c` and study `t`, the frozen anchor weight is

`w_s = 1 / (n_cancers × n_studies_in_cancer_c × n_samples_in_study_tc)`, with `Σ_s w_s = 1`.

For every normalization branch **separately**, define `y_gs = log2(CPM_branch_gs + 1)` and frozen weighted anchor parameters

`mu_g = Σ_{s in frozen anchor} w_s y_gs`

`sd_g = sqrt( Σ_s w_s (y_gs - mu_g)^2 / (1 - Σ_s w_s^2) )`.

Every canonical profile is scaled with `z_gs = (y_gs - mu_g)/sd_g`, retaining the frozen mean and SD when transferring to another dataset; do **not** refit them on recount3. For CPM and TMM: 8,962 anchors. UQ: 8,799 evaluable anchors, excluding 163 canonical UQ fallback profiles; the three branches are **different score spaces**.

Input genes are reordered to match the frozen W base genes. The `input_row_index` in the harmonization table is a processed, 1-based R matrix row number, **not an original HDF5 coordinate**.

The signed scaled gene matrix is converted to nonnegative directional features using `NMF::posneg`, labeled `__pos` and `__neg`. Rows with zero variance are removed and the remaining features intersected with the frozen **500 directional-feature × 10 component** W matrix. The 499 input genes are not 499 independent directional rows.

Exact call verified from Stage 10:

```r
raw_h <- SpatialEcoTyper::DeconvoluteSE(
  dat = scaled_matrix,
  scale = FALSE,
  W = W,
  nsample.per.run = 10000L,
  sum2one = FALSE,
  ncores = 1L
)
```

The script retains nonnegative **raw H coefficients** and derives compositions `p_sk = H_sk / Σ_k H_sk`, with fixed component order `NonSE, SE1,...,SE9`. `NonSE` is a **published model component**, not an unexplained residual. Dominant ecotype is the largest fraction; ties use the first component.

For reconstruction, the *directional input* is `V`, reconstruction `Vhat = W H`, residual `R = V - Vhat`; positive signed residual denotes observed signal stronger than reconstructed. Model adequacy includes normalized generalized KL divergence, Pearson correlation, cosine similarity, RMSE, MAE and `R²=1-||V-Vhat||²/||V-mean(V)||²`. For directional KL, numerator is `Σ[V log(V/Vhat) - V + Vhat]` with zero-value handling and `Vhat` floor; normalize by `ΣV`.

Normalization sensitivity measures per-sample total variation distance `TVD = 0.5 Σ_k |p^CPM_sk - p^alternative_sk|`, maximum absolute component difference, Aitchison distance of CLR compositions (floor 1e-8), and dominant-component agreement. The code distinguishes diagnostic metrics from formally evaluable UQ metrics.

Frozen output verification targets: `03/04/05_pan_cancer_anchor_scaling_*`; `17_primary_CPM_fraction_scores_with_metadata.tsv.gz`; `19_primary_CPM_fixed_W_model_adequacy.tsv.gz`; `22_normalization_sensitivity_by_sample.tsv.gz`; `27_pan_cancer_projection_QC_contract.tsv`; `28_pan_cancer_projection_status.json`.

## 3. Layer 1 adequacy / robustness — Stage 11

**Script:** `11_build_ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1.R`.

**No new deconvolution and no hypothesis testing.** Creates relative, descriptive *within-ARCHS4* evidence tiers and study fit cautions. For the 9,465 samples, calculate four empirical percentile ranks (lower normalized KL is better; higher Pearson, cosine and R² better), then average and rank this composite.

TMM stable = `TVD(CPM,TMM) ≤ 0.10` **and** same dominant component. UQ stable = UQ formally evaluable, `TVD(CPM,UQ) ≤ 0.15`, same dominant component. High sensitivity includes TMM TVD ≥ 0.25, UQ TVD ≥ 0.35, or both alternatives disagreeing in dominance (where UQ evaluable). Extreme scaled-Z warning uses `|z|>8` and severe diagnostic uses `|z|>10` on at least 5 features (check exact comparison in frozen code if reporting boundary equality).

Tier A requires global fit percentile >=0.25, stable TMM, evaluable and stable UQ, no severe Z flag. Tier B requires fit percentile >=0.10, stable TMM, not highly sensitive, acceptable UQ state (or non-evaluable), no severe Z. All others Tier C. **These percentiles and tiers are empirical ARCHS4 policy, not externally calibrated validity probabilities.**

Expected frozen tier counts: A=4,200; B=1,026; C=4,239. It also computes 5 poor-fit cancer-study outliers in the frozen cohort and residual-feature summaries. `NonSE` is never used as the unexplained-expressions fraction.

Key outputs: `01_sample_level_adequacy_robustness_atlas.tsv.gz`, `03_study_level_transferability_summary.tsv`, `12_layer1_adequacy_robustness_QC_contract.tsv`, `13_layer1_adequacy_robustness_status.json`.

## 4. Cancer-specific fixed-anchor projection — Stage 12

**Script:** `12_run_ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3.R`.

For each of 18 cancer classes, build separate within-cancer study-balanced anchor weights:

`w_s^(c) = 1 / (n_studies_c × n_anchor_samples_in_study_cs)`, normalized within cancer.

Cancer-specific anchor statistics and `log2(branch_CPM+1)` transformations follow the same weighted mean and unbiased variance formula **but with cancer-specific anchor samples and weights**, producing distinct cancer-specific coordinate systems. Zero/invalid within-cancer anchor SDs are set to zero for those gene rows in this stage, not treated identically to the pan-cancer branch. Directional feature filtering follows the original code, requiring at least 450 retained directional features and at least 0.90 of each component's retained W weight. The original W is unchanged.

Computes cancer-specific raw H, normalized `NonSE/SE1...SE9` fractions, adequacy and normalization sensitivity, comparison with pan-cancer coordinates and reconstruction residual attenuation/persistence. Do **not** compare cancer-specific score values directly across different cancer classes.

Key outputs: `03_cancer_specific_scaling_parameters_all_branches.tsv.gz`, `17_primary_CPM_cancer_specific_fraction_scores_with_metadata.tsv.gz`, `18_primary_CPM_cancer_specific_fixed_W_model_adequacy.tsv.gz`, `20_normalization_sensitivity_by_sample.tsv.gz`, `23_pan_vs_cancer_specific_coordinate_shift_by_sample.tsv.gz`, `31_cancer_specific_projection_QC_contract.tsv`, `32_cancer_specific_projection_status.json`.

## 5. Release use policy / analysis-ready atlas — Stage 13

**Script:** `13_build_ARCHS4_projection_use_analysis_ready_atlas_v1_0_3.R`.

No new projection. Joins Stage 09–12 frozen outputs into a 9,465-row core; retains **separate** `pan_` and `cs_` ecotype scores, reconstruction results, normalization robustness and fit ranks. Pan-cancer score space is `pan_cancer_shared` and is intended for between-cancer use; `cancer_specific_<cancer>` is for within-cancer use only.

Pan-cancer Level 1 tiers are carried from Stage 11; a separate cancer-specific empirical fit/robustness tier is derived for within-cancer interpretation. Strict `cross_cancer_primary_eligible_strict` means pan Tier A and no pan poor-fit study outlier; strict `within_cancer_primary_eligible_strict` requires cancer-specific Tier A, no cancer-specific poor-fit study outlier, and not the heterogeneous `other_solid_tumour` class. A supported-use manifest includes Tier A or B under corresponding restrictions.

Key outputs: `03_projection_analysis_ready_core.tsv.gz`, `05_pan_cancer_primary_manifest.tsv.gz`, `06_pan_cancer_supported_manifest.tsv.gz`, `07_within_cancer_primary_manifest.tsv.gz`, `08_within_cancer_supported_manifest.tsv.gz`, `12_projection_use_QC_contract.tsv`, `13_projection_use_status.json`.

## Critical adapter restrictions for recount3

1. The validated source data are ARCHS4 **Kallisto pseudocount abundances**; recount3 yields **G026 gene-body base-pair coverage**. The same algebra on different measurements need not transfer ecotype scores. Cross-assay equivalence requires direct matched-specimen **model-output** tests, not only expression Pearson correlations.
2. Reconstruct 499 model features using the original uppercase-symbol collapse and seven alias rules; track all nine different row-ID memberships. Native G026 source rows and shared-ID sensitivity variants must remain explicitly distinguished.
3. Normalize using **all valid gene symbols**, never just 499. The denominator is assay-dependent; preserve branch labels and UQ non-evaluability.
4. Apply **frozen ARCHS4 anchor means/SDs**, rather than estimating recount3 reference means/SDs: otherwise the coordinate system changes. The pan and cancer-specific anchor transformations and the three normalization branches must remain separate.
5. Preserve **raw H** and normalized H fractions, W dimensionality and `DeconvoluteSE` settings. Never equate `NonSE` with reconstruction residuals.
6. Do **not** reuse ARCHS4 empirical percentile/tier thresholds as independently validated for recount3. Compare with the frozen reference and explicitly study calibration/adequacy.
7. No statistical integration, pooled count matrices, training, or new biological claims until the recount3 cohort, patient independence, overlap and transfer gates pass.

## What the companion output validator establishes

`audit_frozen_outputs_READONLY.py --base /.../ARCHS4_V2_6 --output /.../OUT.json` checks Stage 09–13 status/version/QC files and source hashes; expected row counts of original manifests, anchor/score/scaling tables; unique and ordered 9,465 sample IDs; nonnegative sum-to-one components; Stage 13 replay of Stage 10/12 score fields; Stage 11 empirical tier counts; and Stage 13 score-space contract columns.

It **does not** read the huge HDF5 expression matrix, load large RDS directional residuals, rerun the statistical or projection analysis, or assert recount3 equivalence. Missing status or mismatched fields are explicit HOLD conditions. The v1 cluster output validator passed all checks apart from the reviewed Stage 11 source SHA. The source discrepancy has now been adjudicated. Run the v1.1 validator once to issue a clean machine-readable reference-output contract and verify the five study outliers. This is not new ARCHS4 projection or independent recount3 validation.
