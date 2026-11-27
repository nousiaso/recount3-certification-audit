# Technical calibration — current consolidated report
**Status: completed within the available matched controls, 27 November 2026.** Clinical release remains unsigned.

## Method and rationale

ARCHS4 Kallisto pseudocounts and recount3 G026 GENCODE v26 gene-body base-pair coverage are different processing representations. The matched bridge compares the same source specimens, not independent biological cohorts. Recount3 technical runs are summed within exact SRX before authentic Stage09 symbol collapse and normalization; 499 model genes mapped, with the original frozen ARCHS4 9,465-sample normalized expression, anchor scaling, W and archived H used as reference.

Primary CPM; TMM-CPM and UQ-CPM as assay-specific sensitivity branches. Recount3 factors for TMM/UQ are estimated in the matched study; original ARCHS4 factors retain the historical full-atlas context. No anchor/W refit. The directional input has 500 rows, including zero-variance W features. Independent fixed-W nonnegative KL optimizer (L-BFGS-B) uses the identical W with all rows; it is not the original official SpatialEcoTyper solver. Fractions are normalized H coefficients, not physical cell fractions.

The original replay V7H2 reproduced frozen Z/H to floating-point accuracy. V7O demonstrated that dropping four zero directional W rows in the isolated 143-sample official predictor reproduces the observed 500-versus-496 discrepancy. Full500 archive replay maximum TVD was about 0.000152 in lung and 0.000093 in prostate for CPM. Preserve historical failure records and original archived results.

## Matched controls and agreement

| Study | Cancer | Matched specimens | Runs |
|---|---|---:|---:|
| SRP074349 / GSE81089 | Lung and pleural | 143 | 454 |
| SRP157215 / GSE118435 | Prostate | 22 | 22 |
| Total | Two study contexts | 165 | 476 |

| CPM metric | Lung | Prostate | Pooled |
|---|---:|---:|---:|
| Median fraction TVD | 0.021821 | 0.033576 | 0.023125 |
| P95 TVD | 0.031068 | 0.041359 | 0.036448 |
| Maximum TVD | 0.043395 | 0.049216 | 0.049216 |
| Dominant agreement | 140/143 | 21/22 | 161/165 |

| Pooled sensitivity | Median TVD | P95 | Maximum | Dominant agreement |
|---|---:|---:|---:|---:|
| CPM | 0.02312456 | 0.03644799 | 0.04921636 | 161/165 |
| TMM-CPM | 0.02543063 | 0.04137143 | 0.06329934 | 160/165 |
| UQ-CPM | 0.04441662 | 0.07160037 | 0.1010181 | 155/165 |

Pooled values are the actual specimen-level distribution and are lung-weighted, not study-balanced. TVD = half the sum of absolute component-fraction differences.

## Component bias and model adequacy

| Component | Lung CPM median delta (fraction) | Prostate CPM median delta (fraction) |
|---|---:|---:|
| NonSE | -0.0018331184 | -0.0064122259 |
| SE1 | 0.0053848287 | 0.0129366296 |
| SE2 | -0.0055418859 | -0.0112593054 |
| SE3 | -0.0008523797 | -0.0013923293 |
| SE4 | -0.0019987776 | -0.0015351271 |
| SE5 | -0.0038121937 | -0.0016975748 |
| SE6 | 0.0101634990 | 0.0034406421 |
| SE7 | -0.0019571308 | -0.0039189015 |
| SE8 | 0.0052267765 | 0.0115636137 |
| SE9 | -0.0046482536 | 0.0002078233 |

Deltas are recount3 minus ARCHS4. Pearson range lung 0.99021–0.99968; prostate 0.98323–0.99879. Nine of ten CPM bias directions agree, ten of ten TMM, five of ten UQ; components are compositional and not independent tests. Four discordant CPM classifications become SE8; SE1 and SE4 not represented as original dominant classes.

| CPM median reconstruction metric | Lung ARCHS4 | Lung recount3 | Prostate ARCHS4 | Prostate recount3 |
|---|---:|---:|---:|---:|
| Normalized KL | 0.444886 | 0.479580 | 0.939657 | 0.992584 |
| RMSE | 0.552202 | 0.555746 | 0.537440 | 0.558604 |
| R² | 0.186861 | 0.171201 | 0.109266 | 0.091710 |
| Cosine | 0.775193 | 0.759127 | 0.589818 | 0.569579 |

Prostate reconstruction is poorer in the original reference too; source and cancer/study differences cannot be separated by two studies. Stable inferred coefficients do not imply a complete reconstruction of expression. Nine-gene matched substitution did not consistently improve paired agreement; no correction. Repetition and reversed-order sentinel were exact for tested subsets in each study.

## Decision and limitations

Retain CPM primary, full500 technical estimator, continuous fractions, original frozen assets, TMM/UQ sensitivity; no universal empirical source-bias correction. Use within-study/source-stratified effects and independent replication where appropriate; source adjustment is not identifiable if fully confounded with study. No formal prespecified equivalence margin, no pan-cancer universal validation, no independent biological validation. All matched controls remain calibration-only; external cohort model adequacy and clinical gates remain required.

## Evidence and implementations

- [V7H2 original replay](STATUS_20261127_V7H2_ARCHS4_EXACT_NUMERICAL_REPLAY_PASS.md), [V7O feature investigation](STATUS_20261127_V7O_EXECUTED_500_VS_496_CAUSE_CONFIRMED.md)
- [V7P lung results](STATUS_20261127_V7P_V101_EXECUTED_PAIRED_CALIBRATION_PASS.md), [V7R prostate results](STATUS_20261127_V7R_PROSTATE_RESULTS_PASS_V7S_REVIEW_NEXT.md), [V7S two-study closure](STATUS_20261127_V7S_TWO_STUDY_TECHNICAL_CALIBRATION_COMPLETE.md)
- [Study contract](../configs/v7r_matched_calibration_studies_v1.json), [configured preparation](../tools/v7r_prepare_configured_matched_G026_full499_v1.py), [configured projection](../tools/v7r_project_configured_matched_full500_KL_v1.R), [read-only review](../tools/v7s_two_study_checkpoint_review_READONLY_v1.R)

V7P job 13805260 and V7R job 13805281 completed 0:0; V7S review and both audit checksum manifests passed. Original cluster bytes are not all archived remotely. See [workflow](WORKFLOW.md) and [registry](../AUDIT_REGISTRY.tsv) for provenance. Do not rerun unchanged inputs.
