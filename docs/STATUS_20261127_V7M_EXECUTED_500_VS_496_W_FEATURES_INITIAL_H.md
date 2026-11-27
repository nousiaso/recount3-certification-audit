# V7M executed PASS: exact sources of sample-set noninvariance identified

Date: 27 November 2026. Evidence: user-supplied cluster execution stdout of `scripts/run_v7m_frozen_ARCHS4_initializer_READONLY.sh`, including hashed stdout and script (both OK).

## Direct results from the original 9465-sample frozen ARCHS4 CPM Z

| Diagnostic | Original 9465 columns | Matched 143 columns |
|---|---:|---:|
| Model Z | 499 × 9465 | 499 × 143 |
| Original W | 500 × 10 | same W |
| Auto-log input all nonnegative? | FALSE | FALSE |
| DeconvoluteSE auto log2 branch taken? | FALSE | FALSE |
| Minimum Z | -7.4816981 | -2.38740865 |
| Maximum Z | 20.7261945 | 6.72636633 |
| Retained W directional features | **500** | **496** |
| Original W missing features | 0 | 4 |
| Nonvariable/invalid rows in all posneg directional rows | 0 | 26 |

**Exactly four model-relevant directional features are discarded by the original installed `.nmf.predict` variance filter in the 143-column sample context:** `F12__pos`, `JUND__pos`, `PSMA3__neg`, `SURF4__neg`. The 26 total nonvariable rows include rows not in W; **do not call this 26 missing model features**.

## Random initialization replay from original `.nmf.predict`

The user executed only the original-style initialization `set.seed(39); NMF::rnmf(ncol(W), V)`, separately on the frozen 9465-column versus 143-column directional inputs. No optimization was executed by V7M. For the identical 143 GSMs:

- All 143/143 initial H columns differed by more than 1e-12.
- Maximum initial H absolute difference: **20.1471564896**.
- Median initial H absolute difference: **6.78796993464**.
- Initial normalized-H TVD median: **0.318893550241**, maximum: **0.563862557886**.

These measure **starting coefficients, not final coefficients or final ecotype fraction discrepancies**. Do not interpret the initial median TVD as assay disagreement or the percentage attributable to seeding. Original stopped V7J final H discrepancy was 0.0126240960966 max abs, final fraction TVD 0.00959058827769 max; **no direct causal decomposition between dropped rows, initialization and connectivity stopping has yet been established**.

## Conclusions and methodological significance

The original official SpatialEcoTyper v1.0.3 source has at least two directly demonstrated ways the submitted sample set changes the computation: (i) gene-directional variance filtering reduces the 143 sample model from **500 to 496 W features**, changing the KL objective; (ii) random initial H for identical samples changes when the seed is reused for different matrix shapes/orderings. The original optimizer additionally uses `Stop="connectivity"`, whose contribution remains unquantified.

Consequently, **fixed W and fixed anchor scaling do not guarantee sample-wise identical H in this implementation**. The original ARCHS4 9465-sample projection remains numerically exactly reproducible (V7H2) and should not be overwritten. Requiring new 143-column projection H to exactly equal original 9465-column H would be scientifically inappropriate **without** controlling sample context or explicitly specifying a different inference algorithm. Do not relax a failing tolerance silently.

## Next V7N — read-only impact of four missing original W features

Code: `tools/v7n_original_frozen_143_four_feature_KL_READONLY_v1.R`  
Wrapper: `scripts/run_v7n_original_four_feature_KL_READONLY.sh`

The script reuses original frozen ARCHS4 Z/H and the completed V7J 143 matched GSM manifest only (it never loads recount3 gene expression, runs an NMF optimizer or changes frozen inputs). It verifies exactly which four W features are excluded in the 143 context, reports their observed intensity and frozen W loadings, and evaluates their **KL reconstruction contribution at archived frozen H** relative to all 500 original model features. This is not reoptimization, and therefore cannot alone determine whether feature loss or random initialization dominates the V7J final H discrepancy.

Run after pulling latest main:

```bash
bash scripts/run_v7n_original_four_feature_KL_READONLY.sh
```

Capture the output under `diagnostics/20261127_v7n_original_143_four_feature_KL_v1/01_frozen_four_feature_KL_stdout.txt`.

## Remaining transfer design decision

The intended exact clinical/recount3 calibration should distinguish (a) **frozen historical ARCHS4 ecotype scores** produced in the 9465-sample run, from (b) **assay-matched comparisons within an explicitly controlled common query context**. If the original optimizer is used unchanged, a new 143-sample projection is *not the same numerical procedure* as the original 9465-sample fit even though its W and anchor matrices are identical. The subsequent method must use a prespecified, identical directional feature universe and sample/solver context for both assays, or explicitly disclose any batch/context dependence rather than masking it.

Recount3 454→143 raw aggregation remains PASS; recount3 H/fractions/transfer statistics have not been computed; prior V7J remains stopped and preserved. Full Stage11/12 and clinical admissions remain separate, **zero admitted**.
