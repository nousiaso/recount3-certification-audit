# V7L executed: original SpatialEcoTyper fixed-W solver has cohort-context dependencies

Date: 27 November 2026. Evidence: user-supplied actual cluster V7L source-inspection stdout in Pasted markdown(20261127-112102).md, including the full 69-line installed namespace implementation.

## Verified installed implementation

- SpatialEcoTyper v1.0.3, R 4.4.3, namespace-local `.nmf.predict(W,testdat,scale=FALSE,normalize=FALSE)`.
- Function starts with `set.seed(39)`, then converts inputs to matrices and uses original W genes.
- On directional W, constructs `NMF::posneg(testdat)`; **drops every directional row with zero variance across the samples submitted to that invocation** before intersecting with W. Thus the features retained may differ between the complete 9,465-column original ARCHS4 input and the 143-column subset even with identical per-sample expression and W.
- Initializes H with `dummy = NMF::rnmf(ncol(ws), to_predict)` followed by `coef(model) <- dummy@H`. A fixed RNG seed guarantees reproducible draws for a given shape/order, not invariant H initialization for a different sample set or column ordering.
- Uses a custom `NMFStrategy("my-method","brunet",Update=my_method,objective="KL",Stop="connectivity")` and `NMF::nmf(...,nrun=1,method=nmf_method,seed=my.seeding.method,.opt="P1")`. W is injected as the basis; a KL-divergence update adjusts H. Connectivity-based stopping can depend on the sample ensemble, and requires numerical tests before making a definitive causal attribution.
- No `NMFpredict` explicit 10,000-sample batching for either 143 or 9,465 samples.
- V7K exposed a possible auto-log `if(all(dat>=0) && max(dat)>80) log2(dat+1)` in the wrapper; its actual result on original frozen Z remains to be checked numerically.

**This is now concrete source evidence that the package does NOT guarantee per-sample invariance.** It is not yet direct evidence of the relative contributions of random initializer, varying feature set, and connectivity-based stopping to the V7J measured max-H difference (0.0126240960966) and max fraction TVD (0.00959058827769).

## Scientific implications and continuing hold

- V7H2's original full 9,465 ARCHS4 numerical replay remains valid and essentially exact across all CPM/TMM/UQ branches. Original W, anchors, package version and output remain locked.
- V7J recount3 G026 raw extraction and 454→143 run aggregation passed; recount3 model H has **not** been computed, as V7J deliberately stopped at the prior original 143-vs-9,465 frozen-H invariance gate.
- Do not drop the invariance gate silently, alter the reference H, interpret paired recount3 ecotype scores, or treat source identity as external validation. The 143 same-source matched specimens are assay-calibration controls, not independent biological cohorts.

## V7M: small numerical source-of-noninvariance test

Files (versioned; preserve all old results):
- `tools/v7m_original_ARCHS4_9465_vs_143_initialization_READONLY_v1.R`
- `scripts/run_v7m_frozen_ARCHS4_initializer_READONLY.sh`

The R diagnostic reads the already-frozen 499×9,465 CPM Z and 143 matched IDs from the successful V7J preparation only; **no raw G026 extraction or NMF optimization**. It reports:

1. Actual `DeconvoluteSE` auto-log condition for frozen 9,465 versus 143 inputs.
2. Directional-feature variance retention on identical W-ordered model genes in both sample sets, including exactly which W features are absent from the 143 subset (gene-feature names only; no sample IDs).
3. The NMF::rnmf starting H (seed 39) for 9,465 versus 143 using the exact original `.nmf.predict` seeding call, comparing starting H only for the same 143 samples with same identities.

Only after inspecting these outputs can we distinguish a feature-universe change from different initialization. A later tightly controlled experiment may be needed to isolate the effect of connectivity stopping.

Run:

```bash
bash scripts/run_v7m_frozen_ARCHS4_initializer_READONLY.sh
```

This is read-only to frozen/reference/recount3 data and writes only versioned aggregate diagnostic stdout/checksums. No clinical admissions.
