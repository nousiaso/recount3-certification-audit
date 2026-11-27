# V7O executed PASS: 500 versus 496 frozen-W feature exclusion explains original H context difference

**Date:** 27 November 2026. **Evidence:** user-supplied actual cluster `01_500_vs_496_KL_stdout.txt` ending `V7O_FROZEN_143_FIXED_W_KL_FEATURE_AND_OPTIMIZER_DIAGNOSTIC_COMPLETE`; `V7O_CONVEX_KL_OPTIMUM_STABILITY_PASS`. No separate `sacct` record or wrapper SHA256 verification was supplied. This report certifies **observed computation**, not a particular SLURM job ID.

## Exact original reference and model

- Same original frozen 143 ARCHS4 specimens in the 9,465-sample atlas, full W 500 directional rows, 499 base genes and 10 components.
- Four original W rows absent after isolated-143 sample-wise variance filtering: `F12__pos`, `JUND__pos`, `PSMA3__neg`, `SURF4__neg`. All four had zero measured direction-specific values, but nonzero W coefficients.
- Original independent V7O reference `DeconvoluteSE` under 143-column context reproduces prior V7J maximum raw-H difference **0.0126240960966** and maximum fraction TVD **0.00959058827769**.

## Controlled fixed-W KL optimization

V7O used *independent* `stats::optim(method="L-BFGS-B")` to solve 143 per-sample convex nonnegative KL objectives, comparing two identical specimens' optimization domains (full 500 and reduced 496) and two different H starts. This is a **diagnostic, not a silently substituted original solver**.

- Total fits: **572** (143 × 2 feature sets × 2 starts).
- Nonzero optimization convergence codes: **0**.
- Maximum objective relative difference across starts: **1.55564217768e-08**.
- Maximum optimized fraction TVD across starts: **4.1355825256e-05** (500), **5.80871206277e-05** (496).
- Controlled 500-versus-496 optimized fraction TVD: median **0.00758648444323**, 95th percentile **0.00890178602414**, max **0.00959013558684**.
- Official 143-versus-original 9,465 fraction TVD: median **0.00758743421332**, max **0.00959058827769**.
- Archived original 9,465 H versus controlled all-500 solution: median TVD **0**, 95th percentile **7.78737584174e-06**, max **0.000151514166957**.
- Official isolated-143 H versus controlled 496 solution: median TVD **5.77887366732e-06**, 95th percentile **2.36709277774e-05**, max **0.000172357811293**.
- Median full-500 versus reduced-496 optimal KL objective difference **1.93152967133**, max **2.93634324723**.

## Inference and limits

**Primary confirmed finding:** changing the directional-feature universe from original 500 to isolated-query 496 is sufficient to reproduce *virtually the entire magnitude* of the historical cohort-context discrepancy. Independent fixed-W KL solutions are stable to the two starts tested, making large random-seed and stopping-rule contributions unnecessary to explain this specific discrepancy.

This is strong quantitative numerical evidence, **not** a rigorous additive per-specimen causal attribution: median/max TVD values are summary statistics, not a paired identity mapping of TVD vectors, and optimizers may differ slightly. No recount3 H or cross-assay agreement has been computed. Preserve the archived official 9,465 H/fractions, the V7J STOP and the original frozen W and scaling anchors.

## Decision: end the original-ARCHS4 diagnostic sequence

No more original-reference-only micro-audits are needed at this stage. **Next scientific goal:** cross-assay paired ARCHS4-versus-recount3 calibration of the exact same 143 specimens.

Prespecify and clearly label a **fixed 500-feature nonnegative KL projection** for both assay measurements, using the unchanged frozen W and frozen Stage10 per-gene anchor scaling. The controlled optimizer can be used only as an explicitly versioned technical bridge after checking it numerically reproduces the historical original ARCHS4 9,465 H/fractions to a tolerance justified by V7O. Do not claim the independent optimizer is the official SpatialEcoTyper function. Keep **archived original H/fractions as ground truth on the ARCHS4 side**.

Reuse V7J successful recount3 preparation (all 63,856 G026 rows streamed, original 454 SRRs summed to 143 SRXs, 56,937 valid symbol groups, 499/499 model genes). Do not re-extract raw G026 or overwrite its successful `private_data` outputs.

The next versioned assay-calibration stage should include matched GSM alignment, CPM primary and TMM/UQ sensitivities, concordance of ecotype fractions and dominant component, actual adequacy and gene-expression correlation versus source, technical controls, and clear distinction between measurement agreement and biology. Report paired sample-level data only to private_data and only aggregates in diagnostics/GitHub.

**Clinical admissions remain zero.** All 143 matched specimens are same-source technical controls, not independent external biological validation.
