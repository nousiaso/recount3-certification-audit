# V7J v1 — frozen ARCHS4 versus recount3, same-specimen technical transfer

**Created:** 27 November 2026. **Status: implementation committed; computational execution NOT YET VERIFIED; run read-only static preflight first.** 
This supersedes the design of V7G v1 without rewriting, deleting or silently rerunning previous versions. The original ARCHS4 Stage09-12 inputs/results and frozen W are immutable.

## Why V7J is necessary

V7I proved all **143/143** matched ARCHS4 GSMs are already in the *frozen original* 9,465-sample canonical atlas. Therefore the control side of a paired assay comparison MUST use the actual original Stage09/10 stored 499-gene normalized inputs, original Stage10 frozen standardized Z and original H/fractions, not new TMM/UQ normalization estimated from only 143 samples.

The paired recount3 source is the original SRP074349 G026 gene-body base-pair coverage matrix with **668 SRR columns**, all **454** matching runs (143 distinct SRX/BioSample/GSM pairs). Stage2's 144-run/56-SRX selected HDF5 is entirely disjoint and must not be substituted. Recount3 is a different expression measurement (gene-body bp coverage versus ARCHS4 pseudocounts). Thus applying the same W does not automatically imply assay equivalence.

## Verified prerequisites

- Original Stage09 v1.0.1 SHA256 `46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500`.
- Original Stage10 v1.0.4 SHA256 `16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c`.
- Original Stage11 SHA256 `529551d7ecce253b1ab8626956afd076690631bf7d8b51cfcf0c613440a603f5`; Stage12 SHA256 `a63750e5bfae63383a397aefb8af94a3e29b17bf3952b2cf47d8c3fc9b858f69`.
- Original frozen 9,465-sample Stage09/10 numerical replay **PASS** on all 3 branches: frozen-Z max error ≤4.44e-14, raw-H max error ≤4.44e-15, 100% ecotype dominant-component agreement.
- Original G026 source checksum `b056331eb16d6ec612409fd37a266abfe10c4c15329b486dbb9df213612797c5`.
- Frozen bridge checksum `c7d6d7c754e3ca77d62424f26627a40ceaec10fac98a0741c945773f3d3b1b67`; gene mapping `68be266dd6bd30bfa79a2bd8ac79f0e49c6f5d0f19d2d083e2ca57850ef17391`.
- Recount3 mapping resolves 499/499 original W genes (0 aliases); original ARCHS4 gene model mapping resolves 499/499 (7 aliases).

## Algorithm

1. Python reads the exact original bridge, 63,856-row G026 gene map and source gzip; verifies source SHA256s. It checks all 143 original ARCHS4 HDF5 GSM column identities but **does not extract or re-normalize ARCHS4 abundance**.
2. From the full 668-run G026 matrix, Python streams 63,856 gene rows and uses all 454 matched SRR columns. For each gene it **sums raw gene-body base-pair coverage within SRX**, respecting all run multiplicities (84×1, 4×2, 22×3, 1×8, 32×9). It never averages run-normalized values. A separate per-gene summation-conservation audit requires a maximum relative difference no greater than 1e-12 between all 454 original SRR values and the 143 aggregate SRX values.
3. Invoke original immutable Stage09 Python implementations for gene-symbol harmonization/collapse, full-symbol CPM mass, TMM reference/factors/centering, upper-quartile factors, and all 499 model-gene matrices. Here **CPM is the primary technical bridge**; recount3 TMM/UQ are explicitly *143-SRX cohort-specific sensitivities*, not the original ARCHS4 9,465-sample factors.
4. R loads the actual **frozen original 9465-cohort ARCHS4** 499-gene Stage09 input, Stage10 scaled Z, raw H and fraction outputs, subset by exact 143 GSMs for CPM, TMM and UQ. It reconstructs Z from archived normalized inputs and fixed anchor means/SDs; it must agree with the original stored Z (tolerance 1e-8); original archived H and fractions must be consistent.
5. R uses **verbatim original Stage10 v1.0.4 function bodies** for `harmonize_expression_to_W`, `make_directional_bundle`, `project_raw_h`, `score_fractions`, `calculate_reconstruction_adequacy`, `normalization_sensitivity`, `dominant_component`, and `clr_matrix`. Exact original code source SHA256 is verified; the static preflight checks each function body byte-for-byte.
6. For each normalization branch, as an additional batch-dependence sentinel, project the 143 **archived original** ARCHS4 frozen-Z columns alone using the same original Stage10 code and verify that the resulting raw H and fraction vectors agree with the original 9,465-sample H (maximum absolute H and fraction TVD each ≤1e-8). This check does NOT replace the original frozen 9,465-sample values in any comparison.
7. Original W (500 directional features × 10 components NonSE/SE1–SE9), SpatialEcoTyper **1.0.3**, frozen pan-cancer Stage10 scaling, original `DeconvoluteSE(scale=FALSE,sum2one=FALSE,nsample.per.run=10000,ncores=1)`; no W fit, no scaling refit. Require all **500** original directional features retained for the 143-specimen recount3 projection (fail with HOLD if zero/nonvariable missing features).
8. Replay original ARCHS4 **CPM adequacy** for 143 GSMs using the unchanged frozen Z/H/W and original Stage10 calculation; compare with `19_primary_CPM_fixed_W_model_adequacy.tsv.gz` within absolute tolerance 1e-8. Evaluate recount3 adequacy using the same original function.
9. Save paired within-specimen compositional TVD, component-wise fraction and raw-H differences, reconstruction RMSE/R², normalized KL per input mass, gene-wise frozen-Z differences, dominant-component switches, descriptive p90/p95/max, and original normalization sensitivity. The sample/run-linked tables and RDS stay under ignored `private_data/`.
10. CPM repeat and reverse-order sentinel uses original Stage10 thresholds (max repeat difference 1e-12; reversed median fraction TVD 1e-5; reversed max fraction TVD 1e-3). An explicit nine-gene counterfactual replaces recount3 frozen-Z for nine known source-annotation differences with that specimen's original ARCHS4 frozen-Z and reprojects; it is **not** a gene correction, explanation of causality or alternative model fit.
11. Do not compute sample-independent bootstrap confidence intervals: 143 specimens are **not certified as 143 distinct biological patients**. These are technical assay comparisons, not external biological validation or clinical admissions.

## Code and execution gates

Code (versioned, original V7G untouched):
- `tools/v7j_prepare_recount3_454runs_143pairs_frozen_ARCHS4_v1.py`
- `tools/v7j_project_recount3_vs_original_frozen_ARCHS4_143pairs_v1.R`
- `tools/v7j_static_frozen_code_contract_READONLY_v1.py`

Preflight (run first on login; **no expression loaded**, no model projection):

```bash
bash scripts/run_v7j_static_frozen_transfer_preflight_READONLY.sh
```

This must print `V7J_STATIC_FIDELITY_CONTRACT_PASS`, `V7J_R_SOURCE_SYNTAX_PASS`, `V7J_R_PACKAGE_PREFLIGHT_PASS` and `V7J_PYTHON_PACKAGE_PREFLIGHT_PASS`. The R projection file has **not yet been parsed on supekgate** at creation.

**Only after reviewing that preflight PASS**, submit on a compute node:

```bash
sbatch scripts/slurm_v7j_frozen_ARCHS4_vs_recount3_calibration.sbatch
```

Output is strictly versioned:
- sample-linked local-only: `private_data/v7j_143_frozen_archs4_vs_recount3_full499_v1/`
- aggregate diagnostic stdout/checksums: `diagnostics/20261127_v7j_143_frozen_archs4_vs_recount3_full499_v1/`
- static preflight: `diagnostics/20261127_v7j_static_frozen_transfer_preflight_v1/`

Never overwrite an existing output directory. Preserve V7C, V7G and V7H failure evidence. Do not push local-only sample/run files to GitHub.

## Explicit interpretation boundaries

The ARCHS4 control is numerically identical to its frozen original 9,465-sample workflow by design and direct verification, **not** a newly normalized 143-sample atlas. The recount3 arm uses the same fixed model mathematical and scaling framework, but a different original measurement and cohort-specific TMM/UQ size factors. The paired assay data will decide transfer fidelity; neither a high per-gene cross-sample correlation nor a Stage09/10 replay PASS is evidence that recount3 ecotype fractions agree.

Stage11's original model adequacy metrics are reused in their Stage10 implementation but the complete Stage11 atlas robustness/study outlier process, Stage12 cancer-specific fixed anchors, and all clinical patient/specimen admission gates are NOT being silently represented as complete by V7J.

**Clinical admissions: 0. Model and source files: unchanged. Actual V7J run: pending.**
