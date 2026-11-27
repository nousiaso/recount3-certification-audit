# V7P v1 — paired frozen ARCHS4 versus recount3 assay calibration using all 500 W directional features

**Date:** 27 November 2026. **Status:** Implementation committed to GitHub, **NOT YET EXECUTED OR NUMERICALLY CERTIFIED ON CLUSTER**. First run static/read-only reuse preflight; inspect stdout before any full SLURM job.

## Why V7P supersedes stopped V7J, without altering original work

V7J's original extra 143-vs-9,465 H invariance check was not scientifically appropriate for different model feature universes. V7M/V7N/V7O established that the exact official SpatialEcoTyper v1.0.3 predictor drops W features with zero variance in the submitted query batch. For the 143 original ARCHS4 matched specimens, this changes the W feature universe from 500 to 496. V7O independently optimized the same convex fixed-W KL objective and reproduced the magnitude of the official 143-vs-9,465 discrepancy from the **500-vs-496 feature change alone**, with numerical stability across two starts.

The original complete 9,465-sample ARCHS4 projection, frozen W, original 499-gene anchor means/SDs, and archived H/fractions remain immutable. V7P does **not** modify or rerun them. It also preserves the earlier stopped V7J and successful V7J 454-SRR→143-SRX gene-body coverage aggregation; there is **no recount3 raw G026 re-extraction**.

## Exact inputs and mathematical contract

- Frozen original ARCHS4 Stage09 CPM, TMM_CPM and UQ_CPM original 9,465-sample normalized matrices and Stage10 499-gene frozen Z and 9465-sample H/fraction outputs. Restrict by exactly matched 143 GSMs only *after* original cohort-level normalization/anchor scaling.
- Reuse already-prepared V7J 143-matched recount3 model-gene CPM/TMM_CPM/UQ_CPM matrices from `private_data/v7j_143_frozen_archs4_vs_recount3_full499_v1`. They derive from original G026 63,856 genes, 454 SRRs grouped into 143 SRXs, mapped to the original 499 W genes.
- Preserve original Stage10 source helpers and frozen gene harmonization ledger. All model genes and all 500 directional W features must be aligned in original order.
- The V7P estimator is **an independent, explicitly identified method**, not `SpatialEcoTyper::DeconvoluteSE`. It uses the original directional `NMF::posneg` transformation, then **retains all 500 W rows, including sample-constant or zero rows** and minimizes nonnegative KL divergence per specimen:
  `h_i = argmin_{h>=0} Σ_g [(Wh)_g - V_{gi} + V_{gi} log(V_{gi}/(Wh)_g)]`
  with `stats::optim(method="L-BFGS-B")`, fixed W, 2000-iteration limit and 1e-8 projected-gradient tolerance. 4 SLURM CPU workers via `parallel::mclapply`.
- V7P **does not** refit W or original anchors or substitute a new fixed-W optimizer for the original archived H. Archived H is the actual ARCHS4 comparison baseline.
- For each normalization branch, recompute 143 original ARCHS4 H using the **independent** full500 estimator only as a fidelity sentinel (warm-started from archived original H). Require max fraction TVD against the archived original H **≤5e-4**. V7O measured CPM replay max TVD 1.515e-4, but V7P's other branches remain untested and must satisfy this gate independently.
- For recount3 projection, use identical full500 KL estimator, original W and matching original frozen branch anchors, with deterministic non-informative uniform start 0.1 for each component. Unlike warm-started archival-reference fidelity test, recount3 does **not** use the matched ARCHS4 H as starting values.
- Full 500-feature adequacy reconstruction and comparison to original archived H/fractions. CPM reference adequacy must replay original Stage10 (tolerance 1e-8). The original helper `calculate_reconstruction_adequacy` is reused verbatim.
- Paired fraction TVD, dominant ecotype concordance, frozen-Z gene expression differences, reconstruction KL/R2/RMSE, TMM/UQ normalization sensitivities, nine known annotation-difference genes and 12-specimen repeat/reverse-order sentinel; all local per-pair and per-gene tables remain in ignored private_data.
- UQ analyzability reflects each assay's explicitly separate evaluability mask. TMM/UQ recount3 factors are 143-SRX cohort-specific, which is part of the technical bridge and not the original 9465-sample normalization.
- No patient bootstrap confidence intervals, mechanistic assertions, external biological validation or clinical sample admissions.

## Files

- `tools/v7p_frozen_ARCHS4_vs_recount3_143pairs_full500_KL_v1.R`
- `scripts/run_v7p_reused_input_full500_static_preflight_READONLY.sh`
- `scripts/run_v7p_full500_frozen_ARCHS4_recount3_calibration.sh`
- `scripts/slurm_v7p_full500_ARCHS4_vs_recount3_calibration.sbatch`

Only proceed through these gates:

1. Run (after `git pull --ff-only origin main`):
   `bash scripts/run_v7p_reused_input_full500_static_preflight_READONLY.sh`
2. Inspect three required PASS markers:
   - `V7P_R_SYNTAX_PASS`
   - `V7P_R_PACKAGES_PASS`
   - `V7P_V7J_PREPARED_INPUT_REUSE_SHA256_PASS`
3. Only after preflight review, submit:
   `sbatch scripts/slurm_v7p_full500_ARCHS4_vs_recount3_calibration.sbatch`
   (normal_prio, no_limit_account, 4 CPUs, 64GB, 8h).
4. Check audit `diagnostics/20261127_v7p_143_frozen_ARCHS4_vs_recount3_full500_KL_v1/01_paired_full500_KL_stdout.txt` and `aggregate_*.tsv` only if job completed all model and computational gates.

### Frozen archive and privacy discipline

No source, archive, model or completed result is overwritten. Exact 143-GSM and 454-SRR run-linked data, raw H and gene-level matrices stay local under ignored `private_data`; GitHub should contain scripts and aggregate-only audit files. Do not treat the matched 143 as independent biological subjects or promote clinical admissions; **clinical admissions 0**.

**Non-execution disclaimer:** the scripts have been committed and source-inspected in GitHub, but have not been parsed or run on the actual cluster. Any preflight failure should be reviewed using exact stdout before modifying versioned code. No new recount3 fraction concordance has been observed.
