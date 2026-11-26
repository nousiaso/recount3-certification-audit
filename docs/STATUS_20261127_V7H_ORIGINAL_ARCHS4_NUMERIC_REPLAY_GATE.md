# V7H gate: authentic ARCHS4 original provenance and numerical replay

Date: 27 November 2026.

## Provenance: VERIFIED ON USER CLUSTER

The user ran `cmp -s` plus SHA256 on original, frozen release Stage09–12 files and their checked-in copies under `scripts/originals/frozen_source`. **All four are byte-for-byte identical:**

| Original stage | Version | SHA256 |
| --- | --- | --- |
| Stage09 normalization Python | 1.0.1 | `46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500` |
| Stage10 fixed-anchor pan-cancer R | 1.0.4 | `16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c` |
| Stage11 Layer1 adequacy/robustness R | 1.0.2.1 | `529551d7ecce253b1ab8626956afd076690631bf7d8b51cfcf0c613440a603f5` |
| Stage12 cancer-specific fixed-anchor R | 1.0.3 | `a63750e5bfae63383a397aefb8af94a3e29b17bf3952b2cf47d8c3fc9b858f69` |

Identical file SHA256 proves **source provenance**, not numerical reproducibility. No claim of numerical equivalence yet.

## Critical fidelity differences in proposed V7G

V7G's Stage09 functions mirror the original operations, but ARCHS4's frozen TMM/UQ size factors were estimated over **9,465 canonical samples**, whereas V7G independently estimates factors within **143 matched bridge samples** for each assay. Consequently, V7G's 143-ARCHS4 TMM/UQ inputs are not guaranteed identical to the historical normalized ARCHS4 inputs.

Stage10 also has authoritative gene-symbol and `NMF::posneg` directional-feature/variance checks, plus reconstruction/adequacy and normalization robustness definitions. V7G's independently authored R code must be reconciled with these exact original implementations. Cancer-specific Stage12 is separate from the Stage10 pan-cancer projection and is not implemented in V7G.

**Decision:** Keep V7G on HOLD (no claim of exact original-ARCHS4 method equivalence) until the V7H numerical replay and code reconciliation are reviewed. Do not overwrite prior failed V7C artifacts or frozen sources.

## V7H numerical replay

Files:
- `tools/v7h_replay_frozen_ARCHS4_9465_stages09_10_READONLY_v1.R`
- `scripts/run_v7h_ARCHS4_9465_original_numeric_replay.sh`
- `scripts/slurm_v7h_ARCHS4_9465_original_numeric_replay.sbatch`

The job is fully read-only toward original/frozen data. It first runs the unmodified original Stage09 and Stage10 self-tests, then for the **complete 9,465 canonical samples and all three normalization branches**:

1. Load frozen Stage09 original 499-gene pseudocount matrix, three normalized matrices and exact sample-specific normalization denominators; independently recompute CPM, TMM-CPM and UQ-CPM; require small numerical discrepancies.
2. Reconstruct fixed-anchor standardized expression from original archived Stage10 model-gene means/SDs and original `16A_frozen_W_gene_symbol_harmonization.tsv` gene mapping; compare with immutable 499×9,465 stored scaled RDS.
3. Check the original Stage10 raw H and fraction files are internally consistent.
4. Re-run the exact original Stage10 `SpatialEcoTyper::DeconvoluteSE(dat=Z, scale=FALSE, W=W, nsample.per.run=10000L, sum2one=FALSE, ncores=1L)` on the complete original 9,465 samples in every branch, compare replayed raw H and fractions with frozen originals; report maximum errors, median/95th/max compositional TVD, dominant-component concordance. Do **not** label these numerical differences PASS using sample-order sentinel tolerances without scientific justification.

This exercise independently replays Stage09 *arithmetic from frozen per-gene pseudocounts and stored size factors*, but does **not** independently re-extract all 67,186 ARCHS4 HDF5 gene rows nor re-estimate the full cohort TMM/UQ factors. Those processes were handled in original Stage09 execution QC. Scope is explicit; nothing is assumed beyond the checked source and archived outputs.

SLURM command (after pull):
```bash
sbatch scripts/slurm_v7h_ARCHS4_9465_original_numeric_replay.sbatch
```

Output (aggregate only):
`diagnostics/20261127_v7h_ARCHS4_9465_original_numeric_replay_v1/`

Numerical results have **not** run on the cluster at time of this document; review stdout before any V7G modifications or external assay calibration interpretation.

No source files modified; clinical admissions remain zero.
