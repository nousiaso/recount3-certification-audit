# V7H2 — original ARCHS4 Stage09/Stage10 numerical replay certified

Date: 27 November 2026.

Status: **PASS: complete original 9,465-sample archived normalization/scaling/deconvolution numerical replay (all branches).** This statement applies to the original ARCHS4 source and is NOT a recount3 assay-transfer validation.

Evidence is the user's executed `V7H2_ORIGINAL_ARCHS4_FULL_NUMERICAL_REPLAY_COMPLETE` stdout.

### Locked source code

Original Stage09, Stage10, Stage11, Stage12 are byte-identical to GitHub copies (user-confirmed by SHA256 and cmp in original frozen cluster release):

- Stage09 v1.0.1: `46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500`
- Stage10 v1.0.4: `16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c`
- Stage11 v1.0.2.1: `529551d7ecce253b1ab8626956afd076690631bf7d8b51cfcf0c613440a603f5`
- Stage12 v1.0.3: `a63750e5bfae63383a397aefb8af94a3e29b17bf3952b2cf47d8c3fc9b858f69`

### Numerical replay (9,465 samples × 499 W genes × 3 branches)

| Metric | CPM | TMM_CPM | UQ_CPM |
|---|---:|---:|---:|
| Stage09 max normalized expression absolute error | 4.61180316051468e-05 | 4.99966990901157e-05 | 3.42204235494137e-04 |
| Stage09 max normalized expression relative error | 4.99236050696712e-10 | 4.99743879445086e-10 | 4.97865250779322e-10 |
| Stage10 frozen-Z max absolute error | 3.37507799486048e-14 | 4.44089209850063e-14 | 4.44089209850063e-14 |
| Archived H-derived fraction max absolute error | 1.33226762955019e-15 | 1.83186799063151e-15 | 1.77635683940025e-15 |
| Replayed raw H max absolute error | 1.66533453694e-15 | 1.55431223448e-15 | 4.4408920985e-15 |
| Replayed H max relative error | 1.66533453694e-15 | 1.55431223448e-15 | 3.86309091025e-15 |
| Replayed fraction TVD median | 7.4593109467e-16 | 7.19042880792e-16 | 7.47882658581e-16 |
| Replayed fraction TVD p95 | 1.15601972439e-15 | 1.11369247158e-15 | 1.17093834628e-15 |
| Replayed fraction TVD max | 1.75727488116e-15 | 1.63324215263e-15 | 1.7612864292e-15 |
| Dominant-component concordance | 1.0 | 1.0 | 1.0 |

All three archived harmonization-ledger partitions have 499 rows with index range 1..499. SpatialEcoTyper version is **1.0.3** and W has 500 directional rows.

The replay re-derives normalized expression from archived 499-gene pseudocounts and **archived original per-sample size factors**, reconstructs frozen scaling, and directly reruns `DeconvoluteSE` on 9,465 samples per branch. It **does not** independently rerun the original 67,186-gene extraction or re-estimate the original 9,465-sample TMM/UQ factors. The original Stage09 and Stage10 saved QC statuses are PASS; this is an additional numerical validation, not a new end-to-end extraction certification.

### Remaining transfer gate

The 143 SRP074349 bridge experiments and all 454 original G026 gene-sums SRRs have been verified. However **V7G is still not approved for execution**:

1. The V7G v1 source contains two literal bridge and gene-map expected SHA256 values that were not established from cluster evidence. Run V7I to verify the actual hashes; do not weaken source checks.
2. V7G's TMM/UQ normalization factors are estimated on only 143 paired samples, whereas original ARCHS4 factors were estimated on the full 9,465 canonical atlas. Distinguish frozen original replay from novel cohort normalization. Primary CPM is the first direct comparison; TMM/UQ are explicitly cohort-dependent sensitivity branches.
3. Reconcile original Stage10 gene harmonization, directional-feature variance retention, and Stage11 adequacy definitions before declaring external pipeline methodological equivalence. The W and frozen anchors remain untouched, and the nine source-annotation differences remain an explicit sensitivity analysis.
4. All 143 overlapping bridge specimens remain **same-source technical calibration only**, not 143 independent biological validations or newly admitted clinical samples.

### Immediate next step

Run `scripts/run_v7i_post_v7h2_transfer_contract_READONLY.sh` and inspect its aggregate stdout. It hashes the actual frozen and selected source inputs, validates the 143/454 identity, separately resolves all 499 model genes in the ARCHS4 and recount3 gene annotations, and counts matched GSM overlap with original 9,465-sample atlas without loading expression matrices. Only then revise/freeze V7G under a new immutable version and submit actual calibration.

No changes were made to the frozen ARCHS4 atlas or source files. Clinical admissions: 0.
