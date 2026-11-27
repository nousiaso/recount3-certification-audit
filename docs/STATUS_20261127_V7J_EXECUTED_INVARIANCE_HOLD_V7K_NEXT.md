# V7J executed: recount3 extraction PASS; frozen 143-versus-9465 H invariance HOLD

**Date:** 27 November 2026. **Evidence:** user-supplied original cluster stdout from executed V7J compute job.  
**Scientific decision:** preserve the stopped V7J run and **do not** proceed with recount3 projection or compare assay scores until the installed `SpatialEcoTyper::DeconvoluteSE` 1.0.3 code has been inspected for sample-set dependence.

## Recount3 preparation: PASS

- `BRIDGE_143_EXPERIMENTS_454_RUNS_PASS`
- Exact 143 GSM positions in ARCHS4 original H5: PASS.
- ARCHS4 reference retained as original frozen 9,465-sample Stage09/10 matrices; not renormalized.
- Original G026 file: 668 SRR columns; header found at 0-based line 2; all 454 matched SRRs present.
- All 63,856 G026 gene rows streamed; 454 technical runs summed to 143 matching SRX experiments.
- `MAX_RELATIVE_RAW_COVERAGE_AGGREGATION_MASS_ERROR 0`.
- 56,937 valid symbol groups, 499 W model genes, 0 recount3 aliases, 143/143 UQ-evaluable samples.
- `V7J_RECOUNT3_ONLY_143_MATCHED_FULL499_PREPARATION_PASS`.
- No input/frozen/W modification, no clinical admissions.

## Frozen ARCHS4 reference replay and invariance: mixed

- For CPM, frozen `Z` replay max absolute error `2.04281036531e-14`: **PASS**.
- Frozen Stage10 gene harmonization: 499 rows **PASS**.
- Original ARCHS4 archived H-to-fraction self-consistency max absolute error `8.881784197e-16`: **PASS**.
- **The same 143 ARCHS4 frozen-Z columns projected alone with unchanged official W using the same `DeconvoluteSE` call DO NOT exactly reproduce their original coefficients from the 9,465-sample projection:**
  - `MAX_H_ABSOLUTE_ERROR 1.26240960966e-02`.
  - `MAX_FRACTION_TVD 9.59058827769e-03`.
  - Stopped at `Frozen original 143-column projection differs from the archived 9465-sample H; HOLD CPM`.

Interpretation: this is a query-batch/reference-batch numerical discrepancy of the *same measured ARCHS4 specimens*, not a recount3 result, and does not by itself prove any defect in W/anchors/expression. The precise mechanism (software preprocessing, sample-wise vs cohort-wise fitting, initialization, algorithmic convergence etc.) remains **unknown**; do not speculate that it is already proven to be one cause. Since original full 9,465-sample H was reproduced to ≤4.44e-15 in V7H2, this discrepancy must be investigated before external assay comparisons.

The following were **not executed** by V7J: recount3 projection, TMM and UQ frozen-reference replay (stopped in CPM), paired ecotype TVD, adequacy comparisons, normalization sensitivity, nine-gene swap, repeat/order sentinel, aggregate transfer TSV generation. It is expected that `aggregate_*.tsv` are absent. Do not call V7J a transfer calibration PASS.

## Immediate diagnostic: V7K static installed-source inspection

The next step is **read-only code inspection**, not a new full projection.

- `tools/v7k_inspect_installed_DeconvoluteSE_SOURCE_READONLY_v1.R`
- `scripts/run_v7k_installed_DeconvoluteSE_source_READONLY.sh`

The diagnostic prints the **actual installed** `SpatialEcoTyper::DeconvoluteSE` 1.0.3 function body, arguments and direct local helper signatures and checksums the stdout. It does not load expression values or execute any model.

Run after pulling latest `main`:
```bash
bash scripts/run_v7k_installed_DeconvoluteSE_source_READONLY.sh
```

Output:
`diagnostics/20261127_v7k_DeconvoluteSE_installed_source_v1/01_installed_DeconvoluteSE_source_stdout.txt`.

**No changes to V7J code or reruns until the observed installed model code tells us which minimal, evidence-grounded numerical diagnostic to perform next.** Preserve the existing V7J partial data/logs and the failed audit directory. New versions, not patched-in-place results, are required for any eventual calibration rerun.

**Clinical admissions: zero.** All 143 overlap specimens remain technical-calibration controls.
