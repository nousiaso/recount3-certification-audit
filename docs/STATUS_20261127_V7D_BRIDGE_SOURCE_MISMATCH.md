# V7D — SRP074349 cross-assay calibration input membership checkpoint

Date: 27 November 2026.

The **executed cluster stdout** establishes:

- Bridge: 143 source records, 84 single-run and 59 multirun.
- Run multiplicity: 84 × 1, 4 × 2, 22 × 3, 1 × 8, 32 × 9 = **454 run references** (not necessarily 454 distinct runs).
- All 143 records are `CALIBRATION_ONLY_NOT_IN_PRIMARY`.
- Available `stage2_output/studies/SRP074349/selected_run_coverage_G026.h5`: **144 runs from 56 HDF5 library identifiers**.
- Exact library-ID overlap bridge → HDF5: **0**.
- Original Stage5 gene mapping agrees in all 63,856 ordered gene rows; 499/499 model genes resolve. The 0 observed alias-resolved count in this recount3 source must not be conflated with the 7 historical ARCHS4 alias resolutions.
- Frozen Stage10 H and anchor/scaling files exist.
- No normalized recount3 matched inputs, projections or aggregate calibration metrics were produced in the failed V7C attempt.

**Clarification of prior checkpoint:** the earlier diagnostic's final message `FAILURE_REPRODUCED_AT_BRIDGE_CARDINALITY_ONLY` is too narrow. The executed output also shows a complete library-ID mismatch and a larger bridge run universe than the available selected HDF5. Thus the baseline is **NOT YET MATCHED TO EXPRESSION**.

**Do not fix** by assuming 143 bridge rows correspond to 144 HDF5 columns, reordering, merging arbitrary run columns, or dropping specimens. An existing HDF5 may only represent the separately selected non-overlap cohort rather than the calibration-only bridge, but this explanation is **a hypothesis** until exact run membership is measured.

Next diagnostic: `scripts/run_v7d_srp074349_bridge_expression_identity.sh`. It compares exact run-token intersection, run/experiment/BioSample fields in Stage2 and Stage5, and original SRA study metadata without printing sample-linked IDs. It is strictly read-only and preserves failed V7C outputs.

After that evidence, decide whether the 143 matched specimens already have accessibly extracted recount3 coverage elsewhere or whether a separately approved, versioned extraction of bridge run data is necessary. No extraction, normalization, calibration or biological validation is authorized by this diagnostic.
