# V7C initial execution failure — evidence checkpoint

Date: 27 November 2026.

User-observed execution after Git pull to `5354754ee8a38ba491da2002cbf9f7bd6866f707`:

- Python syntax and required Python modules: PASS.
- R packages in official conda environment: PASS.
- Preparation began; stopped at line 102 of `tools/v7c_prepare_matched_143_G026_499.py` with `ValueError: Expected one unambiguous recount3 run per bridge pair`.
- No normalized 499-gene matrices, R projection, matched H metrics, or `aggregate_transfer_metrics.tsv` were produced in this run.
- This is a validation of the *assumption encoded in V7C* (single run per pair), **not evidence that the source or the bridge is wrong**.
- V7B independently found 143 expected bridge specimens and a 63,856 × 144-run SRP074349 recount3 HDF5, so multiple runs for a single bridge library are plausible but **not yet proven**.

Preserve original failed local outputs at:
`diagnostics/20261127_v7c_srp074349_matched_transfer_v1/`
and any `private_data/v7c_srp074349_matched_transfer_143_v1/`.

Do **not** delete/overwrite these folders or rerun V7C v1 into their immutable paths. Any approved correction must use a new versioned script and output directory.

Next: run `scripts/run_v7c_run_mapping_failure_diagnostic.sh` to examine actual 143 bridge records and 144 HDF5 run/library axes, original Stage2 QC run multiplicity, and frozen Stage09/10 mapping availability. Print only aggregate categories. Then design explicit per-library run aggregation or record-join correction if—and only if—evidence supports it.

No model execution, calibration interpretation or clinical admissions occurred.
