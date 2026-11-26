# ARCHS4 identity-input schema and file preflight

Date: 2026-11-26
Status: INPUTS_LOCATED_PARQUET_ENGINE_UNAVAILABLE_IN_SELECTED_ENV

- Stage04, Stage05 and 9,465-GSM atlas inputs found in two locations.
- Stage41 3,602 candidates and 5,221 full-source inputs found once.
- DRP005276 has 103 rows in each Stage41 table; later QC selection is
  83 retained and 20 historical holds, not re-adjudicated here.
- Parquet metadata not opened: pyarrow absent in selected Python.
- No atlas-side accession join, treatment certification or approvals.

Evidence: diagnostics/20261126_v6_archs4_identity_input_schema_v1/
Next: read-only file identity, available reader and membership probe.
