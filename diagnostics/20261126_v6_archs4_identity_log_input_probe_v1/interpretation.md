# ARCHS4 identity-join historical log search

Date: 2026-11-26
Decision: ORIGINAL_EXECUTION_LOG_NOT_FOUND_IN_BOUNDED_SEARCH

The original source-package SHA256 and program SHA256 were
verified before inspecting the historical source definitions.

The bounded search visited 759 text files. No file contained
the original identity-join completion marker.

No files were skipped for size, and the configured cap was
not reached. The search does not cover every possible filename,
archive, location or retained terminal session.

Original input definitions identified:
- ARCHS4 Stage04 sample metadata Parquet.
- ARCHS4 Stage05 adjudication metadata Parquet.
- Frozen Stage14 9,465-GSM analysis-ready atlas.
- Stage41 3,602 candidate integrated evidence.
- Stage41 5,221 full-source experiment identities.

The original atlas-wide identity join was NOT executed here.

Next: inspect actual locations, sizes, schemas and source
coverage of these five existing inputs before selecting a
targeted DRP005276 accession comparison.

Cross-atlas patient independence: UNSIGNED.
Pretreatment status: UNSIGNED.
Clinical admissions: ZERO.
Frozen atlas and model: UNCHANGED.
