# ARCHS4/recount3 historical identity execution locator

Date: 2026-11-26

Decision:
ORIGINAL_PROGRAM_AND_PROVENANCE_IDENTIFIED_EXECUTION_LOG_NOT_YET_LOCATED

## Verified source package

- Original ZIP SHA256:
  20c22c10e07ec8d9bc5847b42c65eae3786614eabdf675afd46ad21b48502c05
- Original program and README SHA256s match the package manifest.
- The source program was inspected as static text, never executed.

## Input provenance

The historical program references:
- Frozen 9,465-GSM ARCHS4 atlas.
- Stage 04 and Stage 05 ARCHS4 GEO `relation` Parquet metadata.
- Stage 41 candidate integrated evidence: 3,602 experiments.
- Stage 41 full-source identity evidence: 5,221 experiments.

Both ARCHS4 `relation` fields derive from the SAME GEO metadata.
Their agreement is an internal consistency check, not independent
source corroboration.

The source excerpt reports reconciliation results through console
printing. No stored crosswalk table was found inside the original ZIP.

Historical 9,465/9,465 accession coverage and overlap counts have
not been independently reproduced by this inspection.

## Remaining work

- Locate any saved original execution stdout.
- Confirm actual Stage04/05 and Stage41 input paths and schemas.
- Use original identity evidence for a targeted DRP005276 join.
- Do not rerun the original atlas-wide program unless necessary.

DRP005276 clinical admissions: ZERO.
Pretreatment status: UNSIGNED.
Cross-study/atlas patient independence: UNSIGNED.

No biological expression data or patient-level records were read.
