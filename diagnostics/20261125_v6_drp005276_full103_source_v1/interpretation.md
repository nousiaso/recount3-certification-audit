# DRP005276 — complete 103-experiment source census (2026-11-25)

Disposition: **FULL_SOURCE_MEMBERSHIP_PASS_QC_HOLDS_IDENTIFIED_CLINICAL_UNSIGNED**.

## Source evidence and frozen membership

- Original DDBJ/SRA metadata: 103 distinct experiments, 103 metadata rows.
- All 103 occur in the frozen V5 5,241-library study inventory.
- 83 selected experiments occur in the frozen V5 2,055-library administrative screen.
- All 20 nonselected experiments occur in V5 inventory but not its primary screen; all 20 have `stage41_partition_NOT_ADMISSION=QC_HOLD_NOT_ADMITTED`.
- No source specimen-metadata conflicts detected across runs.
- Selected and nonselected records are uniformly annotated `tissue=bladder` and `disease=bladder cancer` (83/83 and 20/20, respectively).
- Sample names, SRA sample accessions, original titles and isolate identifiers are unique within each group and share no literal values between groups.
- Study title, abstract, library strategy, source and selection show one distinct populated value apiece at full-study scope.

## Important limits

- The original run-level technical reasons for the 20 historical QC holds **have not yet been retrieved**. Hold membership is established; its scientific justification is not yet audited in this checkpoint.
- Distinct source identifiers are not validated patient identifiers or proof of 103 independent biological specimens.
- The source labels are not specimen-level pathology confirmation, primary-vs-metastatic verification, or pretreatment evidence.
- The 83 selected remain clinically unsigned and must not be automatically admitted. The 20 held remain held; none is reinstated.
- No V5, V6, ARCHS4, W-model, or downstream inference inputs have been changed.

## Reproducibility

The actual executed read-only code and captured stdout are `script.py` and `stdout.txt`; both have recorded SHA256 values in `SHA256SUMS.txt`. The diagnostic verified the hashes of the original SRA metadata, frozen 599-source table, V5 5,241 inventory and V5 2,055 administrative screen. No accession-level or patient-level identifiers are printed in stdout.

Next: check the historical 409-run hold queue and 3,602-library readiness evidence for the **specific QC reasons** behind the 20 held experiments; then establish original study/patient/treatment crosswalk independently.
