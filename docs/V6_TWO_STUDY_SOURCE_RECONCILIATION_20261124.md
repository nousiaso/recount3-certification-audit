# V6 two-study source and patient-crosswalk reconciliation

**Date:** 2026-11-24

**Status:** Evidence reconciliation PASS; clinically UNSIGNED.

## Scope

- 75 historical V5 records across SRP072494 and SRP090849.
- 42 experiments from the previously verified specimen-source review.
- The remaining 33 V5 records were retained for full-context review.

## SRP072494 — glioblastoma

- 36 V5 experiments: 20 baseline, 16 other.
- 20 submitted patient groups, with 16 proposed pairs and 4 singletons.
- Ten selected experiments have SRA primary-diagnosis=yes.
- Ten selected experiments have primary-diagnosis=no.
- The baseline time point is before BEVACIZUMAB, not necessarily
  before previous radiotherapy and temozolomide.
- The original paper states that diagnostic-surgery blocks were used
  where suitable pre-bevacizumab relapse blocks were unavailable.
- Patient grouping is supported by study design but remains unsigned.

Source:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79671

Publication:
https://doi.org/10.1186/s12885-017-3251-3

## SRP090849 — osteosarcoma

- 39 experiments represented in the historical V5 ledger.
- 22 administrative primary-screen experiments.
- Nine technical-QC holds and eight outside Stage 41.
- The 22 selected experiments have 17 submitted primary-tumour-type
  labels and five unknown-type labels.
- The historical full V5 subset includes eight METASTATIC_CONTEXT
  grades, not automatically certified sampled metastatic lesions.
- All 39 curated patient keys are unique, but patient independence
  is not certified.
- The original GEO series lists 52 samples, including tumour tissues,
  normal bone and cell-line controls. The 39 V5 experiments must NOT
  be presented as the complete original GEO source census.
- Exact source-census membership and treatment timing remain open.

Source:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87624

## Technical and scientific limitations

- No individual clinical admission has been granted.
- Submitted specimen-type labels are evidence, not pathology signoff.
- The absence of culture-related keywords is not proof that a
  specimen is patient tissue.
- No run QC, cross-study patient independence or full-W model-transfer
  validation was performed.
- The original V5 and frozen ARCHS4 files are unchanged.

## Provenance

Script: `tools/v6_two_study_full_universe_crosswalk_v1.py`
Script SHA256: `8b786f9c0dc629886492302378fcc06ce01d35e37e6b2898eb374361a7e53360`

Stdout: `diagnostics/20261124_v6_two_study_full_crosswalk/stdout.txt`
Stdout SHA256: `a2f2966df35c617d84454d1b51dbf5f1f35ef93853cc8e3a44dc7116be71b39f`

Aggregate summary: `diagnostics/20261124_v6_two_study_full_crosswalk/summary.json`
Summary SHA256: `9e96a56282fbcce3f6d052637da94c26846739eeb04809633d6c806352a9bee6`

Accession-level crosswalk: `diagnostics/20261124_v6_two_study_full_crosswalk/75_full_study_crosswalk_LOCAL_ONLY.tsv` (cluster only)
Local detail SHA256: `d33b0ac431284a8369140e6e6ff31e43ef3995558cb911179cc4218956168dce`

## Disposition

SOURCE_CROSSWALK_COMPLETE; SPECIMEN_AND_TREATMENT_ADJUDICATION_OPEN.

Do not rerun the 75-record join unless source files or the source
interpretation contract change.
