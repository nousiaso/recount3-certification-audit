# V6 — 599-experiment integrated clinical evidence

Date: 2026-11-24
Status: SOURCE INTEGRATION PASS; CLINICALLY UNSIGNED.

## Scope

All 599 historical-primary-unattested experiments, across 39 studies,
have undergone first-pass original-SRA value inspection.

The remaining 28-study diagnostic completed: 217 experiments.
The earlier 11-study review covered the other 382 experiments.

## Integrated source reconciliation

- 599/599 exact experiment records joined.
- Zero SRA-to-V5 sample-accession conflicts.
- 571 records with recognized SRA specimen attributes.
- 441 with recognized SRA disease attributes.
- 24 with recognized SRA timing/treatment attributes.
- 75 with recognized SRA patient tokens.
- 599 nonempty V5 patient_key values.
- 599 nonempty V5 treatment_timing values.

The last two counts indicate populated curated fields, NOT
599 verified independent patients or treatment-naive specimens.

## Original evidence provenance

Remaining-28 script: `tools/v6_remaining_28_studies_values_READONLY_v1.py`
Script SHA256: `e90ca5ff7c18db5bb3c05026bc81c01fb140a702e88a62bceb2361de4450cc98`
Cluster stdout: `diagnostics/20261124_v6_remaining_28_studies_values/stdout.txt`
Stdout SHA256: `cbedd7ed607ea28240a54ddf1b6e3138dd905a926e31a61f53557f0c7488c13c`

Integrated script: `tools/v6_599_integrated_evidence_review_v1.py`
Script SHA256: `5cba0a1d857cd7cba57cb117c8a8396a5d5b22203d6e475d45318a6bd87f79bf`
Cluster stdout: `diagnostics/20261124_v6_599_integrated_review/stdout.txt`
Stdout SHA256: `e262d30b0d4313d28fc5265323fbfa7671d19f2ca313af656a483c89ccc948f3`

Integrated summary: `diagnostics/20261124_v6_599_integrated_review/summary.json`
Summary SHA256: `ac6d2d42065cc98bdf9a88ad541964ba43621d6b79410cfcb4c673e19b34a248`

39-study summary: `diagnostics/20261124_v6_599_integrated_review/39_study_integrated_coverage_UNSIGNED.tsv`
SHA256: `1c38798e7528875d7e17660807c4c407014fc00f1abd72cac3618ee38b198e99`

Local 599-record matrix: `diagnostics/20261124_v6_599_integrated_review/599_accession_integrated_evidence_LOCAL_ONLY.tsv`
SHA256: `c77599f38ecfcf16be5f8a4648f0386421f65de215a9a3a3c13bcb965f82291a`

## Reviewer-critical limitations

An available SRA field is not a certified clinical diagnosis.
A nonempty V5 patient key is not proof of patient independence.
A nonempty treatment field is not proof of treatment-naive status.

The integrated matrix uses NO for certification fields because
certification has not occurred. This must not be interpreted as
evidence of biological absence or a negative clinical finding.

Historical exclusions, metastatic-context flags, technical holds
and frozen cohort membership are unchanged.

Next step: examine patient-key and treatment-value semantics
against their historical evidence grades before clinical adjudication.

No source files modified. No clinical admissions.
No SpatialEcoTyper modeling or projection performed.
