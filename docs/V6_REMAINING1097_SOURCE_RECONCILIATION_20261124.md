# V6 — Original SRA reconciliation of remaining 1097 records

Date: 2026-11-24
Status: SOURCE_JOIN_PASS_NOT_CLINICALLY_CERTIFIED.

## Original historical-review universe

- 1696 experiments in HISTORICAL_EVIDENCE_REVIEW_PENDING.
- 599 historical-primary-unattested records previously reconciled.
- 1097 remaining records across 28 studies.
- All 1097 exact experiment-level original-SRA matches.
- Zero detected specimen-metadata conflicts or source holds.

## Historical evidence

- 603 surgery-context, not pathology-certified.
- 494 primary-tumour study-source-supported.
- 674 at surgery without pretreatment claim.
- 334 cohort-level pretreatment source-supported.
- 75 pretreatment not established.
- 14 GEO baseline text source-supported.

## Original SRA attribute availability

- Specimen: 946/1097.
- Disease: 167/1097.
- Treatment: 6/1097.
- Patient token: 96/1097.
- No entirely unparsable attribute records.

Attribute presence is not clinical proof.
Study-level pretreatment evidence requires separate provenance.
Surgery context is not proof of primary tumour pathology.

## Provenance

Script: `tools/v6_remaining1097_original_sra_audit_v1.py`
Script SHA256: `b9e29036c261ef04ef4b8c87b741abc02ad0e56aa732863edb09696b36468bde`
Stdout SHA256: `3c0aefada22a3b8b711b149d787369069d9a9c971d636e95d4b5c972b477386c`

- `diagnostics/20261124_v6_remaining1097_source_v1/1097_original_source_evidence_LOCAL_ONLY.tsv` SHA256 `194a214e48b0105e9e6a40b1b72a8d6c5c7cf08abf4a728d274f80575f04bf1c`
- `diagnostics/20261124_v6_remaining1097_source_v1/1097_attribute_values_LOCAL_ONLY.tsv` SHA256 `daea3d4755d28983199e9bedcda479363dbfb8789931ca01c33bfc07836802db`
- `diagnostics/20261124_v6_remaining1097_source_v1/28_study_source_summary_UNSIGNED.tsv` SHA256 `2b107551a8269dc53dc5b2c01adace86ca96db9d87c06571e460cffdc373f82b`
- `diagnostics/20261124_v6_remaining1097_source_v1/historical_grade_crosstab_UNSIGNED.tsv` SHA256 `f5dfb15b1a8e5aa5f97854b37abe0f9447da0cf706121024d5a8a6f645cf81f2`
- `diagnostics/20261124_v6_remaining1097_source_v1/summary.json` SHA256 `004d7757432541b5b54edbaea6f4bb9f7c243887aeffc98a79a9fc7808ecf92b`

No clinical approvals.
No changes to frozen V5/V6 records or SpatialEcoTyper models.
