# V6 — 599-record V5 patient and treatment semantics

Date: 2026-11-24
Status: EVIDENCE REVIEW COMPLETE; CLINICALLY UNSIGNED.

## Patient identity evidence

- 467 curated study patient keys, unverified.
- 90 BioSample proxies, not certified person identifiers.
- 42 Stage41 patient proposals, unverified.
- 591 distinct study-scoped keys, NOT 591 verified patients.
- Eight repeated keys covering 16 experiments.
- Repeated groups: SRP058107 (4), SRP089805 (3),
  SRP131744 (1).

## Treatment evidence

- 552 PRETREATMENT_NOT_ESTABLISHED.
- 24 CURATED_PRETREATMENT_LABEL_REQUIRES_SOURCE_CONFIRMATION.
- 23 GEO_BASELINE_TEXT_SOURCE_SUPPORTED.
- Only two literal V5 treatment-timing values across 599 records,
  with frequencies 575 and 24. These are NOT equivalent to the
  three evidence-grade categories.

## Interpretation

Populated patient keys and treatment-timing fields do not demonstrate
independent patients or treatment-naive specimens.

The eight repeated-key groups require evaluation of whether they
represent repeated biological specimens, technical replicates or
different regions. Keys from distinct studies are not certified as
globally independent.

The 23 GEO-supported baseline records are not automatically
treatment-naive. Baseline must be defined relative to treatment.

All 599 remain in their original unsigned V6 classification.

## Reproducibility

Script: `tools/v6_599_v5_field_semantics_READONLY_v1.py`
Script SHA256: `326ebf7b5586e9f55cee124e148479490b64141823e6f6e403451af4e69fed5c`

Cluster stdout: `diagnostics/20261124_v6_599_v5_field_semantics/stdout.txt`
Stdout SHA256: `f630dfb328a59eaa555ee0ff7c4bcb1717435f01b54c2dc58543e703ed0fd22d`

Frozen integrated input:
`diagnostics/20261124_v6_599_integrated_review/599_accession_integrated_evidence_LOCAL_ONLY.tsv`

Input SHA256:
`c77599f38ecfcf16be5f8a4648f0386421f65de215a9a3a3c13bcb965f82291a`

No original files changed. No clinical approvals.
