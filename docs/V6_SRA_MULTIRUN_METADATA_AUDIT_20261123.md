# V6 original-SRA repeated-row metadata audit

Status: COMPLETE FOR SPECIMEN-FIELD CONCORDANCE; CLINICAL CERTIFICATION OPEN.

Source diagnostic:
`tools/v6_duplicate_sra_rows_READONLY_v1.py`

Original printout:
`diagnostics/20261123_v6_duplicate_sra_rows/stdout.txt`

Printout SHA256:
`715633781421b1055c10ce52dfb249b18ee78a7496c81504d07b941b61006c4e`

## Observed results

599 unique selected experiments across 39 studies were matched exactly.
No experiment was missing.

59 experiments had multiple rows in their original SRA metadata:
- SRP012656: 6 experiments with 2 rows each
- SRP018008: 1 experiment with 2 rows
- SRP040442: 8 experiments with 2 rows each
- SRP090849: 22 experiments with 6 rows each
- SRP094795: 22 experiments with 2 rows each

There were 147 additional rows beyond one row per experiment.
All repeated rows agreed within the four inspected specimen fields:
sample_acc, biosample, sample_title and sample_attributes.

Different full-row signatures correspond to different run-related metadata.
The presence of run_alias values does NOT certify official run accession
identities, technical QC, or aggregation eligibility.

## Scientific interpretation

A specimen-evidence record may represent a unique experiment when its
repeated source rows have concordant specimen metadata. The original
run-specific source records must remain available independently.

This resolves the parsing failure in the earlier v1 source assembly;
it does NOT resolve sampled primary-tumour site, pretreatment, patient
independence, or frozen-model assay transfer.

No clinical approvals. Frozen ARCHS4 unchanged.
