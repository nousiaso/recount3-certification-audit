# V6 — 599 historical-primary-unattested records

Status: SOURCE ASSEMBLY COMPLETE; CLINICAL CERTIFICATION UNSIGNED.

## Verified observations

- 599 unique experiments across 39 studies.
- 599/599 exact original-SRA experiment matches.
- Zero unmatched experiments.
- 59 experiments with multiple original metadata rows.
- 147 additional run-related metadata rows.
- Zero disagreements in the four inspected specimen fields:
  sample_acc, biosample, sample_title, sample_attributes.

The complete source observations were assembled into one experiment-level
evidence record without treating repeated metadata rows as patients.

## Provenance

Script: `tools/v6_primary_unattested_source_review_v1_1.py`
Script SHA256: `c3a1dce451f94032837d320c99d371031e24fcecd7123967023dfac3b2cf7d1a`

Cluster stdout: `diagnostics/20261123_v6_primary_unattested_review_v1_1/stdout.txt`
Stdout SHA256: `453ffe18e0efff02d4e43cd1b3e110532d472acc23757ce800f5c9b41571095b`

Detailed accession-level output (cluster only):
`diagnostics/20261123_v6_primary_unattested_review_v1_1/source_evidence_599_LOCAL_ONLY.tsv`
SHA256: `8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb`

39-study summary:
`diagnostics/20261123_v6_primary_unattested_review_v1_1/study_source_review_summary.tsv`
SHA256: `b8d815fb4b508ebaefa1bec7d9ffefb6885e073ca0a5ef1cd20b5c80582c99f8`

## Interpretation and limitations

The earlier v1 script failed because it assumed one raw metadata row
per experiment. That failure is preserved as historical evidence.

V1.1 accepts repeated experiment rows only when the four selected
specimen fields agree. Run-specific original records were not altered.

Exact source matching establishes provenance, not that the sampled
material is an untreated invasive primary tumour.

The existence of an attribute key does not establish the biological
meaning or completeness of its value.

No patient identities or official run accession mappings are certified
by this assembly. No model projections or biological analyses ran.

## Next scientific action

Review original values for sampled tissue/site, disease, tumour
presence, treatment timing and patient grouping within these same
599 records. Do not rerun source accession matching.

Frozen ARCHS4 unchanged. No clinical admissions.
