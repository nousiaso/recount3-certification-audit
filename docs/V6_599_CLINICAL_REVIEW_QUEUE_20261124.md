# V6 — 599-experiment clinical adjudication queue

Date: 2026-11-24
Status: REVIEW QUEUE ASSEMBLED, NO ADJUDICATIONS.

## Completed evidence

- 599 unique experiments from 39 studies.
- 591 distinct study-scoped V5 patient keys, NOT certified patients.
- Eight reused keys covering 16 experiments.
- 552 pretreatment-not-established evidence grades.
- 24 curated pretreatment labels requiring confirmation.
- 23 GEO baseline-text-supported evidence grades.
- Two literal V5 treatment timing categories.
- The V5 timing token shared by the GEO-baseline group also occurs
  in 552 pretreatment-unestablished records.

## Review requirements

Each accession requires evidence for its actual sampled site,
primary tumour identity, collection-versus-treatment chronology,
patient grouping, technical release, and cross-study independence.

Patient and specimen labels cannot substitute for validated identities.
The queue must not be used as an admission list.

## Original GEO references requiring accession-level follow-up

SRP112751 / GSE101588:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101588

The RNA-seq subseries contains 17 samples; the published design
describes predominantly primary untreated colorectal specimens.
Exact patient, specimen, assay and treatment evidence remain to be
reconciled with these 17 selected experiments.

SRP194626 / GSE130469:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130469

The study uses sequential patient-matched breast cancer biopsies
across different expression-profiling platforms.
The four selected experiments need precise patient/assay grouping.

## Reproducibility

Script: `tools/v6_599_clinical_adjudication_queue_v1.py`
Script SHA256: `68efcf0126404ecf617a449c420cbaaf35ee138c8215c1bff3cde4b418074944`

Original stdout: `diagnostics/20261124_v6_599_adjudication_queue/stdout.txt`
Stdout SHA256: `63a2363470cf99920b89d9aa7c2663ac1038865a2826187d4d6ad4f50c3cbcfa`

Aggregate JSON: `diagnostics/20261124_v6_599_adjudication_queue/summary.json`
SHA256: `b24534da4f8a1cb25a23dd918bad4d62e286b76b3d27eb06bc23fcc4ac643fd2`

Study summary: `diagnostics/20261124_v6_599_adjudication_queue/39_study_gate_summary_UNSIGNED.tsv`
SHA256: `ba92beb8419e8f0388ce1539ecdf8c2323a95a781ae72c5ddd8f2ea0b48102d0`

Local-only accession queue: `diagnostics/20261124_v6_599_adjudication_queue/599_clinical_review_queue_LOCAL_ONLY.tsv`
SHA256: `b5fdf922e4791b4b30c2ffa798c6ffed99b43478a0681c02a1faa2aaa248f3bd`

No clinical approvals, changes to frozen inputs, or model execution.
