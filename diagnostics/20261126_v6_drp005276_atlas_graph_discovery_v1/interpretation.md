# DRP005276 — archived atlas-graph discovery

Date: 2026-11-26
Decision: HISTORICAL_NO_EXACT_OVERLAP_SUPPORTED_GRAPH_PROVENANCE_PENDING

## Findings from previously frozen V5

- Original DRP005276: 103 experiments (83 selected; 20 QC-held).
- All 103 V5 atlas evidence grades read
  `NO_EXACT_OVERLAP_DETECTED_NOT_PATIENT_CERTIFIED`.
- Known atlas exact BioSample and experiment match flags are `FALSE`
  for both groups; no atlas BioSample/experiment conflicts are flagged.
- The 83 retained rows carry 83 distinct recorded BioSample identifiers,
  83 distinct SRA sample accessions and 83 distinct `source_biosamples`.
  These are deposit/technical identifiers, not cross-study patient IDs.
- Blank fields such as `atlas_match_by_biosample` and
  `atlas_exact_same_atlas_samples` have NOT been treated as negatives.

## Located historical crosswalk assets (NOT yet inspected)

- `26_recount3_external_primary_v1/ARCHS4_RECOUNT3_EXACT_IDENTITY_JOIN_READONLY_20261120.zip`
- `27_recount3_final_cert_v1_0_0/draft/ARCHS4_overlap_status.tsv`
- `28_recount3_final_cert_v1_0_1/draft/ARCHS4_overlap_status.tsv`
- Frozen ARCHS4 project tree contains other older overlap tables;
  different purposes and dates require provenance review.

A prior historical report describes 9,465/9,465 atlas-side exact
SRX/BioSample identity coverage, but the discovered archive's contents,
actual schema, evidence hashes and joins have NOT yet been inspected
in this diagnostic. Its reported absence of overlaps does not certify
that DRP005276 participants are globally independent of ARCHS4 or of
other external studies. The earlier zero-GSM alias scan provides no
negative atlas evidence when it recovered zero GSM aliases.

Clinical admissions: ZERO. Pretreatment status: UNSIGNED.
Cross-study/atlas patient independence: UNSIGNED.
Technical QC-held rows: remain excluded.
No frozen V5/V6, ARCHS4 or SpatialEcoTyper assets were modified.

Next: read-only inspect discovered identity archive ZIP central directory,
metadata, exact file schemas and checksums before any sample-level join.
The executed script, stdout and SHA256 manifest are archived here.
