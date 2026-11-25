# DRP005276 — V6 original SRA source-identity review

Date: 2026-11-25

Status: EXACT_SOURCE_IDENTITY_PASS_SPECIMEN_PATIENT_TREATMENT_UNSIGNED

## Evidence

- Original local SRA metadata: 103 distinct experiments.
- Frozen V6 selected subset: 83 experiments.
- Exact original SRA experiment/source joins: 83/83.
- Each selected experiment has one original metadata row.
- SRA sample names: 83 populated, 83 distinct.
- SRA sample accessions: 83 populated, 83 distinct.
- Source isolate values: 83 populated, 83 distinct.
- These three identifier types show no pairwise literal matches.
- Curated V5 patient keys: 83 populated, 83 distinct, NOT verified.
- Explicit source tissue 'bladder': 83/83.
- Explicit source disease 'bladder cancer': 83/83.
- No accession-level primary, metastatic, recurrent or
  pretreatment designation established from inspected fields.
- All 83 frozen source records are PRIMARY_TUMOUR_NOT_ATTESTED,
  CURATED_STUDY_PATIENT_KEY_UNVERIFIED and
  PRETREATMENT_NOT_ESTABLISHED.
- All 83 clinical approvals remain NO_UNSIGNED.

## Interpretation

Exact source identity is resolved for the selected 83.
The 20 other original study experiments have not yet been
classified against the full original source-study design.

Unique isolate, sample and curated-patient tokens do not prove
independent biological specimens or independent individuals.

The common curated treatment-timing token does not establish
presampling treatment status.

No primary-tumour-site, patient-independence or treatment
certification follows from the SRA attributes alone.

## Provenance

The exact executed diagnostic script and captured cluster stdout
are archived in this directory. Their SHA256 values are recorded
in SHA256SUMS.txt.

The script checked the frozen 599-record source evidence,
the frozen clinical review queue and the original compressed
DRP005276 SRA metadata.

Original SRA metadata SHA256:
cba60a28b225d7852e17e180ddaad2680504e7f1cd5cf94de5d73401db12d36c

Original source files and frozen clinical decisions were not
modified. No patient or individual sample identifiers were
printed in the archived stdout.

## Open gates

- Classify the 20 nonselected original study experiments.
- Verify original study design/publication and sample crosswalk.
- Establish primary sampled tumour site and pathology.
- Establish treatment at the time of tissue collection.
- Verify within-study patient grouping and replicates.
- Verify cross-study and ARCHS4 atlas independence.
- Complete the frozen 499-gene model-transfer validation.

Clinical admissions: ZERO.
SpatialEcoTyper execution: NOT AUTHORIZED.
