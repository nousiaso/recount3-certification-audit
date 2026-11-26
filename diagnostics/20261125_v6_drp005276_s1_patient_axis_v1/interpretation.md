# DRP005276 — PLOS S1 patient rows and stage mapping

Date: 2026-11-26

Decision:
EXACT_S1_ISOLATE_POSITION_103_OF_103_STAGE_MAPPING_PASS_CLINICAL_UNSIGNED

## Original evidence

- Original DRA/SRA study: 103 experiments.
- V5 selected: 83; historical QC-held: 20.
- PLOS S1 workbook SHA256:
  d43789f6dab336a047dc10690600f0b57c48de61156c9d72b59ae4e9a1508270
- All 103 original isolate identifiers occur exactly once in S1.
- They occupy 103 distinct worksheet rows in one column.
- A patient-ID-labelled metadata field was detected at column 1.
- Stage was detected at column 3.
- Age, sex, grade and smoking fields were also detected.

The diagnostic stdout does not explicitly print the isolate-column
index. This must be checked against the patient-ID-labelled column
before treating their equivalence as fully verified.

## Tumour stage

| Cohort | NMIBC | MIBC | Total |
|---|---:|---:|---:|
| Complete source cohort | 59 | 44 | 103 |
| Selected 83 | 48 | 35 | 83 |
| Historical QC-held 20 | 11 | 9 | 20 |

The 59/44 full-cohort distribution matches Kawahara et al.,
PLOS ONE (2019), DOI: 10.1371/journal.pone.0220173.

The diagnostic uses explicit stage-code classification rules.
Original field and stage-label semantics require final review.

## Interpretation

The experiment-to-published-patient-row crosswalk is strongly
supported by 103 unique isolate-to-row relationships.

The 83 selected experiments occupy 83 distinct published isolate
rows, but final within-study patient-ID certification requires
direct verification of the isolate column and its header.

QC retention is approximately 81.4% for NMIBC and 79.5% for
MIBC, without a large descriptive stage imbalance.

## Outstanding gates

- Confirm isolate-column versus patient-ID column semantics.
- Establish treatment chronology before tissue collection.
- Assess specimen-level clinical eligibility.
- Confirm independence from other studies and ARCHS4.
- Complete downstream technical/model-transfer validation.

Clinical admissions: ZERO.
Technical reinstatements: ZERO.
Frozen V5, V6, ARCHS4 and model assets: UNCHANGED.

## Reproducibility

The executed script, captured stdout and SHA256SUMS.txt are
preserved in this diagnostic directory.

No individual patient or sample values are included in stdout.

Next: direct patient-header verification and treatment evidence.
