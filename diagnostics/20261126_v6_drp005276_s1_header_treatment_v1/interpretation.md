# DRP005276 — S1 patient-header and treatment-metadata audit

Date: 2026-11-26

Decision: S1_PATIENT_IDENTIFIER_HEADER_LINK_PASS_TREATMENT_NOT_CERTIFIED

## Completed patient-key evidence

- Original study: 103 SRA experiments and 103 distinct isolate values.
- All 103 isolates occur exactly once in published PLOS S1.
- Isolates occupy 103 separate patient rows in one column.
- The new diagnostic confirms that this isolate column equals the
  worksheet's unique patient-ID-labelled column.
- Thus the original experiment -> isolate -> published patient-record
  crosswalk is supported for 103/103.
- The selected 83 correspond to distinct within-study patient records.

Previously verified stage distribution:

| Cohort | NMIBC | MIBC | Total |
|---|---:|---:|---:|
| Source | 59 | 44 | 103 |
| Selected | 48 | 35 | 83 |
| Historical QC holds | 11 | 9 | 20 |

## Treatment-history evidence

- The inspected worksheet contains 18 populated columns.
- Six columns were recognized by the controlled header dictionary.
- Twelve columns remain unclassified.
- No header matched the defined pretreatment, untimed treatment,
  collection-procedure or outcome/follow-up categories.

This is a negative header screen, NOT proof that the patients were
treatment-naive or that clinically relevant data are absent.

The 12 unclassified headings and publication sampling chronology
must be investigated before making a treatment-timing decision.

## Outstanding eligibility gates

- Individual pretreatment status: NOT CERTIFIED.
- Individual clinical eligibility: NOT FULLY CERTIFIED.
- Cross-study/ARCHS4 patient independence: NOT CERTIFIED.
- Historical 20 technical QC holds: REMAIN EXCLUDED.
- Clinical admissions: ZERO.
- Frozen V5, V6, ARCHS4 and fixed-W model assets: UNCHANGED.
- SpatialEcoTyper biological inference: NOT PERFORMED.

## Reproducibility

The exact executed read-only script and captured stdout are
script.py and stdout.txt.

Their SHA256 hashes are recorded in SHA256SUMS.txt.

The diagnostic verifies two frozen local sources and the published
S1 workbook content hash.

No individual patient or sample identifiers appear in stdout.

Next: inspect the 12 unclassified S1 column headers and examine
the original publication for precollection treatment evidence.
