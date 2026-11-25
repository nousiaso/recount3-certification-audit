# V6 SRP218935 — consolidated source-selection and QC reproducibility checkpoint

Date: 2026-11-25. This is a newly executed, checksum-captured integrated read-only diagnostic. It verifies previously reported aggregate observations; it is not a reconstruction of the four former script files.

- ENA PRJNA553096: 95 experiments (76 primary-labelled, 19 normal-labelled); 75 selected experiments have exact ENA experiment/SRA sample matches and ENA numeric material-alias concordance.
- V5 original inventory contains 76 source-primary samples for SRP218935; 75 entered the unsigned administrative primary screen and one remained QC_HOLD_NOT_ADMITTED.
- The held experiment has a recorded historical run hold with top-10 gene coverage 0.689711896, average mapped read length 22.34 bp, and WITHIN_STUDY_OUTLIER. Three recorded QC reasons are for **one** excluded experiment.
- The 75 selected retain historical source-primary support but patient identity and pretreatment are NOT certified. The historical 75 `PRETREATMENT_NOT_ESTABLISHED` and `CURATED_STUDY_PATIENT_KEY_UNVERIFIED` grades remain unchanged.
- Direct clinical accession-to-patient mapping, individual pre-treatment chronology, cross-study/ARCHS4 independence, assay/model-transfer calibration and all-499 frozen-W feature equivalence remain open. No clinical admissions. Frozen cohort and ARCHS4 model unchanged.

## Provenance

Exact executed script and captured stdout are in this directory, with `SHA256SUMS.txt`. Frozen input SHA256 asserted in the script against original V5 inventory, V5 screen, V6 proposed review, 1097-source ledger, 3602-readiness and 409-run hold queue. No original source files were modified. Current ENA source data were retrieved read-only. No patient/sample values were printed. This diagnostic supersedes the earlier error in an overly restrictive sample-title classifier.

Disposition: `SOURCE_SELECTION_AND_HISTORICAL_QC_RECONCILED_CLINICAL_UNSIGNED`. Do not repeat without changed source inputs. Next review: DRP005276.
