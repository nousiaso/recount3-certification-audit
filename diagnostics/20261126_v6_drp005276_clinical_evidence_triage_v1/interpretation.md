# DRP005276 clinical evidence triage

Date: 2026-11-26

Decision: STAGE41_CLINICAL_EVIDENCE_EXHAUSTED_GATES_UNSIGNED

Stage41 source and candidate identity reconciled:
103 total experiments, 83 retained, 20 QC-held.

All 83 retained have:
- Pretreatment status NOT_ESTABLISHED.
- Pretreatment certification NO.
- Primary tumour certification NO.
- Patient independence certification NO.
- Cohort approval NO.

No retained specimen has populated Stage41 treatment, tissue
or patient-group evidence pointers.

The original Kawahara et al. PLOS ONE 2019 publication, DOI
10.1371/journal.pone.0220173, describes primary bladder tumour
tissue, 103 patients, and prospective multicentre collection.
This supports study-level specimen context, not per-specimen
pretreatment chronology.

The Stage41 study_exposure_detected field is not accepted as
a clinical therapy-history assertion.

DRP005276 has no recovered GEO series in the candidate table.
Thus its zero GSE-series intersection with ARCHS4 is uninformative
for patient independence.

Previously verified zero exact SRX/BioSample overlaps remain
valid technical evidence, but do not establish patient identity.

Clinical admissions: ZERO.
Original 20 technical holds unchanged.
No frozen source/model files modified.

Next: inspect original DDBJ/GEA study/sample metadata availability.
