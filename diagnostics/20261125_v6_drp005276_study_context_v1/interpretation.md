# DRP005276 — verified original study identity and cohort-level pathology evidence

Date: 2026-11-25
Disposition: STUDY_PUBLICATION_IDENTITY_MATCH_PRIMARY_BLADDER_COHORT_SUPPORTED_INDIVIDUAL_GATES_UNSIGNED

## Exact local provenance

The cluster's original compressed SRA metadata has one identical populated study title and study abstract across 103 experiments. The original V6 clinical queue supplies 83 selected experiments. Frozen input hashes passed. The executed diagnostic and its captured stdout are in this directory; their SHA256 values are in SHA256SUMS.txt.

The study title exactly matches: Kawahara T et al. (2019), TP53 codon 72 polymorphism is associated with FGFR3 and RAS mutation in non-muscle-invasive bladder cancer, PLOS ONE 14(8): e0220173. DOI: 10.1371/journal.pone.0220173.
Publication: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0220173

## External publication: study-level evidence, NOT accession-level verification

The authors report a prospective multicenter cohort of 144 patients; the included 103 had primary bladder cancer tissue (59 NMIBC, 44 MIBC). The paper describes pathology review of slides and >=10% malignant tumour cell nuclei, and stored FFPE plus frozen tumour samples. The other 41 patients were excluded for upper-tract origin, non-urothelial histology or unavailability of fresh frozen tissue. RNA was extracted from frozen tissue and mRNA libraries sequenced on Illumina NextSeq 500.

The paper's Data Availability statement names DRA008030 and E-GEAD-306. These publication accessions remain to be reconciled directly with the DRP005276 experiment IDs; a matching study title and matching overall count alone do not establish one-to-one patient assignment. The publication's S1 patient metadata and mutation-status table and S2 expression table are potential sources for the crosswalk.

## Important unresolved gates

- The 103 reported patients versus 103 DRP005276 RNA-seq experiments is suggestive, but experiment-to-patient one-to-one correspondence is not yet proved.
- The 83 selected V6 experiments need their actual patient/stage/specimen crosswalk (59/44 is a full-cohort distribution, NOT the 83-selection distribution).
- The retrieved publication methods do not affirm pretreatment sampling for all individuals. Neoadjuvant, prior BCG and other prior therapies are not ruled out by this evidence.
- Independence from ARCHS4 and cross-study duplicated subjects remain unsigned.
- The 20 original run-QC holds remain excluded; no technical reinstatements.
- V5/V6 historical evidence and clinical approval fields are frozen; this is an external evidence addendum, not an approval or a change to a past grade.

Clinical admissions: ZERO. No model projection or downstream SpatialEcoTyper inference.
