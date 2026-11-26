# DRP005276 — exhaustive PLOS S1 header inventory

Date: 2026-11-26
Decision: ALL_18_PUBLISHED_FIELDS_EXPLAINED_TREATMENT_UNSIGNED

## Verified workbook structure

- Original 103 RNA-seq experiments matched 103 unique S1 patient-ID rows.
- Selected after technical QC: 83 (48 NMIBC, 35 MIBC).
- Historically QC-held: 20 (11 NMIBC, 9 MIBC); none reinstated.
- Published S1 contains exactly 18 populated metadata columns.
- Column 1: patient ID.
- Column 2: TP53 codon 72 germline genotype.
- Columns 3–7: clinical stage, age, sex, grade, smoking status.
- Columns 8–18: somatic mutation status for AKT1, ATM, BRAF,
  CTNNB1, FGFR3, HRAS, KRAS, NRAS, PIK3CA, RB1 and TP53.
- All 18 fields are populated for all 103 patient rows according to
  cell-presence checks; the coding and missing-value semantics were
  not independently validated by this diagnostic.

The 12 previously unclassified fields are the germline TP53 genotype
and the 11 somatic mutation columns, not treatment-history fields.
The column-8 heading `somatic mutation status` groups the gene
columns that follow; the automated semantic labels are descriptive.

## Scientific interpretation and remaining gates

- No S1 column records precollection therapy or collection timing.
- Absence of a treatment column does NOT imply treatment-naive tissue.
- Patient-to-experiment identity and NMIBC/MIBC stage are supported.
- Individual pretreatment status and cross-study/ARCHS4 independence
  remain UNSIGNED.
- Somatic mutation/genotype annotations may support later stratified
  biological analyses only after their value coding is independently
  checked and clinical eligibility is established.
- Clinical admissions: ZERO. Technical reinstatements: ZERO.
- Frozen V5/V6, ARCHS4, W-model and SpatialEcoTyper analyses unchanged.

## Provenance

The executed diagnostic source, captured stdout and SHA256SUMS are
preserved beside this interpretation. Frozen local sources and the
public S1 workbook were hash-checked. No individual patient or
sample values are included in the public archive.

Next: obtain source-backed treatment-before-collection evidence;
separately check cross-study identity/overlap. Do not repeat these
finished S1 header or stage audits on unchanged sources.
