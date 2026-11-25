# DRP005276 — public DDBJ and PLOS S1 crosswalk

Date: 2026-11-25

Disposition:
PUBLIC_DRA_EXACT_103_PASS_S1_ISOLATE_VALUES_103_CLINICAL_UNSIGNED

## Verified original deposit

- The public DRA008030 deposit contains 103 experiments.
- All 103 experiment accessions match the original local SRA study.
- No additional or missing experiments were detected.
- All 103 public experiment-to-sample links agree with local SRA sample accessions.
- All 103 experiments reference DRP005276.
- The selected V6 set contains 83 experiments.
- The other 20 remain excluded by reconciled historical QC failures.

## Published supplementary evidence

Original publication:

Kawahara et al., PLOS ONE (2019).
DOI: 10.1371/journal.pone.0220173

The PLOS S1 supplementary workbook was retrieved and parsed.

- The workbook contains one worksheet.
- All 103 source isolate identifiers match S1 cell values.
- No original SRA sample names match S1 cell values.
- No original SRA sample accessions match S1 cell values.

These matches establish a strong identifier bridge.

However, cell-value membership alone does not demonstrate
103 independent patient records or complete patient metadata.

The patient axis, table geometry, repeated identifiers,
clinical fields, tumour stages and treatment information
must be resolved before individual clinical certification.

## Exact external evidence hashes

DRA experiment XML:
11055236a400ee93f6947086524657b943aa9ab74eedeb5037a7cdc7cbd6c1c6

DRA sample XML:
fe3b0df8421b41066477823b86aa30389eb50b4a4d236f4920903018bb28f736

DRA study XML:
fd06afd393918fb49f4a57fabae2349b92a9f98b8164d8a758c1a6f9548dd545

DRA submission XML:
2126a48215b8a92cc067d77def72c8f143d52b83827aceee001437b2968c7fb0

PLOS S1 workbook:
d43789f6dab336a047dc10690600f0b57c48de61156c9d72b59ae4e9a1508270

## Remaining scientific gates

1. Verify exact isolate-to-patient mapping in S1.
2. Determine NMIBC/MIBC composition of selected 83.
3. Test whether QC exclusion differs by tumour stage.
4. Establish specimen-level tumour origin and pathology.
5. Determine treatment history before tissue collection.
6. Establish independence from ARCHS4 and other studies.

Clinical admissions: ZERO.

Technical reinstatements: ZERO.

Frozen V5/V6 and ARCHS4: UNCHANGED.

SpatialEcoTyper model inference: NOT PERFORMED.

## Reproducibility

The exact executed diagnostic and captured stdout are archived
as script.py and stdout.txt. Their byte-level SHA256 values
are recorded in SHA256SUMS.txt.

Individual sample and patient values were not printed.
