# V6 — Original SRA clinical attribute coverage

**Status:** Completed metadata audit; clinically unsigned.

## Original evidence

- 599 unique experiments in 39 studies.
- 571 have a recognized specimen/site attribute key.
- 441 have a recognized disease/diagnosis key.
- 24 have a recognized treatment/timing key.
- 75 have a recognized patient-token key.

Categories overlap and must not be added together.

An available field is not evidence that its value establishes
a primary invasive tumour, pretreatment status or patient identity.

These counts concern the 599 historically unattested records only,
not the complete 2,055-library primary screen.

## Source provenance

Script: `tools/v6_source_attribute_values_READONLY_v1.py`
Script SHA256: `51d919fe9d07328829565b1d967a8fe1a637aaaf89d3f1b87abb9c9529bc2254`

Original stdout: `diagnostics/20261123_v6_source_attribute_values/stdout.txt`
Stdout SHA256: `69863153362d13de7ec07e4372d4b0ed98c668f8399a9ebf1188f29d635868fd`

Study coverage: `diagnostics/20261123_v6_source_attribute_values/study_attribute_coverage_39.tsv`
Coverage SHA256: `4bbfcf2769fb5cc4108e04c3c3443f8805b2c1f21c29bf506128257ff5eb220e`

Full value evidence: `diagnostics/20261123_v6_source_attribute_values/attribute_value_counts_599_LOCAL_ONLY.tsv` (cluster only)
Value table SHA256: `400329705cb6fc82c223b7a2384f2c0b24e02650e192b1d5a5988465f35fd452`

## Scientific limitations

- Missing SRA treatment keys do not establish treatment-naive status.
- Patient tokens do not establish cross-study independence.
- Tissue and diagnosis are distinct from sampled primary tumour site.
- The review does not supersede V5, Stage 41 or signed study evidence.
- The original source-evidence table and frozen ARCHS4 are unchanged.

## Next decision

Review existing attribute values, original study documentation,
patient relationships and sampled-site evidence before changing
any clinical classification.

**Clinical admissions: zero.**
