# DRP005276 to frozen ARCHS4 — exact accession screen

Date: 2026-11-26

Status:
EXACT_SRX_BIOSAMPLE_SCREEN_PASS_PATIENT_INDEPENDENCE_UNSIGNED

- Stage04 and Stage05 source/frozen copies are SHA256-identical.
- 9,465/9,465 atlas GSMs have populated GEO relation fields.
- 9,465/9,465 have experiment and BioSample associations.
- Stage04/05 parsed mappings agree, but share GEO provenance.
- DRP005276 source universe: 103 experiments.
- Technical selection: 83 retained, 20 historical QC-held.
- Selected 83: zero exact experiment or BioSample overlaps.
- Held 20: zero exact experiment or BioSample overlaps.
- All 103: zero detected exact accession overlaps.
- Candidate/full-source BioSample fields agree for all 103.
- Zero exact accession reuse is NOT proof of patient independence.
- Biological patient independence: UNSIGNED.
- Pretreatment status: UNSIGNED.
- Clinical admissions: ZERO.
- No SpatialEcoTyper modeling resumed.

Evidence:
- diagnostics/20261126_v6_archs4_parquet_reader_membership_preflight_v1/
- diagnostics/20261126_v6_archs4_relation_coverage_conda_v2/
- diagnostics/20261126_v6_drp005276_atlas_exact_join_v1/

Next:
Review biological specimen/patient identity and original
precollection treatment evidence as separate clinical gates.
