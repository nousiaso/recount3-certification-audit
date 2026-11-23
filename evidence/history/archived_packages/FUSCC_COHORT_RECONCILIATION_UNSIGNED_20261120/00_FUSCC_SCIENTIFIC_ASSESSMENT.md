# FUSCC TNBC / SRP157974 — accession-level eligibility reconciliation (UNSIGNED)

**Date:** 20 November 2026. This is a new, standalone analysis of already uploaded Stage 40/41 evidence and user-provided live atlas-join results. It does not change earlier frozen files or certify clinical eligibility.

## Conclusions supported by the supplied data

1. Stage 41 contains **360** FUSCC candidate RNA-seq libraries with **360 distinct SRX, SRS, BioSample and `FUSCCTNBC` RNA specimen tokens**; the full-source Stage 41 set contains **448 experiments**, **88** not selected among candidate tumour libraries. Published FUSCC data descriptor describes **360 tumour RNA-seq samples from 360 patients plus 88 adjacent normal RNA-seq samples**. The count match is strong study-membership corroboration; inspect the 88 noncandidate source names before independently assigning every one to normal.
2. Of the 360 candidates, **334** lack a historical run hold and enter the focused primary screen; **26** are held. Critically, the uploaded Stage 41/40 extracts say only `HISTORICAL_RUN_HOLD`, without exact underlying failure reasons. Do **not** equate the 334 nonheld records to passing expression QC.
3. Among focused 334, **279** have independently submitted OncoScan GEO patient tokens matching RNA specimen labels; **55** lack this cross-assay corroboration but carry distinct SRA RNA specimen labels. Because the paper describes one tumour RNA-seq profile per patient, it is reasonable to propose all 334 as distinct study-scoped patient units **for conditional review**. Nonavailability of OncoScan evidence does not itself establish duplication. A 279-only core would select on availability of another assay and could introduce selection bias.
4. The independent atlas-side read-only comparison reported **9,465/9,465** SRX and BioSample mappings recovered, **850** previously known exact candidate overlaps, **zero new**; **zero for all 360 FUSCC candidates** (and the 334 focused). This establishes recovered accession-disjoint and study-alias-held-out status, *not a guarantee of unrelated patients*.
5. The original SRA `library_selection` field shows 318 `cDNA`, 42 `Inverse rRNA`; curated annotation is `Inverse rRNA` for 360. The **published TruSeq Stranded Total RNA Ribo-Zero Gold protocol** and the inspected NCBI SRA experiment protocol support cohort-level rRNA depletion; metadata label discrepancies are not by themselves 318 failures. Assay source-level support is distinct from expression QC in recount3.
6. Published WES sample QC identifies **three** problematic WES-labelled samples, corresponding to FUSCCTNBC044_WBC, FUSCCTNBC140_WBC, and FUSCCTNBC030_TT. FUSCCTNBC032 was reported as the alternative matched tumour for the FUSCCTNBC030 WES pair, **not as a fourth proven mislabeled WES sample**. Stage 41 flags four associated RNA profiles; none establishes an RNA identity swap. Keep these RNA records flagged for WES-linked analyses, but do not exclude them automatically from RNA-only analysis.

## Recommended analysis definition (not yet released)

- Proposed **primary RNA-only FUSCC biological cohort:** the 334 nonheld primary-tumour RNA profiles, provided quantitative expression-compatibility/QC and original hold-rule checks succeed. Study-namespaced RNA specimen token is the provisional patient unit. This avoids restricting to the OncoScan-overlap subset.
- **Patient-token sensitivity:** the 279 GEO cross-assay corroborated samples; test for material differences versus the full 334.
- **WES-caution sensitivity:** omit 044, 140, 030, 032 from RNA-only sensitivity analysis if a phenotype/result depends on cross-assay identity (334 → 330; 279 → 275). Do not treat 032 as independently mislabeled.
- **Hold** the other 26 pending documented historical-run hold reasons, including two names with `.rep` and three with `.2018` suffixes. They do not share a `FUSCCTNBC` token with any focused-334 library in the Stage 40 data.
- **No approval:** no sample is signed as a final strict study-held-out participant by this analysis.

## Important remaining technical evidence

The Stage 41 summary does not expose why the 26 historical holds were created, nor the recount3 matrix-level gene and depth QC of each focused library. Execute the supplied small read-only cluster diagnostic, then apply the prespecified QC thresholds from the locked pipeline. Do not retrain the SpatialEcoTyper model or create a new collection of unrelated stages.

## Sources

- Jiang et al., *Scientific Data* (2022), https://www.nature.com/articles/s41597-022-01681-z (360 RNAseq patients, 360 tumour and 88 normal samples, TruSeq/Ribo-Zero, tissue and WES QC).
- FUSCC cohort treatment history stated in https://pmc.ncbi.nlm.nih.gov/articles/PMC11082668/ (previously untreated primary cancers).
- Original public SRA experiment protocol, https://www.ncbi.nlm.nih.gov/sra/SRX5322041.
- Original Stage 41 v1.0.1 uploaded TSVs plus Stage 40 candidate source names; frozen 9,465-atlas SRX/BioSample live read-only join output supplied by the user.

**Evidence semantics:** accession matches and cross-assay tokens do not constitute genetic identity confirmation. The published WES QC findings must not be interpreted as proof of a problem in the RNA-seq samples.
