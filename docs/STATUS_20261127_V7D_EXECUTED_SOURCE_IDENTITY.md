# V7D completed — 143-specimen bridge is disjoint from Stage2 selected extraction

Status: **SOURCE_IDENTITY_PASS; MATCHED_EXPRESSION_NOT_AVAILABLE_IN SRP074349 SELECTED H5**
Date: 27 November 2026.

Grounded in user-pasted V7D cluster stdout; the numeric results are transcribed from that executed read-only diagnostic.

## Exact identity results

- 143 distinct SRX bridge experiments; each has an identified BioSample; 454 **distinct SRR runs**.
- Multiplicity: 84 experiments have 1 run, 4 have 2, 22 have 3, 1 has 8, and 32 have 9.
- All 143 bridge experiments are classified `CALIBRATION_ONLY_NOT_IN_PRIMARY`.
- The existing `stage2_output/studies/SRP074349/selected_run_coverage_G026.h5` has **144 SRR runs from 56 distinct SRX libraries**.
- Exact run overlap between bridge and that HDF5: **0/454**.
- Exact library overlap: **0/143**; exact BioSample overlap in Stage2 QC: 0.
- Historical Stage2 run QC covers 144 runs, 56 libraries, also 0 bridge run/library matches.
- Stage5 decisions cover the same 56 distinct libraries; Stage41 marks 50 `DIRECT_ATLAS_ACCESSION_EXCLUDE_STRICT` and 6 `STUDY_EXPOSED_DIFFERENT_ALIAS_PENDING_GATES`; 54 nonheld and 2 technical holds.
- Original native SRA metadata includes **454 records referencing all 143 bridge experiments**.
- Frozen Stage5 mapping agrees in order for all 63,856 source gene rows, with 499 model genes resolved. All 9 frozen Stage10 scaling/H assets checked by the prior diagnostic are present.

## Interpretation and next decision

This is **not** a case of wrong HDF5 column order or an unusual run delimiter. The existing SRP074349 selected expression dataset and the 143 bridge experiments are **disjoint**. The presence of a study-level G026 HDF5 cannot be used as evidence that the 454 bridge SRRs have been extracted.

The source linkage for the matched cohort exists in native SRA, but **actual recount3 gene-body coverage for those 454 SRRs has not been identified**. Therefore the earlier V7C calibration was not executed and the previous one-run assumption was incorrect.

**Next:** Run the one-time V7E read-only census `scripts/run_v7e_matched_calibration_coverage_inventory.sh`. It scans existing Stage2 selected HDF5 run-ID axes across all studies, checks the Stage2 QC run inventory, and inventories nearby extraction filenames without reading expression values or downloading data.

If none of the 454 runs exists in the already-extracted files, decide on an **explicit separate, versioned recount3 retrieval/extraction** for the 143 paired specimens, using all 454 technical runs and preserving run-to-experiment grouping. Do not force these into the independent primary-tumour cohort, which is a separate scientific population. Verify source semantics/available public recount3 files and extraction method before implementing.

Frozen ARCHS4/Stage5/V5/V6 and the failed V7C outputs are unchanged. The 2,055 clinical candidate baseline still has zero final admissions. No new ecotype projection or biological inference was made.
