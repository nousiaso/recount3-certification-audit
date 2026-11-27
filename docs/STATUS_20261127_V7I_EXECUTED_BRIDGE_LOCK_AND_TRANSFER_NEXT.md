# V7I executed — exact bridge identity, original model-gene mapping and source SHA256 PASS

**Date:** 27 November 2026  
**Evidence:** user-executed V7I stdout from supekgate, `diagnostics/20261127_v7i_post_v7h2_transfer_contract_v1/01_readonly_preflight_stdout.txt` (cluster-local; output may not yet be committed to GitHub)  
**Scientific status:** frozen ARCHS4 Stage09/Stage10 numerical replay PASS; external recount3 transfer PENDING; clinical admissions **0**.

## V7I exact input fingerprints

| Input | SHA256 (observed on cluster) |
|---|---|
| Stage5 165-row bridge, 143 SRP074349 matches | `c7d6d7c754e3ca77d62424f26627a40ceaec10fac98a0741c945773f3d3b1b67` |
| Stage5 63,856-row G026 gene mapping | `68be266dd6bd30bfa79a2bd8ac79f0e49c6f5d0f19d2d083e2ca57850ef17391` |
| Full original SRP074349 G026 gene-sums gzip | `b056331eb16d6ec612409fd37a266abfe10c4c15329b486dbb9df213612797c5` |
| Original Stage09 v1.0.1 Python code | `46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500` |
| V7G v1 Python preparer (cluster checkout) | `0740fda79e09b2e1e688abaaafda8d984ffaee120b58146c393150be8b4fae30` |

The **previously unverified hardcoded V7G v1 bridge and gene-map hashes now match**, both `True`. The uncertainty concerning those hashes is closed.

## Source identity and mapping

- `EXACT_CALIBRATION_IDENTITY_PASS SRX=143 GSM=143 SRR=454 BIOSAMPLE=143`
- All 143 are **CALIBRATION_ONLY_NOT_IN_PRIMARY** same-specimen ARCHS4–recount3 bridge records.
- `G026_SOURCE_GENE_MAP_63856_PASS`
- `RECOUNT3_W_GENE_MAPPING_PASS RESOLVED 499 ALIASES 0`
- `ARCHS4_GSM_AND_GENE_MAPPING_PASS GSM 143 MODEL_GENES 499 ALIASES 7`
- **All 143 matched GSMs occur in the frozen 9,465-sample ARCHS4 atlas**; 0 are outside. This is the decisive new evidence for referencing original archived ARCHS4 normalization/scaled Z/raw H directly rather than re-estimating TMM/UQ on a 143-GSM subset.
- Frozen Stage10 W-gene harmonization ledger: **499 rows per branch (3 branches)**.
- Executed V7H2 replay log exists and carries completion and all three-branch markers.
- No expression values loaded by V7I, no extraction, no model run and no modification of frozen files.

## Interpretation: revise V7G before any transfer

**Methodological decision:** the ARCHS4 arm of the paired calibration should be the *original frozen* Stage09/Stage10 result for those exact 143 GSMs: all three archived 499-gene normalized expression matrices, fixed-anchor scaled matrices and raw H/fraction outputs, with their historical 9,465-sample normalization context. This avoids a spurious difference caused by the v1 V7G implementation recomputing ARCHS4 TMM/UQ factors on the 143-sample subset.

For the recount3 arm, retain all 454 runs grouped by exact 143 SRX experiment IDs, correctly sum source **gene-body base-pair coverage** before symbol-level normalization, and use the frozen 499-gene W plus frozen pan-cancer anchor parameters and original DeconvoluteSE call. Normalization of an external assay necessarily requires a separately declared cohort/reference strategy. Same formula and W do **not** make recount3 TMM/UQ size factors the same as ARCHS4's historical factors. Define and label any cohort-normalized TMM/UQ branch as **assay-specific sensitivity**, with CPM as the primary transparent matched comparison.

Before running, reconcile W-gene harmonization/alias resolution, all 500 directional-feature variance and missing-feature checks, original Stage10 adequacy and original Stage11 evaluation metrics. Include an independent archived-ARCHS4 baseline equality gate for the 143 GSMs across all three branches, with explicit numerical tolerances learned from V7H2. Keep alternate/ablation analyses distinct from the original 499-gene primary model.

**V7G v1 is held/superseded, not approved or executed** as an assay calibration. Do not reuse its output directory to hide or overwrite past evidence. A new versioned implementation and output directory are required.

## Original replay provenance and limitations

The separate original full-9,465-sample V7H2 replay conclusively reproduced archived Stage09 normalized expression from frozen pseudocounts and frozen per-sample size factors, fixed-anchor Z and DeconvoluteSE H with maximum raw H differences of 1.67e-15 (CPM), 1.55e-15 (TMM) and 4.44e-15 (UQ), and **100%** dominant-component concordance in all three branches. See `docs/STATUS_20261127_V7H2_ARCHS4_EXACT_NUMERICAL_REPLAY_PASS.md` for exact numbers and the key limitation that original full-gene TMM/UQ reference-factor estimation was not rerun.

**Remaining scientific gates:** exact external assay preprocessing and transfer, 499-gene model adequacy in recount3, normalization sensitivity, context-specific Stage12 projection when appropriate, study/patient/specimen clinical certification and true atlas independence. None is certified by this V7I preflight.

Do not push patient- or run-linked matrices or manifests to GitHub; summarize aggregate statuses only.
