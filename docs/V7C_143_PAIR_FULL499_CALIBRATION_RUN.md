# V7C: execute 143-pair, full-499 frozen-W assay calibration

**Purpose:** actual matched-sample computation (not another schema-only audit). Compare recount3 G026 gene-body base-pair coverage against the archived ARCHS4 frozen Stage10 outputs for exactly matched SRP074349 specimens. This is **technical assay transfer**, not an external biological validation cohort.

## Observed preflight V7B (cluster stdout, 27 November 2026)

- Stage5 bridge: 165 rows (143 SRP074349 calibration pairs expected).
- Gene map: 63,856 rows; expected recount3 H5 axes 63,856 × 144 (SRP074349) and 63,856 × 360 (FUSCC).
- Frozen W: 25,410 bytes, SHA256 `39a96b598ec74ccbe67bea49254ac94ae568ca22ec162241ab7849976b52cdcf`.
- Stage09 source SHA256 `46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500`.
- Original Stage10 source SHA256 `16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c`.
- Rscript available from the SpatialEcoTyper v1.0.3 conda environment even though absent from base PATH.

## Execution

After `git pull --ff-only origin main`:

```bash
bash scripts/run_v7c_matched_143_frozen_transfer.sh
```

The launcher checks Python AST/importability and R package availability before starting. It executes:

1. `tools/v7c_prepare_matched_143_G026_499.py`: validates exact bridge-to-H5 and ARCHS4 GSM-to-HDF5 indices, reconstructs full 499 frozen model symbols using the **original Stage09** mapping/alias/collapse/TMM/UQ functions, and creates CPM/TMM/UQ normalized recount3 G026 feature inputs. The all-valid-symbol source matrix supplies denominators, never the 499-gene subset.
2. `tools/v7c_project_matched_143_frozenW.R`: applies **existing frozen Stage10 pan-cancer** anchor means/SDs, unchanged W and `SpatialEcoTyper::DeconvoluteSE(scale=FALSE, sum2one=FALSE)`. Reads archived ARCHS4 same-GSM H and scaled matrices (never re-estimates anchors or W) and calculates 143-pair transfer metrics.

**Strict gates:** exactly 143 unique SRP074349 GSM/library/run pairs; direct ARCHS4 accession/index concordance; exact recount3 gene-order and run alignment; full 499 source/model features; authoritative frozen scaling and H; reference GSM membership; package/W identity. Any failure stops **without silently substituting alternative mappings**, and is for review.

**Limitations to interpret:** recount3 gene-body coverage is not original ARCHS4 read count/pseudocount. The original Stage09 TMM/UQ algorithms are reused, but recount3 factors are estimated across this defined 143-pair technical calibration set, not the original 9,465-sample reference. Differences reflect transfer and source-assay/normalization effects; they are not automatically biological.

**Outputs (private):**
`private_data/v7c_srp074349_matched_transfer_143_v1/` contains accession-linked matrices, per-pair metrics and RDS. Not for Git. Summary aggregate data are copied to `diagnostics/20261127_v7c_srp074349_matched_transfer_v1/aggregate_transfer_metrics.tsv`; aggregate stdout/checksums are captured there too.

**Scientific exit:** Review numerical matched-pair concordance and uncertainty *before* defining calibrated eligibility thresholds and any FUSCC projection. No biological validation claims or clinical releases in this phase. Nine known Ensembl-membership differences are retained in the full 499 input and explicitly flagged in local per-gene output; a full prespecified nine-gene sensitivity remains a downstream requirement.

Frozen ARCHS4 source is never overwritten, and the work does not refit the model.
