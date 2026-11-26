# V7B full-499 assay-transfer experiment: frozen contract

Status: **PRESPECIFIED TECHNICAL CALIBRATION; NOT EXECUTED** (27 November 2026).

## Scientific question

Does recount3 G026 gene-body **base-pair coverage**, after documented source-to-gene harmonization and normalization, reproduce the same-specimen frozen SpatialEcoTyper component outputs obtained from the ARCHS4 Kallisto-derived data?

This must be established before interpreting recount3 ecotype differences as independent biological replication.

## Prespecified calibration cohort

- SRP074349: 143 documented ARCHS4–recount3 same-specimen bridge pairs.
- 144 available recount3 selected runs; **pair inclusion requires exact verified original accession mapping**, not presumed ordering.
- These 143 are overlapping **technical calibration**, not external independent validation.
- No 334-sample FUSCC TNBC inference before transfer/cancer-specific gate assessment.

## Frozen model contract — do not approximate

1. Exact authoritative Stage09 gene symbol collapse (uppercase symbols, merge duplicated symbol rows) and seven documented aliases.
2. Full **499 gene** signed input interface, not the exploratory 465-gene subset; preserve all nine Ensembl-row membership differences (FAM20C, FCAR, ISLR, MAL2, MATR3, PODNL1, PSME2, SLC39A4, TNNT3), and evaluate a sensitivity specification without altering the W.
3. Compute CPM, TMM-CPM and UQ-CPM across **all valid collapsed genes** for each assay, not only the 499 W input genes. Preserve formal UQ evaluability and fallback flags.
4. Transform `log2(normalized_abundance + 1)`, standardize with the unchanged frozen pan-cancer means and SDs, apply original signed directional `NMF::posneg` convention, and call original `SpatialEcoTyper::DeconvoluteSE` with `scale=FALSE`, fixed W, `sum2one=FALSE`. Do not refit W, anchors or switch to NNLS.
5. Report raw H for all ten components (NonSE, SE1-SE9), normalized fractions, dominant ecotype, reconstruction adequacy `V - W %*% H` and within-assay normalization sensitivity. NonSE is a model component, **not** the reconstruction residual.

## Primary matched-pair readouts

- n exact verified pairs; no accession, axis, branch or feature-order leakage.
- Gene-level input agreement and systematic bias across all 499 genes; explicit nine-gene annotation sensitivity.
- Per-pair `H` component differences (including absolute/relative bias), fraction differences, dominant-component concordance, total variation distance, and distribution of score shifts, including 90th/95th percentiles.
- Calibration stratified by normalization branch, UQ evaluability, study outliers and sequencing/coverage summaries, with uncertainty intervals.
- Actual directional reconstruction-fit differences and extremes in frozen Z values.
- Document a numerical transferability decision, with defensible thresholds informed by the frozen reference/replicate behavior. Do not reuse historical within-assay TMM/UQ stability cutoffs as unverified cross-assay equivalence thresholds.
- An uncalibrated gene-body coverage vector must **never** be relabeled as raw RNA-seq read counts.

## Sequencing of analyses

**First** run the V7B asset preflight (presence/schema only). **Next** implement and run the complete paired full-499 frozen transfer experiment using real resolved asset paths and archived original Stage09/10 methods. **Then**, if technically supported, execute clearly labeled FUSCC projection smoke tests and finalize patient/specimen-independent clinically suitable validation subset before inferential analysis.

All sample-linked matrices and accession mappings are local-only. Archive code, aggregate stdout, input/output hashes, numerical calibration summaries and explicit decisions, not personal/source linked details.

Source methodological record: `docs/README_MANUAL_CURATION_UPDATED_20261120.md`; original audit `tools/v6_fuscc_448_source_bridge_evidence_20261124.sh`.
