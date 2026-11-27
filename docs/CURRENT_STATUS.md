# Current project status
**27 November 2026 | Clinical certification UNSIGNED | Clinical admissions: 0**

This is the authoritative human-readable checkpoint for current work. Dated reports retain their historical scope and should not be used as current run instructions.

## Objective

Expand the frozen ARCHS4 SpatialEcoTyper atlas using independent, clinically documented primary tumours from recount3; preserve the original atlas and model, then enable study-aware biological analyses only after a signed cohort lock.

## Cohort accounting

| Universe | Count | Meaning |
|---|---:|---|
| Frozen ARCHS4 | 9,465 samples | Immutable original reference |
| Historical manual-curation universe | 5,241 libraries / 6,337 runs / 176 studies | Historical inventory, not final cohort |
| Stage41 original source identity table | 5,221 experiment records | Distinct source ledger; do not equate to historical library count |
| Stage41 provisional candidate universe | 3,602 libraries / 126 studies | Administrative integrated evidence |
| Primary screen | 2,055 libraries / 77 studies | Unsigned, physical screen unchanged |
| V5 source-supported non-tumour flags | 21 | Proposed exclusion; 2,034 is not a new certified cohort |
| Publication-reviewed primary-screen studies | 29 / 1,755 libraries | Reviewed, not signed |
| Remaining-queue primary-screen studies | 48 / 300 libraries | Prioritize evidence resolution |
| Exact matched technical bridge | 165 specimens / 476 runs / 2 studies | Calibration only, excluded from primary cohort |
| Clinical admissions | 0 | No release |

Sources: [historical counts](STATUS_20261123.md), [clinical triage](STATUS_20261126_V6_77_STUDY_CLINICAL_TRIAGE.md), [priority reconciliation](STATUS_20261126_V6_PRIMARY_STUDY_PRIORITY.md), [V7S closure](STATUS_20261127_V7S_TWO_STUDY_TECHNICAL_CALIBRATION_COMPLETE.md). Counts refer to distinct historical ledgers and units; no new filtering or reconciliation was performed during documentation cleanup.

## Completed and frozen

- Original ARCHS4 Stage09/10 numerical replay: PASS, 499 genes and 500 directional features. Original Stage09–13 outputs remain unchanged.
- Full-W feature investigation: four removed constant directional rows explain the observed isolated-subset discrepancy to numerical accuracy; preserve original implementation separately from independent full500 KL estimator.
- Full available same-specimen technical calibration: 143 lung and 22 prostate; CPM median pooled TVD 0.02312456, dominant agreement 161/165. All component Pearson correlations >0.98 in both studies. Systematic signed differences remain.
- Source-level evidence and historical QC reconciliation are documented, including 850 exact candidate overlaps in the historical screen. No proof of global patient independence from accession nonoverlap.
- [Technical report](TECHNICAL_CALIBRATION.md) contains the current methodological contract. Do not repeat completed calibration.

## Current clinical focus: SRP057793

| Membership | Count |
|---|---:|
| Original source | 149 |
| Integrated candidates | 148 |
| Historical technical holds | 8 |
| Primary-screen candidates | 140 |
| Source-only record | 1 |

Eight exact queue joins: WITHIN_STUDY_OUTLIER, REVIEW_QC, HOLD, P2_STUDY_RELATIVE_OUTLIER, HOLD_REASON_NOT_YET_DECISIVE. Keep the holds; these are not proven permanent failures or acceptable samples. Source-only record remains excluded. All 140 candidates lack populated proposed patient groups and patient/pretreatment evidence pointers. Populated certification flags have not yet been interpreted from their actual values.

**Next diagnostic:** inspect aggregate actual values of clinical status flags and native study source metadata; establish exact specimen and patient identity, pretreatment and assay evidence, and atlas independence. Do not admit or reinstate any record based on cardinalities or outlier ranges. [Latest evidence](STATUS_20261127_SRP057793_HISTORICAL_HOLD_REASONS_RESOLVED.md).

Other documented holds: [DRP005276](STATUS_20261126_DRP005276_CLINICAL_HOLD.md), [SRP185703](STATUS_20261126_SRP185703_CROSSWALK_HOLD.md). SRP060016 specimen classification resolved but treatment and admission unsigned; [evidence-level decisions](../decisions/SPECIMEN_DECISIONS_SRP060016.tsv).

## Release gates (all required)

1. Exact source/run/specimen identity; assay chemistry and compatible expression QC.
2. Explicit primary sampled lesion and material; exclude normal, metastatic, or ambiguous contexts as applicable.
3. Evidence of treatment timing; do not infer from surgery or first-line treatment alone.
4. Patient grouping, repeated regions, technical replicates, and canonical independent unit.
5. Cross-source and frozen-atlas independence beyond simple accession nonoverlap.
6. Frozen model coverage, adequacy and normalization sensitivity for the external candidates.
7. Signed adjudication, immutable cohort manifest, reproducible analysis-ready release.

## Priorities

First archive the remaining verified original aggregate execution records without rerunning. Then resume targeted clinical review, consolidate evidence across the 77 screen studies, resolve the broader unresolved source universe, and build one versioned end-to-end release workflow. An administrative candidate is not an admission. No inference on the external cohort before signoff.
