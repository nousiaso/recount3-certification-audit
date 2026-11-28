# Current project status
**28 November 2026 | Clinical certification UNSIGNED | Clinical admissions: 0**

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
| V5 source-supported non-tumour flags | 21 | Historical V5 proposed exclusions |
| V6 cumulative proposed non-tumour exclusions | 22 | Earlier precedence audit; 2,033 administrative remainder, NOT approved |
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

Eight exact queue joins: WITHIN_STUDY_OUTLIER, REVIEW_QC, HOLD, P2_STUDY_RELATIVE_OUTLIER, HOLD_REASON_NOT_YET_DECISIVE. Keep the holds; these are not proven permanent failures or acceptable samples. Source-only record remains excluded. All 140 candidates lack populated proposed patient groups and patient/pretreatment evidence pointers. All five certification flags explicitly NO, confirmed in later diagnostic.

**Latest native evidence:** actual values confirm all five clinical certification gates are NO for 140/140; no explicit patient IDs or pretreatment proof. All 140 match native SRA and have unique RNA-seq IDs, seven tissue categories and nine disease categories. SRP057793 links to sarcoma BioProject PRJNA282597 and GEO GSE71119/GSE75885; published subset membership now reconciled below. [Evidence and archival inventory](STATUS_20261127_SRP057793_NATIVE_EVIDENCE_AND_ARCHIVAL_INVENTORY.md).

**Native categories verified:** 140 frozen sarcoma tumours, PolyA paired RNA-seq, HiSeq 2000; nine histologies and seven anatomical labels, unique specimen IDs not patient IDs. [Executed evidence](STATUS_20261127_SRP057793_NATIVE_CATEGORIES_VERIFIED.md).

**Exact GEO crosswalk PASS:** 138/140 retained linked (129 GSE71119, 106 GSE75885; 97 overlap), original hashes verified. Treatment protocol absent; extraction protocol present for all linked. [Executed evidence](STATUS_20261127_SRP057793_GEO_EXACT_CROSSWALK_PASS.md).

**Protocol and reciprocal pairing executed PASS:** 138/140 source-supported frozen PolyA chemistry (remaining 2 native labels only); 90/140 reciprocal frozen microarray links with explicit untreated primary tumour wording. 89/90 metastasis annotations concordant, one exception preserved. No patient independence established. [Executed review](STATUS_20261127_SRP057793_PROTOCOL_PAIRING_EXECUTED.md).

**Cohort-wide exact accession screen executed PASS:** all 2055 candidates / 77 studies, zero detected SRX/ERX/DRX or BioSample overlaps and zero missing or discordant BioSample token sets. 138 exact sarcoma GEO aliases also zero atlas hits. Full 9465 atlas graph verified. Dedicated PyArrow environment retained. [Execution and limits](STATUS_20261128_V8_PRIMARY2055_ATLAS_EXACT_SCREEN_PASS.md). Do not repeat unchanged inputs.

**Cohort-wide patient grouping census executed PASS:** 2055 distinct BioSamples, 1969 provisional study-scoped keys, 54 repeated groups spanning 140 libraries across 13 studies; no verified patient count. All five historical certification flags NO. [Executed evidence](STATUS_20261128_V8_PRIMARY2055_PATIENT_GATE_CENSUS_PASS.md).

**Integrated unsigned ledger executed PASS:** 2055 rows, 1696/599 historical source attachments, 138/90 sarcoma evidence attachments. After 22 proposed exclusions: 2033 records, 1965 provisional keys, 36 repeated groups spanning 104 libraries in 11 studies. Private ledger and study aggregate checksums verified. [Executed evidence](STATUS_20261128_V8_PRIMARY2055_INTEGRATED_LEDGER_PASS.md).

**Repeated-group source review executed PASS:** 22 groups/75 repeated libraries have concordant V5 source partitions across NSCLC, gastric and hepatoblastoma (95 total study records); remaining 14 groups/29 records across eight studies unresolved, including one explicit SRP163173 conflict. [Executed evidence](STATUS_20261128_V8_REPEATED_GROUP_SOURCE_REVIEW_PASS.md).

**Grouping addendum executed PASS:** 124 private evidence records, 95 source-supported and 29 unresolved; all five checksums passed. Preserve 14 unresolved groups and SRP163173 conflict; no automatic collapse. [Executed evidence](STATUS_20261128_V8_GROUP_ADDENDUM_PASS_FUSCC_NEXT.md).

**FUSCC release-gate review executed PASS:** 334 nonheld tumour libraries, all one verified run, zero historical hold hits; 26 held (25 within-study outlier only, one additionally top10 concentration >50%). Retain holds. Native protocol absent, source-level Ribo-Zero evidence; quantitative compatibility and formal clinical gates unresolved. [Executed evidence](STATUS_20261128_V8_FUSCC_RELEASE_GATE_REVIEW_PASS.md).

**FUSCC author annotation join executed PASS:** 448/448 base+material, 360 tumour bases, three unresolved REP suffix differences among retained. All 334 subtype/batch annotations at base level; 291 purity estimates (43 missing all batch3), native selection labels perfectly aliased with batch. [Executed evidence](STATUS_20261128_V8_FUSCC_AUTHOR_ANNOTATION_JOIN_PASS.md).

**FUSCC author annotation addendum executed PASS:** 360 private rows, 334 retained/26 held, three suffix exceptions preserved; all three checksum entries passed. FUSCC annotation complete; clinical Table S1, suffix exceptions and quantitative compatibility remain open. [Executed checkpoint](STATUS_20261128_V8_FUSCC_ANNOTATION_ADDENDUM_PASS_SRP158272_NEXT.md).

**SRP158272 membership review executed PASS:** 139 source, 134 candidates, 129 retained, five held and five source-only. All native Oligo-dT labels; construction protocol absent; unique isolate labels not patient evidence. [Executed evidence](STATUS_20261128_V8_SRP158272_MEMBERSHIP_PASS.md).

**SRP158272 clinical-key screen executed PASS, linkage negative:** 0/139 exact RR matches; hold all 129 retained until documented clinical bridge. [Checkpoint and Porto source review](STATUS_20261128_V8_SRP158272_CLINICAL_KEY_NEGATIVE.md).

**Joint clinical evidence priority executed PASS:** 424 leads (334 FUSCC, 90 sarcoma), 1609 without joint evidence under current historical/additive criteria; not eligible patients. [Executed checkpoint](STATUS_20261128_V8_JOINT_EVIDENCE_PRIORITY_PASS.md).

**Next:** [direct GEO evidence for prostate/cervical studies](../tools/v8_two_study_geo_clinical_evidence_v1.py) (37 retained candidates), where original clinical source evidence is promising. Persist accession-level source evidence privately; no release. FUSCC clinical workbook and sarcoma patient identity remain open.

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

The local inventory identified 24 untracked V7 audit directories and two untracked V6 implementations. Archive verified aggregate tables first, review remaining logs and scripts before staging, without rerunning. Then resume targeted clinical review, consolidate evidence across the 77 screen studies, resolve the broader unresolved source universe, and build one versioned end-to-end release workflow. An administrative candidate is not an admission. No inference on the external cohort before signoff.
