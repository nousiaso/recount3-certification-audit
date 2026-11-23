# ARCHS4 × recount3 — unified manual curation workbook (17 November 2026)

## What you are reviewing

- **5,241 library/experiment records** in one Excel sheet: `CURATE_5241`. These correspond to **6,337 sequencing runs** across **176 studies**; a library is *not* a unique tumour, biopsy or human.
- **61 columns** per library: 3 identity headers, 7 review-editable columns, the other 47 frozen source fields, and 4 supplementary study/source-review fields. The original 50-field source library record is preserved without rewriting.
- `STUDIES_176`: 176 study records with **full original curation notes** (not excerpts), protocol information, and review notes.
- `FIELD_GUIDE`: interpretation and cautions.
- `RUN_LINKS_6337.tsv`: all 6,337 run accessions, joinable through `library_id`. A single library can have several runs.

## Steps

1. Open `RECOUNT3_UNIFIED_MANUAL_CURATION_5241.xlsx` with Excel or LibreOffice. Start at `START_HERE` and then go to `CURATE_5241`.
2. Filter by `cancer`, `study`, `evidence_disposition`, `audit_flags`, `library_selection`, `metastatic`, `treatment_timing`, `atlas_exact_same_atlas_samples`, etc. The first three identity fields remain to the left of review controls.
3. **Edit only yellow columns D–J**. In column D (`REVIEW_DECISION`), use one of:
   - `EXCLUDE`: you believe the record should not enter the primary whole-TME cohort.
   - `KEEP_FOR_VERIFICATION`: you believe it is a reasonable candidate, **not** that it is independently certified.
   - `NEEDS_REVIEW`: uncertain, seek study/source evidence.
   - `SEPARATE_CONTEXT`: valid biological material but requires a different stratum (metastatic, post-treatment, enriched compartment, etc.).
   - Blank: not yet reviewed.
4. In `REVIEW_REASON`, `REVIEW_NOTES` and `EVIDENCE_URL`, explain important or contested decisions. Use optional `SPECIMEN_CLASS_OVERRIDE` and `VERIFIED_PATIENT_ID_OVERRIDE` only when supported by actual source evidence. Do not replace blank biological annotations with assumptions.
5. Filter `STUDIES_176` by `study` to read its complete methodological comment, linked DOI, study-level audit and protocol details. The `STUDY_NOTE_PREVIEW` in the master is only a 440-character excerpt.
6. **Strongly prefer `EXCLUDE` instead of deleting physical rows.** This retains a traceable reason and the frozen source data. If you *do* delete rows, I can reconstruct which `library_id` values disappeared, but not your rationale, and they will be flagged for verification rather than silently accepted as scientifically excluded.
7. Save the edited `.xlsx` and upload it back to our chat. I will compare all identifiers with the 5,241-library baseline, check changed metadata, look for potentially incorrect exclusions/inclusions, apply patient and technical caveats, and issue a reconciled review rather than automatically permitting refitting.

## Critical scientific boundaries

- The original provisional statuses are **4,242 potential core**, **548 explicit core holds**, **451 separate-context review**. **Zero** are approved by the final source/patient/expression model release gate.
- 165 exact ARCHS4 BioSample/experiment links are known. The prior GEO crosswalk covers only **740 of 9,465 ARCHS4** records; *absence* of an exact crosswalk match **does not prove independence**.
- 14 experiment identifiers were missing and run identifiers are being used provisionally as `library_id` keys.
- The study-scoped curated `patient_key` and BioSample fallback must not be interpreted as guaranteed global patient independence.
- Histology such as carcinoma in situ or papilloma, 3SEQ, microdissection, target capture, metastatic lesions, post-treatment specimens and study-specific cell-enriched protocols require appropriately differentiated treatment. `KEEP_FOR_VERIFICATION` is not an automatic waiver.

## Deliverables

- `RECOUNT3_UNIFIED_MANUAL_CURATION_5241.xlsx` — **main editable workbook**, filterable library master + full study notes + guide.
- `LIBRARY_MASTER_EDITABLE_5241.tsv` — same 61 columns in a plain text master; optional for command line review/import.
- `STUDY_REVIEW_EDITABLE_176.tsv` — full length study comments and evidence.
- `RUN_LINKS_6337.tsv` — accession crosswalk for sequencing runs.
- `LIBRARY_SOURCE_LOCKED_5241.tsv.gz` — **unaltered source evidence ledger for independent verification**.
- `BASELINE_COUNTS_AND_HASHES.json` — baseline counts and checksum.
- `verify_edited_workbook.py` — compare your returned workbook (or TSV) to frozen source IDs and original metadata. Produces a change report, *not* scientific approval.

**Return-file policy:** Keep original `library_id` values unchanged, and preferably do not reorder/rename column headers. Sorting rows and filtering are fully supported.



================================================================================
UPDATE — 19 NOVEMBER 2026
RECOUNT3 CERTIFICATION, FROZEN ARCHS4 MODEL AND EXECUTION-CONTRACT AUDIT
================================================================================

PURPOSE
-------

This appendix records the additional technical, biological and methodological
verification completed after the original 17 November 2026 manual-curation
workbook was created.

The original workbook remains the 5,241-library, 176-study review universe.
Later recount3 subsets represent successive selection and QC stages; they
must not be confused with the full manual-curation workbook.

SCIENTIFIC POLICY:

1. The original ARCHS4 SpatialEcoTyper analyses are completed and frozen.
2. They must NOT be rerun, refitted, overwritten or silently modified.
3. The original frozen W model, scaling parameters and score conventions
   are authoritative.
4. Recount3 is a separate candidate external cohort requiring its own
   biological, technical, patient-level and model-transfer certification.
5. No new recount3 biological inference is authorized until the relevant
   release gates pass.


1. RECOUNT3 CANDIDATE UNIVERSE AND TECHNICAL REVIEW
---------------------------------------------------

Original manual-curation universe:
    5,241 libraries
    6,337 sequencing runs
    176 studies

Subsequent expression-extraction selection:
    3,639 libraries
    4,166 runs
    129 studies

All 129 selected studies underwent recount3 expression extraction.

The source expression is recount3 G026 gene-body base-pair coverage,
NOT the original ARCHS4 Kallisto pseudocount representation.

G026 gene reference:
    63,856 gene rows
    63,811 distinct Ensembl identifiers
    45 duplicate identifier records

Thirty-seven short-tag DGE libraries were excluded from the provisional
compatible whole-transcriptome candidate set:

    SRP007947: 17
    SRP003901: 10
    SRP024268: 10

Remaining provisional candidate universe:
    3,602 libraries
    126 studies

Additional assay holds/reviews remain, including four Lexogen QuantSeq
libraries in SRP194626.

Technical audit:
    409 runs held for technical review
    350 candidate libraries affected
    341 affected libraries without a retained run
    9 affected libraries with mixed retained and held runs

The 409 held runs include 14 sparse kidney runs in SRP082503.

These figures describe technical evidence, not a final biologically
approved cohort. Run-level evidence must be reconciled with biological
specimen identity, assay type and patient grouping.

Technical source-evidence and adjudication materials include:

    RECOUNT3_SOURCE_EVIDENCE_3602.tar.gz
    RECOUNT3_CERTIFICATION_ADJUDICATION_20261119.zip
    RECOUNT3_CERTIFICATION_V1_0_1_20261118.zip

Important source-evidence script:

    recount3_accession_source_evidence_v1_0_0.py

The 3,602-library source review retained unresolved specimen,
patient and study-context questions. No library should acquire final
approval solely from passing expression or technical QC.


2. ARCHS4 CROSS-ATLAS OVERLAP AND INDEPENDENCE
-----------------------------------------------

Original frozen ARCHS4 analysis-ready atlas:
    9,465 samples
    18 cancer classes

The frozen atlas sample IDs are GEO GSM accessions.

Read-only overlap investigation:

    probe_W34_and_ARCHS4_atlas_READONLY.py

Associated evidence files:

    01_W34_metadata_candidates.tsv
    02_GSM_direct_overlap_candidates.tsv
    03_readonly_probe_summary.json

Findings among 3,602 provisional recount3 libraries:

    2,161 libraries with a recovered GSM accession
    1,441 without a recovered GSM in the inspected source fields

    850 libraries with an exact GSM match to the frozen
    9,465-sample ARCHS4 atlas

    These 850 exact-GSM overlaps span 43 recount3 studies.

DECISION:

    The 850 exact-GSM overlaps cannot count as independent external
    ARCHS4 validation specimens.

    They may be used in explicitly labeled technical calibration where
    methodologically appropriate.

    The remaining 2,752 libraries are NOT thereby certified independent.

The earlier 165 direct ARCHS4 BioSample/experiment links and the later
850 frozen-atlas GSM overlaps reflect different crosswalk methods and
must not be added together without accession-level deduplication.

Independence certification still requires crosswalking:

    GSM / GEO series
    SRX / SRA experiment
    SRR / sequencing run
    BioSample
    patient / donor identifier
    study and source provenance

Repeated specimens, multiple regions, longitudinal biopsies and
multiple runs must not be treated as independent patients.


3. FROZEN W MODEL — COMPLETE GENE-IDENTITY RECONSTRUCTION
---------------------------------------------------------

Official frozen model asset:

    /maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/
    02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds

The original W interface consists of:

    499 distinct input genes
    500 directional model features
    10 output components:
        NonSE
        SE1-SE9

Three normalization branches:

    CPM
    TMM_CPM
    UQ_CPM

Original ARCHS4 gene reference:

    00_download/human_gene_v2.latest.h5

    67,186 HDF5 gene rows

The original Stage 09 gene processing was reconstructed from the
actual source code.

RULE:

    Normalize source symbols to uppercase.
    Sum expression rows sharing the same normalized symbol.
    Apply the seven documented historical-symbol aliases.
    Do not select an arbitrary canonical Ensembl ID.

Documented aliases:

    CXORF36  -> DIPK2B
    FAM71A   -> GARIN4
    GPR1     -> CMKLR2
    H2AFY    -> MACROH2A1
    MRVI1    -> IRAG1
    RARRES3  -> PLAAT4
    WARS     -> WARS1

The frozen W harmonization record:

    16A_frozen_W_gene_symbol_harmonization.tsv

records branch, W feature position, input model gene, processed input
row position and harmonization.

Its input_row_index refers to the processed model matrix, NOT the
original HDF5 gene row.

FINAL FULL-499 ENSEMBL MEMBERSHIP AUDIT:

    490 / 499 genes:
        Identical Ensembl source-row ID multisets between
        recount3 G026 and ARCHS4.

    9 / 499 genes:
        Known differences in annotation row membership.

The nine genes are:

    FAM20C
        ARCHS4-only: ENSG00000288499

    FCAR
        ARCHS4-only:
            ENSG00000283953
            ENSG00000284004
            ENSG00000284245

    ISLR
        ARCHS4-only: ENSG00000288508

    MAL2
        recount3-only: ENSG00000253972

    MATR3
        recount3-only: ENSG00000280987

    PODNL1
        ARCHS4-only: ENSG00000288195

    PSME2
        ARCHS4-only: ENSG00000284889

    SLC39A4
        ARCHS4-only: ENSG00000285243

    TNNT3
        ARCHS4-only: ENSG00000288250

DECISION:

    Gene-identity discovery is complete for the frozen 499-gene
    interface at the inspected annotation-membership level.

    Preserve all nine mismatches in a feature-provenance ledger.

    Evaluate the effect of these mismatches in a prespecified
    model-transfer sensitivity analysis.

    Matching Ensembl identifiers is NOT equivalent to matching
    measured expression abundance across assays.


4. ORIGINAL ARCHS4 SCRIPTS — AUTHORITATIVE EXECUTION CHAIN
----------------------------------------------------------

The following original scripts were identified and reviewed as the
authoritative Stage 09-13 methodology.

Original scripts are stored in:

    ARCHS4_V2_6/
      RELEASES/
        ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/
          PROJECT_TREE/

STAGE 06 — SOURCE EXPRESSION COMPATIBILITY

    06_audit_ARCHS4_primary_expression_compatibility_v1.py

    Validates original ARCHS4 HDF5 axes, identifiers and primary
    expression compatibility.

STAGE 09 — EXPRESSION PREPARATION

    09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py

    Version 1.0.1

    Performs:
        Source gene-symbol cleaning and aggregation
        Seven historical-gene aliases
        Canonical sample and anchor locking
        CPM normalization
        TMM-CPM normalization
        UQ-CPM normalization
        Sample normalization QC

    Uses all valid-symbol gene groups to calculate the normalization
    denominator, not just the 499 frozen W genes.

    For gene g and sample s:

        CPM[g,s] = 1,000,000 * X[g,s] / L[s]

    where L[s] is the total abundance over all valid collapsed genes.

    TMM adjusts the effective library size using trimmed M-values.
    UQ uses an upper-quartile-based effective library size.

    Original ARCHS4 UQ:
        9,302 evaluable samples
        163 fallback / non-evaluable samples

    Fallback UQ projections exist for diagnostics, but are not
    considered formally UQ-evaluable.

STAGE 10 — PAN-CANCER FIXED-ANCHOR PROJECTION

    10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R

    Version 1.0.4

    Performs:
        Frozen study/cancer-balanced anchor scaling
        Gene harmonization with unchanged W
        log2(CPM + 1)
        Frozen gene-wise standardization
        Directional feature construction with NMF::posneg
        SpatialEcoTyper::DeconvoluteSE
        Raw H coefficient extraction
        Normalized ecotype fraction calculation
        Reconstruction adequacy
        Normalization sensitivity

    Scaling:

        Y[g,s] = log2(CPM[g,s] + 1)

        Z[g,s] = (Y[g,s] - frozen_mean[g])
                 / frozen_sd[g]

    Frozen pan-cancer anchors:
        8,962 cross-sectional samples

    The UQ anchor subset excludes the non-evaluable fallback
    profiles.

    Original deconvolution call:

        SpatialEcoTyper::DeconvoluteSE(
            dat = scaled_matrix,
            scale = FALSE,
            W = W,
            nsample.per.run = 10000L,
            sum2one = FALSE,
            ncores = 1L
        )

    Raw nonnegative coefficients H are retained.

    Fractions are calculated afterward:

        fraction[k,s] = H[k,s] / sum_j(H[j,s])

    NonSE is a model component, not the reconstruction residual.

    The actual directional reconstruction residual is:

        R = V - W %*% H

STAGE 11 — ADEQUACY AND ROBUSTNESS

    11_build_ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1.R

    Version 1.0.2.1

    Performs:
        Relative reconstruction-fit ranking
        TMM and UQ normalization robustness
        Extreme frozen-Z diagnostics
        Evidence Tier A/B/C assignment
        Study-level transferability review
        Directional reconstruction-residual summaries

    TMM stability:
        TVD <= 0.10
        Same dominant model component

    UQ stability:
        UQ evaluable
        TVD <= 0.15
        Same dominant model component

    Frozen ARCHS4 evidence tiers:

        Tier A = 4,200
        Tier B = 1,026
        Tier C = 4,239

    These are empirical within-ARCHS4 categories.
    They are NOT externally calibrated validity probabilities.

    Stage 11 identified five poor-fit cancer-study combinations:

        colorectal_cancer / GSE152395
        colorectal_cancer / GSE183984
        colorectal_cancer / GSE299799
        lung_and_pleural_cancer / GSE126044
        lung_and_pleural_cancer / GSE190266

STAGE 12 — CANCER-SPECIFIC FIXED-ANCHOR PROJECTION

    12_run_ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3.R

    Version 1.0.3

    Performs:
        Within-cancer study-balanced anchor scaling
        Original frozen W projection
        Cancer-specific score fractions
        Adequacy and normalization robustness
        Pan-versus-cancer-specific coordinate comparison

    Cancer-specific anchor parameters are separate for each cancer.

    Cancer-specific scores are valid for within-cancer interpretation;
    they are not directly comparable between cancer classes.

STAGE 13 — ANALYSIS-READY ATLAS

    13_build_ARCHS4_projection_use_analysis_ready_atlas_v1_0_3.R

    Version 1.0.3

    Performs:
        Integration of frozen Stage 09-12 results
        Pan-cancer and cancer-specific score-space separation
        Evidence-tier propagation
        Study-level exclusion/caution policies
        Strict primary and supported-use manifest creation

    Stage 13 performs NO new SpatialEcoTyper projection.

    Main result:

        14_projection_use_analysis_ready_atlas_v1_0_3/
          03_projection_analysis_ready_core.tsv.gz

        9,465 analysis-ready ARCHS4 samples


5. FROZEN EXECUTION-CONTRACT AUDIT — VERIFIED PASS
--------------------------------------------------

The original Stage 09-13 implementation was audited against its frozen
source files and output contracts.

Audit package:

    ARCHS4_EXECUTION_CONTRACT_AUDIT_V1_1_20261119.zip

Authoritative cluster audit directory:

    ARCHS4_V2_6/
      31_frozen_execution_contract_audit_v1_1_20261119/

Audit script:

    audit_frozen_outputs_READONLY_v1_1.py

Primary audit output:

    frozen_output_contract_v1_1.json

FINAL RESULT, 19 NOVEMBER 2026:

    PASS_REFERENCE_OUTPUT_CONTRACT

    problems = []

Verified:

    All five Stage 09-13 authoritative source hashes match.
    All five stage status contracts report PASS.
    All five stages report zero QC failures.
    Frozen sample universe: 9,465.
    Cancer classes: 18.
    Cross-sectional anchors: 8,962.
    Stage 10 directional features: 500.
    Pan-cancer scaling tables: 499 rows per branch.
    Cancer-specific scaling: 26,946 records.
    Stage 11 Tier A/B/C: 4,200 / 1,026 / 4,239.
    UQ evaluable: 9,302.
    UQ non-evaluable: 163.
    Five Stage 11 study outliers reproduced exactly.
    Stage 13 pan-cancer score replay: maximum absolute difference 0.
    Stage 13 cancer-specific score replay: maximum absolute difference 0.

STAGE 11 SOURCE ADJUDICATION:

    An earlier reviewed Stage 11 source copy differed from the
    authoritative frozen source at one executable line:

        centre = centre    (reviewed copy; invalid R argument)
        center = centre    (authoritative frozen; valid argument)

    The correct R stats::mad argument is "center".

    The authoritative frozen source hash is:

        529551d7ecce253b1ab8626956afd076690631bf7d8b51cfcf0c613440a603f5

    This difference was resolved in audit v1.1.
    Do not replace or modify the actual frozen Stage 11 source.

AUDIT INTERPRETATION:

    The original Stage 09-13 source and inspected output contracts
    are verified and may serve as the methodological reference.

    The audit did NOT rerun the historical expression processing,
    deconvolution or downstream statistical analyses.

    It does NOT prove historical scheduler-job execution provenance
    solely from the surviving source and output files.

    It does NOT establish recount3 cross-assay equivalence.

NO ORIGINAL ARCHS4 ANALYSES NEED TO BE RESTARTED.


6. RECOUNT3 — NEXT EXPERIMENT AND RELEASE GATES
-----------------------------------------------

NEXT SCIENTIFIC COMPUTATION:

    Complete 499-gene frozen-model assay-transfer validation.

    Use the existing same-source ARCHS4/recount3 bridge specimens
    for matched-specimen technical comparison.

    The previously established bridge includes 143 matched
    lung-study specimens.

    These are a TECHNICAL TRANSFER cohort, not 143 independent
    biological validation specimens.

Required transfer assessment:

    A. Build the 499 recount3 gene features using the verified
       Stage 09 symbol-collapse rule.

    B. Record all nine annotation row-membership differences.

    C. Calculate CPM, TMM_CPM and UQ_CPM using all valid G026
       expression features, while retaining formal UQ evaluability.

    D. Apply the unchanged ARCHS4 frozen anchor means and SDs.

    E. Apply the original directional-feature transformation and
       unchanged SpatialEcoTyper W projection.

    F. Compare matched-specimen outputs:
           Gene-level input values
           Raw ten-component H coefficients
           Normalized ecotype fractions
           Dominant ecotypes
           Total variation distance
           Calibration bias and absolute differences
           Reconstruction adequacy
           Normalization-branch robustness

    G. Quantify the effects of the nine annotation differences.

    H. Issue an explicit transferability decision supported by
       numerical results and uncertainty.

    The previous 465-gene comparison is an exploratory assay bridge,
    NOT a certified complete frozen-W transfer.

CRITICAL REMAINING COHORT RELEASE GATES:

    Biological specimen eligibility
    Primary / metastatic / treated / enriched context
    Assay compatibility
    Run-level QC and run aggregation
    Patient and donor identity
    Repeated biopsy and multiregion grouping
    Cross-study and cross-atlas independence
    Frozen 499-gene model-transfer validation
    Statistical analysis eligibility

Passing the original ARCHS4 execution-contract audit does NOT
automatically approve the recount3 cohort.

No patient-level statistical replication, refitting, pooled
pan-cancer inference or new biological discovery should proceed
until the corresponding recount3 certification decisions are made.


7. SCRIPTS AND EVIDENCE USED IN THESE INVESTIGATIONS
-----------------------------------------------------

ORIGINAL ARCHS4 SOURCE SCRIPTS:

    06_audit_ARCHS4_primary_expression_compatibility_v1.py

    09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py

    10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R

    11_build_ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1.R

    12_run_ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3.R

    13_build_ARCHS4_projection_use_analysis_ready_atlas_v1_0_3.R

RECOUNT3 AND INDEPENDENCE INVESTIGATIONS:

    recount3_accession_source_evidence_v1_0_0.py

    probe_W34_and_ARCHS4_atlas_READONLY.py

    Read-only full-HDF5 Ensembl-axis diagnostic
    Read-only frozen-W feature provenance diagnostic
    Read-only original source-code inventory
    Read-only 34-gene symbol-group comparison
    Read-only full-499 gene membership comparison

    The latter five diagnostics included inline read-only code
    executed on the cluster; they are not necessarily separate
    versioned script files.

FROZEN EXECUTION-CONTRACT VALIDATORS:

    audit_frozen_outputs_READONLY.py          (v1)

    audit_frozen_outputs_READONLY_v1_1.py     (authoritative v1.1)

SUPPORTING PACKAGES:

    ARCHS4_FROZEN_GENE_PROVENANCE_SOURCE.tar.gz

    RECOUNT3_W34_ATLAS_CANDIDATES.tar.gz

    RECOUNT3_W34_OVERLAP_EVIDENCE_REVIEW_20261119.zip

    ARCHS4_EXECUTION_CONTRACT_AUDIT_V1_1_20261119.zip

    ARCHS4_STAGE11_FROZEN_SOURCE.tar.gz

For original script execution, normalization, formulas, file names,
reconstruction metrics and score-space restrictions, consult:

    EXECUTION_CONTRACT_SOURCE_AUDIT_V1_1.md

For the final machine-readable frozen reference checks, consult:

    31_frozen_execution_contract_audit_v1_1_20261119/
      frozen_output_contract_v1_1.json


8. CURRENT STATUS / DO NOT REPEAT COMPLETED WORK
-------------------------------------------------

COMPLETED:

    Original ARCHS4 analyses and frozen atlas.
    Original Stage 09-13 methodological reconstruction.
    Full frozen model feature-identity review.
    Identification of nine G026/ARCHS4 annotation differences.
    Exact-GSM frozen-atlas overlap screening.
    Source hash and output-contract audit v1.1: PASS.

NOT YET COMPLETED:

    Complete 499-gene recount3 model-output transfer.
    Formal calibration and sensitivity analysis.
    Full source- and patient-level recount3 certification.
    Global independence certification for the remaining candidates.
    Final biological release of the recount3 primary cohort.

DO NOT:

    Rerun the original ARCHS4 discovery or statistical pipeline.
    Modify frozen W or anchor scaling.
    Treat Kallisto pseudocounts and gene-body coverage as identical.
    Substitute the 465-gene bridge for the full 499-gene projection.
    Interpret NonSE as reconstruction error.
    Compare cancer-specific ecotype scores across cancer classes.
    Treat exact-GSM overlaps as independent validation.
    Treat remaining unmatched accessions as proven independent.
    Interpret technical PASS as biological cohort approval.
    Automatically apply ARCHS4 empirical adequacy tiers to recount3.
    Approve manual workbook rows merely because they passed QC.

NEXT ACTION:

    Implement and execute ONE versioned, auditable, full-499-feature
    recount3 frozen-model transfer experiment using the verified
    original Stage 09-13 execution contract.

    Maintain the 5,241-library manual curation workbook and the
    original frozen ARCHS4 release unchanged as reference evidence.

================================================================================
END OF 19 NOVEMBER 2026 UPDATE
================================================================================

================================================================================
UPDATE — 20 NOVEMBER 2026
RECOUNT3 IDENTITY, FUSCC QC, FROZEN-W SOURCE AND TRANSFER PREFLIGHT
================================================================================

PURPOSE AND STATUS
------------------
This appendix continues the 19 November record. It consolidates findings from
read-only cluster diagnostics, the Stage 40/41 evidence integration, FUSCC
cohort source review, the original frozen Stage 06/09-13 scripts and the
20 November transfer-contract preflight. It is a provenance/status update,
NOT an approval of recount3 model inputs or a completed fixed-W transfer.

Preserve the original manual-curation universe (5,241 library records,
6,337 runs, 176 studies) and the frozen ARCHS4 release unchanged. The
3,602-library subset below is a later selection; do not overwrite the
original workbook or reinterpret its unsigned review fields.

1. COMPLETE ARCHS4 ACCESSION-LEVEL IDENTITY RECONCILIATION
----------------------------------------------------------
The frozen analysis-ready ARCHS4 atlas has 9,465 canonical GSM samples.
An independent read-only join recovered an exact SRX experiment and BioSample
link for every atlas GSM (9,465/9,465); two independent source-metadata
copies agreed. Joining these identities against recount3 found:

    850 exact frozen-atlas overlaps in the 3,602 recount3 candidates;
    0 additional candidate overlaps from the recovered SRX/BioSample links;
    299 additional same-SRX overlaps among noncandidate source experiments;
    0 exact atlas overlaps among the 360 FUSCC tumour RNA libraries.

The 850 candidate overlaps must remain EXCLUDED from strict independent
ARCHS4 external validation, although some may serve technical calibration.
A negative accession match does not prove distinct human patients. The
additional 299 noncandidate source overlaps must not be added to the 850
candidate-library count or mistaken for newly eligible candidates.

Stage 41 v1.0.1 contains 3,602 candidates across 126 studies, drawn from
5,221 full-source experiments. Its unsigned screening strata were:

    2,055  no detected study exposure and no historical run hold;
      439  different samples from studies with detected atlas exposure;
      258  technical holds outside the exact-overlap category;
      850  exact atlas accession overlaps;
    -----
    3,602  TOTAL.

These are screening strata, not signed patient- or treatment-level approvals.
Stage 41 reviewer decisions were not signed and cohort_approved was NO.

2. FUSCC / SRP157974: SOURCE AND BIOLOGICAL ELIGIBILITY
------------------------------------------------------
The original source contains 448 experiments; 360 selected candidate records
are distinct tumour RNA-seq experiment, BioSample and RNA-specimen tokens.
The remaining 88 experiments numerically match the publication's adjacent-
normal count, but their individual normal-tissue labels were not all
independently verified in this audit.

Published FUSCC cohort descriptions support 360 primary TNBC tumours from
360 previously untreated patients and 88 adjacent normal RNA-seq specimens;
the source names label the selected 360 as tumour tissue. The reported
pre-treatment and TruSeq Stranded Total RNA / Ribo-Zero Gold properties are
cohort-level source evidence, not independently reviewed medical records for
each RNA library.

    360 source tumour RNA-seq libraries;
    334 libraries without historical technical holds;
     26 libraries with historical technical holds.

Among the 334, 279 have patient-token corroboration from the independent GEO
OncoScan assay; 55 have distinct RNA isolate tokens without that extra assay
corroboration. Use 334 as the PROVISIONAL primary RNA-only biological screen;
279 is a prespecified cross-assay-token sensitivity subset, not an automatic
restriction that would select patients by availability of OncoScan data.
Distinct specimen labels and the one-tumour-per-patient study description
support study-scoped patient grouping but are not genetic proof of distinct
individuals. The original SRA selection labels (318 cDNA, 42 Inverse rRNA)
and curated rRNA-depletion labels differ; do not infer per-library assay
failure from generic SRA metadata without the published protocol evidence.

Keep WES-context warnings for RNA profiles associated with FUSCCTNBC044,
140, 030 and 032. The published WES identity concerns relate to 044_WBC,
140_WBC and 030_TT; 032 is a comparison tumour, not a fourth proven WES
mislabel. These warnings do not establish an RNA sample swap. If clinically
relevant associations depend on cross-assay identity, include an explicit
four-record sensitivity analysis rather than silently dropping them.

3. FUSCC ORIGINAL TECHNICAL HOLDS: ROOT CAUSES RESOLVED
------------------------------------------------------
Read-only joins against the existing Stage 2, 3, 5 and run-hold queue confirm:

    All 360: 1 expected sequencing run, 1 verified expression run,
             run_count_status = PRESENT_ALL.
    All 360: stage3 library_aggregation_state = NOT_PERFORMED,
             as expected for one run per library.
    334 nonheld: zero low-detection, top-10-concentration,
             within-study-outlier or mapped-length-below-40-bp run flags;
             empty Stage 5 run_qc_review_reasons.
     25 held: WITHIN_STUDY_OUTLIER only.
      1 held: WITHIN_STUDY_OUTLIER plus top-10-gene concentration;
             experiment SRX5321794.
    All 360: zero low-detection and zero short-mapped-length flags.

Detected genes per run:
    NONHELD n=334: median 39,703.5; range 34,274–44,026.
    HELD    n=26: median 38,033;   range 33,896–50,128.

The 26 holds are expression-distribution outlier decisions, NOT missing
extraction, low gene detection or absent runs. Retain them outside the
primary QC screen; assess separately for biological versus technical
outliers. No multirun correlation is applicable to this one-run cohort.

These findings support technical SCREENING of the 334, not final model-input
approval, clinical certification or completed normalization compatibility.

4. MODEL SOURCE: EXACT FROZEN PROJECTION CONTRACT RECOVERED
----------------------------------------------------------
Original authoritative code locations are below:

    ARCHS4_V2_6/RELEASES/
      ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE/

Original authoritative scripts:
    06_audit_ARCHS4_primary_expression_compatibility_v1.py
    09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py
    10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R
    11_build_ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1.R
    12_run_ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3.R
    13_build_ARCHS4_projection_use_analysis_ready_atlas_v1_0_3.R

The frozen model asset in RELEASE/02_MODEL/Bulk_SE_Recovery_W.rds has SHA256:
    39a96b598ec74ccbe67bea49254ac94ae568ca22ec162241ab7849976b52cdcf

W dimensions: 500 directional feature rows x 10 components; 499 distinct
model genes; one gene, PPIB, has both positive and negative features.
Component columns are NonSE and SE1 through SE9. NonSE is a genuine model
component and MUST NOT be treated as the reconstruction residual.

Original Stage 09 input harmonization:
    - clean/uppercase symbols and SUM source rows by normalized symbol;
    - preserve all 499 distinct frozen model genes;
    - apply seven historical aliases:
         CXORF36 -> DIPK2B; FAM71A -> GARIN4; GPR1 -> CMKLR2;
         H2AFY -> MACROH2A1; MRVI1 -> IRAG1;
         RARRES3 -> PLAAT4; WARS -> WARS1;
    - compute normalization over ALL valid source gene groups, not only W;
    - use CPM, TMM_CPM, UQ_CPM with original evaluability conventions.

Pan-cancer Stage 10 calculation:
    Y[g,s] = log2(CPM[g,s] + 1)
    Z[g,s] = (Y[g,s] - FROZEN_ANCHOR_MEAN[g]) / FROZEN_ANCHOR_SD[g]

The original signed, ordered 499-gene Z matrix is passed to:
    SpatialEcoTyper::DeconvoluteSE(
        dat = scaled_matrix,
        scale = FALSE,
        W = W,
        nsample.per.run = 10000L,
        sum2one = FALSE,
        ncores = 1L
    )

The model procedure constructs the correct 500 directional features using
the original NMF::posneg convention. Preserve raw nonnegative H coefficients
before calculating fractions[k,s] = H[k,s] / sum_j H[j,s]. Evaluate the
actual directional reconstruction V - W %*% H, reconstruction adequacy and
normalization-branch robustness. Do not replace the original estimator with
an unvalidated NNLS or hand-assembled 500-feature approximation.

The frozen breast_cancer Stage 12 anchors must be used for final FUSCC
cancer-specific fractions; pan-cancer anchors are the common coordinate
system for initial lung-study same-specimen calibration. Cross-cancer
comparison of cancer-specific score coordinates is not permitted.

5. GENE IDENTITY: CORRECTION TO THE EXPLORATORY 465-GENE BRIDGE
--------------------------------------------------------------
Stage 5's strict one-to-one mapping has 1,497 branch records: 499 genes x
3 branches; 1,428 unambiguous and 69 ambiguous records, the latter covering
23 distinct symbols. An older exploratory bridge used 465 genes, excluding
23 recount3 symbol ambiguities plus 11 ARCHS4 alias/annotation cases.
These 34 are NOT shown to be missing from the original W input.

The original Stage 09 symbol-COLLAPSE and seven-alias rules resolve the
full 499-gene interface. The separate full-499 Ensembl-membership audit,
as documented in the 19 November README, found 490/499 identical source-
row-ID multisets and NINE differences:

    FAM20C, FCAR, ISLR, MAL2, MATR3,
    PODNL1, PSME2, SLC39A4, TNNT3.

Keep all nine as matched gene-symbol features in the intended primary
projection; record differing contributing Ensembl rows and quantify an
explicit nine-gene sensitivity. Do not arbitrarily select one Ensembl row,
zero a gene, drop it, or refit W. The membership-audit finding was
previously established; the source-code archive independently confirmed
the original collapse/alias method, not a new rerun of both complete axes.

The nine genes represent approximately 4.1% of NonSE W loading, 3.85% of
SE5 loading and 2.42% of SE2 loading (lower in other components). These
loadings do not predict the numerical fraction error.

6. CALIBRATION AND EXPRESSION INPUT READINESS — 20 NOVEMBER
-----------------------------------------------------------
A read-only cluster preflight verified:

    frozen W SHA256: MATCH;
    pan-cancer CPM, TMM_CPM, UQ_CPM scaling: 499 rows per branch;
    breast_cancer cancer-specific scaling rows: PRESENT;
    recount3 gene-order SHA256:
      3e4bde5146353005a4cdcbe7d5ec3dcc95dbc6f9abdefe8d303014ad8a5c7d69;
    SRP074349 calibration HDF5: 63,856 genes x 144 selected runs;
    SRP157974 FUSCC HDF5:     63,856 genes x 360 selected runs;
    HDF5 datasets: coverage_bp, gene_id, library_id, run_id;
    HDF5 orientation: genes x runs;
    SRP074349 documented matched bridge: 143 specimens;
    SRP157215: 22 candidate bridge rows, but its expression source file
      was not present in the prior bridge-availability check.

Relevant recount3 input paths (relative to ARCHS4_V2_6):
    26_recount3_external_primary_v1/stage2_output/studies/
      SRP074349/selected_run_coverage_G026.h5
      SRP157974/selected_run_coverage_G026.h5
    26_recount3_external_primary_v1/stage5_output/
      03_bridge_165_availability.tsv
      00_gene_row_mapping.tsv
      04_library_preanalysis_decision_3639.tsv

The G026 measure is gene-body BASE-PAIR COVERAGE, not gene counts. The
ARCHS4 source is Kallisto-derived pseudocounts. An estimator that reads
coverage as if it were counts is a transfer hypothesis requiring direct
matched-specimen tests and, if appropriate, prespecified gene-length/
coverage sensitivity; it is NOT established as assay-equivalent by shared
gene IDs or a matching W SHA256.

An earlier exploratory 465-gene comparison across 143 matched specimens
reported median per-gene Pearson r about 0.999 and median within-sample
across-gene r about 0.984, but several genes had major absolute shifts.
It did NOT compare complete 499-gene raw H, ecotype fractions, dominance
or calibrated reconstruction metrics. The matched 143 are technical
TRANSFER pairs, not 143 independent external biological validations.

At the end of this update:
    Stage 5 status: TECHNICAL_EVIDENCE_COMPLETE_MODEL_INPUT_NOT_APPROVED;
    FUSCC model_input_approved: False for all 360;
    patient_independence_certified: False for all 360;
    ARCHS4_independence_certified: False for all 360;
    recount3_to_ARCHS4_calibration_passed: False for all 360.

These flags are historic unsigned gate records; no new calibration score,
projection result or certification was written by the read-only preflight.

7. NEXT COMPUTATION — ONE FROZEN FULL-499 TRANSFER EXPERIMENT
-------------------------------------------------------------
Before using FUSCC ecotype outcomes, run a separately versioned calibration
ONLY on the 143 matched SRP074349 specimens, with no writes to the frozen
ARCHS4 release and no refitting W, anchors or model components:

    A. Verify 143 exact ARCHS4 GSM/index to recount3 SRX/run pairs.
    B. Reconstruct all 499 source gene symbols using frozen Stage 09 rules.
    C. Generate CPM/TMM/UQ branches with full-gene denominators and original
       frozen standardization; preserve UQ evaluability and fallback flags.
    D. Use exact original SpatialEcoTyper DeconvoluteSE calls; retain H,
       fractions and directional reconstruction diagnostics.
    E. Compare to original ARCHS4 paired outputs: H and fractions,
       fraction total-variation distance, component rankings/dominance,
       systematic bias, gene-level shifts, fit and normalization robustness.
    F. Quantify uncertainty and the specific effect of the nine annotation
       membership differences and coverage-vs-pseudocount representation.
    G. Prespecify PASS / FAIL / INDETERMINATE transfer rules; do not select
       thresholds retrospectively from FUSCC ecotype or clinical outcomes.

Only with a defensible technical transfer decision should the 334 nonheld
FUSCC libraries proceed to FROZEN projection, using the unchanged breast-
cancer-specific anchors for breast-only interpretation. Keep the 26 QC-held
outliers as a separate sensitivity stratum. Biological/patient identity,
assay and source documentation must then receive explicit sign-off before
external outcome associations or patient-level validation claims.

Frozen Stage 09-13 reference and atlas: DO NOT MODIFY / RERUN.
The unified residual discovery configuration is a DIFFERENT pipeline and
must not be used as a substitute for frozen-W projection implementation.

8. RELEVANT REPRODUCIBILITY ARTIFACTS
------------------------------------
Cluster project:
    /maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6

Frozen release:
    RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0

Recount3 project:
    26_recount3_external_primary_v1

Frozen execution audit:
    31_frozen_execution_contract_audit_v1_1_20261119/
      frozen_output_contract_v1_1.json

Uploaded/reviewed supporting packages:
    RECOUNT3_STAGE41_V101_REVIEW_RESULTS.tar.gz
    FUSCC_COHORT_RECONCILIATION_UNSIGNED_20261120.zip
    FUSCC_MODEL_COMPATIBILITY_EVIDENCE.tar.gz
    FUSCC_FROZEN_W_PROJECTION_METHODS.tar.gz
    FUSCC_ORIGINAL_PROJECTION_SOURCE.tar.gz
    FUSCC_FIXED_W_SOURCE_REVIEW_AND_PREFLIGHT_20261120.zip

An archive name denotes an evidence package, NOT successful execution.
Use the exact frozen source files and live cluster QC tables as the
computational authorities; review reports preserve interpretation and
outstanding caveats.

================================================================================
END OF 20 NOVEMBER 2026 UPDATE — TRANSFER CALIBRATION NOT YET EXECUTED
================================================================================
