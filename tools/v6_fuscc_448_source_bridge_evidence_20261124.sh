#!/usr/bin/env bash

# =============================================================================
# V6 FUSCC TNBC — COMPLETE RNA-SEQ SOURCE-PROVENANCE CERTIFICATION
# =============================================================================
#
# DATE: 2026-11-24
# STUDY: SRP157974
#
# FINAL SOURCE-IDENTIFICATION STATUS:
#   EXACT_SRA448_TO_FIGSHARE448_CASE_NORMALIZED_PASS
#
# CLINICAL STATUS:
#   NO_UNSIGNED
#
# PURPOSE:
#   Preserve the completed original-SRA-to-published-RNA-seq reconciliation.
#   This file records already executed cluster diagnostics.
#   It does not rerun them, alter frozen metadata or certify patients.
#
# PROVENANCE LIMITATION:
#   Findings below are transcribed from actual cluster outputs supplied
#   during the V6 investigation. This documentation script is not the
#   original executable diagnostic code.
#
# =============================================================================
# 1. ORIGINAL RNA-SEQ STUDY UNIVERSE
# =============================================================================
#
# Original SRA study:
#   448 distinct RNA-Seq experiments
#   448 metadata rows
#   All 448 have library_strategy = RNA-Seq
#
# Frozen V5 tumour-candidate cohort:
#   360 experiments
#
# Stage41 nonheld tumour candidates:
#   334 experiments
#
# Technical QC holds:
#   26 experiments
#   All retain HOLD_TECHNICAL_QC
#
# Outside V5:
#   88 additional SRA RNA-seq experiments
#
# No technical QC hold was released.
# No frozen cohort membership was altered.
#
# =============================================================================
# 2. EXACT ORIGINAL-SRA ISOLATE EVIDENCE
# =============================================================================
#
# All 360 V5 experiments:
#   Have one unambiguous original-SRA isolate identifier
#   Have 360 distinct isolate identifiers
#   Have no isolate normalization collisions
#
# The 334 nonheld experiments:
#   Have 334 distinct isolate identifiers
#   Follow the FUSCCTNBC sample naming convention
#
# The 26 QC-held experiments:
#   Have distinct isolates outside the 334 selected isolates
#
# The original SRA sample titles do NOT establish patient identity.
# SRA isolate uniqueness alone does NOT establish patient independence.
#
# =============================================================================
# 3. PUBLISHED FIGSHARE EXPRESSION MATRIX
# =============================================================================
#
# Dataset:
#   Figshare article 19783498, version 5
#
# File:
#   FUSCCTNBC_Expression_RNAseqCount_hg38.tsv
#
# Figshare file ID:
#   35150692
#
# Published matrix:
#   448 expression columns
#   360 distinct FUSCCTNBC sample bases
#   No duplicated column names
#
# Observed published suffix counts:
#   NO_SUFFIX       331
#   _REP             29
#   _PT              83
#   _PT_REP           5
#   TOTAL           448
#
# The original study reports:
#   360 tumour RNA-seq samples
#    88 matched adjacent-normal RNA-seq samples
#
# =============================================================================
# 4. SOURCE LIBRARY STRUCTURE
# =============================================================================
#
# ORIGINAL SRA LIBRARY NAMES:
#
# V5 tumour candidates (360):
#   No suffix       331
#   .REP             29
#
# Selected tumour candidates (334):
#   No suffix       307
#   .REP             27
#
# Technical QC holds (26):
#   No suffix        24
#   .REP              2
#
# Outside-V5 libraries (88):
#   .PT              83
#   .PT.REP           5
#
# All 88 outside-V5 PT-labelled library bases have a matching
# tumour-associated base among the 360 V5 specimens.
#
# One paired PT base corresponds to a technically held tumour candidate.
#
# Every tumour base has exactly one non-PT expression column.
# Every paired PT base has exactly one PT expression column.
#
# REP-labelled columns are NOT counted as additional patients.
#
# =============================================================================
# 5. INITIAL 448-TO-448 COMPARISON
# =============================================================================
#
# Original SRA library names recognized:
#   448 / 448
#
# Original published expression columns:
#   448
#
# Initial exact matches:
#   414 / 448
#
# Initial unmatched source-derived names:
#   34
#
# Initial unmatched matrix names:
#   34
#
# This produced:
#   HOLD_SRA_FIGSHARE_SPECIMEN_MISMATCH
#
# The HOLD was retained pending investigation.
# No sample names were patched in frozen metadata.
#
# =============================================================================
# 6. ROOT CAUSE OF THE 34 UNMATCHED NAMES
# =============================================================================
#
# Original SRA-derived suffix  Published matrix suffix  Count
# --------------------------  -----------------------  -----
# REP                         rep                         29
# PT_REP                      PT_rep                       5
#
# Total mismatched name pairs:
#   34
#
# Distinct affected sample bases:
#   31
#
# Partition membership of 34 affected libraries:
#   Selected tumour candidates: 27
#   Technical QC holds:          2
#   Outside-V5 PT libraries:     5
#
# All 34 source isolate bases agree with the corresponding
# original SRA library bases.
#
# All 360 sample bases occur in both sources.
#
# ROOT CAUSE:
#   Suffix capitalization difference.
#
# NOT observed:
#   Biological specimen mismatch
#   Missing SRA experiment
#   Missing published sample base
#   Source identifier collision
#
# =============================================================================
# 7. FINAL CASE-NORMALIZED RECONCILIATION
# =============================================================================
#
# After explicitly documented normalization of suffix case:
#
# Original SRA experiments:       448
# Published expression columns:   448
#
# Exact matches before correction:
#   414
#
# Exact matches after correction:
#   448
#
# Residual unmatched:
#   0
#
# Normalized identifier collisions:
#   0
#
# Final partitions:
#   Selected tumour experiments: 334
#   QC-held tumour experiments:   26
#   Additional paired PT:         88
#
# Distinct tumour bases:
#   360
#
# Distinct PT bases:
#   88
#
# PT bases linked to tumour bases:
#   88
#
# FINAL STATUS:
#   EXACT_SRA448_TO_FIGSHARE448_CASE_NORMALIZED_PASS
#
# The earlier 34-name HOLD is superseded by this documented correction.
# It must NOT be deleted from the historical audit trail.
#
# =============================================================================
# 8. FROZEN SOURCE PROVENANCE
# =============================================================================
#
# Original SRA source:
#   26_recount3_external_primary_v1/raw_studies/SRP157974/
#   metadata/sra.sra.SRP157974.MD.gz
#
# Original SRA source SHA256:
#   c2e1c7b5b1169e13059447e14bdc71ed3ffd4166dd90131b2791360ebab28155
#
# Published Figshare RNA-seq matrix header SHA256:
#   8b91e34904251263f02ba89be9eec2aa1cefa009caa3ba4423110e16224fc254
#
# Figshare file ID:
#   35150692
#
# Frozen V5 5,241-row source SHA256:
#   ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a
#
# These hashes are transcribed from the verified cluster diagnostics.
# This comment-only documentation does not independently rehash the files.
#
# =============================================================================
# 9. PATIENT IDENTITY — NOT CERTIFIED
# =============================================================================
#
# The selected 334 tumour experiments have:
#
#   279 SOURCE_CROSS_ASSAY_PATIENT_TOKEN_PROPOSAL
#    55 STUDY_UNIQUE_PATIENT_COHORT_WITH_RNA_ISOLATE_ONLY
#
# Both groups:
#   Match published RNA-seq sample identifiers
#   Have distinct SRA isolate identifiers
#   Lack populated V5 independent patient-crosswalk evidence fields
#
# The original-SRA isolates are embedded in the Stage41 patient proposals,
# but not equal to the frozen V5 patient keys.
#
# No independently validated RNA-to-WES patient mapping was established.
#
# Numerical correspondence of 279 proposals to the published 279 WES
# patients is not proof of cross-assay patient identity.
#
# Original-study WES sample-label discrepancies must remain a separate
# cross-assay identity gate.
#
# REQUIRED:
#   Authoritative SRX/isolate-to-clinical-patient mapping
#   Verified cross-assay sample membership where claimed
#   Patient grouping and independence certification
#
# =============================================================================
# 10. TREATMENT AND CLINICAL CERTIFICATION
# =============================================================================
#
# Published cohort-level treatment evidence exists.
#
# However:
#   Exact 334-library pretreatment status is not individually certified.
#   Treatment timing remains a separate clinical gate.
#   Patient identity remains a separate clinical gate.
#   Technical QC holds remain unchanged.
#   No final clinical sample approvals exist.
#
# FINAL_CLINICAL_APPROVAL:
#   NO_UNSIGNED
#
# CLINICAL APPROVALS:
#   0
#
# No SpatialEcoTyper modelling has resumed.
#
# =============================================================================
# 11. FINAL SCIENTIFIC DECISION
# =============================================================================
#
# CLOSE:
#   FUSCC original-SRA accession identity
#   Original SRA isolate membership
#   Original SRA library-name reconciliation
#   Published RNA-seq matrix sample membership
#   Full 448-to-448 case-normalized source mapping
#   360-to-334 selection accounting
#   26 technical QC hold accounting
#   88 additional PT library accounting
#
# KEEP OPEN:
#   Patient independence
#   279 proposed cross-assay identities
#   55 RNA-isolate-based patient identities
#   Individual specimen pretreatment evidence
#   Final clinical approval
#
# DECISION:
#   SOURCE_IDENTIFICATION_PASS_NOT_CLINICALLY_CERTIFIED
#
# REPEAT POLICY:
#   DO_NOT_RERUN_WITH_UNCHANGED_SOURCE_INPUTS
#
# =============================================================================
# END OF AUDIT
# =============================================================================

exit 0
