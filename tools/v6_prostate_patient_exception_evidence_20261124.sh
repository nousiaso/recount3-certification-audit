#!/usr/bin/env bash
#
# V6 PROSTATE — ORIGINAL GEO PATIENT-ID EXCEPTION
# DATE: 2026-11-24
# STUDY: SRP163173
#
# STATUS:
#   HOLD_PATIENT_ID_SOURCE_CONFLICT
#
# GEO study: GSE120741
#
# Exact original GEO-to-SRX relations:
#   92 / 92 PASS
#
# V5 selected experiments:
#   85
#
# Technical QC holds:
#   7
#
# Patient annotation review:
#   91 of 92 GEO patient IDs agree with their corresponding
#   sample-title identifiers.
#
#   One original GEO record has a patient-ID versus
#   sample-title discrepancy.
#
#   The discrepancy affects a repeated V5 patient-key group
#   containing two experiments.
#
#   The two records have distinct literal GEO patient IDs.
#   Patient independence is NOT established.
#
#   The patient's identity must not be inferred from the
#   shared title base or the provisional V5 grouping.
#
# Exact exception:
#   1 source annotation discrepancy
#   2 affected experiments
#
# LOCAL EVIDENCE:
# diagnostics/20261124_v6_prostate_patient_exception_v1/
#
# Local exception TSV SHA256:
# 1007d95372fce3dee69102127b90a20386953e255204d800d18fd9f1a4fe146f
#
# Summary JSON SHA256:
# 307296075e21f6be8ccc3e2b852d01dc186452acaedf7d534f29b7ec64895900
#
# The source conflict is localized to the affected patient group.
# The original GEO/SRX joins remain valid.
#
# No frozen patient keys changed.
# Technical QC holds retained.
# No clinical approvals.
#
# This file records previously completed cluster results.
# It is NOT the original diagnostic executable.
# Patient-linked local evidence remains outside Git.
#
# DECISION:
# HOLD_PATIENT_ID_SOURCE_CONFLICT_NOT_CERTIFIED
#
# REPEAT POLICY:
# DO_NOT_RERUN_UNLESS_NEW_INDEPENDENT_IDENTITY_EVIDENCE
#
exit 0
