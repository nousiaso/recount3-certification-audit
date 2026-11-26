# SRP057793 exclusion and BioSample audit

Date: 2026-11-26
Decision: NONHELD_PRIMARY_EXACT_MATCH

149 original experiments.
148 Stage41 candidates.
140 nonheld candidates.
8 held candidates.
140 exact V5 primary-screen members.

The 140 nonheld Stage41 experiment IDs
exactly equal the V5 primary-screen IDs.

One original experiment is absent from Stage41.
Its recorded biosamples_all_sources value is empty.

The original full-source table has 148 populated,
distinct recorded BioSample values. These are not
proof of 148 independent biological patients.

Stage41 does not expose explicit technical QC
reason fields in the examined candidate schema.

Next: inspect original-only identity provenance
and locate the authoritative eight hold reasons.

Pretreatment: UNSIGNED.
Patient independence: UNSIGNED.
Clinical admissions: ZERO.
