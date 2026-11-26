# SRP057793 authoritative QC schema and identity join

Date: 2026-11-26
Decision: EIGHT_HELD_JOIN_TO_ORIGINAL_RUN_HOLD_QUEUE

The source has 149 experiments.
Stage41 has 148 candidates, including eight technical holds.
V5 primary selection contains 140 nonheld candidates.

Original Stage2 library and run QC tables cover all
eight technical holds.

The 409-row historical run-hold queue contains exactly
the eight held SRP057793 experiments, with populated
reason_codes in all eight matched rows.

The source-only experiment is absent from these Stage2
QC records. That does not establish its exact exclusion cause.

Next: extract aggregate original hold reason codes and
inspect source-only run/BioSample evidence.

Pretreatment: UNSIGNED.
Patient independence: UNSIGNED.
Clinical admissions: ZERO.
