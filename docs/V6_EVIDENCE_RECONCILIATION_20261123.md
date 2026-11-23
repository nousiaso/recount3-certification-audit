# V6 cumulative specimen-evidence reconciliation

Status: EVIDENCE ONLY — UNSIGNED — NO CLINICAL ADMISSIONS.

Input: V5 administrative primary screen, 2,055 experiments in 77 studies.

Verified categories:
- 22 cumulative non-tumour evidence restrictions, including SRX106046 from historical Stage 41.
- 5 SRP060016 primary colorectal specimen identities established from raw SRA attributes.
- 52 sarcoma clinical metastasis flags with sampled lesion site unproven.
- 280 other V5 tumour/context evidence records.
- 1,696 records requiring historical evidence reconciliation.

Administrative remainder after non-tumour restrictions: 2,033,
NOT a clinically certified tumour population.

Cumulative non-tumour evidence includes 21 V5 flags plus one
additional Stage 41 normal-colon annotation. Exact source and
patient relationship for SRX106046 require targeted confirmation.

The five categories were prespecified in the diagnostic script.
Their counts passing assertions is a consistency check, not
independent clinical validation.

Complete output remains on the cluster:
diagnostics/20261123_v6_evidence_precedence_dryrun/proposed_evidence_2055_UNSIGNED.tsv

Output SHA256:
08efe2c8f469edcbdfb98231dc87f705d801434e78dcbea85f57c4b94da88eb0

Frozen ARCHS4 unchanged. Full-499 model transfer not certified.
