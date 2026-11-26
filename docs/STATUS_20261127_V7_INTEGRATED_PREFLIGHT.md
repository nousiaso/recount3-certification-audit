# V7 integrated 2,055-library gate preflight — 27 November 2026

Status: **FROZEN_EVIDENCE_BASELINE_RECONCILED; CLINICAL_RELEASE_UNSIGNED**

Source: user-provided cluster stdout from `scripts/run_v7_integrated_preflight.sh`. The values below reproduce the stdout; the exact cluster-generated stdout and files are retained on supekgate, not independently fetched by this GitHub write.

- Frozen V5 primary screen: 2,055 experiments, 77 studies.
- Source non-tumour exclusions: 22.
- Historic primary-specimen-specific evidence, all other gates open: 5.
- Sampled-site uncertainty (sarcoma): 52.
- V5 tumour/context evidence, gates open: 280.
- Historic evidence review pending: 1,696.
- Held pending final gates: 2,033.
- Clinical admissions: **0**.
- SpatialEcoTyper projection executed in this step: **NO**.

Local-only outputs (do not commit or disclose sample-linked rows):
`private_data/v7_integrated_candidate_gates_v1/2055_gate_matrix_LOCAL_ONLY.tsv`
`private_data/v7_integrated_candidate_gates_v1/77_study_triage_LOCAL_ONLY.tsv`
`private_data/v7_integrated_candidate_gates_v1/summary.json`

Output hashes reported by cluster:
- 2,055-gate matrix: `1ab9dd3f2c829d03e8b6f56ebbb45bb135ffd107ce54e6c2e884db5d7b4ca295`
- 77-study triage: `2dd3530915a40163592ae6f35185383931abf53cb39483f0114c551c288a75d3`

Next computational milestone is **full 499-gene frozen model assay-transfer calibration** on the existing 143 paired SRP074349 ARCHS4/recount3 specimens. These 143 pairs are assay-calibration controls, **not independent clinical validation**. Recount3 G026 is gene-body base-pair coverage, not Kallisto pseudocounts. Do not treat assays as equivalent without matched H/fraction/adequacy and normalization-branch comparisons. Only after technical transfer is supported should projection on a clinically certified external candidate subset be used for inferential validation.

Study SRP157974 (FUSCC TNBC) has 334 historically non-held tumour libraries and cohort-level pretreatment support, but patient/pretreatment/atlas/model-transfer certification still requires evidence. All clinical approval flags remain unsigned.

No frozen dataset, W, anchor, cohort selection, historical QC hold or source evidence was changed.
