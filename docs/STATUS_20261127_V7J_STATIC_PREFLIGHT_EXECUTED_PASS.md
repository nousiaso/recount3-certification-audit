# V7J static original-code fidelity preflight: EXECUTED PASS

**Date:** 27 November 2026  
**Source:** user-supplied actual supekgate stdout from `scripts/run_v7j_static_frozen_transfer_preflight_READONLY.sh`  
**Status:** STATIC SOURCE/DEPENDENCY PREFLIGHT PASS; full 454-SRR G026 calibration **NOT YET EXECUTED**; clinical admissions **0**.

## Cluster-confirmed exact input fingerprints

| Input | SHA256 |
|---|---|
| Original frozen Stage09 | `46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500` |
| Original frozen Stage10 | `16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c` |
| Stage5 exact bridge | `c7d6d7c754e3ca77d62424f26627a40ceaec10fac98a0741c945773f3d3b1b67` |
| Stage5 63,856-row gene map | `68be266dd6bd30bfa79a2bd8ac79f0e49c6f5d0f19d2d083e2ca57850ef17391` |
| Original 668-SRR G026 source | `b056331eb16d6ec612409fd37a266abfe10c4c15329b486dbb9df213612797c5` |
| V7J recount3-only Python preparer | `7539bf696b8de14eafe7d6b709c6d8e49d73cf0ece5ac936ea0dffaff211c0f4` |
| V7J original-baseline R projection | `ff0d70ae72386bd8294acc099a9d2fbc16d0322442dde284117dfd085e5f6d0c` |

All eight original frozen Stage10 helper function bodies were confirmed byte-exact in V7J:
`score_fractions`, `clr_matrix`, `dominant_component`, `harmonize_expression_to_W`,
`make_directional_bundle`, `project_raw_h`, `calculate_reconstruction_adequacy`,
`normalization_sensitivity`.

Executed PASS markers:
- `V7J_PREPARATION_PYTHON_AST_PASS`
- `V7J_PREPARATION_ARCHS4_UNCHANGED_REFERENCE_PASS`
- `STAGE10_AUTHENTIC_FUNCTIONS_IDENTICAL 8`
- `V7J_STATIC_FIDELITY_CONTRACT_PASS`
- `V7J_R_SOURCE_SYNTAX_PASS`
- `V7J_R_PACKAGE_PREFLIGHT_PASS`
- `V7J_PYTHON_PACKAGE_PREFLIGHT_PASS`
- `V7J_READONLY_STATIC_PREFLIGHT_CAPTURED_PENDING_REVIEW`

Original cluster preflight hash-check lines for the three source scripts and four preflight output files were all `OK`.

## What this proves and what it does not

The V7J full executable is syntactically valid in the specified cluster environments and contains faithful copies of eight authentic Stage10 helpers. Code/data checksums match and ARCHS4 remains frozen. The separate V7H2 full-original 9,465-sample numerical replay passed all three normalization branches at near floating-point precision.

This **does not** demonstrate that the 454-SRR recount3 source has been numerically extracted/aggregated by V7J, that the 143-query batch has invariant frozen ARCHS4 H, that all 500 directional features will be retained, or that recount3 ecotype fractions agree with ARCHS4. Those are mandatory **runtime and scientific gates** in V7J.

Stage09 QC field `upper_quartile_fraction` and Stage10 original constants `CLR_FLOOR=1e-8`, `RECONSTRUCTION_EPS=1e-12`, `nsample.per.run=10000L`, `ncores=1L` have also been confirmed in the original archived scripts.

## Next action

Submit the already-versioned, preflight-gated compute-node job after `git pull --ff-only origin main`:

```bash
sbatch scripts/slurm_v7j_frozen_ARCHS4_vs_recount3_calibration.sbatch
```

Its runner rechecks the static source/input contract before extraction; preserves sample/run-linked outputs strictly under ignored `private_data/`; refuses an existing result/audit directory; and captures aggregate stdout, branch-level calibration results, normalization sensitivity, nine-gene sensitivity and determinism checks under `diagnostics/20261127_v7j_143_frozen_archs4_vs_recount3_full499_v1/`.

**Do not commit, push or share raw patient/run linked matrices.** Capture aggregate stdout and TSVs for scientific interpretation, and preserve any failure without overwriting or running the same output version a second time.

External transfer equivalence is **pending**; Stage11/Stage12 full revalidation and all clinical admission gates remain separate. No clinical admissions.
