# V7N executed PASS; four lost original W features impose nonzero KL penalty

Date: 27 November 2026. Evidence: user's actual V7N stdout, SHA-verified on supekgate.

## Original frozen 143-GSM ARCHS4 evidence

- Original model W: 500 directional rows × 10 components; no W missing/nonfinite entries.
- Original 143 specimens from frozen 9,465-sample atlas, archived H for all 143 × 10 components.
- The 499 base genes produce 998 directional rows. Of these, **26 are nonvariable** across the selected 143 samples, but precisely **four belong to frozen W**; 496 W rows remain.
- The four removed W features are `F12__pos`, `JUND__pos`, `PSMA3__neg`, `SURF4__neg`. Each has exact measured directional expression **0 in all 143** original archived specimens. Each retains positive original W weights across the 10 components (the V7N program tested strictly nonzero elements; do not interpret tiny weights as biologically meaningful on their own).

### KL reconstruction penalty at original archived H

| W feature | Median predicted value at archived H | Median KL penalty |
|---|---:|---:|
| F12__pos | 0.480912342329 | 0.480912342329 |
| JUND__pos | 0.411350838614 | 0.411350838614 |
| PSMA3__neg | 0.454448740125 | 0.454448740125 |
| SURF4__neg | 0.538323166113 | 0.538323166113 |

Since measured V is zero, the KL term for each feature is **(WH)**. Removing a zero-valued row with nonzero W weights still removes the **positive penalty on inferred H**.

- Total KL across all 500 features at archived H: median **137.629504592**, p95 **164.032712815**.
- Four removed features' KL: median **1.89607796836**, p95 **2.553614655**, max **2.88179599801**.
- Their per-specimen contribution to the 500-feature KL total: median **1.39694156012%**, p95 **1.82165657945%**, max **2.46600315342%**.

These percentages are KL-objective fractions **at the archived H**, not effect sizes for H/fraction scores. They cannot be used to attribute 1.4% or any other fixed fraction of the V7J 0.009590588 maximum TVD to four feature omissions.

## Decision: V7O controlled fixed-W KL diagnostic

We now have concrete evidence that a reduced 496-feature projection solves **a different KL objective** than the original 500-feature projection; random initial H also differs and original solver has connectivity-based stopping. Quantifying their separate contributions requires optimization under controlled conditions.

V7O is a new immutable read-only-to-source numerical experiment **on original ARCHS4 only**, not a recount3 calibration:

- `tools/v7o_143_frozen_ARCHS4_controlled_500_vs_496_KL_READONLY_v1.R`
- `scripts/run_v7o_143_frozen_ARCHS4_500_vs_496_KL_READONLY.sh`
- `scripts/slurm_v7o_143_frozen_ARCHS4_500_vs_496_KL.sbatch`

It extracts only the original frozen 499×143 Z and archived 143×10 H. It reconstructs the exact 500 versus 496 feature universes; separately executes original SpatialEcoTyper v1.0.3 on the 143 original samples to confirm historical V7J discrepancy. It then solves **the fixed-W, nonnegative KL objective sample-by-sample** by a *separate*, explicitly labeled deterministic `stats::optim(L-BFGS-B)` diagnostic, with BOTH the original archived-H starting vector and a uniform starting vector for each specimen and feature universe.

The fixed-W KL problem is convex in H for nonnegative W, allowing separate, data-independent optimization. V7O checks solver convergence codes, relative objective stability across starts and H-derived fraction consistency across starts. It reports:

1. Original official 143-versus-9465 H / TVD difference.
2. Controlled 500-versus-496 KL-optimized fraction TVD, isolating feature-universe changes without altering sample set or initial H.
3. Original archived H versus controlled 500-feature KL optimum.
4. Original official 143-context H versus controlled 496-feature KL optimum.

**This alternate diagnostic optimizer is NOT the original SpatialEcoTyper implementation and must not be used to silently replace its archived results.** Discrepancies are nonadditive. If convergence checks fail, numerical conclusions remain HOLD pending a focused investigation.

The SLURM job is reserved for a compute node, preserves V7J outputs and source freezes, and writes only aggregate stdout and SHA256s in `diagnostics/20261127_v7o_143_frozen_ARCHS4_500_vs_496_KL_v1/`. No raw recount3 counts or patient-linked tables are written.

Submit:

```bash
sbatch scripts/slurm_v7o_143_frozen_ARCHS4_500_vs_496_KL.sbatch
```

Review `01_500_vs_496_KL_stdout.txt`. Until the controlled results are available, external assay transfer and clinical admission are HOLD, 0 new admitted specimens.
