# V7H original ARCHS4 replay: v1 stop and exact v2 correction

Date: 27 November 2026.

## Observed v1 cluster output (user-supplied)

- Stage09 frozen status: PASS.
- Stage10 frozen status: PASS.
- Replayed reference samples: 9,465.
- Model genes: 499; directional W features: 500.
- Installed SpatialEcoTyper 1.0.3.
- Stage09 CPM normalized expression arithmetic: maximum absolute difference `4.61180316051468e-05`; maximum relative difference `4.99236050696712e-10`. This passed the prespecified `1e-9` relative numeric check.
- Next checkpoint stopped with `Original W-gene ledger order mismatch`. **No frozen Z nor H replay was completed**, and TMM/UQ were not reached.

## Root cause identified in V7H v1 code

The v1 R statement was:

```r
map <- mapping[as.character(mapping[["branch"]]) == branch]
```

Since `mapping` is a `data.table` with a column named `branch`, the bare RHS `branch` is looked up in the table's column scope, so all three branches can be retained. This is **our diagnostic code bug**, not evidence that the original model's gene order is wrong.

## Corrected immutable V7H2

The corrected code computes the row indices **outside** `data.table[...]`:

```r
branch_indices <- which(as.character(mapping[["branch"]]) == branch)
map <- mapping[branch_indices]
setorder(map, W_gene_index)
```

It also verifies the 3 frozen branch partitions each have 499 rows; then verifies unique sequential W gene indices and exact original W base-gene order for the selected branch. All other replay calculations, frozen sources and numeric thresholds remain as v1.

Versioned files:
- `tools/v7h2_replay_frozen_ARCHS4_9465_stages09_10_READONLY_v2.R`
- `scripts/run_v7h2_ARCHS4_9465_original_numeric_replay.sh`
- `scripts/slurm_v7h2_ARCHS4_9465_original_numeric_replay.sbatch`

**Do not rerun v1 or erase the failed v1 audit:**
`diagnostics/20261127_v7h_ARCHS4_9465_original_numeric_replay_v1/`.

To run v2 once on a compute node:
```bash
sbatch scripts/slurm_v7h2_ARCHS4_9465_original_numeric_replay.sbatch
```

The v2 logs go to a separate path:
`diagnostics/20261127_v7h2_ARCHS4_9465_original_numeric_replay_v2/`.

No raw data, frozen source, original Stage09-12 code or model W was modified. **V7G remains on HOLD** until v2 numeric results (all three normalization branches, frozen Z, full 9,465-sample H and fraction replay) are reviewed.
