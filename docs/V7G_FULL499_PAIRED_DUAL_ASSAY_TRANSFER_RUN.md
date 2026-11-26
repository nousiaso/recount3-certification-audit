# V7G: full 499-gene paired ARCHS4–recount3 SpatialEcoTyper assay transfer

## Why this version supersedes V7C

V7F established the **original full SRP074349 G026 gene-sums gzip** contains all 454 exact SRR runs associated with 143 previously verified ARCHS4 GSMs. Its matrix header is the third decompressed line (0-based line 2), with 668 run columns. The 144-run Stage2 selected HDF5 is a **different, disjoint** subset and must not be used.

V7C v1 stopped before normalization, because it assumed one run per specimen and used the wrong selected HDF5. Its output and log must remain intact.

## V7G computation (separate, new output)

1. Verify SHA256s for the exact raw G026 gene-sums source, Stage5 run bridge, 63,856-row gene map and original frozen Stage09 Python code.
2. Verify all 143 GSM-to-zero-based ARCHS4 H5 column indices and all 454 uniquely grouped bridge run identifiers.
3. Discover the actual raw G026 matrix header by 668 exact run IDs; require all 454 matched runs. Stream all 63,856 gene rows, verify same ordered Ensembl identities as frozen Stage5 mapping, and SUM original **gene-body base-pair coverage** across technical SRRs **within the same original SRX experiment** (84×1, 4×2, 22×3, 1×8 and 32×9 runs). This is not read-count conversion.
4. Separately extract the **same 143 GSM samples** directly from original ARCHS4 HDF5 `data/expression`, preserving sample ordering and verifying original GEO column identities. These are frozen source pseudocount values, not input recalculated from the recount3 measurement.
5. Reuse the authoritative unchanged Stage09 `build_gene_maps`, `collapse_batch`, `quantile75_fraction`, `select_tmm_reference`, `tmm_factor`, and `geometric_mean_positive` implementations. For **each** assay, calculate CPM/TMM_CPM/UQ_CPM with denominators across ALL valid collapsed symbol groups. The 499 W genes are not the denominator. TMM and UQ factors are estimated **within these 143 matched specimens separately for each assay**, rather than reusing historical 9,465-sample TMM/UQ factors; note when interpreting comparison with archived Stage10 H.
6. Apply the original, unchanged pan-cancer Stage10 mean/SD scaling to each assay and branch. Use unchanged official 500×10 W with `SpatialEcoTyper::DeconvoluteSE(scale=FALSE, sum2one=FALSE, ncores=1)`. No W refit, no anchor recomputation, no other deconvolution solver.
7. Compare same-specimen raw H, compositional fractions, dominant component, total variation distance, reconstructed directional fit, frozen-Z gene bias, and within-assay normalization branch robustness.
8. Keep all **499 model genes** in the primary analysis. In an explicit sensitivity for the nine known source-annotation differences, set the nine gene frozen-Z inputs to zero **in both assays**, and re-project with W unchanged. This is a transparent feature-neutralization ablation (not rewriting the raw gene map or dropping weights).
9. Optionally compare projected ARCHS4 primary CPM against archived Stage10 raw H where same GSMs are present in frozen canonical atlas. This is a reference replay, not a requirement that all 143 be in the 9,465 frozen atlas.
10. Leave every patient, study, model and clinical-release status untouched. All outputs are descriptive technical calibration, NOT independent biological validation.

## Execution on compute node

After pulling GitHub `main`, submit:

```bash
sbatch scripts/slurm_v7g_matched_143_dual_assay_full499.sbatch
```

The SLURM job uses `normal_prio`, `no_limit_account`, 8 CPUs and 96 GB RAM. It performs Python AST/import checks and R syntax/package checks before reading data; all source files are read-only. **Do not rerun with the same output folder**. If it fails, preserve logs and diagnose the exact error rather than overwriting.

Private sample-linked outputs:
`private_data/v7g_143_matched_dual_assay_full499_v1/`

Versioned aggregate diagnostics:
`diagnostics/20261127_v7g_143_matched_dual_assay_full499_v1/`

Expected aggregate files if both stages finish:
- `aggregate_transfer_metrics.tsv`
- `aggregate_normalization_sensitivity.tsv`
- `aggregate_nine_gene_sensitivity.tsv`
- `input_provenance_aggregate.json`
- `01_python_preparation_stdout.txt`
- `02_R_projection_stdout.txt`
- `SHA256SUMS.txt`

Please share only aggregate metrics and stdout for interpretation. **Do not Git-push specimen/run-specific TSVs, input matrices or RDS objects.**

## Interpretation and limitations

- In-sample paired study technical controls (143 GSMs; 454 SRRs), not 143 independent patients.
- TMM/UQ size-factor reference cohort differs from the archived 9,465 frozen Stage09 cohorts; assess CPM first.
- Correlation of same-gene expression across specimens, even near 1, does **not** establish within-specimen multi-gene similarity, calibrated component H, ecotype fractions, or reconstruction performance.
- Avoid declaring equivalence from one average number: assess 90th/95th percentiles, dominant-class switches, per-component effects, sequencing depth dependence and the nine-gene sensitivity.
- This phase does not certify SRP157974/FUSCC patient identity, pretreatment status, direct ARCHS4 independence or clinical eligibility. No clinical release is permitted on the basis of V7G alone.

Historical frozen model assets and all V7 stages remain untouched.
