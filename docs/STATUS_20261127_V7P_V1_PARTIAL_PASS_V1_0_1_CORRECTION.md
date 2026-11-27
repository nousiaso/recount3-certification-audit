# V7P v1 partial PASS / HOLD; v1.0.1 correction

2026-11-27. Evidence: user-supplied actual cluster stdout and synthetic error reproduction. V7O job 13805107 COMPLETED 0:0, elapsed 00:01:08, stdout SHA256 77aac6b3d763a29f3e5ffb9707ec3b9a0a06e9312f9d042e92512ebaac3d5712. V7P static preflight PASS at 32a0fad05d7f083b63dead013b1fdb797070c7e7; five V7J input files verified.

V7P v1 all projection branches completed then failed in normalization_sensitivity. Failed stdout SHA256 f1c7bf2bc77a5a1c0142217316c3fea947e214855222dc6bcc0d1797c4246f97. Private output directory empty, audit contains stdout only. No recoverable checkpoint or aggregate files. Preserve v1 evidence.

| Branch | N | Median TVD | P95 TVD | Dominant concordance | Archive replay max TVD |
|---|---:|---:|---:|---:|---:|
| CPM | 143 | 0.02182098 | 0.03106810 | 0.979021 | 0.000151514166957 |
| TMM_CPM | 143 | 0.02394995 | 0.03682443 | 0.972028 | 0.000011302194964 |
| UQ_CPM | 143 | 0.04492710 | 0.06820524 | 0.958042 | 0.000102842724753 |

All 858 primary fits convergence code 0; original CPM adequacy replay max abs 7.21644966006e-16; 499 genes harmonized and 500 features retained. No nine-gene or order sentinel results. Technical agreement not biological validation or clinical equivalence.

Root cause: stored arch vs requested ARCHS4; incorrect lookup returns NULL and original helper fails mask check. Synthetic diagnostic reproduced exact error and corrected lookup passed. v1.0.1 resolves public labels to internal keys explicitly, adds dimension checks and private checkpoint before post-processing, and preflight regression with the actual helper for both assay labels and both alternatives including masked UQ. No automatic checkpoint reuse. Original optimizer, parameters, W, anchors, source helpers, inputs and thresholds unchanged. New versioned output paths. No raw extraction.

Next: fast-forward pull, run scripts/run_v7p_reused_input_full500_static_preflight_v1_0_1_READONLY.sh, review PASS markers, then submit scripts/slurm_v7p_full500_ARCHS4_vs_recount3_calibration_v1_0_1.sbatch. Not executed on cluster. Clinical admissions 0. No patient-linked data committed.
