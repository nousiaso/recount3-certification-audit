# V7P v1.0.1 wrapper failure; shellfix1

2026-11-27. User-supplied SLURM accounting: job 13805257 FAILED 2:0, elapsed 00:00:01. All five checksum entries in the v1.0.1 static preflight manifest verified in SLURM stdout. Stderr: unexpected EOF while looking for matching single quote. R did not start; both v1.0.1 output directories absent. Original V7P v1 job 13805232 FAILED 1:0 after all three projection branches, as documented separately.

At commit a97bddf3fd62acbe6e788f66df70aee9bc13af80, the v1.0.1 execution wrapper contains a truncated grep marker line and duplicate commands. This was a wrapper-generation error, not a model or R failure. Preserve old wrappers and SLURM logs.

Shellfix1 creates two new files, reconstructing the runner from the original v1 runner with v1.0.1 paths, exact regression-marker grep, explicit no-existing-output guard and private checkpoint existence check. New SLURM wrapper runs bash -n on runner before invoking it, with distinct SLURM logs. The R script, v1.0.1 output paths (observed absent), preflight, and frozen input assets are unchanged. No automatic resumption or raw extraction.

Before submitting, run bash -n on both new shell scripts and verify the existing v1.0.1 preflight manifest and regression marker. Submit scripts/slurm_v7p_full500_ARCHS4_vs_recount3_calibration_v1_0_1_shellfix1.sbatch only after those checks and no active duplicate job. No scientific claims beyond prior partial PASS. Clinical admissions 0. Shellfix1 not yet executed on cluster.
