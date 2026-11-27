# V7R prostate result generation PASS; V7S two-study review next

2026-11-27. User supplied successful static preflight (all 11 entries), private checkpoint 586K, all four final aggregate tables and all 12 final audit checksum entries OK. Actual SLURM accounting and final stdout markers not included in this excerpt; formal execution record still to verify. Original output bytes not imported by this documentation commit.

SRP157215/GSE118435: 22 prostate specimens, 22 runs from 41 raw source runs, 499 genes and 500 directional features; original frozen 9465 ARCHS4 reference and independent full500 KL. All branches jointly evaluable n=22.

| Branch | median TVD | P90 | P95 | max | dominant concordance | median Z MAE | archive replay max TVD |
|---|---:|---:|---:|---:|---:|---:|---:|
| CPM | 0.0335756971920466 | 0.0394960538835071 | 0.0413593644021742 | 0.0492163593193023 | 0.954545454545455 | 0.106784777706212 | 0.0000927824861273879 |
| TMM_CPM | 0.0383658421798043 | 0.0530921872214561 | 0.0546824454059611 | 0.0632993419918979 | 0.954545454545455 | 0.154109611213169 | 0.0000203081163748309 |
| UQ_CPM | 0.0402346290250821 | 0.0724454405011703 | 0.0959688089252167 | 0.101018088062328 | 0.818181818181818 | 0.128735512087332 | 0.0000182837778555501 |

Normalized KL medians ARCHS4/recount3 CPM 0.939656702931675/0.99258394990939, TMM 0.876422029316781/1.03386115656186, UQ 0.997857774410337/1.03401925234331. RMSE medians CPM 0.537440063994395/0.558604261250872, TMM 0.55040381579955/0.543389820477467, UQ 0.487038811345076/0.516371510700758. Higher absolute normalized KL in prostate even in original archive, so not a recount3-only artifact. Cannot isolate cancer from study; paired delta not equal to difference of medians. All archived replay errors below original thresholds.

Within-assay CPM vs TMM median TVD ARCHS4 0.0212253919465015, recount3 0.0215991893224576; CPM vs UQ 0.0302835650210746 and 0.0298972980990778. UQ P95 0.155195651399755 and 0.0833386200027277. Nine-gene substitution original median 0.0335756971920466, swapped median 0.0326323768257548, median within-pair delta +0.000585486562584605, P95 absolute delta 0.00695952007642034; no consistent evidence of correction benefit. Order sentinel 0/0/0/1.

CPM combined 161/165 dominant agreement, 97.58%, descriptive only and dominated by lung n=143; 476 runs in total, all bridge records calibration-only not primary. No formal equivalence or independent biological validation. CPM primary remains prespecified; TMM/UQ sensitivities. V7S read-only script tools/v7s_two_study_checkpoint_review_READONLY_v1.R uses both private checkpoints to compare component-specific biases, dominant margins, reconstruction metrics and actual pooled TVDs (not averages of medians). No projection rerun or empirical bias correction. Clinical admissions 0. Return to clinical certification after bounded review, archive aggregate-only local artifacts separately.
