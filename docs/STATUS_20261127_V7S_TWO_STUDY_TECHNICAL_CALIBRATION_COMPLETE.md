# V7S completed: technical calibration closed, clinical certification HOLD
2026-11-27. Evidence: user-provided executed V7S stdout; V7R job 13805281 COMPLETED 0:0 00:02:10, final R marker present; both V7P and V7R checksum manifests verified. V7S syntax and completion markers PASS. Original cluster audit files not yet imported byte-for-byte. No model fitting or files modified in V7S. Clinical admissions 0.

All 165 exact matched controls (476 runs, two studies: SRP074349/GSE81089 lung n=143, SRP157215/GSE118435 prostate n=22) calibrated against the unchanged 9465 ARCHS4 atlas. All 499 genes and 500 directional W features retained; independent fixed-W KL, not official solver. CPM primary.

| Branch | N | Pooled median TVD | P95 | Maximum | Dominant agreement |
|---|---:|---:|---:|---:|---:|
| CPM | 165 | 0.02312456 | 0.03644799 | 0.04921636 | 161/165 |
| TMM_CPM | 165 | 0.02543063 | 0.04137143 | 0.06329934 | 160/165 |
| UQ_CPM | 165 | 0.04441662 | 0.07160037 | 0.1010181 | 155/165 |

CPM component median signed differences, recount3 minus ARCHS4 (fraction units):
| Component | Lung | Prostate |
|---|---:|---:|
| NonSE | -0.0018331184 | -0.0064122259 |
| SE1 | 0.0053848287 | 0.0129366296 |
| SE2 | -0.0055418859 | -0.0112593054 |
| SE3 | -0.0008523797 | -0.0013923293 |
| SE4 | -0.0019987776 | -0.0015351271 |
| SE5 | -0.0038121937 | -0.0016975748 |
| SE6 | 0.0101634990 | 0.0034406421 |
| SE7 | -0.0019571308 | -0.0039189015 |
| SE8 | 0.0052267765 | 0.0115636137 |
| SE9 | -0.0046482536 | 0.0002078233 |

Same direction 9/10 CPM, 10/10 TMM, 5/10 UQ; compositional dependence precludes treating components as independent sign tests. Prostate CPM Pearson range 0.9832264–0.9987916, Spearman 0.9649915–0.9988707. Lung Pearson 0.9902106–0.9996816, Spearman 0.9876150–0.9994008. Four CPM discordant labels all become SE8 (lung SE5 one and SE7 two, prostate NonSE one). Prostate discordant archived margin 0.0113216 vs concordant median 0.0360678. Original dominant classes across studies: NonSE 19, SE1 0, SE2 12, SE3 1, SE4 0, SE5 35, SE6 3, SE7 91, SE8 3, SE9 1. Two dominant classes unrepresented.

CPM prostate paired adequacy median delta / Spearman: normalized KL +0.04799218 / 0.9706381, RMSE +0.01846609 / 0.9491813, R2 -0.01540935 / 0.9616036, cosine -0.02228005 / 0.9706381. Original prostate R2 median 0.1092658 vs lung 0.1868612; prostate normalized KL 0.9396567 vs lung 0.4448863. Model fit differences present in original archive; cannot separate cancer and study. UQ prostate normalized KL median paired delta -0.002460946 despite higher marginal recount3 median; paired differences not marginal differences.

Decision: close available matched technical calibration, no rerun; preserve original official ARCHS4 and explicit independent full500 estimator. CPM primary, TMM/UQ sensitivities, continuous fractions primary, no empirical component or nine-gene correction. Strong descriptive transfer evidence, not formal equivalence, not universal transfer, not independent biological validation. Use source-stratified/within-study effects and external validation, avoid unidentifiable source adjustment when confounded with study; future external projection requires QC, clinical status, patient independence, nonoverlap and model adequacy. All matched controls excluded from primary cohort. No clinical release. Next resume clinical certification at SRP057793 (149 source, 148 candidates, 8 held, 140 primary-screen, one source-only). Last exact QC join verified eight reason_codes, reason interpretation pending; see docs/STATUS_20261126_SRP057793_AUTHORITATIVE_QC_SCHEMA.md.
