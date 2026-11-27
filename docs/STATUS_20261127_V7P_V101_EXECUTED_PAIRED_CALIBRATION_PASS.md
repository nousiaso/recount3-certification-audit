# V7P v1.0.1 shellfix1 executed: computational PASS; technical transfer review

2026-11-27. Evidence: user-supplied actual cluster terminal output including all four aggregate TSVs, SLURM accounting, and final SHA256 verification. Job 13805260 v7p_143_assay500_v101_sf1 COMPLETED 0:0, elapsed 00:01:32; empty stderr as displayed. Both previous failed jobs retained (13805232 and 13805257). The original local stdout/TSVs/checksum manifest have not yet been imported byte-for-byte into GitHub by this documentation commit.

## Exact matched technical calibration

143 identical specimens, 454 runs aggregated into 143 SRX experiments; all 499 model genes and all 500 W directional rows retained. Original frozen 9465 ARCHS4 normalized Z/H/fractions were not re-estimated. Independent full500 nonnegative KL estimator, not official SpatialEcoTyper solver, projected recount3. Three branch replay fidelity gates passed (max TVD below 0.0005); original CPM adequacy max replay absolute error 7.21644966006e-16. All 1037 optimization fits (858 branch fits + 143 nine-gene swap + 36 sentinel fits) returned convergence code 0. Checkpoint saved (3.5 MB, private). Final four aggregate tables and all seven checksum entries verified.

| Branch | N | Median TVD | P90 | P95 | Max | Dominant concordance | Median Z MAE | Archive replay max TVD |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CPM | 143 | 0.0218209751485707 | 0.0298222534293889 | 0.031068097208178 | 0.0433952129428884 | 0.979020979020979 | 0.10965707960928 | 0.000151514166956803 |
| TMM_CPM | 143 | 0.0239499530956963 | 0.032429951053127 | 0.0368244321044379 | 0.0476242414002724 | 0.972027972027972 | 0.134860824403943 | 1.1302194963991e-05 |
| UQ_CPM | 143 | 0.044927100065551 | 0.0616406546753595 | 0.0682052364126112 | 0.0820273982007886 | 0.958041958041958 | 0.266124928455061 | 0.000102842724753109 |

Median normalized KL (ARCHS4/recount3): CPM 0.444886262047484/0.479580011879346, TMM 0.439901327959142/0.491517478433806, UQ 0.397001053824083/0.517342368912432. Median RMSE: CPM 0.552202333210651/0.555746474153234; TMM 0.561474755207551/0.560729423503162; UQ 0.575773738627419/0.539018926292179. Higher normalized KL in recount3 for all branches; similar CPM/TMM RMSE and lower UQ RMSE do not establish equal reconstruction adequacy. Distinct metrics measure distinct aspects.

Within-assay normalization median TVD / P95 / dominant agreement: ARCHS4 CPM–TMM 0.0167822086546303 / 0.0361327278762504 / 0.972027972027972; ARCHS4 CPM–UQ 0.0365708496957587 / 0.0746559147435043 / 0.958041958041958; recount3 CPM–TMM 0.0153319195881858 / 0.0245378071962468 / 0.965034965034965; recount3 CPM–UQ 0.0206909132332691 / 0.047918706259405 / 0.951048951048951.

Nine-gene matched ARCHS4 Z substitution: original median TVD 0.0218209751485707, swapped median 0.0227029805053621, median paired delta +0.00113138324769963, P95 absolute change 0.00273398999022391. No support for this swap as a correction; not a causal decomposition or proof of no gene annotation effects. Repeat max fraction difference 0; reversed-order median/max TVD 0; dominant concordance 1 in 12-specimen CPM sentinel.

## Decision

Computational FULL PASS, strong descriptive same-specimen technical agreement, primary CPM retained as prespecified, TMM/UQ sensitivity. No equivalence margin for cross-assay agreement was prespecified; do not claim formal equivalence, universality, independent biological validation, or clinical release. The scores are model fractions not physical cell fractions. Next bounded review uses existing private checkpoint only: component-specific paired bias/agreement, dominant-class margins and confusion, paired reconstruction adequacy distributions, source study/cancer coverage (exact metadata joins, no assumptions). Then return to clinical/independence certification and scalable release workflow; no repeated extraction or model refitting. Clinical admissions 0. Archive verified aggregate-only original outputs from cluster separately after reviewing git diff; private checkpoint and linked data stay local.
