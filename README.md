# ARCHS4 × recount3 | clinical certification

Reproducible evidence and code for extending the frozen 9,465-sample ARCHS4 SpatialEcoTyper atlas with independently certified primary-tumour RNA-seq from recount3.

**Current status (27 November 2026): technical calibration complete within the available matched bridge; clinical certification UNSIGNED. Zero clinical admissions. No external-cohort biological inference authorized.**

## Start here

| Document | Purpose |
|---|---|
| [Current project status](docs/CURRENT_STATUS.md) | Cohort counts, completed work, unresolved gates, next clinical task |
| [Technical calibration](docs/TECHNICAL_CALIBRATION.md) | Frozen model, 165 matched specimens, results, decisions and limitations |
| [Working procedures](docs/WORKFLOW.md) | Cluster synchronization, safe diagnostic and archival workflow, executable entry points |
| [Historical index](docs/README.md) | Dated reports, original implementations, prior failures and audit evidence |
| [Audit registry](AUDIT_REGISTRY.tsv) | Machine-readable evidence provenance and rerun policy |
| [Provenance policy](docs/PROVENANCE_POLICY.md) | Definitions, admission rules and data governance |

## Status at a glance

| Workstream | State |
|---|---|
| Frozen ARCHS4 normalization/projection replay | PASS |
| Matched ARCHS4–recount3 technical calibration (165 specimens, 2 studies) | CLOSED |
| Clinical review of 77 primary-screen studies (2,055 candidate libraries) | UNSIGNED |
| Patient/specimen independence and final cohort lock | HOLD |
| Combined biological inference | HOLD |

Next: SRP057793, 140 primary-screen candidates, eight study-relative QC holds, one source-only experiment. Resolve actual clinical evidence and independence, not another projection rerun.

## Repository map

| Path | Role |
|---|---|
| `docs/` | Current guides and indexed historical reports |
| `configs/` | Explicit study contracts |
| `scripts/` | Versioned execution wrappers and frozen originals |
| `tools/` | Versioned diagnostics and analysis implementations |
| `diagnostics/`, `evidence/` | Immutable execution evidence and historical artifacts |
| `decisions/` | Evidence-level decisions, not final admissions |

Data and specimen-linked outputs stay on the cluster. Do not change frozen W, anchors, Stage09–13 or historical execution records. Do not stage private data or blindly rerun superseded scripts.

For a new conversation, use the [handoff](docs/NEW_CHAT_HANDOFF.md). The long historical README is preserved at [the previous revision](https://github.com/nousiaso/recount3-certification-audit/blob/77680d130af21a6f81827c14dbf90f1107dbd585/README.md); it is not an active runbook.
