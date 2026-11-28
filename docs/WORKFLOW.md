# Working procedures and reproducibility
**Current active task: clinical certification; see [current status](CURRENT_STATUS.md).** This is the operational guide. Historical commands in dated reports are not current instructions.

## 1. Cluster synchronization and inventory

```bash
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
ENV=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human
cd "$REPO"
git status --short --untracked-files=no
git pull --ff-only origin main
git rev-parse HEAD
git ls-files --others --exclude-standard -- diagnostics/
git ls-files --others --exclude-standard -- scripts/ tools/ configs/ docs/
git ls-files -- private_data/
```

Inspect actual results. Do not use git add ., git clean, reset --hard, or force-push. If there are local modifications or conflicts, stop and reconcile without discarding evidence.

## 2. Diagnostic protocol

When uncertain about paths, schema, mapping, model behavior or job failure, start with a short, read-only, evidence-driven diagnostic. Inspect actual output before a patch or rerun. Use exact paths, versioned source and strict identity/shape/hash gates. Save the original diagnostic script before execution; capture stdout/stderr and checksums under a new diagnostic ID. If already run, import verified existing source/log bytes without re-execution where possible. Never claim reconstructed code or transcribed output is byte-identical.

Use the existing [record_existing.sh](../tools/record_existing.sh) only for already executed Python script/log pairs that have been reviewed as safe for the hosted repository. It creates an immutable diagnostic record; review its interpretation before committing. Record the decision and provenance in [AUDIT_REGISTRY.tsv](../AUDIT_REGISTRY.tsv). For mixed/private logs preserve the originals locally and archive only approved aggregates, explicitly marked as derived. Keep the repository private per institutional policy.

## 3. Current executable entry points

| Purpose | Entry point | Status |
|---|---|---|
| Current SRP057793 historical QC schema/reason reference | [v8_srp057793_qc_reason_and_clinical_schema_READONLY_v1.py](../tools/v8_srp057793_qc_reason_and_clinical_schema_READONLY_v1.py) | Reconstructed aggregate-only variant; historical result already reviewed; not a new clinical diagnostic |
| V7R configured matched-study preflight | [run_v7r_configured_static_preflight_v1.sh](../scripts/run_v7r_configured_static_preflight_v1.sh) | Verified for prostate; do not repeat |
| V7R configured technical calibration | [slurm_v7r_configured_matched_calibration_v1.sbatch](../scripts/slurm_v7r_configured_matched_calibration_v1.sbatch) | Prostate executed PASS; original lung V7P preserved |
| V7S two-study review | [v7s_two_study_checkpoint_review_READONLY_v1.R](../tools/v7s_two_study_checkpoint_review_READONLY_v1.R) | Executed PASS; no repeat |
| Archive import | [record_existing.sh](../tools/record_existing.sh) | Use only with explicit reviewed input paths |

Next: explicit investigator approval of [FUSCC331 proposed release policy](V8_FUSCC331_PROPOSED_RELEASE_POLICY.md); project-key nonoverlap [executed PASS](STATUS_20261128_V8_FUSCC_PROJECT_KEY_NONOVERLAP_PASS.md). Then immutable private signed cohort manifest, no expression rerun.

For new computation, validate Bash syntax, Python/R parsing, dependencies, source fingerprints and regression tests before SLURM submission. Preserve distinct job logs and output directories; no overwrites. Original frozen W, anchors and Stage09–13 outputs are immutable.

## Dedicated metadata reader environment

The successful 26 November temporary reader used Python 3.12.16 and PyArrow 26.0.1; the environment was subsequently removed. The 28 November primary2055 preflight stopped at ModuleNotFoundError: pyarrow, before creating its execution audit. Do not install packages into the frozen SpatialEcoTyper environment.

Use the retained dedicated prefix `$REPO/private_data/environments/atlas_metadata_reader_v1` on the project filesystem, created from [the pinned environment specification](../configs/atlas_metadata_reader_v1.yml). Use `conda env create --prefix "$META_ENV" --file "$REPO/configs/atlas_metadata_reader_v1.yml"` only if the prefix does not exist; if it exists, verify its exact Python and PyArrow versions instead of overwriting or reinstalling. Keep it for subsequent metadata audits. Store explicit package export locally after successful provisioning; the YAML pins the two direct requirements, not every transitive package. No environment binary files or package caches should be committed. Use its bin/python for Parquet diagnostics; continue using the original SpatialEcoTyper environment for frozen model computations.

Dedicated environment provisioned and verified on 28 November 2026; Python 3.12.16 / PyArrow 26.0.1, explicit package export SHA256 abc130c76768c74dc6d63928b53b2500d9e2b7ad643924d6f733821531b733cf. Primary2055 scientific script executed PASS with all four checksum entries verified. Retain environment; no repeated provisioning.

## 4. Outstanding archival task

GitHub contains documented results and implementations; raw local audit bytes for many V7 steps have not yet been imported. The registry marks this distinction. Start from the cluster inventory rather than rerunning. For each candidate archive:
- Confirm original path, existence, execution status, checksum manifest and safety of every file.
- Select only explicit aggregate-only TSVs and reviewed aggregate logs, no specimen/run/patient linked records or private checkpoints. A manifest referencing local private inputs may be retained locally; do not mistake it for a portable fully verifiable remote manifest.
- Preserve original bytes where safe, otherwise label derived/redacted summaries and maintain local original hashes.
- Stage explicit file paths, inspect staged diff and check private data is not tracked; commit/push intentionally.
- Update the registry with actual archival provenance; no retrospective claim of executed bytes.

The original documentation policy remains [PROVENANCE_POLICY.md](PROVENANCE_POLICY.md). Current human-readable status lives in CURRENT_STATUS.md; update it and the corresponding scientific report for material changes. Dated evidence records are history, not additional active READMEs.

## 5. Repository layout decisions

Executable and evidence paths are stable because cluster wrappers and checksum manifests reference them directly. The documentation-only reorganization consolidates the active narrative, keeps original dated documents byte-identical, and indexes them by topic. No deletion, relocation, history rewriting, or scientific change. A physical archive migration would require a separate dependency/path audit and synchronization of the local untracked evidence; do not do that casually.
