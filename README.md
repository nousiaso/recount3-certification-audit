# ARCHS4 × recount3 — clinical certification & reproducible audit log

**Status: UNSIGNED. Do not run SpatialEcoTyper biological inference from the recount3 cohort yet.**

This is a **small provenance repository**, not a copy of the cluster datasets or frozen ARCHS4 release.
It tracks the exact diagnostic scripts, outputs, checksums, interpretation and outstanding gates.
It was initialized from user-supplied November 2026 evidence; historical scripts not physically
provided have **not** been invented or declared archived.

## Latest V6 checkpoint (2026-11-26)

- `docs/STATUS_20261126_SRP185703_MEMBERSHIP_PREFLIGHT.md`: SRP185703 75 Stage41 records; primary-screen check pending.
- `docs/STATUS_20261126_V6_NEXT_STUDY_QUEUE.md`: 147-study unresolved queue; priority ranking pending.
- `docs/STATUS_20261126_DRP005276_CLINICAL_HOLD.md`: DRP005276 treatment/assay hold; 83 retained, zero clinical admissions.
- `docs/STATUS_20261126_DRP005276_CLINICAL_SOURCE_CONTENT.md`: V2/V3 treatment review text found; source interpretation pending.
- `docs/STATUS_20261126_DRP005276_NATIVE_METADATA_MANIFEST.md`: native metadata and historical clinical ZIP inventory.
- `docs/STATUS_20261126_DRP005276_ORIGINAL_METADATA_LOCATOR.md`: historical metadata files located; contents not yet assessed.
- `docs/STATUS_20261126_DRP005276_CLINICAL_EVIDENCE_TRIAGE.md`: Stage41 clinical pointers absent; treatment and patient independence unresolved.
- `docs/STATUS_20261126_DRP005276_ARCHS4_EXACT_ACCESSION_SCREEN.md`: 9,465/9,465 atlas accession mapping; DRP005276 zero exact overlaps in 83 retained; patient independence UNSIGNED.
- Membership, GEO coverage and exact accession-join diagnostics archived with original stdout checksums.
- `docs/STATUS_20261126_ARCHS4_IDENTITY_INPUT_SCHEMA.md`: source inputs found; Stage41 103 versus QC retained 83; Parquet reader unavailable.
- `diagnostics/20261126_v6_archs4_identity_input_schema_v1/`: executed source input availability/schema audit and hashes.
- `docs/STATUS_20261126_ARCHS4_IDENTITY_LOG_SEARCH.md`: 759-file historical log search; no original execution marker recovered.
- `diagnostics/20261126_v6_archs4_identity_log_input_probe_v1/`: executed bounded search and hashes.
- `docs/STATUS_20261126_ARCHS4_IDENTITY_EXECUTION_LOCATOR.md`: original identity script verified; same-GEO metadata is not independent evidence; original stdout pending.
- `diagnostics/20261126_v6_archs4_identity_execution_locator_v1/`: executed static source locator and checksums.
- `docs/STATUS_20261126_DRP005276_IDENTITY_ZIP_SCHEMA.md`: original ZIP contains code, README and checksums only; identity results must be recovered.
- `diagnostics/20261126_v6_drp005276_identity_archive_schema_v1/`: executed ZIP and old-contract schema diagnostic with checksums.
- `docs/STATUS_20261126_DRP005276_ATLAS_GRAPH_DISCOVERY.md`: decoded DRP005276 V5 overlap flags and located historical identity archive; independence unsigned.
- `diagnostics/20261126_v6_drp005276_atlas_graph_discovery_v1/`: executed file-discovery diagnostic, stdout, hashes and interpretation.
- `docs/STATUS_20261126_DRP005276_ATLAS_IDENTITY.md`: no direct accession collisions; zero GSM coverage, atlas independence unsigned.
- `diagnostics/20261126_v6_drp005276_atlas_identity_v1/`: exact executed identity screen, stdout, checksums and interpretation.
- `docs/STATUS_20261126_DRP005276_S1_COMPLETE_HEADERS.md`: all 18 PLOS S1 headings resolved; patient mapping confirmed, treatment timing unsigned.
- `diagnostics/20261126_v6_drp005276_s1_remaining_headers_v1/`: executed header diagnostic, stdout, hashes and interpretation.
- `docs/STATUS_20261126_DRP005276_S1_HEADER_TREATMENT.md`: S1 patient-ID header link PASS; treatment timing unsigned, 12 unclassified fields.
- `diagnostics/20261126_v6_drp005276_s1_header_treatment_v1/`: executed diagnostic, stdout, interpretation and checksums.
- `docs/STATUS_20261126_DRP005276_S1_PATIENT_AXIS.md`: 103 unique S1 isolate rows and 59/44 full-cohort stage mapping; clinical unsigned.
- `diagnostics/20261125_v6_drp005276_s1_patient_axis_v1/`: executed script, stdout, checksums and interpretation.

## Latest V6 checkpoint (2026-11-25)

- `docs/STATUS_20261125_DRP005276_PUBLIC_CROSSWALK.md`: public 103/103 DRA accession crosswalk, 103 S1 isolate matches; clinical unsigned.
- `diagnostics/20261125_v6_drp005276_public_crosswalk_v1/`: executed diagnostic, stdout, interpretation and checksums.
- `docs/STATUS_20261125_DRP005276_STUDY_CONTEXT.md`: exact publication identified; study-level primary-bladder cohort support, patient/treatment unsigned.
- `diagnostics/20261125_v6_drp005276_study_context_v1/`: executed study-context diagnostic, stdout and SHA256 evidence.
- `docs/STATUS_20261125_DRP005276_RUNHOLD_QC.md`: original 20/20 QC exclusions reconciled; clinical gates unsigned.
- `diagnostics/20261125_v6_drp005276_runhold_qc_v1/`: executed run-hold diagnostic, stdout, checksums and interpretation.
- `docs/STATUS_20261125_DRP005276_FULL103.md`: 103 source experiments = 83 screened + 20 historical QC holds; clinical gates unsigned.
- `diagnostics/20261125_v6_drp005276_full103_source_v1/`: executed full-study census and hashes.
- `docs/STATUS_20261125_DRP005276.md`: 83 exact SRA joins; patient, primary-site and treatment gates unsigned.
- `diagnostics/20261125_v6_drp005276_source_identity_v1/`: executed diagnostic, stdout and SHA256 evidence.
- `docs/STATUS_20261125_SRP218935.md`: completed SRP218935 source + historical-QC reconciliation (no clinical admissions).
- `diagnostics/20261125_v6_srp218935_reproduced_v1/`: exact script, captured stdout and SHA256s.
- Default policy: archive and push each completed diagnostic before proceeding, using only redacted public-safe evidence.

## Read first in a new conversation

1. `docs/STATUS_20261123.md` — latest verified status and open gates.
2. `AUDIT_REGISTRY.tsv` — diagnostics for which code and output are in this bootstrap.
   `docs/COMPLETED_AUDITS_INDEX.tsv` — prior completed audits documented in the historical README,
   clearly marked when original script bytes are not yet imported.
3. `decisions/SPECIMEN_DECISIONS_SRP060016.tsv` — evidence-only accession mapping (NOT admissions).
4. `docs/README_MANUAL_CURATION_UPDATED_20261120.md` — historical methods and earlier audits.
5. `docs/PROVENANCE_POLICY.md` — meaning of verified, reviewed, and pending.

## Files and directories

- `diagnostics/<date>_<topic>/script.py`: exact script **as supplied in chat**; independent cluster-source identity not yet checked.
- `diagnostics/<date>_<topic>/stdout.txt`: original cluster printout if provided, otherwise an explicitly marked reconstruction.
- `diagnostics/<date>_<topic>/interpretation.md`: supported conclusions and limitations.
- `evidence/V5/`: QA and checksum manifest only, not the large clinical tables.
- `tools/record_existing.sh`: import an already-run script + stdout safely into a new diagnostic folder.
- `scripts/originals/`: reserved for verified copies of original cluster scripts; empty until imported.

## Rules

- Never rewrite historical outputs or decision history. Add a new ID/version when evidence changes.
- Do not copy 5,241-row source tables, run-level matrices, HDF5s, raw patient-level clinical data, credentials or HPC secrets to a hosted repository.
- Use a **private GitHub repository**, subject to institutional rules, for this small audit log.
- Do not interpret a successful hash/technical test as clinical admission.
- Original ARCHS4 W, anchors and Stage 09–13 outputs are frozen and must not be modified.
- `SRP060016` specimen classification is resolved, but pretreatment and all release gates remain unresolved.

## Daily workflow

For a new inline diagnostic, first save the entire script to `diagnostics/YYYYMMDD_TOPIC/script.py`.
Then execute with your selected cluster interpreter, while capturing stdout and stderr:

```bash
set -o pipefail
REPO=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python
ID=YYYYMMDD_TOPIC
mkdir -p "$REPO/diagnostics/$ID"
"$PYTHON" "$REPO/diagnostics/$ID/script.py" 2>&1 \
  | tee "$REPO/diagnostics/$ID/stdout.txt"
sha256sum "$REPO/diagnostics/$ID/script.py" \
          "$REPO/diagnostics/$ID/stdout.txt" \
  > "$REPO/diagnostics/$ID/SHA256SUMS.txt"
```

For diagnostics already run where both script and output files exist, import without rerunning:

```bash
"$REPO/tools/record_existing.sh" 20261123_audit_name /absolute/code.py /absolute/printout.txt
```

Record evidence in `AUDIT_REGISTRY.tsv`, describe findings in an interpretation markdown file,
review the diff, then `git add`, `git commit`, and `git push` intentionally. Never automate git push
from cluster jobs. 

## GitHub

This repository is initially only a portable local Git repository. A hosted GitHub remote
must be created and authorized separately. ChatGPT can access it **only if** its GitHub
connection is enabled for the account and the repository is accessible to that connection.
In a future chat, give the repository link and ask to read the status/registry before new diagnostics.
