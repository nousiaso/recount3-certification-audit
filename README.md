# ARCHS4 × recount3 — clinical certification & reproducible audit log

**Status: UNSIGNED. Do not run SpatialEcoTyper biological inference from the recount3 cohort yet.**

This is a **small provenance repository**, not a copy of the cluster datasets or frozen ARCHS4 release.
It tracks the exact diagnostic scripts, outputs, checksums, interpretation and outstanding gates.
It was initialized from user-supplied November 2026 evidence; historical scripts not physically
provided have **not** been invented or declared archived.

## Current critical hold — executed V7J CPM invariance failure (2026-11-27)

**Actual V7J source preparation: PASS. Projection calibration: HOLD — not an assay-transfer result.**

- Original recount3 G026 source: all **454 SRRs** from **143 SRXs** processed across 63,856 genes, grouped without raw coverage mass loss (max relative error **0**). 56,937 valid gene-symbol groups, all 499 W genes, UQ-evaluable 143/143.
- Original frozen ARCHS4 CPM 143-GSM Z replay **PASS** (max abs **2.04e-14**), exact 499-gene harmonization PASS and archived H/fraction consistency PASS (**8.88e-16**).
- **New result:** running unchanged frozen W on *only the 143 original ARCHS4 frozen-Z columns* differs from the H originally obtained within the **9,465-sample** projection: max absolute H difference **0.0126241**, maximum ecotype fraction TVD **0.00959059**. This triggered the deliberate stop **before recount3 projection**.
- Therefore there are **no paired assay ecotype scores, no aggregate transfer TSVs and no clinical admissions**. Do not relax the invariance threshold or overwrite the failed V7J audit.
- [Executed V7J hold and V7K next diagnostic](docs/STATUS_20261127_V7J_EXECUTED_INVARIANCE_HOLD_V7K_NEXT.md). The immediate next step is the **read-only installed `DeconvoluteSE` 1.0.3 source inspection**, not another full matrix run:
  `bash scripts/run_v7k_installed_DeconvoluteSE_source_READONLY.sh`.
- The mechanism of the 143-vs-9465 projection discrepancy is **unverified** pending original installed solver/source inspection. Frozen 9,465-sample V7H2 numerical replay remains PASS.

## Next execution gate — V7J frozen-reference matched assay calibration (2026-11-27)

**V7J static preflight EXECUTED PASS; full 454-SRR recount3 calibration NOT YET EXECUTED.** Eight original Stage10 functions were confirmed byte-identical, Python/R source and installed dependencies passed, and original input hashes matched. The earlier V7G remains immutable and on HOLD. See [executed static preflight evidence](docs/STATUS_20261127_V7J_STATIC_PREFLIGHT_EXECUTED_PASS.md).

- [V7J implementation and scientific protocol](docs/V7J_FROZEN_ARCHS4_143_RECOUNT3_TRANSFER_IMPLEMENTATION.md): direct comparison of the **actual original frozen 9,465-cohort ARCHS4 Stage09/Stage10 normalized inputs, Z, H and fractions for the 143 matching GSMs** against new recount3 G026 projection; 143 SRX samples from all 454 original SRR runs.
- Original Stage10 R functions (harmonization, directional features, original DeconvoluteSE call, raw-H fractions, adequacy, normalization robustness) are copied **verbatim**, with a static byte-fidelity check against original verified Stage10 source.
- Three branches CPM (primary), TMM-CPM and UQ-CPM (external-assay cohort-specific sensitivities), unchanged official W and frozen pan-cancer anchors. Strict 500-feature retention, original archived ARCHS4 Z/H/adequacy equality checks, repeat/order sentinel, nine-gene counterfactual and local-only per-specimen outputs.
- **Read-only preflight complete:** `bash scripts/run_v7j_static_frozen_transfer_preflight_READONLY.sh` passed all four markers. Preserve the original preflight logs; do not rerun with the same audit path.
- **Next executable step:** `sbatch scripts/slurm_v7j_frozen_ARCHS4_vs_recount3_calibration.sbatch`. The runner refuses existing versioned output directories. Technical assay transfer results and clinical patient/sample independence remain **unverified**; clinical admissions **0**.

## Current checkpoint — V7I PASS, recount3 calibration not yet run (2026-11-27)

**Frozen-source and numerical reproducibility: PASS. Clinical validation/admission: HOLD (0 admissions).**

- **V7H2 original ARCHS4 replay PASS:** 9,465 original atlas samples, 499 W genes, 500 directional features, and CPM/TMM-CPM/UQ-CPM. The largest frozen-Z discrepancy was `4.44e-14`, the largest raw-H discrepancy `4.44e-15`, and dominant-component concordance was **100%** in every branch. See [executed replay results](docs/STATUS_20261127_V7H2_ARCHS4_EXACT_NUMERICAL_REPLAY_PASS.md).
- **V7I exact-input preflight PASS:** 143 matching SRX experiments, 143 GSMs, 143 BioSamples, 454 SRR runs. Original G026 matrix contains these runs; G026 source map has 63,856 genes, all **499 model genes** resolve in both assays (0 recount3 alias resolutions, 7 ARCHS4). The bridge, map and original G026 SHA256s match V7G's exact source expectations. See [executed V7I audit](docs/STATUS_20261127_V7I_EXECUTED_BRIDGE_LOCK_AND_TRANSFER_NEXT.md).
- **Crucial update: all 143/143 bridge GSMs are in the original 9,465-sample frozen atlas.** The exact archived ARCHS4 Stage09 normalization, Stage10 scaled-Z/raw-H and scores must be the comparison baseline. Do **not** re-estimate ARCHS4 TMM/UQ on 143 specimens.
- **V7G v1 remains HELD, superseded in design, NOT an executed or approved assay-transfer result.** A new version must reconcile authentic Stage10 W-gene harmonization, directional-feature variance gating, adequacy and frozen-ARCHS4 baseline equality. recount3 G026 is gene-body base-pair coverage, not ARCHS4 Kallisto pseudocounts; paired assay equivalence remains unproven.
- Prior V7C/H failures and their versioned audits are preserved; the frozen ARCHS4 source and cohort certifications have not been modified. Zero clinical admissions and no FUSCC release.

**Next scientific step:** implement and check a new immutable matched-assay transfer version, explicitly anchored against the original frozen 143-GSM ARCHS4 outputs. Keep clinical specimen review separate.

## Latest V6 checkpoint (2026-11-26)

- `docs/STATUS_20261126_SRP057793_AUTHORITATIVE_QC_SCHEMA.md`: SRP057793 eight held experiments exactly join to the original QC reason queue.
- `docs/STATUS_20261126_SRP057793_EXCLUSION_BIOSAMPLE.md`: SRP057793 140 nonheld exact V5 match; one source-only record and eight hold reasons pending.
- `docs/STATUS_20261126_SRP057793_COHORT_MEMBERSHIP.md`: SRP057793 149→148→140 membership; exclusions pending.
- `docs/STATUS_20261126_V6_77_STUDY_CLINICAL_TRIAGE.md`: all 77 studies clinically unsigned; next targeted study SRP057793.
- `docs/STATUS_20261126_V6_PRIMARY_STUDY_PRIORITY.md`: 48 queued primary studies contribute 300/2055 libraries.
- `docs/STATUS_20261126_SRP185703_CROSSWALK_HOLD.md`: 75 GIST experiments; published clinical subset not accession-linked; no admissions.
- `docs/STATUS_20261126_SRP185703_JCI_SUPPLEMENT_SCHEMA.md`: JCI supplement verified; clinical crosswalk pending.
- `docs/STATUS_20261126_SRP185703_UNTREATED_PRIMARY_EVIDENCE.md`: SRP185703 published UPG subgroup; per-experiment crosswalk pending.
- `docs/STATUS_20261126_SRP185703_PRIMARY2055_RECONCILED.md`: SRP185703 75/75 primary-screen membership; clinical gates unsigned.
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
