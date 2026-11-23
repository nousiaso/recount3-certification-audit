RECOUNT3 CANCER SPECIMEN AND PATIENT EVIDENCE ADJUDICATION
2026-11-19 | v1.0.0 | REVIEW-ONLY, DO NOT IMPORT AS FINAL DECISIONS

SCOPE AND INPUTS
- 3602 candidate source-linked SRA libraries from 126 studies.
- Unchanged previous user-provided source/triage/run/patient/W review data.
- SHA256 of each immutable input recorded in REVIEW_OUTPUT/09_adjudication_summary.json.
- Scripts tested on locally uploaded records; cluster gene-axis probe NOT yet executed.
- Original frozen ARCHS4 9465-sample analysis and v1.0.0 evidence lock untouched.

NEW FINDINGS FROM THE FULL ACCESSION REVIEW
- 920 accessions have direct submitted primary-origin language in a sample-specific attribute.
  The 920 occur across 19 studies. This supports PRIMARY ORIGIN, not untreated status,
  patient independence, whole-tumour bulk composition, or model compatibility.
- Another 214 have direct tumour descriptions without sufficient primary-origin proof.
- 11 have explicit non-tumour source evidence; review and exclude from primary-tumour set.
- 4 Lexogen QuantSeq 3-prime breast libraries need assay compatibility exclusion review.
- 94 contain explicit treatment context requiring dedicated review (74 SRP114518 and
  20 SRP072494); baseline/pre-guadecitabine does not mean untreated disease.
- 45 are multi-region lung specimens, not 45 independent patients.
- 63 need whole-tumour composition/compartment review (42 salivary Epi/Myo, 21 gastric epithelial).
- Remaining 2251 lack sufficient direct evidence for the intended PRIMARY cohort.
- 281 accessions in 13 studies include explicit patient/subject token fields.
- 0 source patient tokens conflict with curated patient groupings; one curated group
  has two labels, P265T and P265T replicate (two library entries).
- 25 additional plausible relationship groups: 17 salivary Epi/Myo numeric roots and
  8 breast RR- specimen roots. These are NOT source-confirmed patient equivalence.
- Prior unresolved queues remain 409 run holds, 109 repeated curated-patient groups,
  34 W genes, incomplete 9465-sample ARCHS4 alias crosswalk and missing model-transfer test.

PRECAUTIONS
- All 3602 entries remain clinical_specimen_release=HOLD. Reviewer approval is not fabricated.
- No row may be marked APPROVED from this output alone. Exact accession and source citations
  are supplied for downstream signoff in cert_v101 review tables.
- No automatic exclusion from general keyword matches such as 'normal/tumor = tumor' or
  'cell line = NA'. Disease history != specimen identity.
- Any future confirmed non-primary specimens should be excluded from the primary-tumour
  cohort, but can be preserved in an explicitly separate assay/biological-use manifest.
- Model transfer requires all 499 W features, true same-specimen model output comparison,
  cross-study tests, and full frozen atlas overlap audit.
- The original code cert_v101 has an explicit HOLD model gate; do not claim PASS by
  supplying descriptive scores alone.

OUTPUTS
REVIEW_OUTPUT/01_accession_specimen_review_3602.tsv: per-SRX evidence and prioritization.
REVIEW_OUTPUT/02_study_adjudication_queue_126.tsv: review workload by study.
REVIEW_OUTPUT/03_direct_patient_token_evidence_281.tsv: accession-level patient tokens.
REVIEW_OUTPUT/04_additional_possible_patient_relationships.tsv: 25 candidates.
REVIEW_OUTPUT/05_run_hold_queue_409.tsv: original triage evidence, no decisions changed.
REVIEW_OUTPUT/06_held_W_feature_gates_34.tsv: unresolved model features.
REVIEW_OUTPUT/07_public_source_registry_REVIEW_ONLY.tsv: four verified public sources.
REVIEW_OUTPUT/08_repeated_curated_subject_groups_109.tsv: existing group queue.
REVIEW_OUTPUT/09_adjudication_summary.json: mutually exclusive triage totals and hashes.
REVIEW_OUTPUT/SHA256SUMS.txt: generated evidence checksums.

READ-ONLY CLUSTER TARGETED W34 + ATLAS PROBE
From 26_recount3_external_primary_v1:
ROOT="$PWD"
WORK="$ROOT/../29_recount3_adjudication_20261119"
P="$WORK/RECOUNT3_CERTIFICATION_ADJUDICATION_20261119"
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python

# Ensure the ZIP has been extracted to $WORK (not into existing v1.0.1).
(cd "$P" && sha256sum -c PACKAGE_SHA256SUMS.txt)
"$PYTHON" "$P/test_reconciled_review.py"
"$PYTHON" "$P/probe_W34_and_ARCHS4_atlas_READONLY.py" \
  --root "$ROOT" \
  --source-evidence "$ROOT/../28_recount3_final_cert_v1_0_1/source_evidence_v1_0_0/accession_source_evidence_3602.tsv" \
  --out "$WORK/W34_AND_ATLAS_READONLY_PROBE"
cat "$WORK/W34_AND_ATLAS_READONLY_PROBE/03_readonly_probe_summary.json"

The probe examines only G026 63856 gene-row metadata, 1497 W mapping entries,
original ARCHS4 HDF5 gene symbol/Ensembl metadata (67186 rows), the frozen 9465
atlas sample-ID table and optionally 3602 SRA source-evidence rows.
NO large HDF5 expression matrix is read and NO frozen model is modified.
Output contains candidate W gene/row identity pairs and DIRECT GSM-overlap candidates,
NOT approved mappings or certified independence.
If real source schema differs, the command stops rather than guessing.

NEXT SCIENTIFIC ACTION
Review source evidence for every selected accession in the 126-study queue,
starting with the 920 direct-primary candidates and 11+4 exclusion/assay cases.
Complete patient and run decisions and the targeted 34-gene/atlas identity probe.
Do not start the new recount3-based statistical analysis until both biological
cohort and original frozen-W transfer input contracts have been independently locked.
