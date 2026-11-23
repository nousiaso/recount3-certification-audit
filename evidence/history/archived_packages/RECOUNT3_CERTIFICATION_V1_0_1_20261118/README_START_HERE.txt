RECOUNT3 FINAL CERTIFICATION v1.0.1 — SOURCE AND SPECIMEN EVIDENCE GATES
Date: 2026-11-18
Status: SOFTWARE CORRECTION VALIDATED OFFLINE; SCIENTIFIC REVIEW NOT COMPLETE

PURPOSE
-------
Patch the two important release-gate omissions in v1.0.0 without modifying it:
1. Study source_verified=YES does NOT pass without documented sample types,
   patient-key provenance, source reference, and named reviewer.
2. No primary tumour enters the released cohort without a reviewed per-accession
   specimen record and explicit primary-tumour eligibility decision.

Additional safeguards:
- Curated patient keys require per-library patient-evidence and reviewer, not just
  an inherited patient string. The 581 BioSample proxies need independent patient
  evidence; BioSample is never accepted as a patient identifier.
- Repeat-subject groups must be reviewed and remain clustered by patient.
- Different libraries with identical reviewed specimen IDs are held to prevent
  biological pseudoreplication. Shared BioSample IDs also force review.
- Run QC uses v1.0.0's 4,166-run evidence and non-overridable DGE exclusions;
  manually changing 409 HOLDs requires reviewer and evidence.
- Cross-atlas independence and all-499-feature W transfer remain separate gates.
- No direct matrix pooling or fixed-W-499 model equivalence is approved by this code.

SOURCE CODE CONTRACT
--------------------
cert_v101/certify.py                      v1.0.1 controller
cert_v101/_v100_reference_readonly.py    PINNED exact v1.0.0 gate code
cert_v101/tests/test_v101.py             scientific gate tests

The controller verifies that installed final_cert/certify.py is byte-identical to
its pinned v1.0.0 code. It reads v1.0.0's sealed evidence lock and all relevant
Stage 2-5 existing manifests. Do not edit either version of source code.

USAGE ON CLUSTER — from 26_recount3_external_primary_v1
--------------------------------------------------------
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/26_recount3_external_primary_v1
PREV="$ROOT/../27_recount3_final_cert_v1_0_0"
WORK="$ROOT/../28_recount3_final_cert_v1_0_1"
PYTHON=/home/tpn890/orestis/conda_envs/spatialecotyper_v1_0_3_human/bin/python

# Once the ZIP is copied to $ROOT:
unzip -n RECOUNT3_CERTIFICATION_V1_0_1_20261118.zip -d "$ROOT"
sha256sum -c CERT_V101_PACKAGE_SHA256SUMS.txt
"$PYTHON" cert_v101/tests/test_v101.py

# READ ONLY (checks the original 129/4166/3639 inputs and 154-file lock)
"$PYTHON" cert_v101/certify.py --root "$ROOT" --previous "$PREV" preflight

# CREATES A NEW SIBLING work directory, NEVER OVERWRITES previous work
"$PYTHON" cert_v101/certify.py --root "$ROOT" --previous "$PREV" --work "$WORK" init

# First conservative assessment (safe to run with no decisions yet)
"$PYTHON" cert_v101/certify.py --root "$ROOT" --previous "$PREV" --work "$WORK" assess
cat "$WORK/assessment_v101/release_gates.json"

EXPECTED INITIAL ASSESSMENT
---------------------------
0 source-attested studies; 0 accession-reviewed specimens; 0 approved primary
libraries; 409 QC run holds; 581 patient proxies; 109 repeated-patient groups;
34 unresolved W features; 37 protocol-excluded libraries.

EDIT ONLY REVIEW FILES IN THE NEW $WORK:
---------------------------------------
review/study_source_attestations.tsv      126 study-level declarations
review/specimen_eligibility_3602.tsv     3602 accession-level classifications
review/patient_identity.tsv              3602 reviewed patient assignments
review/repeated_subject_relationships.tsv 109 repeated-subject review groups
review/run_decisions.tsv                 4166 individual runs
review/held_W_feature_resolutions.tsv    34 held genes (original verified mapping)
review/model_transfer_review.json        provenance notes; NOT self-certifying
review/ARCHS4_crosswalk_review.json      provenance notes; NOT self-certifying

Specimen template schema: (immutable): library_id, study, biosample, cancer,
original library selection, preservation, cancer variant, source row SHA256.

Editable accession-level fields:
- reviewed_specimen_id: source-grounded specimen identifier, not just library ID.
- specimen_class: PRIMARY_TUMOUR | RECURRENT_TUMOUR | METASTATIC_TUMOUR |
  MALIGNANT_ASCITES | ADJACENT_NORMAL | NORMAL_TISSUE | CELL_LINE | OTHER |
  UNRESOLVED
- treatment_status: PRETREATMENT | TREATMENT_NAIVE_CONFIRMED | ON_TREATMENT |
  POST_TREATMENT | UNKNOWN | NOT_APPLICABLE
- disease_context, collection_timepoint, anatomic_site, tumour_region.
- primary_cohort_decision: INCLUDE_PRIMARY | EXCLUDE_NONPRIMARY | HOLD.
- accession_source_reference: primary GEO/SRA/ENA/publication source or
  provenance-backed on-cluster sample map (record file SHA256 where available).
- accession_evidence_record: specific source record/row, mentioning the EXACT
  SRX/DRX/ERX library accession or BioSample accession; avoid generic article-only
  claims. A citation to a paper with no accession mapping is NOT sufficient.
- reviewer: recorded researcher who inspected the accession-level evidence.

REVIEW SEMANTICS
----------------
INCLUDE_PRIMARY requires PRIMARY_TUMOUR and specimen source/reviewer.
EXCLUDE_NONPRIMARY allows evidence-backed rejection of normal, ascites,
recurrent/metastatic or cell-line material from the PRIMARY cohort.
TREATMENT UNKNOWN is recorded as unknown, not promoted to treatment-naive.
Any ambiguous specimen remains HOLD. Repeated regions and longitudinal
samples are patient-clustered and not considered independent patients.

STUDY ATTESTATION
-----------------
Set source_verified=YES ONLY after confirming its source reference, source
sample types, patient-key provenance and accession-level sample map; reviewer
must be recorded. A typed YES by itself cannot certify any library.

PATIENT EVIDENCE
----------------
All 3602 patient_identity rows require reviewed_patient_id,
source_evidence_reference and reviewer. The 3021 curated keys are candidates,
NOT proven by being populated. The 581 BioSample proxies cannot be promoted to
patient IDs without source-based subject mapping. Repeated-subject relationships
require one supported classification and reviewer; they remain statistically
clustered even after classification.

RUN INCLUSION
-------------
The v1.0.0 decision ledger remains copied intact. Do not bulk-change 395 general
HOLDs. The 14 sparse kidney runs have strong exclusion evidence, but still need
reviewer/source confirmation in the run ledger. At least one run must remain and
all pairs of retained runs must satisfy the conservative correlation check.

REASSESS AFTER ENTERING REAL EVIDENCE
-------------------------------------
"$PYTHON" cert_v101/certify.py --root "$ROOT" --previous "$PREV" --work "$WORK" assess

Review outputs:
assessment_v101/release_gates.json
assessment_v101/library_certification_3639.tsv
assessment_v101/SHA256SUMS.txt

Only source-attested, accession-reviewed, patient-reviewed, QC-approved primary
specimens become CERTIFIED_PRIMARY_WITHIN_STUDY_PATIENT_AWARE_ONLY.
This does NOT claim distinct patients across studies, ARCHS4 independence,
model-input equivalence, or permission to pool raw expression.

OPTIONAL RAW COVERAGE AGGREGATION AFTER SUBSET CERTIFICATION
------------------------------------------------------------
"$PYTHON" cert_v101/certify.py --root "$ROOT" --previous "$PREV" \
  --work "$WORK" aggregate-study --study SRPXXXXX

This refuses to aggregate without a certified subset, verifies assessment
SHA256SUMS, and sums raw, unscaled G026 coverage for reviewed retained runs.
The output is NOT approved frozen-W model input. It does not silently collapse
separate biological libraries or convert technical runs into patient replicates.

INDEPENDENT CROSSWALK AND MODEL VALIDATION (OPTIONAL EVIDENCE)
-------------------------------------------------------------
The assess command also accepts:
  --atlas-crosswalk PATH    complete 9465-sample source-backed identity map
  --projection-pairs PATH   ACTUAL same-specimen 143x10 full-model fractions

The crosswalk must match exact frozen sample IDs and source accessions; an
incomplete crosswalk never proves external independence. The paired fraction
file must be generated by a reproducible actual 499-feature frozen-model
projection, NOT simulated data or a 465-gene approximation. Score comparisons
are descriptive and cannot alone pass the independent-study transfer gate.
For the 34 features, use review/held_W_feature_resolutions.tsv with original
HDF5 and G026 rows, Ensembl IDs, source evidence and named reviewer. The code
independently re-checks actual HDF5 gene labels before recording any rescue.

FOR THE PAPER
-------------
Technical evidence remains locked as v1.0.0 with 154 source-file hashes.
This v1.0.1 is the corrected reviewer-attestation and specimen-level candidate
cohort gate. The statistical release remains on HOLD until all necessary
specimen/patient decisions and model-transfer work are supported.

LIMITATIONS
-----------
Software verifies completeness, accession linkage, internal consistency and
hashes; it cannot independently authenticate a human's stated pathology or
publication content. The reviewer must check original sample records.
No clinical facts are filled from study summaries, and no model scores are
invented. All tests use synthetic attestations, not real source verification.
