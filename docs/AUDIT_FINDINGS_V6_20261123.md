# Recount3 certification — cumulative diagnostic findings

Date: 2026-11-23

**Scientific status: UNSIGNED. No recount3 sample admitted.**

Each record identifies its original script and cluster printout.
A checksum proves file identity, not clinical correctness.
Failed and superseded investigations remain part of the history.

## Diagnostic record

### 20261123_HISTORICAL_INVENTORY

- Status: `COMPLETED`
- Original script: `tools/inventory_historical_artifacts_READONLY_v1.py`
- Script SHA256: `9ae8fdf8a7b2f7f56955d03cd864fa8e4e1e4ff24e5511f5c00154ef2e82910b`
- Original cluster printout: `diagnostics/20261123_historical_recovery_inventory/stdout.txt`
- Printout SHA256: `4cc3610c050857a612683f9fd876f8c189e70002e5aa4948c1b37eb2391adf5c`
- Finding: 30/34 historical artifact names located in inspected directories; subsequent archive inspection recovered additional missing components.

### 20261123_ARCHIVE_MEMBER_INVENTORY

- Status: `COMPLETED`
- Original script: `tools/inventory_archive_members_READONLY_v1.py`
- Script SHA256: `ce9aeaa0e5e5ab3ccefa10db9d907734cc40d931356bf4fe26a64bf7f86439ee`
- Original cluster printout: `diagnostics/20261123_archive_members/stdout.txt`
- Printout SHA256: `6a1c31b5da2148820415980c5794dd542ea61fe49955088209c9f291f686cde8`
- Finding: Existing archive contents catalogued without extraction or scientific execution.

### 20261123_ARCHIVE_RECOVERY

- Status: `COMPLETED`
- Original script: `tools/recover_archive_evidence_v1.py`
- Script SHA256: `35febf399f907550d3f8cb1634ba4bd2f4ba47149561ee8acf6e132103a7fde0`
- Original cluster printout: `diagnostics/20261123_archive_recovery/stdout.txt`
- Printout SHA256: `432958846466bf1a8be268bf95aefa18053f56e01e8a2255fc0a17dcb030dc65`
- Finding: Nine archives verified and 28 historical source/report files recovered; no scientific reruns.

### 20261123_V6_SPECIMEN_PRECEDENCE

- Status: `COMPLETED_EVIDENCE_ONLY`
- Original script: `tools/v6_evidence_precedence_READONLY_v1.py`
- Script SHA256: `3d670c1b6fc71d0897674098e955adc4fe9e6cda2a3b678e58e1afd5b77081d9`
- Original cluster printout: `diagnostics/20261123_v6_evidence_precedence_dryrun/stdout.txt`
- Printout SHA256: `149c950f19eefc915b1fbbe43ccf7ecff7ac46ddf36440844460ebc571da8e97`
- Finding: 2055 screen records; 22 cumulative non-tumour restrictions, five exact primary CRC specimens, 52 sarcoma site holds, 280 V5 context records, 1696 historical reviews. No approvals.

### 20261123_V6_TARGETED_REVIEW

- Status: `COMPLETED_EVIDENCE_ONLY`
- Original script: `tools/v6_targeted_evidence_review_READONLY_v1.py`
- Script SHA256: `48308d71f4c499a0c8b55b165a00e1690f7c4bee9f03774d9bbb5a09b6836be3`
- Original cluster printout: `diagnostics/20261123_v6_targeted_evidence_review/stdout.txt`
- Printout SHA256: `ef40384000085eea7ab486a7d467b7db8079e2a3e2013c7a2d9372a2855cc7b3`
- Finding: 1696 historical records divided into 599 unattested, 603 surgery context, 494 study-source supported. SRX106046 normal; SRX106048 cancerous but primary site and patient relation unverified.

### 20261123_V6_STUDY_MATRIX_V1

- Status: `SUPERSEDED_REPORTING_BUG`
- Original script: `tools/v6_study_review_matrix_v1.py`
- Script SHA256: `96a4011ab91bccc22d28957fb19388f172d49a6011576826d47419f9f10ad6e5`
- Original cluster printout: `diagnostics/20261123_v6_study_review_matrix/stdout.txt`
- Printout SHA256: `2659b24c5a522bf343cac409a9b7d9735f87c1c2b3ce6cbbd47f5b6b6c794fcb`
- Finding: Incorrectly reported 176 screen contributors due to defaultdict key insertion; study-count interpretation invalid. Superseded by v1.1.

### 20261123_V6_STUDY_MATRIX_V1_1

- Status: `COMPLETED_CORRECTED`
- Original script: `tools/v6_study_review_matrix_v1_1.py`
- Script SHA256: `4df92d97fdd4f03279e5074782d64e1f915cf8739f9e0eb6a6bc5df5be4ba678`
- Original cluster printout: `diagnostics/20261123_v6_study_review_matrix_v1_1/stdout.txt`
- Printout SHA256: `d1b02ebea48603e33c37d8add28e23d05746e3dc3bd502f8b288fb47f63a97e8`
- Finding: 176 inventory studies; 77 screen contributors; 99 outside screen. Historical-only 599 unattested, 603 surgery, 494 study-source-supported.

### 20261123_V6_599_SOURCE_ASSEMBLY

- Status: `FAILED_PARTIAL_EVIDENCE`
- Original script: `tools/v6_primary_unattested_source_review_v1.py`
- Script SHA256: `d6b20c145b14f90a70b95c220d7e178d74d282781e574d5c79dc49e916435c37`
- Original cluster printout: `diagnostics/20261123_v6_primary_unattested_review/stdout.txt`
- Printout SHA256: `c67468abbd03516d9b96e0c3a442d80721942db89a51acbf0b8ccb53ed721092`
- Finding: 83 DRP005276, 17 ERP010930, and 25 ERP023222 experiments matched before duplicate SRX145877 in SRP012656 stopped execution. No completed 599-row output.

### 20261123_V6_DUPLICATE_SRA_AUDIT

- Status: `COMPLETED_SPECIMEN_METADATA_CONCORDANCE`
- Original script: `tools/v6_duplicate_sra_rows_READONLY_v1.py`
- Script SHA256: `6eb23b4b43ca0f6ff58bb78fd5c0f3d829a16d08e94686a39be2d46a5b827865`
- Original cluster printout: `diagnostics/20261123_v6_duplicate_sra_rows/stdout.txt`
- Printout SHA256: `715633781421b1055c10ce52dfb249b18ee78a7496c81504d07b941b61006c4e`
- Finding: Inspects all 599 targeted experiment identities for repeated source rows, differing specimen metadata and missing exact matches. No automatic deduplication.

## Historical evidence

Previous diagnostics and scientific conclusions remain indexed in:
- `docs/COMPLETED_AUDITS_INDEX.tsv`
- `docs/HISTORICAL_SOURCE_RECOVERY_MANIFEST.tsv`
- `docs/HISTORICAL_ARCHIVE_RECOVERY_MANIFEST.tsv`
- `evidence/history/`
- `docs/README_MANUAL_CURATION_UPDATED_20261120.md`

The original historical record must not be overwritten by later findings.

## Outstanding clinical gates

1. Resolve primary tumour/site identity at experiment level.
2. Establish treatment timing relative to sampled tissue.
3. Resolve within-study and cross-study patient independence.
4. Reconcile run QC, assay type, and historical technical holds.
5. Establish atlas independence and full-499 model transfer.
6. Obtain explicit clinical and technical approval before cohort lock.

## Data governance

Source printouts containing detailed sample attributes remain
on the cluster unless explicitly cleared for repository upload.
Their paths and checksums preserve traceability.

#### Interpretation after inspecting original cluster output

- 599/599 exact experiments matched across 39 studies; none missing.
- 59 experiments have multiple raw metadata rows, in five studies.
- Repeated-row distribution: SRP012656 6, SRP018008 1,
  SRP040442 8, SRP090849 22, SRP094795 22.
- 147 additional raw rows beyond one per experiment.
- Zero conflicts within the checked fields: sample_acc, biosample,
  sample_title, sample_attributes.
- Distinct run_alias values and sequencing totals explain row differences,
  but run_alias is not a verified official SRA run accession.
- No patient splitting, specimen exclusion, run aggregation or clinical
  admission follows from this diagnostic.
- The earlier 599-source-assembly v1 failure is explained. Preserve that
  failure; create a separate version with specimen-field concordance checks.
