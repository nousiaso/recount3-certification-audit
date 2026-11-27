# SRP057793: exact historical QC hold reasons reviewed
2026-11-27. Evidence: user-supplied executed inline diagnostic stdout. Original executed script and stdout bytes were not captured as files in that invocation. The related read-only script tools/v8_srp057793_qc_reason_and_clinical_schema_READONLY_v1.py is explicitly a reconstruction with aggregate metric ranges instead of per-experiment rows; not claimed to be the byte-identical executed script.

Prior exact membership: 149 source experiments, 148 integrated candidates, 8 technical holds, 140 retained primary-screen candidates, 1 source-only experiment. All 8 held experiments exactly join to historical run-hold queue.

## Observed historical decisions
- reason_codes: WITHIN_STUDY_OUTLIER 8/8
- automatic_proposal: REVIEW_QC 8/8
- decision: HOLD 8/8
- triage_priority: P2_STUDY_RELATIVE_OUTLIER 8/8
- proposed_action_not_applied: HOLD_REASON_NOT_YET_DECISIVE 8/8

Observed ranges among the eight (not individual records):
- gene_detection_fraction 0.465062015–0.569218241
- top10_coverage_fraction 0.122681807–0.295514339
- average_mapped_length_bp 150.4–150.79
- gene_sum_coverage_bp 3589937303–5858830891

These are historical study-relative QC review holds, not demonstrated irreversible assay failures. Absolute values do not establish that they are acceptable either; original outlier reference distribution and rule must be recovered before any readmission. Preserve all eight holds. The 140 remain primary-screen candidates, not clinically certified admissions.

## Clinical evidence gap
For all 140 retained records, proposed_patient_group_NOT_CERTIFIED, patient_group_evidence_pointers and pretreatment_evidence_pointers have zero populated values. patient_group_source_type, pretreatment_source_text_status, stage37_protocol_review_category, study_detected_atlas_overlap and the five certification columns are populated with one distinct value each; actual values were not printed, so no certification inferred from population. No patient independence, treatment, material or protocol conclusion can be drawn from cardinality alone.

Next: recover aggregate actual decision values and original native study metadata / publication evidence, resolve specimen/patient/pretreatment/assay/nonoverlap independently; do not rerun V7 technical calibration or change cohort membership. Clinical admissions 0. Source-only record remains excluded pending original evidence. No frozen files modified.
