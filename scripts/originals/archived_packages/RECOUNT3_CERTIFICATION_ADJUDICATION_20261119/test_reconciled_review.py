#!/usr/bin/env python3
"""Test the immutable source-derived review pack (does not test model transfer)."""
import csv, json, pathlib, collections, hashlib, sys
base=pathlib.Path(__file__).resolve().parent/'REVIEW_OUTPUT'
def read(name):
 with open(base/name,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
for line in (base/'SHA256SUMS.txt').read_text().splitlines():
 expected,name=line.split(None,1); name=name.strip()
 actual=hashlib.sha256((base/name).read_bytes()).hexdigest()
 assert expected==actual,name
records=read('01_accession_specimen_review_3602.tsv')
assert len(records)==3602
assert len({(r['study'],r['library_id']) for r in records})==3602
assert all(r['clinical_specimen_release']=='HOLD' and r['patient_group_verified']=='NO' and r['ARCHS4_independence']=='NOT_CERTIFIED' and r['frozen_W_transfer']=='NOT_CERTIFIED' for r in records)
c=collections.Counter(r['evidence_category_NOT_CERTIFIED'] for r in records)
assert c == {'SPECIMEN_ELIGIBILITY_UNRESOLVED':2251,'DIRECT_PRIMARY_ORIGIN_METADATA':920,'SOURCE_EXPLICIT_NON_TUMOUR':11,'DIRECT_TUMOUR_METADATA_PRIMARY_UNPROVEN':214,'WHOLE_TUMOUR_COMPOSITION_REVIEW':63,'KNOWN_TREATMENT_CONTEXT_REVIEW':94,'MULTIREGION_DESIGN_PATIENT_AWARE':45,'PROTOCOL_COMPATIBILITY_HOLD':4}
assert len(read('02_study_adjudication_queue_126.tsv'))==126
pat=read('03_direct_patient_token_evidence_281.tsv')
assert len(pat)==281
assert sum(p['token_conflicts_with_curated_patient_group']=='YES' for p in pat)==0
assert sum(p['curated_group_has_multiple_direct_tokens']=='YES' for p in pat)==2
assert len(read('04_additional_possible_patient_relationships.tsv'))==25
assert len(read('05_run_hold_queue_409.tsv'))==409
assert len(read('06_held_W_feature_gates_34.tsv'))==34
assert len(read('08_repeated_curated_subject_groups_109.tsv'))==109
meta=json.loads((base/'09_adjudication_summary.json').read_text())
assert meta['n_certified_libraries']==0
assert meta['n_source_primary_origin_phrases']==920
print('PASS full 3602-library adjudication review integrity, 126 studies, 281 patient tokens, 409 runs, 34 W holds; ZERO unauthorized approvals')
