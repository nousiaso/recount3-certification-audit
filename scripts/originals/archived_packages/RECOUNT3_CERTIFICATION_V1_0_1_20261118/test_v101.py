#!/usr/bin/env python3
"""Offline scientific-gate tests, including real 3,639-library review-bundle dry run.
No cluster reads or writes and no user release modifications.
"""
import csv, os, sys, tempfile, unittest, json, shutil
from pathlib import Path
from unittest import mock
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import certify as c
import _v100_reference_readonly as old

class Tests(unittest.TestCase):
 def test_study_source_not_boolean(self):
  x={'source_verified':'YES','verified_sample_types':'PRIMARY_TUMOUR',
     'patient_key_provenance':'SRA sample source mapping','publication_or_SRA_reference':'','reviewer':'A'}
  self.assertFalse(c.reviewed_study(x))
  x['publication_or_SRA_reference']='https://www.ncbi.nlm.nih.gov/sra'
  self.assertFalse(c.reviewed_study(x))
  x['reviewer']='DrA'
  self.assertTrue(c.reviewed_study(x))
 def test_accession_specific_evidence(self):
  x={'library_id':'SRX123456','biosample':'SAMN999999','primary_cohort_decision':'INCLUDE_PRIMARY',
     'specimen_class':'PRIMARY_TUMOUR','reviewed_specimen_id':'PRIMARY123',
     'accession_source_reference':'https://www.ncbi.nlm.nih.gov/sra',
     'accession_evidence_record':'See source study only', 'reviewer':'ABC'}
  self.assertFalse(c.reviewed_specimen(x))
  x['accession_evidence_record']='SRX123456 tumour from primary resection source table row 3'
  self.assertTrue(c.reviewed_specimen(x))
 def test_unsupported_patient_proxy(self):
  x={'reviewed_patient_id':'SAMN999','biosample':'SAMN999','patient_key_from_curation':'SRP:P1',
     'proposed_patient_group':'UNKNOWN', 'patient_identity_state':'BIOSAMPLE_PROXY_NO_PATIENT_ID',
     'source_evidence_reference':'https://sra.ncbi.nlm.nih.gov','reviewer':'ABC'}
  self.assertFalse(c.reviewed_patient(x,True,True,False))
  x['reviewed_patient_id']='SRP:SUBJECT_BIOLOGICAL_1'
  self.assertTrue(c.reviewed_patient(x,True,True,False))
  x['source_evidence_reference']=''
  self.assertFalse(c.reviewed_patient(x,True,True,False))
 def test_specimen_incompatible(self):
  base={'library_id':'SRX123','study':'SRP123','biosample':'SAMN111','cancer':'lung',
   'original_library_selection':'','original_preservation':'','original_cancer_variant':'',
   'original_source_row_sha256':'','specimen_class':'PRIMARY_TUMOUR','treatment_status':'UNKNOWN',
   'primary_cohort_decision':'INCLUDE_PRIMARY'}
  c.validated_specimen([base],{'SRX123':base})
  x=dict(base,specimen_class='MALIGNANT_ASCITES')
  with self.assertRaisesRegex(ValueError,'non-primary'):c.validated_specimen([x],{'SRX123':base})
  x=dict(base,biosample='CHANGED')
  with self.assertRaisesRegex(ValueError,'Immutable'):c.validated_specimen([x],{'SRX123':base})
 def test_real_review_bundle_all_hold(self):
  source=Path('/mnt/data/recount3_review_work')
  if not source.is_dir():self.skipTest('Original review bundle fixture not mounted')
  with tempfile.TemporaryDirectory() as d:
   root=Path(d)/'root'; root.mkdir(); (root/'final_cert').mkdir()
   shutil.copy2(Path(c.__file__).resolve().parent/'_v100_reference_readonly.py',root/'final_cert/certify.py')
   for n in ('manifest/selected_libraries.tsv','manifest/studies.tsv','manifest/selected_runs.tsv'):
    (root/n).parent.mkdir(exist_ok=True,parents=True);shutil.copy2(source/n,root/n)
   prev=Path(d)/'old';prev.mkdir();shutil.copytree(source/'draft',prev/'draft');shutil.copytree(source/'review',prev/'review')
   for name in ('run_decisions','patient_identity','repeated_subject_relationships','study_source_attestations'):
    shutil.copy2(prev/'review'/f'{name}.tsv',prev/'draft'/f'{name}.tsv')
   (prev/'draft/evidence_sha256.tsv').write_text('role\tpath\tsha256\n')
   (prev/'TECHNICAL_EVIDENCE_LOCK_v1_0_0').mkdir()
   lock=prev/'TECHNICAL_EVIDENCE_LOCK_v1_0_0/LOCK_MANIFEST.json';lock.write_text('{}')
   work=Path(d)/'new'
   fakepaths={'libs':root/'manifest/selected_libraries.tsv'}
   with mock.patch.object(c,'preflight',return_value=({'libs':old.readtsv(fakepaths['libs'])},fakepaths)):
    c.init(root,prev,work)
   spec=old.readtsv(work/'review/specimen_eligibility_3602.tsv')
   self.assertEqual(len(spec),3602)
   self.assertTrue(all(s['primary_cohort_decision']=='HOLD' for s in spec))
   self.assertEqual(sum(c.reviewed_specimen(x) for x in spec),0)
   self.assertEqual(len(old.readtsv(work/'review/run_decisions.tsv')),4166)
   self.assertEqual(len(old.readtsv(work/'review/study_source_attestations.tsv')),126)
   self.assertEqual(len(old.readtsv(work/'review/held_W_feature_resolutions.tsv')),34)
   reviewdata={
    'libs':old.readtsv(root/'manifest/selected_libraries.tsv'),
    'runs':old.readtsv(root/'manifest/selected_runs.tsv'),
    'pairs':old.readtsv(source/'stage5_output/04_multirun_concordance_pairs.tsv'),
    'prior_link':old.readtsv(source/'stage3_output/00_ARCHS4_prior_linkage_flags_NOT_FINAL.tsv')}
   with mock.patch.object(c,'check_existing_lock'),mock.patch.object(c,'verify_snapshot'),\
      mock.patch.object(old,'inspect',return_value=(reviewdata,{},None,{})),\
      mock.patch.object(old,'verify_W',return_value=([],{})),\
      mock.patch.object(old,'verify_feature_resolutions',return_value={'n_held':34,'status':'34_UNRESOLVED'}):
    c.assess(root,prev,work)
   result=json.loads((work/'assessment_v101/release_gates.json').read_text())
   self.assertEqual(result['n_candidate_libraries'],3602)
   self.assertEqual(result['n_primary_within_study_certified'],0)
   self.assertEqual(result['n_run_holds'],409)
   self.assertEqual(result['source_attestations_incomplete'],126)
   self.assertEqual(result['n_features_still_held'],34)
   self.assertEqual(result['n_library_statuses']['EXCLUDE_PROTOCOL'],37)
   self.assertEqual(result['n_library_statuses']['HOLD_REVIEW'],3602)
   # Source- and accession-backed approval of one nonflagged curated patient:
   # all other records remain held, demonstrating the two cohort gates separate.
   target='DRX159123'; study='DRP005276'
   def edit_tsv(name,lookup,updates):
    path=work/'review'/name;rows=old.readtsv(path)
    for row in rows:
     if lookup(row):row.update(updates)
    old.writetsv(path,rows,list(rows[0]))
   edit_tsv('study_source_attestations.tsv',lambda r:r['study']==study,{
     'source_verified':'YES','verified_sample_types':'PRIMARY_TUMOUR',
     'patient_key_provenance':'Original sample-to-subject map accession-level cohort',
     'publication_or_SRA_reference':'https://www.ncbi.nlm.nih.gov/sra',
     'reviewer':'Reviewer_A'})
   edit_tsv('patient_identity.tsv',lambda r:r['library_id']==target,{
     'source_evidence_reference':'https://www.ncbi.nlm.nih.gov/sra',
     'reviewer':'Reviewer_A'})
   edit_tsv('specimen_eligibility_3602.tsv',lambda r:r['library_id']==target,{
     'reviewed_specimen_id':'SPECIMEN_X_1','specimen_class':'PRIMARY_TUMOUR',
     'treatment_status':'UNKNOWN','primary_cohort_decision':'INCLUDE_PRIMARY',
     'accession_source_reference':'https://www.ncbi.nlm.nih.gov/sra',
     'accession_evidence_record':'DRX159123 confirmed primary tumour in accession metadata table',
     'reviewer':'Reviewer_A'})
   with mock.patch.object(c,'check_existing_lock'),mock.patch.object(c,'verify_snapshot'),\
      mock.patch.object(old,'inspect',return_value=(reviewdata,{},None,{})),\
      mock.patch.object(old,'verify_W',return_value=([],{})),\
      mock.patch.object(old,'verify_feature_resolutions',return_value={'n_held':34,'status':'34_UNRESOLVED'}):
    c.assess(root,prev,work)
   result=json.loads((work/'assessment_v101/release_gates.json').read_text())
   self.assertEqual(result['n_primary_within_study_certified'],1)
   self.assertEqual(result['n_library_statuses']['HOLD_REVIEW'],3601)
   self.assertEqual(result['n_library_statuses']['EXCLUDE_PROTOCOL'],37)
   self.assertEqual(result['model_transfer_gate'].startswith('HOLD'),True)

if __name__=='__main__':unittest.main(verbosity=2)
