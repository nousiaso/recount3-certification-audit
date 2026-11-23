#!/usr/bin/env python3
"""Deterministic, source-scoped specimen/patient crosswalk; never certifies primary eligibility.

Consumes bundled V4 clinical ledger and immutable Stage40 SRA/GEO evidence snapshots.
Creates one additive 5,241-row evidence layer and targeted exception files.
No network, no original file modifications, no model fitting.
"""
import argparse, csv, hashlib, json, re, shutil
from pathlib import Path
from collections import Counter, defaultdict

HERE=Path(__file__).resolve().parents[1]
INPUT=HERE/'inputs'
BASE=INPUT/'V4_5241_CLINICAL_SOURCE_UNSIGNED.tsv'
SRA=INPUT/'STAGE40_3602_SRA_SAMPLE_NAMES.tsv'
GEO=INPUT/'STAGE40_GEO_CLINICAL_SOURCE_EVIDENCE.tsv'
EXTRA=INPUT/'V4_SOURCE_ONLY_1123_HELD.tsv'
STUDIES=INPUT/'V4_STUDIES_176_SOURCE_LEDGER.tsv'
EXPECTED={
 'V4_5241_CLINICAL_SOURCE_UNSIGNED.tsv':'2c780cbbb0344f7f02882632ea655005f15a03309be3c103064188d3186aaad2',
 'STAGE40_3602_SRA_SAMPLE_NAMES.tsv':'50d3640ece5ea036ecc62a4c93a17f596b0e1cbad13f621390037e5d59f4879f',
 'STAGE40_GEO_CLINICAL_SOURCE_EVIDENCE.tsv':'2513fd25f1beee632f2b04a3b805fe6b3f12c9991a2bd73c27bb4cd6d05e3b41',
 'V4_SOURCE_ONLY_1123_HELD.tsv':'4158097b0ac67ddd338d74781d2fcc370a699309f895f780adb77e763c6adf61',
 'V4_STUDIES_176_SOURCE_LEDGER.tsv':'b9e1ae8d91307cc47c5411f1b8b201db2a6cfe92a8db0d0f7d79917abaef2db7',
}
URL={
 'SRP064894':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5337622/',
 'SRP124458':'https://www.ncbi.nlm.nih.gov/bioproject/PRJNA416439',
 'SRP067524':'https://doi.org/10.1002/cncr.29959',
 'SRP119917':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6555818/',
 'SRP139591':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112996',
 'SRP152822':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116782',
 'SRP119676':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104766',
 'SRP055438':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66207',
 'SRP057793':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75885',
}
NEW=['v5_original_sra_sample_name','v5_original_sra_sample_title','v5_specimen_evidence_class',
 'v5_specimen_evidence_scope','v5_tumour_subtype','v5_source_metastasis_annotation',
 'v5_patient_group_within_study','v5_patient_group_evidence_grade','v5_study_group_size',
 'v5_specimen_region_or_replicate','v5_source_evidence_url',
 'v5_proposed_clinical_action_NOT_APPROVAL','v5_remaining_gates',
 'v5_source_graph_match','v5_final_primary_approval']

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1<<20),b''): h.update(block)
 return h.hexdigest()

def read(p):
 with p.open(newline='',encoding='utf-8') as f:
  dr=csv.DictReader(f,delimiter='\t');return dr.fieldnames,list(dr)

def save(p,fields,rows):
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)

def fail_unless(v,msg):
 if not v: raise AssertionError(msg)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 fail_unless(not a.out.exists() or (a.out.is_dir() and not any(a.out.iterdir())),f'Output not empty: {a.out}')
 for p in [BASE,SRA,GEO,EXTRA,STUDIES]:fail_unless(digest(p)==EXPECTED[p.name],f'Input SHA mismatch {p}')
 orig_fields,base=read(BASE);_,sra=read(SRA);_,geo=read(GEO);_,extra=read(EXTRA);_,studies=read(STUDIES)
 fail_unless(len(extra)==1123 and len(studies)==176,'Preserved source-only or study evidence snapshot row count')
 fail_unless(len(base)==5241 and len({(r['study'],r['library_id']) for r in base})==5241,'V4 base invariant')
 fail_unless(len({r['study'] for r in base})==176,'Study invariant')
 fail_unless(len(sra)==3602 and len({(r['study'],r['library_id']) for r in sra})==3602,'Stage40 graph invariant')
 sramap={(r['study'],r['library_id']):r for r in sra}
 geo_map=defaultdict(lambda:defaultdict(set))
 for g in geo:geo_map[(g['study'],g['experiment_acc'])][g['metadata_field']].add(g['metadata_value'])
 indexed={(r['study'],r['library_id']):r for r in base}
 fail_unless(set(sramap)<=set(indexed),'Stage40 candidate outside base')
 output=[];affected=[];flags=[];counts=Counter();per_study=defaultdict(Counter)
 for old in base:
  r=dict(old);key=(r['study'],r['library_id']);s=r['study'];lib=r['library_id']
  sr=sramap.get(key);gs=geo_map.get(key,{});n=(sr.get('source_sample_title') or sr.get('source_sample_name') or '').strip() if sr else ''
  r.update({f:'' for f in NEW})
  r['v5_original_sra_sample_name']=sr.get('source_sample_name','') if sr else ''
  r['v5_original_sra_sample_title']=sr.get('source_sample_title','') if sr else ''
  r['v5_source_graph_match']='YES_STAGE40_SRX' if sr else 'NOT_IN_STAGE40_3602_GRAPH'
  r['v5_final_primary_approval']='NO_UNCHANGED'
  r['v5_remaining_gates']='ORIGINAL_GATES_UNCHANGED;NO_FINAL_PATIENT_INDEPENDENCE_OR_EXPRESSION_APPROVAL'
  if s in URL and (sr or (s=='SRP057793' and gs)):
   r['v5_source_evidence_url']=URL[s]
   r['v5_specimen_evidence_scope']='EXACT_SRX_LINKED_TO_ORIGINAL_SRA_OR_GEO_SAMPLE_LABEL_AND_COHORT_PUBLICATION'
   grade='STUDY_LABEL_ONLY_NOT_EXTERNAL_PERSON_ID'
   cls='UNRESOLVED_CLASS';action='HOLD_NO_NEW_SAMPLE_LEVEL_CONCLUSION';group='';region='';sub='';meta='NOT_ASSESSED'
   if s=='SRP064894':
    m=re.fullmatch(r'(\d+)T',n);fail_unless(m is not None,f'Unexpected esophageal sample {lib}: {n}')
    cls='TUMOUR_ESCC_SOURCE_LABEL';group=f'{s}:patient_token:{m[1]}'
    grade='WITHIN_STUDY_SRA_T_LABEL_AND_15_PAIRED_PATIENTS';action='TUMOUR_CANDIDATE_TREATMENT_UNKNOWN'
   elif s=='SRP124458':
    m=re.fullmatch(r'(HB\d+)-T',n);fail_unless(m is not None,f'Unexpected hepatoblastoma label {lib}: {n}')
    cls='HEPATOBLASTOMA_TUMOUR_SOURCE_LABEL_ONLY';sub='HEPATOBLASTOMA'
    grade='SOURCE_HB_SPECIMEN_LABEL_NOT_CERTIFIED_PATIENT_ID';group='' # Do not call HB specimen label a person ID.
    action='TUMOUR_LABEL_CANDIDATE_TREATMENT_AND_PATIENT_UNKNOWN'
   elif s=='SRP067524':
    m=re.fullmatch(r'(Epi|Myo|N)_(\d+)',n);fail_unless(m is not None,f'Unexpected salivary label {lib}: {n}')
    sub={'Epi':'EPITHELIAL_DOMINANT_WHOLE_ACC','Myo':'MYOEPITHELIAL_DOMINANT_WHOLE_ACC','N':'NORMAL_SALIVARY'}[m[1]]
    cls='NORMAL_SALIVARY_SOURCE_LABEL' if m[1]=='N' else 'PRIMARY_SALIVARY_ACC_SOURCE_LABEL'
    action='EXCLUDE_STRICT_TUMOUR_RETAIN_NORMAL_COMPARATOR' if m[1]=='N' else 'TUMOUR_CANDIDATE_TREATMENT_UNKNOWN'
    grade='STUDY_MORPHOLOGIC_SUBTYPE_LABEL_NOT_DISSOCIATED_CELL_POPULATION'
    # Epi_n and Myo_n numbers not proven comparable between groups or globally unique patients.
   elif s=='SRP119917':
    m=re.fullmatch(r'([HML]\d+)([TU])',n);fail_unless(m is not None,f'Unexpected prostate label {lib}: {n}')
    group=f'{s}:patient_token:{m[1]}';region=m[2]
    cls='PROSTATE_TUMOUR_PUNCH_LABEL' if m[2]=='T' else 'BENIGN_UNINVOLVED_PROSTATE_PUNCH_LABEL'
    action='RETAIN_TUMOUR_PUNCH_GROUPED_BY_PATIENT_OTHER_GATES_PENDING' if m[2]=='T' else 'EXCLUDE_STRICT_TUMOUR_RETAIN_BENIGN_COMPARATOR'
    grade='SOURCE_PAIRING_TOKEN_MATCHES_PUBLISHED_TUMOUR_AND_BENIGN_PUNCH_DESIGN'
   elif s=='SRP139591':
    m=re.fullmatch(r'(\d+)T-(\d+)',n);fail_unless(m is not None,f'Unexpected lung label {lib}: {n}')
    vals=gs.get('characteristics_ch1:patient',set());fail_unless(len(vals)==1,f'GEO patient ambiguity: {lib} {vals}')
    group=f'{s}:GEO_patient:{next(iter(vals))}';region=m[2]
    cls='NSCLC_MULTI_REGION_TUMOUR_SAMPLE';grade='DIRECT_GEO_PATIENT_FIELD_EXACT_SRX';action='RETAIN_MULTI_REGION_ONLY_WITH_PATIENT_GROUPING_TREATMENT_UNKNOWN'
   elif s=='SRP152822':
    m=re.fullmatch(r'patient([1-7])-([ABC])',n);fail_unless(m is not None,f'Unexpected gastric label {lib}: {n}')
    vals=gs.get('characteristics_ch1:diagnosis',set());fail_unless(vals=={'gastric cancer'},f'GEO diagnosis mismatch {lib} {vals}')
    group=f'{s}:GEO_patient:{m[1]}';region=m[2]
    cls='GASTRIC_TUMOUR_REGION_PURITY_VARIABLE';grade='DIRECT_GEO_SAMPLE_TITLE_AND_COHORT_PATIENT_DESIGN'
    action='RETAIN_MULTI_REGION_ONLY_WITH_PATIENT_GROUPING_TREATMENT_UNKNOWN'
   elif s=='SRP119676':
    m=re.fullmatch(r'(P\d+)(?:_([abc]))?_T(?:_(rep\d+))?',n);fail_unless(m is not None,f'Unexpected hepatoblastoma sample {lib}: {n}')
    vals=gs.get('characteristics_ch1:cell type',set());fail_unless(vals=={'tumour'},f'GEO cell type mismatch {lib} {vals}')
    group=f'{s}:GEO_patient_token:{m[1]}';region=((m[2] or '')+' '+(m[3] or '')).strip()
    cls='HEPATOBLASTOMA_TUMOUR_GEO_CONFIRMED';grade='DIRECT_GEO_TUMOUR_TISSUE_AND_PATIENT_PREFIX_REPLICATES'
    action='RETAIN_TUMOUR_ONLY_WITH_PATIENT_AND_REPLICATE_GROUPING_TREATMENT_UNKNOWN'
   elif s=='SRP055438':
    fail_unless('Normal colon acquired' in n and 'neuroendocrine tumor' in n,f'Unexpected normal colon label {lib}: {n}')
    vals=gs.get('characteristics_ch1:disease',set());fail_unless('Normal, SI Neuroendocrine tumor' in vals,f'SRP055438 disease mismatch {lib} {vals}')
    cls='NORMAL_COLON_NO_CANCER_IN_SAMPLED_TISSUE';grade='DIRECT_GEO_SAMPLE_TITLE_AND_DISEASE';action='EXCLUDE_STRICT_TUMOUR_RETAIN_NORMAL_COMPARATOR'
   elif s=='SRP057793':
    vals=gs.get('characteristics_ch1:metastasis',set());sv=gs.get('source_name',set())
    if vals:
     fail_unless(len(vals)==1 and vals <= {'Yes','No'},f'Metastasis tag inconsistent {lib}: {vals}')
     meta=next(iter(vals));r['v5_source_metastasis_annotation']='GEO_METASTASIS_'+meta.upper()+'_LESION_SITE_NOT_PROVEN'
     cls='SARCOMA_GEO_TUMOUR_CONTEXT_METASTASIS_CLINICAL_FLAG';grade='DIRECT_GEO_METASTASIS_FIELD_NOT_METASTATIC_SPECIMEN_SITE'
     action='HOLD_PRIMARY_SITE_METASTATIC_DISEASE_ANNOTATED' if meta=='Yes' else 'TUMOUR_CANDIDATE_METASTASIS_FLAG_NO_NOT_SITE_CERTIFIED'
    else: cls='SARCOMA_SOURCE_LABEL_NO_GEO_METASTASIS_FIELD';action='HOLD_MISSING_GEO_SITE_EVIDENCE';grade='SRA_NAME_ONLY'
    if n.startswith('N') and re.fullmatch(r'N\d{3}',n):
     sub=('N_PREFIX_GEO_SARCOMA_SOURCE_CONFIRMED_NOT_NORMAL' if sv and all('sarcoma' in v.lower() for v in sv) else 'N_PREFIX_AMBIGUOUS_NEVER_AUTO_NORMAL')
   r['v5_specimen_evidence_class']=cls;r['v5_proposed_clinical_action_NOT_APPROVAL']=action
   r['v5_tumour_subtype']=sub;r['v5_patient_group_within_study']=group
   r['v5_patient_group_evidence_grade']=grade;r['v5_specimen_region_or_replicate']=region
   if meta=='NOT_ASSESSED' and not r['v5_source_metastasis_annotation']:
    r['v5_source_metastasis_annotation']='NOT_ESTABLISHED'
   r['v5_remaining_gates']='TREATMENT_IF_NOT_EXPLICITLY_VERIFIED;TUMOUR_SITE_WHERE_AMBIGUOUS;GLOBAL_PATIENT_LINKAGE;TECHNICAL_QC;STUDY_OVERLAP;MODEL_TRANSFER;NO_FINAL_APPROVAL'
   if action.startswith('EXCLUDE_STRICT_TUMOUR'):
    r['v5_remaining_gates']='DEFINITIVE_NON_TUMOUR_FOR_STRICT_TUMOUR_SELECTION;RETAIN_AS_COMPARATOR_WITH_ASSAY_QC_AND_PATIENT_FLAGS'
   affected.append(r)
   per_study[s][cls]+=1
   counts[cls]+=1
  output.append(r)
 # Compute group sizes strictly from source labels in this review; never interpret as global patient identity.
 sizes=Counter(r['v5_patient_group_within_study'] for r in output if r['v5_patient_group_within_study'])
 for r in output:
  if r['v5_patient_group_within_study']: r['v5_study_group_size']=str(sizes[r['v5_patient_group_within_study']])
 primary=[r for r in output if r['stage41_partition_NOT_ADMISSION']=='NO_DETECTED_STUDY_ALIAS_OVERLAP_PENDING_GATES']
 non_tumour=[r for r in primary if r['v5_proposed_clinical_action_NOT_APPROVAL'].startswith('EXCLUDE_STRICT_TUMOUR')]
 flagged_metastatic=[r for r in primary if r['v5_proposed_clinical_action_NOT_APPROVAL']=='HOLD_PRIMARY_SITE_METASTATIC_DISEASE_ANNOTATED']
 # Assurance checks for known source-paper-vs-study members. Every test covers all study records in V4.
 expected_by_study={
   'SRP064894':Counter({'TUMOUR_ESCC_SOURCE_LABEL':15}),
   'SRP124458':Counter({'HEPATOBLASTOMA_TUMOUR_SOURCE_LABEL_ONLY':24}),
   'SRP067524':Counter({'PRIMARY_SALIVARY_ACC_SOURCE_LABEL':42,'NORMAL_SALIVARY_SOURCE_LABEL':5}),
   'SRP119917':Counter({'PROSTATE_TUMOUR_PUNCH_LABEL':24,'BENIGN_UNINVOLVED_PROSTATE_PUNCH_LABEL':20}),
   'SRP139591':Counter({'NSCLC_MULTI_REGION_TUMOUR_SAMPLE':45}),
   'SRP152822':Counter({'GASTRIC_TUMOUR_REGION_PURITY_VARIABLE':21}),
   'SRP119676':Counter({'HEPATOBLASTOMA_TUMOUR_GEO_CONFIRMED':30}),
   'SRP055438':Counter({'NORMAL_COLON_NO_CANCER_IN_SAMPLED_TISSUE':1}),
 }
 for s,wanted in expected_by_study.items():
  got=Counter(r['v5_specimen_evidence_class'] for r in output if r['study']==s)
  fail_unless(got==wanted,f'{s} sample-count mismatch: got {got}, expected {wanted}')
 sarcoma=Counter(r['v5_source_metastasis_annotation'] for r in output if r['study']=='SRP057793')
 fail_unless(sarcoma['GEO_METASTASIS_YES_LESION_SITE_NOT_PROVEN']==54 and sarcoma['GEO_METASTASIS_NO_LESION_SITE_NOT_PROVEN']==92,
             f'Sarcoma GEO evidence inconsistency {sarcoma}')
 fail_unless(len([r for r in output if r['study']=='SRP057793'])==149,'Sarcoma V4 149-row check')
 fail_unless(len({r['v5_patient_group_within_study'] for r in output if r['study']=='SRP139591'})==12,'Lung GEO 12 patients')
 fail_unless(len({r['v5_patient_group_within_study'] for r in output if r['study']=='SRP152822'})==7,'Gastric 7 patients')
 fail_unless(len({r['v5_patient_group_within_study'] for r in output if r['study']=='SRP119917'})==24,'Prostate 24 people')
 fail_unless(len({r['v5_patient_group_within_study'] for r in output if r['study']=='SRP119676'})==24,'HB 24 source grouping tokens')
 fail_unless(len(primary)==2055,'primary administrative universe must remain untouched')
 fail_unless(len(non_tumour)==21,f'Non tumour count changed: {len(non_tumour)}')
 fail_unless(len(flagged_metastatic)==52,f'Sarcoma flag count changed: {len(flagged_metastatic)}')
 fail_unless(all(not r['v5_patient_group_within_study'] for r in output if r['study']=='SRP124458'),'HB specimen id must not become patient')
 fail_unless(all(r['FINAL_CLINICAL_APPROVAL']!='APPROVED' for r in output),'Source clinically approved unexpectedly')
 fail_unless(all(r['v5_final_primary_approval']=='NO_UNCHANGED' for r in output),'V5 must not approve')
 # Original 107 V4 fields must remain strictly unchanged by the supplement.
 for old,new in zip(base,output):fail_unless(all(old[f]==new[f] for f in orig_fields),'V4 original field modified')
 # Flag other ambiguous names; do not reinterpret globally.
 for r in primary:
  if not r['v5_specimen_evidence_class'] and r['v5_original_sra_sample_title']=='' and r['v5_original_sra_sample_name']=='':continue
  n=r['v5_original_sra_sample_title'] or r['v5_original_sra_sample_name']
  if r['study'] not in URL:
   if re.search(r'(?i)(^|[_\-\s])(normal|adjacent|benign|non.?tumou?r|metastasis|post.?treatment|relapse)(?=$|[_\-\s])',n):
    flags.append({'study':r['study'],'library_id':r['library_id'],'source_sample_text':n,'review_reason':'SUSPICIOUS_TEXT_REQUIRE_STUDY_SPECIFIC_INTERPRETATION','stage41_partition':r['stage41_partition_NOT_ADMISSION']})
 # Output
 a.out.mkdir(parents=True,exist_ok=True)
 save(a.out/'01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_CROSSWALK_UNSIGNED.tsv',orig_fields+NEW,output)
 save(a.out/'02_376_FOCUSED_ACCESSIONS_SPECIMEN_PATIENT_EVIDENCE.tsv',orig_fields+NEW,affected)
 save(a.out/'03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv',orig_fields+NEW,primary)
 save(a.out/'04_21_SOURCE_SUPPORTED_NON_TUMOUR_IN_PRIMARY_SCREEN.tsv',orig_fields+NEW,non_tumour)
 save(a.out/'05_52_SARCOMA_GEO_METASTATIC_DISEASE_FLAGS_IN_PRIMARY_SCREEN.tsv',orig_fields+NEW,flagged_metastatic)
 save(a.out/'06_OTHER_SAMPLE_NAME_REVIEW_QUEUE.tsv',['study','library_id','source_sample_text','review_reason','stage41_partition'],flags)
 study_summary=[]
 for s in sorted(URL):
  seen=[r for r in output if r['study']==s];src=[r for r in seen if r['v5_specimen_evidence_class']]
  study_summary.append({'study':s,'n_original':len(seen),'n_newly_crosswalked':len(src),'n_primary_administrative':sum(r['stage41_partition_NOT_ADMISSION']=='NO_DETECTED_STUDY_ALIAS_OVERLAP_PENDING_GATES' for r in seen),
   'source_url':URL[s],'specimen_evidence_classes':json.dumps(dict(per_study[s]),sort_keys=True),
   'n_exact_patient_group_tokens':len({r['v5_patient_group_within_study'] for r in src if r['v5_patient_group_within_study']}),
   'scientific_limit':'SOURCE_SPECIMEN_AND_WITHIN_STUDY_GROUPS_ONLY;NO_GLOBAL_PATIENT_INDEPENDENCE_OR_PRETREATMENT_CERTIFICATION'})
 save(a.out/'07_NINE_STUDY_SOURCE_CROSSWALK_SUMMARY.tsv',list(study_summary[0]),study_summary)
 shutil.copyfile(EXTRA,a.out/'10_1123_SOURCE_ONLY_EXPERIMENTS_UNCHANGED_FROM_V4.tsv')
 shutil.copyfile(STUDIES,a.out/'11_176_STUDY_SOURCE_LEDGER_UNCHANGED_FROM_V4.tsv')
 report={'status':'V5_SOURCE_SPECIMEN_AND_PATIENT_GROUP_CROSSWALK_NOT_FINAL_COHORT_LOCK',
  'n_original':len(base),'n_source_only_held':len(extra),'n_studies':len(set(r['study'] for r in output)),'n_stage40_graph':len(sra),'n_affected_source_crosswalk':len(affected),
  'source_backed_specimen_class_counts':dict(counts),'primary_screen_original':len(primary),
  'primary_screen_non_tumour_found':len(non_tumour),'primary_screen_non_tumour_by_study':dict(Counter(r['study'] for r in non_tumour)),
  'remaining_primary_screen_after_source_non_tumour_restriction_NOT_APPROVED':len(primary)-len(non_tumour),
  'sarcoma_GEO_metastatic_disease_flags_primary_screen_NOT_PROVEN_METASTATIC_SAMPLED_LESIONS':len(flagged_metastatic),
  'salivary_normal_vs_ACC':'5 normals, 42 ACC tumours by Epi/Myo morphology, not purified cells',
  'prostate_study_tissue':'24 T vs 20 U, only 24 patient tokens across 44 libraries',
  'lung_patient_groups':'45 multiregion libraries in 12 GEO patient groups',
  'gastric_patient_groups':'21 gastric tumour locations in 7 GEO patient groups',
  'hepatoblastoma_patient_groups':'30 selected tumour libraries in 24 patient tokens, some regions/technical replicates',
  'source_gene_clinical_reviewer_queue_flags':len(flags),'remaining_gates':'STUDY_LEVEL_REVIEW_ALL_176;PATIENT_INDEPENDENCE;TREATMENT;TUMOUR_SAMPLE_SITE;TECHNICAL_QC;MODEL_TRANSFER',
  'input_sha256':EXPECTED}
 (a.out/'08_V5_QA.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 with (a.out/'09_OUTPUT_SHA256SUMS.txt').open('w') as f:
  for path in sorted(a.out.iterdir()):
   if path.name=='09_OUTPUT_SHA256SUMS.txt':continue
   f.write(f'{digest(path)}  {path.name}\n')
 print(json.dumps(report,indent=2,sort_keys=True))
 print('NO CLINICAL APPROVAL. ORIGINAL V4 COLUMNS INTACT. NO SOURCE DATA MODIFIED.')

if __name__=='__main__':main()
