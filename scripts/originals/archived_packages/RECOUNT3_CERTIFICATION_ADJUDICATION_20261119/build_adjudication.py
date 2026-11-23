#!/usr/bin/env python3
"""Build a NON-DECISIONAL, accession-grounded recount3 certification review.

No existing manifest or model is modified. Every output remains HOLD pending an
identified reviewer and accession-level evidence. No expression analysis is run.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path

VERSION = 'recount3_accession_evidence_adjudication_v1_0_0'

def rd(p):
    with open(p, newline='', encoding='utf-8') as f: return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(p, rows, fields):
    with open(p, 'w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, delimiter='\t', fieldnames=fields, lineterminator='\n', extrasaction='ignore')
        w.writeheader();w.writerows(rows)

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4<<20), b''):h.update(b)
    return h.hexdigest()

def cleaned(x):return ' '.join(str(x or '').split())

def norm(x):return cleaned(x).casefold()

PRIMARY_TOKENS=[
    'primary breast tumor tissue','primary invasive ovarian carcinoma',
    'primary ovarian tumor tissue','primary tumor tissue from high grade serous ovarian tumors',
    'primary prostate cancer','primary melanoma','primary breast tumor',
    'primary invasive breast carcinoma','frozen primary tumor tissue',
    'primary tumor','primary crc total rna','primary hcc tissue',
    'primary colon cancer','primary gallbladder cancer tissue',
    'melanoma primary tumor','primary prostate tumour',
    'nasopharyngeal carcinoma – primary tumor',
]

def primary_source_fields(attrs):
    eligible=[]
    for k, values in attrs.items():
        if k in ('history','disease stage','diagnosis','study','study title','study_title','study_description'):
            continue
        if k in ('isolate','source_name','tissue','tissue nature','rna source','sample_type','tumor type','description','biopsy_type'):
            for val in values:
                low=norm(val)
                if any(token in low for token in PRIMARY_TOKENS) or (k=='tumor type' and low=='primary'):
                    eligible.append(f'{k}={cleaned(val)}')
    return sorted(set(eligible))

def patient_source(row):
    obj=json.loads(row['source_patient_token_candidates_json'])
    return obj

def create(source,triage,runqueue,groups,featuremap,out):
    src=rd(source);tr=rd(triage);runs=rd(runqueue);repeat=rd(groups);w=rd(featuremap)
    if len(src)!=3602 or len(tr)!=3602 or len(runs)!=409 or len(w)!=34:
        raise ValueError('Unexpected counts; REFUSE to adjudicate mismatched input evidence')
    ids={(x['study'],x['library_id']) for x in src}
    if len(ids)!=3602:raise ValueError('Duplicate source library IDs')
    tmap={(x['study'],x['library_id']):x for x in tr}
    if set(tmap)!=ids:raise ValueError('Source and prior triage do not identify the same libraries')
    if len({x['run_id'] for x in runs}) !=409:raise ValueError('Run review duplicates')
    if len({x['gene_symbol'] for x in w}) !=34:raise ValueError('W-gene review duplicates')
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    study=collections.defaultdict(lambda:collections.Counter())
    record=[];patient=[];reviewpriority=[]
    aliases=collections.defaultdict(list)
    direct_map=collections.defaultdict(set);reverse_map=collections.defaultdict(set)
    patients_perstudy=collections.defaultdict(set)
    for a in sorted(src,key=lambda z:(z['study'],z['library_id'])):
        key=(a['study'],a['library_id']); b=tmap[key]
        tags=set(filter(None,b['source_review_tags'].split('|')))
        attrs=json.loads(a['parsed_sample_attributes_json']);direct=patient_source(a)
        primary=primary_source_fields(attrs)
        # This is specimen-specific source-submitted metadata, never text about
        # the wider study. Disease-subtype terminology like "primary GBM" in
        # the history field is not sufficient for untreated/primary-site inference.
        non_tumour=('SOURCE_EXPLICIT_NON_TUMOUR' in tags)
        quantseq='ASSAY_3PRIME_QUANTSEQ_INCOMPATIBILITY_REVIEW' in tags
        temporal='TREATMENT_EXPOSURE_CONFOUNDED' in tags or 'PRIOR_THERAPY_BEFORE_BEVA' in tags
        multiregion='MULTIREGION_SAME_PATIENT' in tags
        special='EPITHELIAL_SOURCE_WHOLE_TUMOUR_REVIEW' in tags or 'SALIVARY_EPI_MYO_COMPARTMENT_CONTEXT_REVIEW' in tags
        primary_flag=bool(primary)
        if non_tumour:category='SOURCE_EXPLICIT_NON_TUMOUR'
        elif quantseq:category='PROTOCOL_COMPATIBILITY_HOLD'
        elif temporal:category='KNOWN_TREATMENT_CONTEXT_REVIEW'
        elif multiregion:category='MULTIREGION_DESIGN_PATIENT_AWARE'
        elif special:category='WHOLE_TUMOUR_COMPOSITION_REVIEW'
        elif primary_flag:category='DIRECT_PRIMARY_ORIGIN_METADATA'
        elif any(t in tags for t in ['HEPATOBLASTOMA_TUMOUR_SOURCE_DOCUMENTED','SARCOMA_FROZEN_TUMOUR_SOURCE_DOCUMENTED','PROSTATE_TUMOUR_SOURCE']):
            category='DIRECT_TUMOUR_METADATA_PRIMARY_UNPROVEN'
        else:category='SPECIMEN_ELIGIBILITY_UNRESOLVED'
        # Findings are evidence-based proposals only; no clinician-level claims
        # of treatment naivety or independent participants are inferred.
        proposed=(
          'CONFIRM_EXCLUDE_FROM_PRIMARY_TUMOUR' if non_tumour else
          'ASSAY_PROTOCOL_REVIEW_BEFORE_ANY_MODEL_INPUT' if quantseq else
          'REVIEW_TREATED_RECURRENT_CONTEXT' if temporal else
          'REVIEW_PATIENT_AWARE_MULTIREGION' if multiregion else
          'REVIEW_COMPARTMENT_OR_ENRICHMENT' if special else
          'VERIFY_PRIMARY_TUMOUR_AND_TREATMENT_STATUS' if primary_flag else
          'VERIFY_TUMOUR_CLASS_AND_PATIENT' )
        sample=str(a['source_sample_acc'])
        proof='; '.join(primary) if primary else (b['evidence_fields'] or '')
        row={
            'study':a['study'],'library_id':a['library_id'],'source_SRA_sample':sample,
            'existing_cancer_label':a['cancer_from_existing_curation'],
            'curated_patient_key':a['existing_patient_key'],
            'curated_patient_basis':a['existing_patient_key_basis'],
            'evidence_category_NOT_CERTIFIED':category,
            'source_explicit_primary_phrase':('; '.join(primary)),
            'source_review_tags':b['source_review_tags'],
            'source_sample_evidence_excerpt':cleaned(proof)[:650],
            'source_metadata_path':a['metadata_source_path'],
            'source_metadata_sha256':a['metadata_source_sha256'],
            'source_accession_evidence_pointer':f"{a['metadata_source_path']}#experiment_acc={a['library_id']}&sample_acc={sample}",
            'review_priority':('CRITICAL' if non_tumour or quantseq else 'HIGH' if temporal or multiregion or special else 'NORMAL'),
            'review_action_PROPOSAL_ONLY':proposed,
            'study_source_attestation':'PENDING_REVIEWER',
            'specimen_reviewer':'','clinical_specimen_release':'HOLD',
            'patient_group_verified':'NO',
            'ARCHS4_independence':'NOT_CERTIFIED',
            'frozen_W_transfer':'NOT_CERTIFIED'
        }
        record.append(row)
        s=study[a['study']];s['libraries']+=1;s[category]+=1
        if primary_flag:s['primary_metadata']+=1
        if not a['source_sample_attributes_raw'].strip():s['missing_sample_attrs']+=1
        if 'BioSample proxy' in a['existing_patient_key_basis']:s['proxy_libraries']+=1
        if b['source_review_tags']:s['prior_triage_alerts']+=1
        # Direct patient tokens are preserved, never interpreted as global IDs.
        if direct:
            if len(direct)!=1:
                raise ValueError('Unexpected multi-field patient token candidates: '+a['library_id'])
            field, values=next(iter(direct.items()))
            if len(values)!=1:raise ValueError('Unexpected patient token multiplicity: '+a['library_id'])
            token=cleaned(values[0]);tok=norm(token)
            direct_map[(a['study'],field,tok)].add(a['existing_patient_key'])
            reverse_map[(a['study'],a['existing_patient_key'],field)].add(tok)
            patients_perstudy[a['study']].add(tok)
            patient.append({
                'study':a['study'],'library_id':a['library_id'],'source_SRA_sample':sample,
                'attribute_name':field,'source_patient_token':token,
                'curated_patient_key':a['existing_patient_key'],
                'source_pointer':row['source_accession_evidence_pointer'],
                'source_sha256':row['source_metadata_sha256'],
                'token_conflicts_with_curated_patient_group':'PENDING_COMPUTATION',
                'curated_group_has_multiple_direct_tokens':'PENDING_COMPUTATION',
                'patient_id_reviewer':'','patient_independence_decision':'HOLD'
            })
        # Additional candidate relationships from source sample labels
        isolate=(attrs.get('isolate') or [''])[0]
        if a['study']=='SRP067524':
            m=re.fullmatch(r'(Epi|Myo)_(\d+)',isolate,flags=re.I)
            if m:aliases[(a['study'],'EPI_MYO_SUFFIX',m.group(2))].append((a['library_id'],isolate))
        if a['study']=='SRP158272':
            m=re.match(r'^(RR-\d+)',isolate,flags=re.I)
            if m:aliases[(a['study'],'RR_PREFIX',norm(m.group(1)))].append((a['library_id'],isolate))
    for p in patient:
        k=(p['study'],p['attribute_name'],norm(p['source_patient_token']))
        k2=(p['study'],p['curated_patient_key'],p['attribute_name'])
        p['token_conflicts_with_curated_patient_group']='YES' if len(direct_map[k])>1 else 'NO'
        p['curated_group_has_multiple_direct_tokens']='YES' if len(reverse_map[k2])>1 else 'NO'
    cross_token_conflicts=sum(len(v)>1 for v in direct_map.values())
    reverse_token_conflicts=sum(len(v)>1 for v in reverse_map.values())
    candidates=[]
    for (s,kind,root),rr in sorted(aliases.items()):
        if len(rr)>1:
            candidates.append({'study':s,'possible_relation_type':kind,'source_token':root,
              'n_libraries':len(rr),'library_ids':';'.join(x[0] for x in rr),
              'source_specimen_tokens':';'.join(x[1] for x in rr),
              'claim':'POSSIBLE_RELATED_SAMPLES_NOT_PATIENT_IDENTITY',
              'review_status':'HOLD_REQUIRES_ORIGINAL_STUDY_SAMPLE_MAP'})
    stud=[]
    for s,c in sorted(study.items(),key=lambda kv:(-(kv[1]['SOURCE_EXPLICIT_NON_TUMOUR']*100+kv[1]['PROTOCOL_COMPATIBILITY_HOLD']*100+kv[1]['KNOWN_TREATMENT_CONTEXT_REVIEW']*2+kv[1]['MULTIREGION_DESIGN_PATIENT_AWARE']*2+kv[1]['WHOLE_TUMOUR_COMPOSITION_REVIEW']*2+kv[1]['libraries']/5),kv[0])):
        stud.append({'study':s,'n_candidate_libraries':c['libraries'],
            'source_explicit_non_tumour':c['SOURCE_EXPLICIT_NON_TUMOUR'],
            'nonstandard_assay_holds':c['PROTOCOL_COMPATIBILITY_HOLD'],
            'treatment_context_review':c['KNOWN_TREATMENT_CONTEXT_REVIEW'],
            'multiregion_patient_aware_review':c['MULTIREGION_DESIGN_PATIENT_AWARE'],
            'whole_tumour_composition_review':c['WHOLE_TUMOUR_COMPOSITION_REVIEW'],
            'direct_primary_origin_source_phrases':c['primary_metadata'],
            'no_direct_primary_origin_phrase':c['libraries']-c['primary_metadata'],
            'biosample_proxy_libraries':c['proxy_libraries'],
            'direct_source_patient_tokens':sum(x['study']==s for x in patient),
            'source_attribute_missing':c['missing_sample_attrs'],
            'study_attestation_status':'PENDING_ORIGINAL_SOURCE_AND_REVIEWER',
            'next_step':'REVIEW_ACCESSION_SOURCE_AND_ORIGINAL_STUDY_SAMPLES'})
    wr=[{
        'gene_symbol':r['gene_symbol'],'hold_reason':r['hold_reason'],
        'G026_candidate_resolution':'NOT_RESOLVED',
        'ARCHS4_H5_candidate_resolution':'NOT_RESOLVED',
        'required_evidence':'GENCODE_G026_GENE_ROWS_AND_ARCHS4_H5_GENE_IDENTITY;FROZEN_W_FEATURE_DEFINITION',
        'strict_prohibition':'NO_DROP_NO_ZERO_FILL_NO_WEIGHT_CHANGE',
        'reviewer':'','decision':'HOLD'
    } for r in w]
    ref=[
      {'study':'SRP157974','external_source':'https://www.ncbi.nlm.nih.gov/sra/SRX5322041','source_kind':'NCBI_SRA_EXPERIMENT_AND_PROJECT','establishes':'FUSCC primary TNBC project; source-backed context, not automatic all-sample attestation','status':'PUBLIC_CONTEXT_NOT_ATTESTED'},
      {'study':'SRP114518','external_source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102118','source_kind':'NCBI_GEO_STUDY','establishes':'recurrent platinum-resistant ovarian trial; tumour biopsy versus ascites and timepoints need accession review','status':'PUBLIC_CONTEXT_NOT_ATTESTED'},
      {'study':'SRP139591','external_source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112996','source_kind':'NCBI_GEO_STUDY','establishes':'multi-locus lung tumour sampling; individual patient and region map in SRA attributes','status':'PUBLIC_CONTEXT_NOT_ATTESTED'},
      {'study':'SRP194626','external_source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130469','source_kind':'NCBI_GEO_STUDY','establishes':'method comparison includes Lexogen QuantSeq; gene-body compatibility hold','status':'PUBLIC_CONTEXT_NOT_ATTESTED'},
    ]
    write_tsv(out/'01_accession_specimen_review_3602.tsv',record,list(record[0]))
    write_tsv(out/'02_study_adjudication_queue_126.tsv',stud,list(stud[0]))
    write_tsv(out/'03_direct_patient_token_evidence_281.tsv',patient,list(patient[0]))
    write_tsv(out/'04_additional_possible_patient_relationships.tsv',candidates,list(candidates[0]))
    write_tsv(out/'05_run_hold_queue_409.tsv',runs,list(runs[0]))
    write_tsv(out/'06_held_W_feature_gates_34.tsv',wr,list(wr[0]))
    write_tsv(out/'07_public_source_registry_REVIEW_ONLY.tsv',ref,list(ref[0]))
    write_tsv(out/'08_repeated_curated_subject_groups_109.tsv',repeat,list(repeat[0]))
    ctr=collections.Counter(x['evidence_category_NOT_CERTIFIED'] for x in record)
    summary={
      'version':VERSION,
      'status':'EVIDENCE_ADJUDICATION_PREPARED_NO_RELEASE_DECISIONS',
      'n_libraries':len(record),'n_studies':len(stud),
      'n_source_primary_origin_phrases':sum(bool(x['source_explicit_primary_phrase']) for x in record),
      'evidence_category_distribution':dict(ctr),
      'n_source_explicit_nontumour':ctr['SOURCE_EXPLICIT_NON_TUMOUR'],
      'n_assay_incompatibility_holds':ctr['PROTOCOL_COMPATIBILITY_HOLD'],
      'n_direct_source_patient_tokens':len(patient),
      'n_studies_with_direct_patient_tokens':len(patients_perstudy),
      'n_source_patient_tokens_mapping_to_multiple_curated_groups':cross_token_conflicts,
      'n_curated_groups_with_multiple_direct_source_tokens':reverse_token_conflicts,
      'n_possible_additional_relationship_groups':len(candidates),
      'n_run_holds':len(runs),
      'n_repeated_curated_patient_groups':len(repeat),
      'n_W_feature_holds':len(w),
      'crossatlas_alias_coverage':'740_OF_9465_ONLY_FROM_PREVIOUS_AUDIT',
      'frozen_model_score_transfer':'NOT_PERFORMED',
      'n_certified_libraries':0,
      'scope':'source evidence triage; no execution of release assessor or expression projection',
      'source_input_sha256':{str(p.name):sha256(p) for p in (source,triage,runqueue,groups,featuremap)},
    }
    (out/'09_adjudication_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    readme='''RECOUNT3 ACCESSION EVIDENCE / PATIENT ADJUDICATION V1.0.0

THIS IS AN AUDITED REVIEW QUEUE, NOT A COHORT LOCK OR MODEL PASS.
Compiled from the user-supplied exact SRA source evidence and previous scientific review.
No existing cluster release, run decision, patient decision, or frozen W feature has been edited.

RELEASE LOGIC
1. Study-level source design is NOT evidence for every accession.
2. Direct source attribute explicitly saying "primary [cancer/tumour]" supports primary anatomical origin,
   not treatment naivety, whole tumour composition, a unique patient or compatibility with frozen W.
3. Sample source 'normal tissue' and QuantSeq flags merit review for exclusion; they are not
   automatically imported into cert_v101 or assigned an invented reviewer.
4. Disease-history term 'primary glioblastoma' may distinguish de novo from secondary disease;
   it is not evidence that every specimen is untreated, primary-site bulk tumour.
5. Only a source attribute explicitly naming the patient is treated as a direct patient token;
   arbitrary "isolate", BioSample and accession IDs are never assumed globally independent patients.
6. Strong source patient token consistency within a study does not prove cross-study independence.
7. Two fields point to distinct blocks: 34 gene-mapping holds and incomplete ARCHS4 alias coverage.
   Neither can be passed using the 465-gene same-source correlation diagnostic.

OUTPUTS
01: exactly 3602 accession-level evidence-based review proposals (all HOLD).
02: 126-study priority table with evidence coverage and review burden.
03: 281 direct source patient tokens, with curation concordance flags.
04: candidate Epi/Myo and RR-root relationships, explicitly not called patient identities.
05: 409 run QC holds copied unaltered from the triage ledger.
06: 34 frozen model features requiring symbol/Ensembl/H5 row resolution.
07: small public study reference registry (not formal attestations).
08: 109 repeated curated patient groups, source verification pending.
09: deterministic accounting and input SHA256 hashes.

WORKFLOW BEFORE ANY NEW STATISTICS
(a) Confirm source annotation and original study publication / GEO supplement for each
    library; reviewer must explicitly sign study and accession-level decisions.
(b) Resolve or exclude unknown patients without inferring identity from BioSample proxies.
(c) Resolve all run holds and final aggregation (only accepted RAW G026 coverage summed).
(d) Build complete 9465-atlas alias overlap audit, do not infer independence from missing aliases.
(e) Resolve all 34 W genes and compare actual full frozen-model output on paired specimens.
(f) Only then use cert_v101 assessor and establish a controlled release. Its current
    model-transfer gate is explicitly not a full automated PASS implementation.

NO FILE IN THIS PACKAGE SHOULD BE COPIED OVER cert_v101/review/*.tsv AS AN APPROVAL.
'''
    (out/'README_REVIEW_FIRST.txt').write_text(readme)
    names=[p for p in sorted(out.iterdir()) if p.is_file()]
    with (out/'SHA256SUMS.txt').open('w') as f:
        for p in names:f.write(f'{sha256(p)}  {p.name}\n')
    return summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('source','triage','runqueue','groups','featuremap','out'):
        p.add_argument('--'+k,required=True,type=Path)
    a=p.parse_args()
    print(json.dumps(create(a.source,a.triage,a.runqueue,a.groups,a.featuremap,a.out),indent=2))
if __name__=='__main__':main()
