#!/usr/bin/env python3
"""Recount3 certification v1.0.1: explicit accession-level eligibility and strict release gates.

Standalone read-only import of exact v1.0.0 validation implementation. NEVER changes
v1.0.0 source, work, evidence lock, frozen model, or extraction data.
Commands: preflight, init, assess. `init` creates a NEW work directory.
The source-review attestations are evidence annotations by a researcher, not an API
verification of external provenance. Missing information always produces HOLD.
"""
from __future__ import annotations
import argparse, collections, csv, datetime, hashlib, json, math, os, shutil, sys
from pathlib import Path
import _v100_reference_readonly as v100

VERSION='recount3_final_certification_v1_0_1'
SPECIMEN_CLASSES={'PRIMARY_TUMOUR','RECURRENT_TUMOUR','METASTATIC_TUMOUR',
                  'MALIGNANT_ASCITES','ADJACENT_NORMAL','NORMAL_TISSUE',
                  'CELL_LINE','OTHER','UNRESOLVED'}
TREATMENT={'PRETREATMENT','TREATMENT_NAIVE_CONFIRMED','ON_TREATMENT',
           'POST_TREATMENT','UNKNOWN','NOT_APPLICABLE'}
ELIGIBILITY={'INCLUDE_PRIMARY','EXCLUDE_NONPRIMARY','HOLD'}
PATIENT_RELATIONS={'REPEATED_BIOLOGICAL_SPECIMENS','LONGITUDINAL','MULTIREGION',
                   'MULTIPLE_TUMOURS','TECHNICAL_DUPLICATE_RESOLVED'}
SPEC_FIELDS=['library_id','study','biosample','cancer','original_library_selection',
             'original_preservation','original_cancer_variant','original_source_row_sha256',
             'reviewed_specimen_id','specimen_class','treatment_status','disease_context',
             'collection_timepoint','anatomic_site','tumour_region','primary_cohort_decision',
             'accession_source_reference','accession_evidence_record','source_data_sha256',
             'reviewer','review_note']
IMMUTABLE_SPEC_FIELDS=SPEC_FIELDS[:8]

def canonical(path):return Path(path).expanduser().resolve()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def nonempty(row,*fields):return all(str(row.get(x,'')).strip() for x in fields)
def sha256(p):return v100.digest(p)
def safe_reference(s):
    """Reject obvious placeholders, not a substitute for reviewing source content."""
    s=str(s).strip();bad={'NA','N/A','UNKNOWN','NONE','TODO','TBD','UNVERIFIED','NO','YES','SAME','-'}
    return len(s)>=8 and s.upper() not in bad

def check_program(root):
    user_code=Path(root)/'final_cert/certify.py'
    packaged=Path(__file__).resolve().parent/'_v100_reference_readonly.py'
    if not user_code.is_file():raise ValueError('Original v1.0.0 certify.py not found: '+str(user_code))
    if sha256(user_code)!=sha256(packaged):
        raise ValueError('Installed v1.0.0 code differs from pinned, reviewed reference; investigate before certification')

def check_existing_lock(previous):
    lock=Path(previous)/'TECHNICAL_EVIDENCE_LOCK_v1_0_0'
    needed=['LOCK_MANIFEST.json','SHA256SUMS.txt','source_evidence_sha256.tsv']
    if any(not (lock/n).is_file() for n in needed):raise ValueError('v1.0.0 technical evidence lock missing or incomplete')
    for line in (lock/'SHA256SUMS.txt').read_text().splitlines():
        h,n=line.split(None,1)
        if not (lock/n.strip()).is_file() or sha256(lock/n.strip())!=h:
            raise ValueError('v1.0.0 technical evidence lock checksum mismatch: '+n)
    state=json.loads((lock/'LOCK_MANIFEST.json').read_text())
    if state.get('status')!='EVIDENCE_HASH_LOCK_COMPLETE_NOT_BIOLOGICAL_COHORT_CERTIFICATION':
        raise ValueError('Unexpected v1.0.0 lock status')
    return lock

def preflight(root,previous):
    root=canonical(root);previous=canonical(previous)
    check_program(root);lock=check_existing_lock(previous)
    data,paths,reports,summary=v100.inspect(root,quiet=True)
    mandatory=['run_decisions','patient_identity','repeated_subject_relationships','study_source_attestations',
               'held_W_feature_resolutions']
    for name in mandatory:
        if not (previous/'review'/f'{name}.tsv').is_file():raise ValueError('Missing v1.0.0 reviewed input '+name)
    v100.check_reviews(data,previous)
    libs=v100.readtsv(paths['libs'])
    if len(libs)!=3639 or len([r for r in libs if r['study'] not in v100.DGE])!=3602:raise ValueError('Library selection changed')
    print(json.dumps({'status':'V101_PREFLIGHT_PASS_REVIEW_ONLY',
       'previous_evidence_lock':str(lock),'n_studies':129,'n_runs':4166,
       'n_protocol_exclusions':37,'n_primary_candidates_requiring_specimen_review':3602,
       'existing_archs4_independence_certified':False,'model_transfer_certified':False},indent=2))
    return data,paths

def init(root,previous,work):
    root,previous,work=map(canonical,(root,previous,work))
    if work.exists():raise ValueError('Refuse to overwrite existing work: '+str(work))
    if any(work==x or work in x.parents or x in work.parents for x in (root,previous)):
        raise ValueError('v1.0.1 work must be a separate sibling directory')
    data,paths=preflight(root,previous)
    libs=v100.readtsv(paths['libs']);source={r['library_id']:r for r in libs if r['study'] not in v100.DGE}
    candidate=v100.readtsv(previous/'draft/library_candidates.tsv')
    if {r['library_id'] for r in candidate}!=set(source):raise ValueError('Candidate IDs differ between v1.0.0 and manifests')
    temp=work.with_name(work.name+'.partial')
    if temp.exists():raise ValueError('Stale temporary work directory: '+str(temp))
    temp.mkdir(parents=True)
    try:
        shutil.copytree(previous/'draft',temp/'draft')
        shutil.copytree(previous/'review',temp/'review')
        # Do not replicate or edit v1.0.0 sealed TECHNICAL_EVIDENCE_LOCK.
        specimen=[]
        for lid,r in sorted(source.items()):
            specimen.append({'library_id':lid,'study':r['study'],'biosample':r['biosample'],
              'cancer':r['cancer'],'original_library_selection':r.get('library_selection',''),
              'original_preservation':r.get('preservation',''),
              'original_cancer_variant':r.get('cancer_variant',''),
              'original_source_row_sha256':r.get('SOURCE_ROW_SHA256',''),
              'reviewed_specimen_id':'','specimen_class':'UNRESOLVED',
              'treatment_status':'UNKNOWN','disease_context':'','collection_timepoint':'',
              'anatomic_site':'','tumour_region':'','primary_cohort_decision':'HOLD',
              'accession_source_reference':'','accession_evidence_record':'',
              'source_data_sha256':'','reviewer':'','review_note':''})
        v100.writetsv(temp/'draft/specimen_eligibility_3602.tsv',specimen,SPEC_FIELDS)
        v100.writetsv(temp/'review/specimen_eligibility_3602.tsv',specimen,SPEC_FIELDS)
        # Audit linked source evidence and review initial state, but no auto decisions.
        source_entries=[{'role':'V100_SOURCE_EVIDENCE_LOCK','path':str(previous/'TECHNICAL_EVIDENCE_LOCK_v1_0_0/LOCK_MANIFEST.json'),
                         'sha256':sha256(previous/'TECHNICAL_EVIDENCE_LOCK_v1_0_0/LOCK_MANIFEST.json')},
                        {'role':'V100_REVIEW_SEED','path':str(previous/'draft/evidence_sha256.tsv'),
                         'sha256':sha256(previous/'draft/evidence_sha256.tsv')},
                        {'role':'V100_CERTIFICATION_CODE','path':str(root/'final_cert/certify.py'),
                         'sha256':sha256(root/'final_cert/certify.py')}]
        v100.writetsv(temp/'draft/v101_input_provenance.tsv',source_entries,['role','path','sha256'])
        v100.savejson(temp/'draft/v101_work_manifest.json',{'version':VERSION,'created_utc':utc(),
          'previous_work':str(previous),'source_root':str(root),
          'status':'UNREVIEWED_HOLD','n_specimen_templates':3602,
          'release_policy':'No study-level YES, sample BioSample, or unvalidated W mapping alone certifies an analysis-ready cohort'})
        os.replace(temp,work)
    except Exception:
        if temp.exists():shutil.rmtree(temp)
        raise
    print(json.dumps({'status':'V101_REVIEW_TEMPLATES_CREATED_NOT_CERTIFIED',
      'new_work':str(work),'n_specimen_eligibility_rows':3602,
      'edit':'review/specimen_eligibility_3602.tsv, review/study_source_attestations.tsv, review/patient_identity.tsv, review/repeated_subject_relationships.tsv and review/run_decisions.tsv',
      'original_evidence_lock_unchanged':True},indent=2))

def validated_specimen(rows,immutable):
    if len(rows)!=len(immutable) or len({r['library_id'] for r in rows})!=len(immutable):
        raise ValueError('Specimen review must contain exactly 3602 unique libraries')
    if set(x['library_id'] for x in rows)!=set(immutable):raise ValueError('Specimen library ID set altered')
    for r in rows:
        old=immutable[r['library_id']]
        for c in IMMUTABLE_SPEC_FIELDS:
            if r[c]!=old[c]:raise ValueError('Immutable specimen field altered: '+r['library_id']+'/'+c)
        if r['specimen_class'] not in SPECIMEN_CLASSES:raise ValueError('Invalid specimen_class: '+str(r['specimen_class']))
        if r['treatment_status'] not in TREATMENT:raise ValueError('Invalid treatment_status: '+str(r['treatment_status']))
        if r['primary_cohort_decision'] not in ELIGIBILITY:raise ValueError('Invalid primary_cohort_decision: '+str(r['primary_cohort_decision']))
        if r['primary_cohort_decision']=='INCLUDE_PRIMARY' and r['specimen_class']!='PRIMARY_TUMOUR':
            raise ValueError('Cannot INCLUDE_PRIMARY a non-primary specimen: '+r['library_id'])
        if r['primary_cohort_decision']=='EXCLUDE_NONPRIMARY' and r['specimen_class']=='PRIMARY_TUMOUR':
            raise ValueError('Cannot EXCLUDE_NONPRIMARY an attested PRIMARY_TUMOUR: '+r['library_id'])

def reviewed_study(row):
    return (row['source_verified']=='YES' and all(safe_reference(row.get(x)) for x in
       ('verified_sample_types','patient_key_provenance','publication_or_SRA_reference'))
       and len(row.get('reviewer','').strip())>=3)

def reviewed_specimen(row):
    if row['primary_cohort_decision']=='HOLD':return False
    if row['specimen_class']=='UNRESOLVED':return False
    if not nonempty(row,'reviewed_specimen_id','accession_source_reference','accession_evidence_record','reviewer'):
        return False
    if not all(safe_reference(row.get(k)) for k in ('accession_source_reference','accession_evidence_record')) or len(row['reviewer'].strip())<3:
        return False
    # Exact accession linkage must be documented; record free text must contain
    # either the library accession or its BioSample ID, preventing generic paper-only attestations.
    record=row['accession_evidence_record'].upper()
    if row['library_id'].upper() not in record and row['biosample'].upper() not in record:
        return False
    return True

def reviewed_patient(row,study_good,group_ok,need_group):
    if not study_good:return False
    if not nonempty(row,'reviewed_patient_id','source_evidence_reference','reviewer'):
        return False
    if not safe_reference(row['source_evidence_reference']) or len(row['reviewer'].strip())<3:return False
    # Curated patient labels must be linked back to per-library source evidence,
    # not merely inherited from v1.0.0 defaults.
    if row['patient_identity_state']=='BIOSAMPLE_PROXY_NO_PATIENT_ID':
        pid=row['reviewed_patient_id'].strip()
        if pid in {row['biosample'],row['patient_key_from_curation'],row['proposed_patient_group']}:
            return False
    if need_group and not group_ok:return False
    return True

def validate_review_integrity(data,work):
    # This checks evidence immutable fields and manual run override provenance.
    rr,pr,groups,att,holds=v100.check_reviews(data,work)
    immutable={r['library_id']:r for r in v100.readtsv(work/'draft/specimen_eligibility_3602.tsv')}
    specs=v100.readtsv(work/'review/specimen_eligibility_3602.tsv')
    validated_specimen(specs,immutable)
    return rr,pr,groups,att,holds,specs

def verify_snapshot(work,previous):
    wm=json.loads((work/'draft/v101_work_manifest.json').read_text())
    if canonical(wm['previous_work'])!=canonical(previous):raise ValueError('Work belongs to different previous evidence snapshot')
    for r in v100.readtsv(work/'draft/v101_input_provenance.tsv'):
        if not Path(r['path']).is_file() or sha256(r['path'])!=r['sha256']:
            raise ValueError('Input seed/evidence changed: '+r['role'])
    for r in v100.readtsv(work/'draft/evidence_sha256.tsv'):
        if not Path(r['path']).is_file() or sha256(r['path'])!=r['sha256']:
            raise ValueError('Technical source evidence hash changed: '+r['path'])

def score_projections(projection_pairs,bridge):
    if not projection_pairs:return {'status':'NOT_PROVIDED'}
    return v100.evaluate_projection_pairs(projection_pairs,bridge)

def assess(root,previous,work,atlas_crosswalk=None,projection_pairs=None):
    root,previous,work=map(canonical,(root,previous,work))
    if not work.is_dir():raise ValueError('Run v1.0.1 init first')
    check_program(root);check_existing_lock(previous);verify_snapshot(work,previous)
    data,paths,reports,summary=v100.inspect(root,quiet=True)
    run,patient,groups,att,hold,specs=validate_review_integrity(data,work)
    mapp,W_summary=v100.verify_W(data,root)
    W_rescue=v100.verify_feature_resolutions(root,work,mapp)
    runlib=collections.defaultdict(list)
    for r in run:runlib[r['library_id']].append(r)
    pairs=collections.defaultdict(list)
    for r in data['pairs']:pairs[r['library_id']].append(r)
    patientby={r['library_id']:r for r in patient}
    studies={r['study']:r for r in att}
    groupby={(r['study'],r['patient_key']):r for r in groups}
    specby={r['library_id']:r for r in specs}
    # A biological specimen may have more than one prepared library. Such records
    # must not be auto-counted as independent even when patient labels differ.
    spec_groups=collections.defaultdict(list)
    for row in specs:
        if row['reviewed_specimen_id'].strip():
            spec_groups[(row['study'],row['reviewed_specimen_id'].strip())].append(row['library_id'])
    duplicated_specimens={lid for ids in spec_groups.values() if len(ids)>1 for lid in ids}
    bios_groups=collections.defaultdict(list)
    for lib in data['libs']:
        if lib['study'] not in v100.DGE and lib['biosample']:
            bios_groups[(lib['study'],lib['biosample'])].append(lib['library_id'])
    duplicated_biosamples={lid for ids in bios_groups.values() if len(ids)>1 for lid in ids}
    overlap_rows=None
    overlap={'status':'HOLD_INCOMPLETE_9465_ATLAS_ALIAS_CROSSWALK'}
    if atlas_crosswalk:
        overlap_rows,overlap=v100.verify_complete_atlas_crosswalk(data,root,canonical(atlas_crosswalk))
    overlapby={r['library_id']:r for r in overlap_rows} if overlap_rows else {}
    prior={r['library_id']:r for r in data['prior_link']}
    held_run_count=sum(r['decision']=='HOLD' for r in run)
    ledger=[]
    for lib in data['libs']:
        lid=lib['library_id'];is_dge=lib['study'] in v100.DGE
        if is_dge:
            ledger.append({'library_id':lid,'study':lib['study'],'cancer':lib['cancer'],
             'biosample':lib['biosample'],'specimen_class':'','primary_cohort_decision':'EXCLUDE_PROTOCOL',
             'reviewed_patient_id':'','independence_unit':'','included_run_ids':'','n_included_runs':0,
             'status':'EXCLUDE_PROTOCOL','hold_reasons':'','ARCHS4_overlap':'NOT_ASSESSED','source_verified':'NOT_APPLICABLE'})
            continue
        spec=specby[lid];pat=patientby[lid];sr=studies[lib['study']]
        good_study=reviewed_study(sr)
        group=groupby.get((lib['study'],lib['patient_key']))
        good_group=(group is None or (group['reviewed_relationship'] in PATIENT_RELATIONS
          and group['reviewed_relationship']!='TECHNICAL_DUPLICATE_RESOLVED'
          and safe_reference(group.get('source_evidence_reference')) and len(group.get('reviewer','').strip())>=3))
        good_patient=reviewed_patient(pat,good_study,good_group,group is not None)
        good_spec=reviewed_specimen(spec)
        rr=runlib[lid];kept=[r['run_id'] for r in rr if r['decision']=='RETAIN']
        pending=[r['run_id'] for r in rr if r['decision']=='HOLD']
        involved=[p for p in pairs[lid] if p['run_id_a'] in kept and p['run_id_b'] in kept]
        concordant=all(v100.positive(p['log2CPM_Pearson_r']) and float(p['log2CPM_Pearson_r'])>=0.9 for p in involved)
        reasons=[]
        if pending:reasons.append('RUN_QC_UNRESOLVED')
        if not kept:reasons.append('NO_RETAINED_RUN')
        if not concordant:reasons.append('RETAINED_RUNS_DISCORDANT')
        if not good_study:reasons.append('SOURCE_ATTESTATION_INCOMPLETE')
        if not good_spec:reasons.append('ACCESSION_SPECIMEN_ELIGIBILITY_UNVERIFIED')
        if not good_patient:reasons.append('PATIENT_IDENTITY_UNVERIFIED')
        if not good_group:reasons.append('REPEATED_PATIENT_RELATION_UNRESOLVED')
        if lid in duplicated_specimens:reasons.append('MULTIPLE_LIBRARIES_ONE_SPECIMEN_NOT_RESOLVED')
        if lid in duplicated_biosamples:reasons.append('SHARED_BIOSAMPLE_NOT_RESOLVED')
        if good_spec and spec['primary_cohort_decision']=='EXCLUDE_NONPRIMARY':
            # Exclusion can be confirmed even if other gates are incomplete.
            status='EXCLUDE_SOURCE_VERIFIED_NONPRIMARY'
        elif not reasons and spec['primary_cohort_decision']=='INCLUDE_PRIMARY':
            status='CERTIFIED_PRIMARY_WITHIN_STUDY_PATIENT_AWARE_ONLY'
        else:status='HOLD_REVIEW'
        prior_match=(v100.truth(prior[lid].get('is_known_atlas_exact_biosample_match','False')) or
                     v100.truth(prior[lid].get('is_known_atlas_exact_experiment_match','False')))
        overlapstatus=('KNOWN_OVERLAP_IN_PRIOR_PARTIAL_CROSSWALK' if prior_match else 'INCOMPLETE_CROSSWALK_NO_INDEPENDENCE_CLAIM')
        if overlap_rows:overlapstatus=overlapby[lid]['crossatlas_status']
        ledger.append({'library_id':lid,'study':lib['study'],'cancer':lib['cancer'],
         'biosample':lib['biosample'],'specimen_class':spec['specimen_class'],
         'primary_cohort_decision':spec['primary_cohort_decision'],
         'reviewed_patient_id':pat['reviewed_patient_id'],
         'independence_unit':(lib['study']+'|'+pat['reviewed_patient_id']) if good_patient else '',
         'included_run_ids':';'.join(kept),'n_included_runs':len(kept),
         'status':status,'hold_reasons':';'.join(reasons),
         'ARCHS4_overlap':overlapstatus,'source_verified':str(good_study)})
    certified=[r for r in ledger if r['status']=='CERTIFIED_PRIMARY_WITHIN_STUDY_PATIENT_AWARE_ONLY']
    independent_possible=[r for r in certified if r['ARCHS4_overlap']=='NO_MATCH_IN_COMPLETE_VERIFIED_CROSSWALK']
    if projection_pairs:
        proj=score_projections(canonical(projection_pairs),data['bridge_sample'])
    else:proj={'status':'NOT_PROVIDED'}
    # SCORE COMPARISON IS DESCRIPTIVE; no automatic equivalence PASS.
    counts=dict(collections.Counter(x['status'] for x in ledger))
    holdreasons=dict(collections.Counter(r for row in ledger for r in row['hold_reasons'].split(';') if r and row['status']=='HOLD_REVIEW'))
    independent_units=len({x['independence_unit'] for x in certified})
    result={'version':VERSION,'created_utc':utc(),
       'technical_evidence_lock':'PASS_V100_HASH_SNAPSHOT_VERIFIED',
       'source_attestations_complete':sum(reviewed_study(r) for r in att),
       'source_attestations_incomplete':sum(not reviewed_study(r) for r in att),
       'specimen_accession_records_reviewed':sum(reviewed_specimen(x) for x in specs),
       'n_candidate_libraries':3602,'n_primary_within_study_certified':len(certified),
       'n_within_study_patient_units_among_certified':independent_units,
       'n_eligible_for_external_overlap_review_if_full_crosswalk':len(independent_possible),
       'n_run_holds':held_run_count,'n_specimen_duplicate_library_holds':len(duplicated_specimens),
       'n_biosample_duplicate_library_holds':len(duplicated_biosamples),'n_proxy_patient_review_holds':sum(
            x['patient_identity_state']=='BIOSAMPLE_PROXY_NO_PATIENT_ID' and not nonempty(x,'reviewed_patient_id','source_evidence_reference','reviewer') for x in patient),
       'n_repeated_patient_group_unresolved':sum(not (x['reviewed_relationship'] in PATIENT_RELATIONS
           and x['reviewed_relationship']!='TECHNICAL_DUPLICATE_RESOLVED' and
           safe_reference(x.get('source_evidence_reference')) and len(x.get('reviewer','').strip())>=3) for x in groups),
       'n_features_still_held':W_rescue['n_held'],
       'W_feature_rescue':W_rescue,
       'model_projection_comparison':proj,
       'ARCHS4_crosswalk_gate':overlap['status'],
       'model_transfer_gate':'HOLD_NOT_EQUIVALENCE_CERTIFIED_REQUIRE_FULL_499_AND_INDEPENDENT_STUDIES',
       'pooled_ARCHS4_recount3_gate':'HOLD',
       'standalone_primary_cohort_gate':('ELIGIBLE_REVIEWED_SUBSET_WITHIN_STUDY_PATIENT_AWARE_ONLY' if certified else 'HOLD'),
       'n_library_statuses':counts,'library_hold_reason_counts':holdreasons,
       'scientific_limits':['Source classification is reviewer-attested and cannot be established by typed YES alone',
          'Within-study patient grouping is not cross-study global patient linkage',
          'Same-source 143 paired assays are not biological replication',
          'Full W-499 frozen-model score transfer requires executable proof, not a 465-feature approximation',
          'No automatic declaration of cross-platform pooled expression equivalence']}
    out=work/'assessment_v101';out.mkdir(parents=True,exist_ok=True)
    if overlap_rows:v100.writetsv(out/'complete_crossatlas_overlap_decisions.tsv',overlap_rows,list(overlap_rows[0]))
    v100.writetsv(out/'library_certification_3639.tsv',ledger,list(ledger[0]))
    v100.savejson(out/'release_gates.json',result)
    outputs=[]
    for d in (work/'review',out):
        for file in sorted(d.iterdir()):
            if file.is_file() and file.name!='SHA256SUMS.txt':
                outputs.append(f'{sha256(file)}  {file.relative_to(work)}')
    (out/'SHA256SUMS.txt').write_text('\n'.join(outputs)+'\n')
    print(json.dumps(result,indent=2))


def aggregate_study(root,previous,work,study):
    """Only after v1.0.1 assessor has approved a within-study cohort subset.

    This produces a stand-alone RAW G026 library matrix, not frozen-W model input.
    Only unscaled gene-body coverage is added across reviewed retained technical runs.
    Never aggregates two libraries into one independent biological observation.
    """
    import h5py, numpy as np
    root,previous,work=map(canonical,(root,previous,work))
    check_program(root);check_existing_lock(previous);verify_snapshot(work,previous)
    outputdir=work/'assessment_v101'
    ledgerpath=outputdir/'library_certification_3639.tsv'
    sums=outputdir/'SHA256SUMS.txt'
    if not ledgerpath.is_file() or not sums.is_file():raise ValueError('Run v1.0.1 assess first')
    # Exact check of reviewed evidence and assessor outputs. Any change forces reassessment.
    for line in sums.read_text().splitlines():
        sha,name=line.split(None,1)
        p=work/name.strip()
        if not p.is_file() or sha256(p)!=sha:raise ValueError('Assessment or review modified; reassess before aggregation: '+str(p))
    gate=json.loads((outputdir/'release_gates.json').read_text())
    if gate['technical_evidence_lock']!='PASS_V100_HASH_SNAPSHOT_VERIFIED':
        raise ValueError('Evidence lock not verified')
    eligible=[r for r in v100.readtsv(ledgerpath) if r['study']==study and
      r['status']=='CERTIFIED_PRIMARY_WITHIN_STUDY_PATIENT_AWARE_ONLY']
    if not eligible:raise ValueError('No source-verified, patient-aware certified primary specimens in '+study)
    # Verify same actual source manifest run membership before any matrix write.
    src=root/'stage2_output/studies'/study/'selected_run_coverage_G026.h5'
    if not src.is_file():raise ValueError('Source expression HDF5 not found: '+str(src))
    outfile=work/'aggregates_v101'/f'{study}.reviewed_raw_gene_body_coverage_G026.h5'
    outfile.parent.mkdir(parents=True,exist_ok=True)
    if outfile.exists():raise ValueError('Refusing to overwrite existing aggregate: '+str(outfile))
    tmp=outfile.with_suffix(outfile.suffix+'.partial')
    with h5py.File(src,'r') as f:
        runids=list(f['run_id'].asstr()[:]);libids=list(f['library_id'].asstr()[:]);genes=list(f['gene_id'].asstr()[:]);x=f['coverage_bp']
        if len(set(runids))!=len(runids) or x.shape!=(len(genes),len(runids)):
            raise ValueError('Source HDF5 run/gene axis invalid')
        idx={rid:i for i,rid in enumerate(runids)};selection=[]
        for r in sorted(eligible,key=lambda x:x['library_id']):
            keep=r['included_run_ids'].split(';')
            if not keep or not all(k in idx and libids[idx[k]]==r['library_id'] for k in keep):
                raise ValueError('Source run membership did not validate for '+r['library_id'])
            selection.append((r['library_id'],[idx[k] for k in keep],r['reviewed_patient_id']))
        try:
            with h5py.File(tmp,'w') as h:
                s=h5py.string_dtype('utf-8')
                h.create_dataset('gene_id',data=np.array(genes,dtype=object),dtype=s)
                h.create_dataset('library_id',data=np.array([a[0] for a in selection],dtype=object),dtype=s)
                h.create_dataset('patient_id',data=np.array([a[2] for a in selection],dtype=object),dtype=s)
                dst=h.create_dataset('coverage_bp',shape=(len(genes),len(selection)),dtype='float64',
                     chunks=(min(2000,len(genes)),1),compression='gzip',compression_opts=4)
                for i in range(0,len(genes),2000):
                    j=min(i+2000,len(genes))
                    for k,(_,run_index,_) in enumerate(selection):
                        block=np.asarray(x[i:j,sorted(run_index)],dtype=np.float64)
                        if not np.isfinite(block).all() or (block<0).any():
                            raise ValueError('Invalid unscaled coverage block')
                        dst[i:j,k]=block.sum(axis=1)
                h.attrs['representation']='UNSCALED_G026_GENE_BODY_COVERAGE_BP_ACCEPTED_RUN_SUM'
                h.attrs['model_input_approved']=False
                h.attrs['external_validation_independence_certified']=False
                h.attrs['source_study']=study
                h.attrs['v101_assessment_sha256']=sha256(ledgerpath)
                if 'gene_order_sha256' in f.attrs:h.attrs['gene_order_sha256']=f.attrs['gene_order_sha256']
            os.replace(tmp,outfile)
        finally:
            if tmp.exists():tmp.unlink()
    print(json.dumps({'status':'RAW_COVERAGE_AGGREGATED_NOT_W_MODEL_INPUT','study':study,
      'n_approved_primary_libraries':len(selection),'output':str(outfile)},indent=2))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',required=True,type=Path)
    ap.add_argument('--previous',type=Path,help='sealed v1.0.0 certification work directory')
    ap.add_argument('--work',type=Path,help='NEW sibling v1.0.1 work directory')
    sub=ap.add_subparsers(dest='command',required=True)
    sub.add_parser('preflight');sub.add_parser('init')
    a=sub.add_parser('assess');a.add_argument('--atlas-crosswalk',type=Path)
    a.add_argument('--projection-pairs',type=Path)
    agg=sub.add_parser('aggregate-study');agg.add_argument('--study',required=True)
    args=ap.parse_args()
    root=canonical(args.root);prev=canonical(args.previous or root.parent/'27_recount3_final_cert_v1_0_0')
    work=canonical(args.work or root.parent/'28_recount3_final_cert_v1_0_1')
    if args.command=='preflight':preflight(root,prev)
    elif args.command=='init':init(root,prev,work)
    elif args.command=='aggregate-study':aggregate_study(root,prev,work,args.study)
    else:assess(root,prev,work,args.atlas_crosswalk,args.projection_pairs)
if __name__=='__main__':
    try:main()
    except Exception as e:
        print('CERTIFICATION_V101_HOLD',type(e).__name__,str(e),file=sys.stderr)
        sys.exit(2)
