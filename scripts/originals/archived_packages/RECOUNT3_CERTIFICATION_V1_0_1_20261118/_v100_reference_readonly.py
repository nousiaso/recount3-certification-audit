#!/usr/bin/env python3
"""Recount3 final two-gate certification, v1.0.0.

Single entrypoint with inspect, draft, assess and aggregate-study subcommands.
All source files and frozen ARCHS4 release are read-only. All scientific holds are
explicit; this code NEVER infers patient independence from distinct BioSamples or
approves a 465/499-gene approximation as the original frozen 499-gene model.
"""
from __future__ import annotations
import argparse, collections, csv, datetime, hashlib, json, math, os, shutil, sys
from pathlib import Path
VERSION='recount3_final_certification_v1_0_0'
DGE={'SRP007947','SRP003901','SRP024268'}
KIDNEY='SRP082503'
KIDNEY_BAD={
 'SRR4046679','SRR4046680','SRR4046687','SRR4046688',
 'SRR4046695','SRR4046696','SRR4046703','SRR4046704',
 'SRR4046711','SRR4046712','SRR4046719','SRR4046720',
 'SRR4046727','SRR4046728'}
EXPECTED={'studies':129,'runs':4166,'libraries':3639,'dge_libraries':37,'remaining':3602,
          'w_unique':499,'w_usable':465,'bridge':143}

def readtsv(path):
    path=Path(path)
    with path.open('r',encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))
def writetsv(path,rows,cols):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.partial')
    with tmp.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader();w.writerows(rows)
    os.replace(tmp,path)
def savejson(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.partial');temp.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n');os.replace(temp,path)
def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def clean(x):return x.decode('utf-8','replace') if isinstance(x,bytes) else str(x)
def base(x):return clean(x).split('.')[0]
def truth(s):return str(s).lower() in {'true','1','yes'}
def flag(x):return str(x).lower() in {'true','1','yes'}
def positive(x):
    try:return math.isfinite(float(x)) and float(x)>0
    except (ValueError,TypeError):return False

def paths(root):
    root=Path(root)
    return {
      'runs':root/'manifest/selected_runs.tsv',
      'libs':root/'manifest/selected_libraries.tsv',
      'studies':root/'manifest/studies.tsv',
      'lock':root/'manifest/intake_lock.json',
      'stage2':root/'stage2_output/expression_intake_summary.json',
      'stage3':root/'stage3_output/02_descriptive_qc_summary.json',
      'stage4':root/'stage4_contract_output/00_contract_collection_summary.json',
      'stage5':root/'stage5_output/04_final_gate_summary.json',
      'run_qc':root/'stage5_output/04_run_coverage_conversion_qc.tsv',
      'library_qc':root/'stage5_output/04_library_preanalysis_decision_3639.tsv',
      'pairs':root/'stage5_output/04_multirun_concordance_pairs.tsv',
      'outliers':root/'stage3_output/02_within_study_outlier_review.tsv',
      'prior_link':root/'stage3_output/00_ARCHS4_prior_linkage_flags_NOT_FINAL.tsv',
      'ref_contract':root/'stage4_contract_output/03_frozen_reference_file_contract.tsv',
      'w_map':root/'stage5_output/00_frozen_W_to_G026_candidate_mapping.tsv',
      'bridge':root/'stage5_output/bridge_SRP074349/bridge_report.json',
      'bridge_gene':root/'stage5_output/bridge_SRP074349/per_W_gene_across_specimen_transfer.tsv',
      'bridge_sample':root/'stage5_output/bridge_SRP074349/per_specimen_same_source_gene_profile.tsv',
    }
def inspect(root, quiet=False):
    p=paths(root)
    missing=[k for k,path in p.items() if not path.is_file()]
    if missing:raise ValueError('Required existing files not found: '+str([(x,str(p[x])) for x in missing]))
    reports={k:json.loads(p[k].read_text()) for k in ('lock','stage2','stage3','stage4','stage5','bridge')}
    for k,n in [('studies',EXPECTED['studies']),('runs',EXPECTED['runs']),('libs',EXPECTED['libraries'])]:
        rows=readtsv(p[k]);
        if len(rows)!=n:raise ValueError(f'{k} has {len(rows)} rows, expected {n}')
    s2,s3,s4,s5,br=(reports[x] for x in ('stage2','stage3','stage4','stage5','bridge'))
    if any(s2.get(k)!=n for k,n in [('studies_extracted',129),('selected_runs_extracted',4166),('selected_libraries_with_all_runs',3639)]):
        raise ValueError('Stage2 contract failed')
    if any(s3.get(k)!=n for k,n in [('n_studies_completed',129),('n_runs',4166),('n_libraries',3639)]):
        raise ValueError('Stage3 contract failed')
    if s4.get('manifest_run_sha256')!=digest(p['runs']):raise ValueError('Manifest run digest mismatch vs Stage4')
    if s4.get('source_intake_xlsx_sha256')!=reports['lock'].get('xlsx_sha256'):raise ValueError('Source workbook fingerprint mismatch')
    if any(s5.get(k)!=n for k,n in [('n_studies',129),('n_runs',4166),('n_libraries',3639),('n_DGE_excluded',37),('n_remaining_provisional_libraries',3602)]):
        raise ValueError('Stage5 finalized cohort accounting mismatch')
    if br.get('status')!='ASSAY_BRIDGE_EXPRESSION_DIAGNOSTIC_COMPLETE_NOT_RELEASE_CERTIFIED' or br.get('n_bridge_libraries_target')!=143 or br.get('n_common_unambiguous_W_genes')!=465 or br.get('bridge_exact_GSM_verification')!=143:
        raise ValueError('Bridge report not the verified 143-library 465-gene diagnostic')
    data={k:readtsv(p[k]) for k in ('runs','libs','studies','run_qc','library_qc','pairs','outliers','prior_link','ref_contract','w_map','bridge_gene','bridge_sample')}
    if any(len(data[x])!=n for x,n in [('run_qc',4166),('library_qc',3639),('prior_link',3639),('w_map',1497),('bridge_gene',465),('bridge_sample',143)]):
        raise ValueError('Input table row-count mismatch')
    if len({r['run_id'] for r in data['run_qc']})!=4166 or len({r['library_id'] for r in data['library_qc']})!=3639:
        raise ValueError('Duplicate run or library IDs in QC')
    if set(r['run_id'] for r in data['runs']) != set(r['run_id'] for r in data['run_qc']):raise ValueError('Run ID set mismatch')
    if len({x['library_id'] for x in data['libs']})!=3639:raise ValueError('Duplicate library IDs in source manifest')
    dge={x['library_id'] for x in data['libs'] if x['study'] in DGE}
    if len(dge)!=37:raise ValueError('DGE mismatch: '+str(len(dge)))
    bylib=collections.Counter(x['library_id'] for x in data['runs'])
    if sum(n*(n-1)//2 for n in bylib.values())!=len(data['pairs']):raise ValueError('Pair count mismatch')
    if sum(1 for r in data['run_qc'] if r['run_id'] in KIDNEY_BAD)!=14:raise ValueError('Kidney 14-run set missing')
    study_bad=[]
    for row in data['studies']:
        file=Path(root)/'stage5_output/studies'/row['study']/'report.json'
        if not file.is_file() or json.loads(file.read_text()).get('status')!='TECHNICAL_DIAGNOSTIC_COMPLETE':study_bad.append(row['study'])
    if study_bad:raise ValueError('Stage5 study reports missing/not completed: '+str(study_bad))
    qc_summary={'status':'TECHNICAL_INPUTS_VERIFIED_DIAGNOSTIC_ONLY','version':VERSION,'n_studies':129,
       'n_runs':4166,'n_libraries':3639,'n_primary_candidates_after_DGE':3602,
       'n_pairs':len(data['pairs']),'n_kidney_sparse_runs':14,
       'bridge_samples':143,'W_gene_unambiguous_for_diagnostic':465,
       'study_report_failures':len(study_bad),
       'source_archive_frozen_unchanged':True,'full_atlas_overlap_independence_certified':False,
       'full_frozen_W_499_transfer_certified':False}
    if not quiet:print(json.dumps(qc_summary,indent=2))
    return data,p,reports,qc_summary

def verify_W(data,root, skip_h5=False):
    """Proof by exact symbol+Ensembl IDs and explicit G026 source row, not W ordinal."""
    wrows=data['w_map'];by_symbol=collections.defaultdict(list)
    for r in wrows:by_symbol[r['W_model_gene']].append(r)
    if len(by_symbol)!=499 or any(len(v)!=3 or {r['branch'] for r in v}!={'CPM','TMM_CPM','UQ_CPM'} for v in by_symbol.values()):
        raise ValueError('W branch contract changed')
    if skip_h5:return [],{'status':'H5_GENE_AXIS_NOT_VERIFIED','n_W_genes':499}
    import h5py
    ref={r['role']:r for r in data['ref_contract']}
    h5=Path(ref['ARCHS4_source_H5_EXTERNAL']['path'])
    if not h5.is_file():raise FileNotFoundError('Original ARCHS4 H5 not found: '+str(h5))
    wfile=Path(ref['frozen_W_harmonization']['path'])
    if not wfile.is_file() or digest(wfile)!=ref['frozen_W_harmonization']['sha256_small_file']:
        raise ValueError('Frozen W manifest file changed')
    first=data['studies'][0]['study']
    source=Path(root)/'stage2_output/studies'/first/'selected_run_coverage_G026.h5'
    with h5py.File(source,'r') as f:src=[base(x) for x in f['gene_id'].asstr()[:]]
    if len(src)!=63856:raise ValueError('Unexpected G026 row dimension')
    with h5py.File(h5,'r') as f:
        if f['data/expression'].shape[0]!=67186:raise ValueError('ARCHS4 gene-axis length changed')
        sym=[clean(x) for x in f['meta/genes/symbol'][:]]
        ens=[base(x) for x in f['meta/genes/ensembl_gene'][:]]
    axis=collections.defaultdict(list)
    for i,pair in enumerate(zip(sym,ens)):axis[pair].append(i)
    joined=[];reason=collections.Counter()
    for g,entries in sorted(by_symbol.items()):
        record={'gene_symbol':g,'branches':'CPM;TMM_CPM;UQ_CPM','G026_gene_id':'','G026_source_row_zero_based':'',
                'ARCHS4_H5_row_zero_based':'','map_state':'HOLD','hold_reason':''}
        if any(r['mapping_status']!='UNAMBIGUOUS_ONE_TO_ONE' for r in entries):
            record['hold_reason']='G026_MAPPING_AMBIGUOUS'
        else:
            cand={(base(r['matched_source_gene_id']),r['source_row_zero_based']) for r in entries}
            if len(cand)!=1:record['hold_reason']='G026_BRANCH_INCONSISTENT'
            else:
                gene_id,pos=next(iter(cand));j=int(pos)
                if j<0 or j>=len(src) or src[j]!=gene_id:
                    raise ValueError('Incorrect G026 identity for '+g)
                hits=axis.get((g,gene_id),[])
                if len(hits)!=1:record['hold_reason']='H5_SYMBOL_ENSEMBL_NOT_UNIQUE'
                else:
                    record.update({'G026_gene_id':gene_id,'G026_source_row_zero_based':j,
                       'ARCHS4_H5_row_zero_based':hits[0],'map_state':'EXACT','hold_reason':''})
        if record['hold_reason']:reason[record['hold_reason']]+=1
        joined.append(record)
    count=sum(x['map_state']=='EXACT' for x in joined)
    if count!=465 or len(joined)!=499 or reason['G026_MAPPING_AMBIGUOUS']!=23 or reason['H5_SYMBOL_ENSEMBL_NOT_UNIQUE']!=11:
        raise ValueError('Unexpected gene-mapping change: '+str(count)+' '+str(reason))
    return joined,{'status':'VERIFIED_465_OF_499_DIAGNOSTIC_ONLY','n_W_genes':499,'n_verified':count,
                   'n_held':499-count,'hold_reasons':dict(reason),'frozen_W_499_projection_allowed':False}

def build_proposals(data):
    qc={r['run_id']:r for r in data['run_qc']}
    stage={r['run_id']:r for r in data['runs']}
    outliers={r['run_id'] for r in data['outliers']}
    perlib=collections.defaultdict(list)
    for r in data['runs']:perlib[r['library_id']].append(r['run_id'])
    runrows=[]
    for r in data['runs']:
        rid=r['run_id'];q=qc[rid];why=[];auto='RETAIN';dec='RETAIN';evidence='AUTOMATED_DESCRIPTIVE_QC_ONLY'
        if r['study'] in DGE:
            auto=dec='EXCLUDE';why.append('CONFIRMED_INCOMPATIBLE_DGE_TAG_ASSAY');evidence='STAGE3_PROTOCOL_REVIEW_V1'
        else:
            if rid in KIDNEY_BAD:
                why.append('KIDNEY_SPARSE_14_RUNS_RETAIN_6_OF_8_PROPOSAL')
            if float(q['gene_detection_fraction'])<0.10:why.append('DETECTION_LT_0_10')
            if float(q['top10_coverage_fraction'])>0.50:why.append('TOP10_COVERAGE_GT_0_50')
            if float(q['avg_mapped_len'])<40:why.append('MAPPED_LEN_LT_40BP')
            if rid in outliers:why.append('WITHIN_STUDY_OUTLIER')
            if why:
                auto='EXCLUDE_SPARSE_KIDNEY' if rid in KIDNEY_BAD else 'REVIEW_QC'
                dec='HOLD';evidence=''
        runrows.append({'study':r['study'],'cancer':r['cancer'],'library_id':r['library_id'],
         'run_id':rid,'biosample':r['biosample'],'gene_detection_fraction':q['gene_detection_fraction'],
         'top10_coverage_fraction':q['top10_coverage_fraction'],'average_mapped_length_bp':q['avg_mapped_len'],
         'gene_sum_coverage_bp':q['gene_sum_coverage_bp'],'automatic_proposal':auto,
         'decision':dec,'reason_codes':';'.join(why),'evidence_reference':evidence,
         'reviewer':'','review_note':''})
    prior={r['library_id']:r for r in data['prior_link']}
    n_by_patient=collections.Counter((r['study'],r['patient_key']) for r in data['libs'] if r['study'] not in DGE and r['patient_key_basis']=='curated subject ID')
    libs=[];patient=[];patientgroups={};studies={r['study']:r['cancer'] for r in data['libs'] if r['study'] not in DGE}
    for r in data['libs']:
        if r['study'] in DGE:continue
        lid=r['library_id'];p=prior[lid]
        if r['patient_key_basis']=='curated subject ID':
            identity='CURATED_WITHIN_STUDY';review_pid=r['patient_key'];context='REPEATED_SUBJECT_GROUP_REVIEW' if n_by_patient[(r['study'],r['patient_key'])]>1 else 'SINGLE_LIBRARY_SUBJECT'
            if context=='REPEATED_SUBJECT_GROUP_REVIEW':
                key=(r['study'],r['patient_key']);patientgroups.setdefault(key,[]).append(lid)
        else:identity='BIOSAMPLE_PROXY_NO_PATIENT_ID';review_pid='';context='PATIENT_ID_UNRESOLVED'
        patient.append({'library_id':lid,'study':r['study'],'biosample':r['biosample'],
             'patient_key_from_curation':r['patient_key'],'patient_key_basis':r['patient_key_basis'],
             'proposed_patient_group':review_pid or 'UNKNOWN',
             'reviewed_patient_id':review_pid,'patient_identity_state':identity,
             'source_evidence_reference':'','reviewer':'',
             'specimen_relationship':'UNVERIFIED' if context!='SINGLE_LIBRARY_SUBJECT' else 'ONE_LIBRARY_SUBJECT',
             'review_note':''})
        match=truth(p.get('is_known_atlas_exact_biosample_match','False')) or truth(p.get('is_known_atlas_exact_experiment_match','False'))
        libs.append({'library_id':lid,'study':r['study'],'cancer':r['cancer'],'biosample':r['biosample'],
          'n_runs':int(r['n_runs']),'patient_key':r['patient_key'],'patient_identity_basis':identity,
          'patient_context':context,'known_ARCHS4_match_from_partial_crosswalk':str(match),
          'ARCHS4_full_independence':'NOT_VERIFIED','all_final_gates_passed':'False'})
    groups=[]
    for (study,pid),ll in sorted(patientgroups.items()):
        groups.append({'study':study,'patient_key':pid,'n_libraries':len(ll),'library_ids':';'.join(ll),
          'reviewed_relationship':'UNREVIEWED','source_evidence_reference':'','reviewer':'','review_note':''})
    attest=[]
    for study in sorted(studies):
        attest.append({'study':study,'cancer_example':studies[study],
          'source_verified':'NO','verified_sample_types':'',
          'patient_key_provenance':'','publication_or_SRA_reference':'','reviewer':'','review_note':''})
    return runrows,libs,patient,groups,attest

def check_reviews(data,work):
    runrows=readtsv(work/'review/run_decisions.tsv');proposed=readtsv(work/'draft/run_decisions.tsv')
    patient=readtsv(work/'review/patient_identity.tsv')
    groups=readtsv(work/'review/repeated_subject_relationships.tsv')
    studies=readtsv(work/'review/study_source_attestations.tsv')
    original_patients={r['library_id']:r for r in readtsv(work/'draft/patient_identity.tsv')}
    original_groups={(r['study'],r['patient_key']):r for r in readtsv(work/'draft/repeated_subject_relationships.tsv')}
    original_studies={r['study']:r for r in readtsv(work/'draft/study_source_attestations.tsv')}
    if {r['library_id'] for r in patient}!=set(original_patients):raise ValueError('Patient review library IDs changed')
    if {(r['study'],r['patient_key']) for r in groups}!=set(original_groups):raise ValueError('Patient group review keys changed')
    if {r['study'] for r in studies}!=set(original_studies):raise ValueError('Study attestation IDs changed')
    for r in patient:
        for k in ('study','biosample','patient_key_from_curation','patient_key_basis','patient_identity_state'):
            if r[k]!=original_patients[r['library_id']][k]:raise ValueError('Immutable patient evidence altered: '+r['library_id']+' '+k)
    for r in groups:
        for k in ('n_libraries','library_ids'):
            if r[k]!=original_groups[r['study'],r['patient_key']][k]:raise ValueError('Immutable group evidence altered')
    if len(runrows)!=4166 or len(patient)!=3602 or len(studies)!=len({r['study'] for r in data['libs'] if r['study'] not in DGE}):
        raise ValueError('Review table sizes changed')
    orig={r['run_id']:r for r in proposed}
    if len(orig)!=4166 or {r['run_id'] for r in runrows}!=set(orig):raise ValueError('Run review IDs changed')
    rbad=[]
    for r in runrows:
        base_row=orig[r['run_id']]
        for k in ('library_id','study','biosample','automatic_proposal','reason_codes'):
            if r[k]!=base_row[k]:raise ValueError(f'Review has altered evidence field {k} for {r["run_id"]}')
        decision=r['decision']
        if decision not in {'RETAIN','EXCLUDE','HOLD'}:raise ValueError('Unknown run decision '+decision)
        if base_row['automatic_proposal']=='EXCLUDE' and decision!='EXCLUDE':raise ValueError('Cannot override DGE exclusion')
        if base_row['automatic_proposal'] in {'REVIEW_QC','EXCLUDE_SPARSE_KIDNEY'}:
            if decision=='HOLD':rbad.append(r['run_id'])
            elif not r['evidence_reference'].strip() or not r['reviewer'].strip():
                raise ValueError('Manual QC override missing evidence/reviewer: '+r['run_id'])
    patients_by={x['library_id']:x for x in patient}
    if len(patients_by)!=3602:raise ValueError('Duplicate patient review IDs')
    studiesby={r['study']:r for r in studies}
    studies_hold=[r['study'] for r in studies if r['source_verified']!='YES' or not r['publication_or_SRA_reference'].strip() or not r['reviewer'].strip()]
    proxy_hold=[r['library_id'] for r in patient if r['patient_identity_state']=='BIOSAMPLE_PROXY_NO_PATIENT_ID' and (not r['reviewed_patient_id'].strip() or not r['source_evidence_reference'].strip() or not r['reviewer'].strip())]
    group_hold=[(r['study'],r['patient_key']) for r in groups if r['reviewed_relationship'] not in {'REPEATED_BIOLOGICAL_SPECIMENS','LONGITUDINAL','MULTIREGION','TECHNICAL_DUPLICATE_RESOLVED'} or not r['source_evidence_reference'].strip() or not r['reviewer'].strip()]
    for r in groups:
        if r['reviewed_relationship']=='TECHNICAL_DUPLICATE_RESOLVED':
            group_hold.append((r['study'],r['patient_key'])) # requires explicit representative map, not yet supported
    return runrows,patient,groups,studies,{'runs_held':rbad,'n_runs_held':len(rbad),'studies_unverified':studies_hold,
      'n_studies_unverified':len(studies_hold),'proxy_id_held':proxy_hold,'n_proxy_id_held':len(proxy_hold),
      'repeated_groups_held':group_hold,'n_repeated_groups_held':len(group_hold)}

def coverage_duplicate_check(data,root,work):
    # No absence-as-independence: 740/9465 historical crosswalk is partial.
    prior={r['library_id']:r for r in data['prior_link']}
    row=[];flagged=0
    for l in data['libs']:
        if l['study'] in DGE:continue
        p=prior[l['library_id']]
        exact=truth(p.get('is_known_atlas_exact_biosample_match','False')) or truth(p.get('is_known_atlas_exact_experiment_match','False'))
        flagged+=int(exact)
        row.append({'library_id':l['library_id'],'biosample':l['biosample'],
           'known_exact_overlap_in_partial_crosswalk':str(exact),
           'crosswalk_coverage':'740_OF_9465_ARCHS4_GEO_SAMPLES',
           'certified_independent_of_ARCHS4':'NO',
           'disposition':'KNOWN_OVERLAP_HOLD' if exact else 'INCOMPLETE_CROSSWALK_HOLD'})
    writetsv(work/'draft/ARCHS4_overlap_status.tsv',row,list(row[0]))
    return {'status':'HOLD_ARCHS4_INDEPENDENCE_NOT_ESTABLISHED','prior_exact_flags':flagged,
            'GEO_crosswalk_evidence_coverage':740,'frozen_atlas_samples':9465,
            'reason':'No complete source-verifiable accession crosswalk supplied; absence of known overlap is not evidence of independence'}

def model_transfer(data,W,summary):
    bridge=data['bridge_gene'];cor=[float(x['across_specimens_pearson_r']) for x in bridge]
    biases=[float(x['median_log2CPM_difference_recount3_minus_ARCHS4']) for x in bridge]
    import statistics
    return {'status':'HOLD_FULL_FROZEN_W_TRANSFER_NOT_VALIDATED',
       'n_model_W_unique':499,'n_W_diagnostic_exact':465,'n_W_unresolved':34,
       'n_lung_same_source_bridge':143,'gene_median_across_specimen_r':statistics.median(cor),
       'gene_median_abs_log2CPM_bias':statistics.median(abs(x) for x in biases),
       'n_genes_abs_bias_gt_1':sum(abs(x)>1 for x in biases),
       'evidence_scope':'TECHNICAL_REPRODUCIBILITY_SAME_SEQUENCING_EXPERIMENTS',
       'required':'full 499-W transfer pathway or justified, separately validated 465-W model; actual frozen programme/ecotype output comparison; independent studies'}

def evidence_paths(p, root, work):
    files=[]
    for name in ('runs','libs','studies','lock','stage2','stage3','stage4','stage5','run_qc','library_qc','pairs','outliers','prior_link','ref_contract','w_map','bridge','bridge_gene','bridge_sample'):
        files.append((name,p[name]))
    for x in sorted((Path(root)/'stage5_output/studies').glob('*/report.json')):files.append(('stage5_study_report',x))
    for x in sorted((Path(work)/'draft').glob('*.tsv')):files.append(('draft',x))
    return [{'role':kind,'path':str(path),'sha256':digest(path),'bytes':path.stat().st_size} for kind,path in files]

def draft(root,work):
    data,p,reports,summary=inspect(root,quiet=True)
    W,w_summary=verify_W(data,root)
    run,libs,patient,groups,studies=build_proposals(data)
    work=Path(work)
    if (work/'review').exists():raise ValueError('Review directory already exists: refusing to overwrite existing review decisions. Use a new --work directory')
    (work/'draft').mkdir(parents=True,exist_ok=True)
    (work/'review').mkdir(parents=True,exist_ok=False)
    for name,rows in [('run_decisions',run),('library_candidates',libs),('patient_identity',patient),('repeated_subject_relationships',groups),('study_source_attestations',studies),('frozen_W_mapping',W)]:
        if rows:writetsv(work/'draft'/f'{name}.tsv',rows,list(rows[0]))
    for name,rows in [('run_decisions',run),('patient_identity',patient),('repeated_subject_relationships',groups),('study_source_attestations',studies)]:
        if rows:writetsv(work/'review'/f'{name}.tsv',rows,list(rows[0]))
    # Provide a separate, evidence-requiring template for the 34 unresolved model genes.
    held=[{'gene_symbol':x['gene_symbol'],'hold_reason':x['hold_reason'],
           'decision':'HOLD','reviewed_G026_gene_id':'','reviewed_G026_row_zero_based':'',
           'reviewed_ARCHS4_H5_row_zero_based':'','source_evidence_reference':'',
           'reviewer':'','review_note':''} for x in W if x['map_state']=='HOLD']
    if len(held)!=34:raise ValueError('Expected 34 W holds')
    writetsv(work/'review/held_W_feature_resolutions.tsv',held,list(held[0]))
    savejson(work/'review/model_transfer_review.json',{'status':'HOLD',
        'full_frozen_499_feature_projection_validated':False,
        'paired_projection_comparison':'NOT_PERFORMED',
        'independent_study_transfer':'NOT_PERFORMED',
        'reviewer':'','evidence_reference':'',
        'warning':'Editing this file alone cannot approve the model: scientific gate remains HOLD without executable full-feature proof.'})
    savejson(work/'review/ARCHS4_crosswalk_review.json',{'status':'INCOMPLETE',
        'frozen_atlas_expected':9465,'previously_fetched_GEO_crosswalk_rows':740,
        'source_verified_complete_crosswalk':'NOT_AVAILABLE',
        'warning':'No duplicate-free claim allowed from incomplete crosswalk; reviewer attestation alone cannot certify this.'})
    overlap=coverage_duplicate_check(data,root,work)
    model=model_transfer(data,W,w_summary)
    savejson(work/'draft/gene_mapping_contract.json',w_summary)
    savejson(work/'draft/ARCHS4_overlap_contract.json',overlap)
    savejson(work/'draft/model_transfer_contract.json',model)
    savejson(work/'draft/technical_evidence_summary.json',summary)
    evidence=evidence_paths(p,root,work)
    writetsv(work/'draft/evidence_sha256.tsv',evidence,list(evidence[0]))
    report={'version':VERSION,'status':'DRAFT_REVIEW_NOT_RELEASED','input_contract':summary,
       'W_mapping':w_summary,'ARCHS4_overlap':overlap,'transfer':model,
       'n_library_candidates':3602,'n_runs_proposed_hold':sum(r['decision']=='HOLD' for r in run),
       'n_kidney_sparse_runs_held_for_source_review':sum(r['run_id'] in KIDNEY_BAD and r['decision']=='HOLD' for r in run),
       'n_proxy_patient_identity_review':sum(r['patient_identity_state']=='BIOSAMPLE_PROXY_NO_PATIENT_ID' for r in patient),
       'n_repeated_subject_groups':len(groups),'all_approval_gates':'HOLD'}
    savejson(work/'draft/draft_summary.json',report)
    print(json.dumps({'status':report['status'],'n_run_review_holds':report['n_runs_proposed_hold'],
      'n_repeated_patient_groups':len(groups),'n_proxy_patient_ids':report['n_proxy_patient_identity_review'],
      'W_exact':465,'W_hold':34,'out':str(work)},indent=2))

def verify_complete_atlas_crosswalk(data,root,crosswalk_path):
    """A complete, evidence-backed alias crosswalk is necessary for an external-independence claim.
    All 9465 atlas sample IDs must match the frozen atlas exactly. Absence of
    a known duplicate in an incomplete crosswalk is never called independent.
    """
    import gzip
    crosswalk_path=Path(crosswalk_path)
    if not crosswalk_path.is_file():raise FileNotFoundError(crosswalk_path)
    rows=readtsv(crosswalk_path)
    expected={'sample_id','biosample','sra_experiment','sra_runs','source_verified','source_reference'}
    if not rows or not expected.issubset(rows[0]):raise ValueError('Crosswalk requires fields: '+str(sorted(expected)))
    if len(rows)!=9465 or len({r['sample_id'] for r in rows})!=9465:raise ValueError('Complete crosswalk must have 9465 unique atlas sample IDs')
    ref={r['role']:r for r in data['ref_contract']}
    atlas=Path(ref['frozen_9465_projection_atlas']['path'])
    if not atlas.is_file():raise FileNotFoundError(atlas)
    opener=gzip.open if atlas.name.endswith('.gz') else open
    with opener(atlas,'rt',encoding='utf-8',newline='') as f:
        expected_ids={r['sample_id'] for r in csv.DictReader(f,delimiter='\t')}
    if len(expected_ids)!=9465 or {r['sample_id'] for r in rows}!=expected_ids:
        raise ValueError('Crosswalk does not match the exact frozen 9465 sample IDs')
    if any(r['source_verified']!='YES' or not r['source_reference'] for r in rows):
        raise ValueError('All 9465 crosswalk aliases need source verification and citation')
    if any(not (r['biosample'] or r['sra_experiment'] or r['sra_runs']) for r in rows):
        raise ValueError('Cannot claim complete alias coverage: some frozen samples lack all SRA/BioSample accessions')
    bs=collections.defaultdict(set);ex=collections.defaultdict(set);run=collections.defaultdict(set)
    for r in rows:
        if r['biosample']:bs[r['biosample']].add(r['sample_id'])
        if r['sra_experiment']:ex[r['sra_experiment']].add(r['sample_id'])
        for t in r['sra_runs'].split(';'):
            if t:run[t].add(r['sample_id'])
    runby=collections.defaultdict(list)
    for r in data['runs']:runby[r['library_id']].append(r['run_id'])
    out=[];n=0
    for r in data['libs']:
        if r['study'] in DGE:continue
        linked=bs[r['biosample']] | ex[r['library_id']]
        for t in runby[r['library_id']]:linked=linked|run[t]
        n+=bool(linked)
        out.append({'library_id':r['library_id'],'matched_ARCHS4_GSMs':';'.join(sorted(linked)),
                    'crossatlas_status':'EXACT_MATCH_HOLD_OUT_OF_EXTERNAL_VALIDATION' if linked else 'NO_MATCH_IN_COMPLETE_VERIFIED_CROSSWALK'})
    return out,{'status':'FULL_9465_ALIAS_CROSSWALK_VERIFIED','crosswalk_sha256':digest(crosswalk_path),
       'n_frozen_atlas_samples':9465,'n_recount3_libraries_with_exact_overlap':n,
       'scope':'exact accession-level overlap; patient aliasing across separately submitted datasets can remain unknown'}

def verify_feature_resolutions(root,work,mapped):
    """Review fields cannot overrule a mismatch in actual HDF5 and G026 gene labels."""
    review=readtsv(Path(work)/'review/held_W_feature_resolutions.tsv')
    held={r['gene_symbol']:r for r in mapped if r['map_state']=='HOLD'}
    if len(review)!=34 or {r['gene_symbol'] for r in review}!=set(held):
        raise ValueError('Held-feature review manifest no longer matches immutable W ledger')
    resolved=[r for r in review if r['decision']=='MAP']
    if not resolved:return {'n_rescued':0,'n_held':34,'status':'34_UNRESOLVED'}
    import h5py
    mapping=readtsv(Path(work)/'draft/frozen_W_mapping.tsv')
    if [(r['gene_symbol'],r['map_state']) for r in mapping]!=[(r['gene_symbol'],r['map_state']) for r in mapped]:
        raise ValueError('Frozen W evidence drift')
    original=Path(root)
    refs={r['role']:r for r in readtsv(original/'stage4_contract_output/03_frozen_reference_file_contract.tsv')}
    source=original/'stage2_output/studies'/readtsv(original/'manifest/studies.tsv')[0]['study']/'selected_run_coverage_G026.h5'
    with h5py.File(source,'r') as f:sgene=[base(x) for x in f['gene_id'].asstr()[:]]
    with h5py.File(refs['ARCHS4_source_H5_EXTERNAL']['path'],'r') as f:
        ss=[clean(x) for x in f['meta/genes/symbol'][:]]
        es=[base(x) for x in f['meta/genes/ensembl_gene'][:]]
    used_g026={int(r['G026_source_row_zero_based']) for r in mapped if r['map_state']=='EXACT'}
    used_h5={int(r['ARCHS4_H5_row_zero_based']) for r in mapped if r['map_state']=='EXACT'}
    for r in resolved:
        if not r['reviewer'].strip() or not r['source_evidence_reference'].strip():raise ValueError('Missing evidence for feature rescue: '+r['gene_symbol'])
        a=int(r['reviewed_G026_row_zero_based']);b=int(r['reviewed_ARCHS4_H5_row_zero_based']);eid=base(r['reviewed_G026_gene_id'])
        if not 0<=a<len(sgene) or not 0<=b<len(ss) or sgene[a]!=eid or ss[b]!=r['gene_symbol'] or es[b]!=eid:
            raise ValueError('Feature rescue source identities did not verify: '+r['gene_symbol'])
        if a in used_g026 or b in used_h5:raise ValueError('Feature rescue reuses row assigned to another W gene: '+r['gene_symbol'])
        used_g026.add(a);used_h5.add(b)
    return {'status':'499_FEATURE_IDENTITIES_VERIFIED_NOT_MODEL_TRANSFER' if len(resolved)==34 else 'FEATURE_RESCUE_INCOMPLETE',
        'n_rescued':len(resolved),'n_held':34-len(resolved),
        'warning':'This proves gene coordinates only, not equivalence of the frozen projection or normalization'}

def evaluate_projection_pairs(path,bridge):
    """Descriptive model-output transfer, never an automatic assay equivalence PASS."""
    rows=readtsv(path)
    expected={'sample_id','component','ARCHS4_fraction','recount3_fraction'}
    if not rows or not expected.issubset(rows[0]):raise ValueError('Projection pairs need '+str(sorted(expected)))
    components={'NonSE'}|{f'SE{i}' for i in range(1,10)}
    by=collections.defaultdict(dict)
    for r in rows:
        sid,comp=r['sample_id'],r['component']
        if comp not in components or comp in by[sid]:raise ValueError('Unexpected or repeated component: '+str((sid,comp)))
        a=float(r['ARCHS4_fraction']);b=float(r['recount3_fraction'])
        if not all(math.isfinite(x) and 0<=x<=1 for x in (a,b)):raise ValueError('Fractions must be finite and in [0,1]')
        by[sid][comp]=(a,b)
    expected_samples={r['ARCHS4_GSM'] for r in bridge}
    if len(by)!=143 or set(by)!=expected_samples or len(rows)!=1430 or any(set(x)!=components for x in by.values()):
        raise ValueError('Need 10 complete components for each of 143 same-specimen samples')
    import statistics
    tvd=[]
    for sid,comp in by.items():
        if abs(sum(x[0] for x in comp.values())-1)>1e-3 or abs(sum(x[1] for x in comp.values())-1)>1e-3:
            raise ValueError('Fraction columns do not sum to 1 for '+sid)
        tvd.append(.5*sum(abs(a-b) for a,b in comp.values()))
    return {'status':'DESCRIPTIVE_143_SAME_SOURCE_PROJECTED_SCORES_NOT_GENERALIZABLE',
      'n_pairs':143,'median_total_variation_distance':statistics.median(tvd),
      'p95_total_variation_distance':sorted(tvd)[int(.95*(len(tvd)-1))],
      'max_total_variation_distance':max(tvd),'input_sha256':digest(path),
      'approval':'HOLD_REQUIRES_SOURCE_CODE_AUDIT_FULL_499_FEATURES_AND_INDEPENDENT_STUDIES'}

def assess(root,work,atlas_crosswalk=None,projection_pairs=None):
    data,p,rep,summary=inspect(root,quiet=True)
    W,w_summary=verify_W(data,root)
    work=Path(work)
    if not (work/'draft/evidence_sha256.tsv').is_file():raise ValueError('Must run draft first')
    for row in readtsv(work/'draft/evidence_sha256.tsv'):
        if not Path(row['path']).is_file() or digest(row['path'])!=row['sha256']:
            raise ValueError('SOURCE OR DRAFT EVIDENCE HASH CHANGED: '+row['path'])
    run,patient,groups,attest,hold=check_reviews(data,work)
    run_by_lib=collections.defaultdict(list)
    for r in run:run_by_lib[r['library_id']].append(r)
    pair_by_lib=collections.defaultdict(list)
    for r in data['pairs']:pair_by_lib[r['library_id']].append(r)
    person_by={r['library_id']:r for r in patient}
    subject_counts=collections.Counter((r['study'],r['reviewed_patient_id']) for r in patient if r['reviewed_patient_id'])
    attest_by={r['study']:r for r in attest}
    group_by={(r['study'],r['patient_key']):r for r in groups}
    out=[]
    for lib in data['libs']:
        lid=lib['library_id'];rr=run_by_lib[lid]
        if lib['study'] in DGE:
            state='EXCLUDE_PROTOCOL';keep=[];why=['DGE_SHORT_TAG_INCOMPATIBLE']
        else:
            keep=[r['run_id'] for r in rr if r['decision']=='RETAIN'];held=[r for r in rr if r['decision']=='HOLD']
            why=[]
            if held:why.append('UNRESOLVED_RUN_QC_REVIEW')
            if not keep:why.append('NO_RUNS_RETAINED')
            related=[p for p in pair_by_lib[lid] if p['run_id_a'] in keep and p['run_id_b'] in keep]
            low_or_missing=[p for p in related if not positive(p['log2CPM_Pearson_r']) or float(p['log2CPM_Pearson_r'])<0.90]
            if low_or_missing:why.append('RETAINED_MULTIRUN_PAIR_DISCORDANCE')
            ident=person_by[lid]
            if not ident['reviewed_patient_id'].strip():why.append('PATIENT_ID_UNRESOLVED')
            if ident['patient_identity_state']=='BIOSAMPLE_PROXY_NO_PATIENT_ID' and not ident['source_evidence_reference'].strip():why.append('NO_PATIENT_SOURCE_EVIDENCE')
            if attest_by[lib['study']]['source_verified']!='YES':why.append('STUDY_SOURCE_NOT_VERIFIED')
            gr=group_by.get((lib['study'],lib['patient_key']))
            if gr and (gr['study'],gr['patient_key']) in hold['repeated_groups_held']:why.append('REPEATED_SUBJECT_RELATION_UNRESOLVED')
            # Cross-atlas overlap is deliberately a separate gate, not an in-dataset QC exclusion.
            state='TECHNICAL_PATIENT_READY_FOR_SEPARATE_ANALYSIS' if not why else 'HOLD_REVIEW'
        out.append({'library_id':lid,'study':lib['study'],'cancer':lib['cancer'],
          'biosample':lib['biosample'],'patient_id_reviewed':person_by[lid]['reviewed_patient_id'] if lid in person_by else '',
          'included_run_ids':';'.join(keep),'n_included_runs':len(keep),'status':state,'hold_reasons':';'.join(why),
          'independent_of_ARCHS4':'NOT_CERTIFIED'})
    ready=sum(x['status']=='TECHNICAL_PATIENT_READY_FOR_SEPARATE_ANALYSIS' for x in out)
    counts=dict(collections.Counter(x['status'] for x in out))
    overlap=json.loads((work/'draft/ARCHS4_overlap_contract.json').read_text())
    model=model_transfer(data,W,w_summary)
    rescue=verify_feature_resolutions(root,work,W)
    overlap_rows=None
    if atlas_crosswalk:
        overlap_rows,overlap=verify_complete_atlas_crosswalk(data,root,atlas_crosswalk)
    projection=evaluate_projection_pairs(projection_pairs,data['bridge_sample']) if projection_pairs else {'status':'NOT_PROVIDED'}
    release={
      'version':VERSION,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'technical_evidence_lock':'PASS_INPUT_HASHES_VERIFIED','biological_cohort_separate_analysis':('ELIGIBLE_SUBSET_ONLY_NOT_FULL_COHORT_LOCK' if ready>0 else 'HOLD'),
      'n_candidate_libraries_after_protocol':3602,'n_technical_and_patient_ready_standalone':ready,
      'n_library_statuses':counts,'n_run_review_holds':hold['n_runs_held'],
      'n_patient_proxy_holds':hold['n_proxy_id_held'],'n_repeated_patient_group_holds':hold['n_repeated_groups_held'],
      'n_study_source_attestation_holds':hold['n_studies_unverified'],
      'n_frozen_W_features_held':rescue['n_held'],
      'held_W_rescue_status':rescue['status'],
      'paired_frozen_projection_evidence':projection,
      'model_transfer_review_status':'HOLD_REQUIRES_FULL_MODEL_OUTPUTS_AND_INDEPENDENT_STUDIES',
      'ARCHS4_independent_external_validation_gate':overlap['status'],
      'frozen_W_model_transfer_gate':model['status'],
      'pooled_ARCHS4_recount3_analysis_gate':'HOLD_NOT_APPROVED',
      'scientific_restrictions':[
         'No independent-patient claims from BioSample proxies',
         'No same-source bridge specimens counted as external validation',
         'No model projection with a silent 465-for-499 feature substitution',
         'No pooled recount3 + ARCHS4 expression matrix without independent cross-platform score calibration',
         'Do not call descriptive QC completion an analysis-ready cohort lock'],
    }
    outdir=work/'assessment';outdir.mkdir(parents=True,exist_ok=True)
    if overlap_rows:
        writetsv(outdir/'complete_crossatlas_overlap_decisions.tsv',overlap_rows,list(overlap_rows[0]))
    writetsv(outdir/'library_release_ledger.tsv',out,list(out[0]))
    savejson(outdir/'release_gates.json',release)
    # Hash all generated assessment and reviewed files without mutating frozen inputs.
    hashes=[]
    for dir in (work/'draft',work/'review',outdir):
        for path in sorted(dir.iterdir()):
            if path.is_file() and path.name!='SHA256SUMS.txt':hashes.append({'path':str(path.relative_to(work)),'sha256':digest(path),'bytes':path.stat().st_size})
    writetsv(outdir/'evidence_outputs_sha256.tsv',hashes,list(hashes[0]))
    print(json.dumps(release,indent=2))

def seal(root,work):
    """Content-addressed, non-overwriting technical evidence snapshot (not a cohort PASS)."""
    work=Path(work)
    inputlist=work/'draft/evidence_sha256.tsv'
    if not inputlist.is_file():raise ValueError('Run draft before seal')
    inputs=readtsv(inputlist)
    if not inputs:raise ValueError('No evidence to lock')
    for r in inputs:
        p=Path(r['path'])
        if not p.is_file() or digest(p)!=r['sha256']:raise ValueError('Input file changed: '+str(p))
    locked=work/'TECHNICAL_EVIDENCE_LOCK_v1_0_0'
    if locked.exists():raise ValueError('Refusing overwrite existing lock '+str(locked))
    temp=work/'TECHNICAL_EVIDENCE_LOCK_v1_0_0.partial'
    if temp.exists():raise ValueError('Stale partial lock exists; inspect manually '+str(temp))
    temp.mkdir(parents=True)
    try:
        shutil.copy2(inputlist,temp/'source_evidence_sha256.tsv')
        shutil.copy2(work/'draft/draft_summary.json',temp/'draft_summary.json')
        shutil.copy2(Path(__file__),temp/'certify_v1_0_0.py')
        source_root=Path(root)
        manifest={'version':VERSION,'release_type':'TECHNICAL_EVIDENCE_ONLY',
                  'status':'EVIDENCE_HASH_LOCK_COMPLETE_NOT_BIOLOGICAL_COHORT_CERTIFICATION',
                  'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  'n_source_files_fingerprinted':len(inputs),
                  'frozen_ARCHS4_release_mutated':False,
                  'warning':'Content-addressed evidence snapshot, not filesystem write-once storage. Preserve in a read-only repository for true immutability.',
                  'analysis_ready':False,'external_validation_certified':False,
                  'full_frozen_W_projection_certified':False}
        savejson(temp/'LOCK_MANIFEST.json',manifest)
        names=[]
        for x in sorted(temp.iterdir()):
            if x.is_file():names.append(f'{digest(x)}  {x.name}')
        (temp/'SHA256SUMS.txt').write_text('\n'.join(names)+'\n')
        os.replace(temp,locked)
    except Exception:
        if temp.exists():shutil.rmtree(temp)
        raise
    print(json.dumps({'lock_status':manifest['status'],'lock_directory':str(locked),
        'source_files_fingerprinted':len(inputs)},indent=2))

def aggregate_study(root,work,study):
    """Add RAW G026 coverage across RETAINED runs, then preserve a raw per-library matrix.
    No AUC-normalized values are ever summed. Requires approved per-library state
    from assess, and confirms report/gene/run axes before writing a new HDF5.
    """
    import h5py,numpy as np
    work=Path(work); root=Path(root)
    assessment=work/'assessment/library_release_ledger.tsv'
    if not assessment.is_file():raise ValueError('Run assess before aggregation')
    gatefile=work/'assessment/release_gates.json'
    evidencefile=work/'assessment/evidence_outputs_sha256.tsv'
    if not gatefile.is_file() or not evidencefile.is_file():raise ValueError('Missing assessed release gates/checksums')
    gate=json.loads(gatefile.read_text())
    if gate.get('technical_evidence_lock')!='PASS_INPUT_HASHES_VERIFIED':raise ValueError('Technical evidence lock not verified')
    entry=[r for r in readtsv(evidencefile) if r['path']=='assessment/library_release_ledger.tsv']
    if len(entry)!=1 or digest(assessment)!=entry[0]['sha256']:
        raise ValueError('Assessed library-release ledger changed since checksums were written')
    led=[r for r in readtsv(assessment) if r['study']==study and r['status']=='TECHNICAL_PATIENT_READY_FOR_SEPARATE_ANALYSIS']
    if not led:raise ValueError('No certified-for-standalone library in '+study)
    source=root/'stage2_output/studies'/study/'selected_run_coverage_G026.h5'
    if not source.is_file():raise FileNotFoundError(source)
    bylib={r['library_id']:r for r in led}
    out=work/'aggregates'/f'{study}.accepted_raw_coverage_G026.h5'
    out.parent.mkdir(parents=True,exist_ok=True)
    if out.exists():raise ValueError('Refusing to overwrite existing aggregate '+str(out))
    with h5py.File(source,'r') as f:
        runids=list(f['run_id'].asstr()[:]);libids=list(f['library_id'].asstr()[:]);geneids=list(f['gene_id'].asstr()[:]);X=f['coverage_bp']
        if X.shape!=(len(geneids),len(runids)):raise ValueError('Bad source dimensions')
        index={r:i for i,r in enumerate(runids)}
        chosen=[]
        for lid,row in sorted(bylib.items()):
            runs=row['included_run_ids'].split(';')
            if not all(r in index and libids[index[r]]==lid for r in runs):raise ValueError('Run/library identity mismatch '+lid)
            chosen.append((lid,[index[r] for r in runs]))
        tmp=out.with_suffix(out.suffix+'.partial')
        try:
            with h5py.File(tmp,'w') as h:
                strings=h5py.string_dtype('utf-8')
                h.create_dataset('gene_id',data=np.array(geneids,dtype=object),dtype=strings)
                h.create_dataset('library_id',data=np.array([x[0] for x in chosen],dtype=object),dtype=strings)
                result=h.create_dataset('coverage_bp',shape=(len(geneids),len(chosen)),dtype='float64',
                    chunks=(min(2000,len(geneids)),1),compression='gzip',compression_opts=4)
                for i in range(0,len(geneids),2000):
                    j=min(len(geneids),i+2000)
                    for col,(_,indices) in enumerate(chosen):
                        block=np.asarray(X[i:j,sorted(indices)],dtype=np.float64)
                        if not np.isfinite(block).all() or (block<0).any():raise ValueError('Invalid source coverage')
                        result[i:j,col]=np.sum(block,axis=1)
                h.attrs['representation']='UNSCALED_GENE_BODY_COVERAGE_BP_SUM_ACROSS_ACCEPTED_RUNS'
                h.attrs['source_study']=study
                h.attrs['n_libraries']=len(chosen)
                h.attrs['not_a_model_input']=True
                h.attrs['gene_order_sha256']=f.attrs['gene_order_sha256']
            os.replace(tmp,out)
        finally:
            if tmp.exists():tmp.unlink()
    print('AGGREGATED',study,'libraries',len(chosen),'genes',len(geneids),'output',out)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,required=True,help='Existing 26_recount3_external_primary_v1 directory')
    ap.add_argument('--work',type=Path,help='Independent output e.g. ../27_recount3_final_cert_v1_0_0')
    sub=ap.add_subparsers(dest='cmd',required=True)
    sub.add_parser('inspect');sub.add_parser('draft');sub.add_parser('seal');review=sub.add_parser('assess')
    review.add_argument('--atlas-crosswalk',type=Path,help='optional 9465-sample source-backed alias crosswalk')
    review.add_argument('--projection-pairs',type=Path,help='optional full-model 143x10 paired fractions table')
    agg=sub.add_parser('aggregate-study');agg.add_argument('--study',required=True)
    args=ap.parse_args();root=args.root.expanduser().resolve();work=(args.work or root.parent/'27_recount3_final_cert_v1_0_0').expanduser().resolve()
    if work==root or root in work.parents or work in root.parents:raise ValueError('--work must be a separate sibling directory, not inside the source root')
    try:
        if args.cmd=='inspect':inspect(root)
        elif args.cmd=='draft':draft(root,work)
        elif args.cmd=='assess':assess(root,work,args.atlas_crosswalk,args.projection_pairs)
        elif args.cmd=='seal':seal(root,work)
        elif args.cmd=='aggregate-study':aggregate_study(root,work,args.study)
    except Exception as exc:
        print('CERTIFICATION HOLD',type(exc).__name__,str(exc),file=sys.stderr)
        sys.exit(2)
if __name__=='__main__':main()
