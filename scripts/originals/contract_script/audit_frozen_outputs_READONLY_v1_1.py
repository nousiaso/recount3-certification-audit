#!/usr/bin/env python3
"""Read-only ARCHS4 stages 09--13 frozen-output contract audit.

Streams only status JSON and modest metadata/score/parameter tables.
Does NOT read large expression, reconstruct coefficients, run SpatialEcoTyper,
or modify source data. Writes an optional small report only outside the frozen release.
"""
import argparse, csv, gzip, hashlib, json, math
from collections import Counter
from pathlib import Path

STAGES={
 '09': ('10_preprojection_normalization_v1_0_1','18_preprojection_normalization_status.json','1.0.1'),
 '10': ('11_pan_cancer_fixed_anchor_projection_v1_0_4','28_pan_cancer_projection_status.json','1.0.4'),
 '11': ('12_layer1_adequacy_robustness_atlas_v1_0_2_1','13_layer1_adequacy_robustness_status.json','1.0.2.1'),
 '12': ('13_cancer_specific_fixed_anchor_projection_v1_0_3','32_cancer_specific_projection_status.json','1.0.3'),
 '13': ('14_projection_use_analysis_ready_atlas_v1_0_3','13_projection_use_status.json','1.0.3'),
}
SRC={
 '09':'09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py',
 '10':'10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R',
 '11':'11_build_ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1.R',
 '12':'12_run_ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3.R',
 '13':'13_build_ARCHS4_projection_use_analysis_ready_atlas_v1_0_3.R',
}
COMPS=['NonSE']+[f'SE{i}' for i in range(1,10)]
AUTHORITATIVE_FROZEN_SOURCE_SHA256={
 '09':'46d6ad82433869386227c30027966a22ee61c918abd4f3142319bbdb2e73e500',
 '10':'16a3ad13fec0a53f3aad6ef4499d235dae8b1a00137d1cb353063aafaaeeb09c',
 '11':'529551d7ecce253b1ab8626956afd076690631bf7d8b51cfcf0c613440a603f5',
 '12':'a63750e5bfae63383a397aefb8af94a3e29b17bf3952b2cf47d8c3fc9b858f69',
 '13':'18ec4761bd1374e9ca6fc63d52048f43e7680cc418d0f61d1fb013026c561907',
}

def stream(path):
    with (gzip.open(path,'rt',newline='') if path.suffix=='.gz' else path.open('r',newline='')) as h:
        yield from csv.DictReader(h,delimiter='\t')

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def check(rows, pred, label, problems):
    if not pred: problems.append(label)
    return bool(pred)

def scan_table(path, expected_n, label, problems, id_column=None, score_prefix='', expected_branch=None):
    if not path.is_file():
        problems.append(label+':FILE_MISSING:'+str(path));return None
    total=0; ids=[]; columns=None; sum_bad=0; nonfinite=0; negative=0; cancer=Counter(); tier=Counter(); uq=Counter()
    for r in stream(path):
        total+=1
        if columns is None:
            columns=list(r)
            if score_prefix != 'DISABLE':
                missing_scores=[score_prefix+c for c in COMPS if score_prefix+c not in r]
                if missing_scores: problems.append(label+':MISSING_COMPONENTS:'+','.join(missing_scores))
        if id_column:
            ids.append(r.get(id_column,''))
        if 'analysis_cancer_class' in r:
            cancer[r['analysis_cancer_class']]+=1
        if 'projection_evidence_tier' in r:
            tier[r['projection_evidence_tier']]+=1
        if 'UQ_evaluable' in r:
            uq[r['UQ_evaluable']]+=1
        if score_prefix is not None and score_prefix!='DISABLE':
            keys=[score_prefix+c for c in COMPS]
            if all(k in r for k in keys):
                try:
                    scores=[float(r[k]) for k in keys]
                    if any(not math.isfinite(x) for x in scores):nonfinite+=1
                    if any(x < -1e-8 for x in scores):negative+=1
                    if abs(sum(scores)-1)>1e-6: sum_bad+=1
                except (ValueError,TypeError):nonfinite+=1
    if expected_n is not None:check(total,total==expected_n,label+f':rows={total},expected={expected_n}',problems)
    if id_column:
        check(total,len(set(ids))==total,label+':nonunique_or_empty_ids',problems)
        check(total,all(ids),label+':empty_ids',problems)
    check(total,sum_bad==0,label+f':fraction_sum_failures={sum_bad}',problems)
    check(total,nonfinite==0,label+f':nonfinite_scores={nonfinite}',problems)
    check(total,negative==0,label+f':negative_scores={negative}',problems)
    return {'path':str(path),'rows':total,'columns':columns,'sample_ids':ids if id_column else None,
            'cancer_counts':dict(cancer),'tier_counts':dict(tier),'uq_counts':dict(uq),
            'sum_bad':sum_bad,'negative_scores':negative,'nonfinite_scores':nonfinite}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--base',required=True);ap.add_argument('--output',default='')
    a=ap.parse_args();base=Path(a.base).resolve()
    frozen=base/'RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE'
    problems=[]; status={}; locations={}; sources={}
    for stage,(rel,fn,version) in STAGES.items():
        d=base/rel; locations[stage]=d;p=d/fn
        if not p.is_file():problems.append(f'Stage {stage}: status missing {p}');continue
        try: s=json.loads(p.read_text())
        except Exception as e:problems.append(f'Stage {stage}:invalid status {e}');continue
        status[stage]=s
        check(None,s.get('status')=='PASS',f'Stage {stage} status={s.get("status")}',problems)
        check(None,s.get('script_version')==version,f'Stage {stage} version={s.get("script_version")} expected={version}',problems)
        check(None,s.get('qc_failures',0)==0,f'Stage {stage} QC failures={s.get("qc_failures")}',problems)
        fp=frozen/SRC[stage]
        if fp.is_file():
            h=sha(fp)
            sources[stage]={'path':str(fp),'sha256':h,'bytes':fp.stat().st_size,
                            'sha256_matches_authoritative_frozen_source':h==AUTHORITATIVE_FROZEN_SOURCE_SHA256[stage]}
            check(None,h==AUTHORITATIVE_FROZEN_SOURCE_SHA256[stage],
                  f'Stage {stage} frozen script SHA256 differs from authoritative verified source: {h}',problems)
        else: problems.append(f'Stage {stage}: missing frozen source {fp}')

    w=base/'../02_model_assets/official_v1.0.3/Bulk_SE_Recovery_W.rds'
    # The installed user's model may be elsewhere; this path is inspected only.
    reference_tables={
      'canonical': ('09','01_canonical_library_manifest.tsv.gz',9465,'sample_id','DISABLE'),
      'anchor': ('09','02_cross_sectional_anchor_manifest.tsv.gz',8962,'sample_id','DISABLE'),
      'pan_scores': ('10','17_primary_CPM_fraction_scores_with_metadata.tsv.gz',9465,'sample_id',''),
      'layer1_sample': ('11','01_sample_level_adequacy_robustness_atlas.tsv.gz',9465,'sample_id','DISABLE'),
      'cs_scores': ('12','17_primary_CPM_cancer_specific_fraction_scores_with_metadata.tsv.gz',9465,'sample_id',''),
      'atlas_core': ('13','03_projection_analysis_ready_core.tsv.gz',9465,'sample_id','DISABLE'),
      'pan_scale_cpm':('10','03_pan_cancer_anchor_scaling_CPM.tsv.gz',499,None,'DISABLE'),
      'pan_scale_tmm':('10','04_pan_cancer_anchor_scaling_TMM_CPM.tsv.gz',499,None,'DISABLE'),
      'pan_scale_uq':('10','05_pan_cancer_anchor_scaling_UQ_CPM.tsv.gz',499,None,'DISABLE'),
      'cs_scale_all':('12','03_cancer_specific_scaling_parameters_all_branches.tsv.gz',499*18*3,None,'DISABLE'),
    }
    table={}
    for key,(stage,name,count,idcol,score_prefix) in reference_tables.items():
        x=scan_table(locations[stage]/name,count,key,problems,id_column=idcol,score_prefix=score_prefix)
        if x:table[key]=x

    for key in ('pan_scores','layer1_sample','cs_scores','atlas_core'):
        if key in table and 'canonical' in table:
            check(None,table[key]['sample_ids']==table['canonical']['sample_ids'],key+':sample_order_mismatch',problems)

    if all(k in table for k in ('pan_scores','cs_scores','atlas_core')):
        # Compare the published Stage 13 score fractions to upstream Stage 10/12.
        ps=locations['10']/'17_primary_CPM_fraction_scores_with_metadata.tsv.gz'
        cs=locations['12']/'17_primary_CPM_cancer_specific_fraction_scores_with_metadata.tsv.gz'
        ac=locations['13']/'03_projection_analysis_ready_core.tsv.gz'
        deltas={'pan':0.0,'cs':0.0}; failed=0
        for x,y,z in zip(stream(ps),stream(cs),stream(ac)):
            if not (x['sample_id']==y['sample_id']==z['sample_id']):
                failed+=1;continue
            for c in COMPS:
                for label,r in (('pan',x),('cs',y)):
                    try:d=abs(float(r[c])-float(z[label+'_'+c]));deltas[label]=max(d,deltas[label])
                    except (KeyError,ValueError,TypeError):failed+=1
        check(None,failed==0,'Stage13:score_comparison_missing_or_sample_mismatch:'+str(failed),problems)
        for label,d in deltas.items():
            check(None,d<=1e-8,'Stage13:'+label+'_score_mismatch_maxabs='+str(d),problems)
        score_replay=deltas
    else:score_replay=None

    if 'anchor' in table and 'canonical' in table:
        check(None,set(table['anchor']['sample_ids']).issubset(set(table['canonical']['sample_ids'])),
              'anchor_not_subset_of_canonical',problems)

    if 'layer1_sample' in table:
        lc=table['layer1_sample']['tier_counts']
        for t,n in [('Tier_A_fit_and_both_normalizations',4200),('Tier_B_fit_and_TMM_supported',1026),('Tier_C_exploratory_or_sensitive',4239)]:
            check(None,lc.get(t)==n,f'layer1_tier_{t}_count={lc.get(t)}_expected={n}',problems)
        uc=table['layer1_sample']['uq_counts']
        if uc:check(None,uc.get('TRUE',0)+uc.get('True',0)+uc.get('1',0)==9302,
                    f'layer1_UQ_true={uc} expected9302',problems)

    # Exact frozen Stage 11 study-level outlier provenance: this is the only downstream
    # computation involving the corrected stats::mad(center=...) argument.
    expected_poor_fit={
        ('lung_and_pleural_cancer','GSE190266'),
        ('lung_and_pleural_cancer','GSE126044'),
        ('colorectal_cancer','GSE299799'),
        ('colorectal_cancer','GSE183984'),
        ('colorectal_cancer','GSE152395'),
    }
    study_path=locations['11']/'03_study_level_transferability_summary.tsv'
    outlier_check={'path':str(study_path),'expected_count':5,'expected_studies':sorted([list(x) for x in expected_poor_fit])}
    if study_path.is_file():
        study_rows=list(stream(study_path))
        check(None,len(study_rows)==489,'Stage11:expected_489_cancer_study_rows_got_'+str(len(study_rows)),problems)
        keys=('analysis_cancer_class','primary_gse_id','poor_fit_study_outlier','adequacy_robust_z_within_cancer','n_samples')
        check(None,bool(study_rows) and all(k in study_rows[0] for k in keys),
              'Stage11:missing_study_outlier_fields',problems)
        if study_rows and all(k in study_rows[0] for k in keys):
            positives=[r for r in study_rows if r['poor_fit_study_outlier'].strip().upper() in ('TRUE','T','1')]
            found={(r['analysis_cancer_class'],r['primary_gse_id']) for r in positives}
            outlier_check.update({'rows':len(study_rows),'found_count':len(positives),'found_studies':sorted([list(x) for x in found])})
            check(None,found==expected_poor_fit and len(positives)==5,
                  'Stage11:poor_fit_outliers_disagree_with_frozen_contract',problems)
            for r in positives:
                try:
                    z=float(r['adequacy_robust_z_within_cancer']); ns=int(r['n_samples'])
                    check(None,math.isfinite(z) and z<=-2.5 and ns>=5,
                          'Stage11:invalid_outlier_statistics:'+r['primary_gse_id'],problems)
                except (ValueError,TypeError):
                    problems.append('Stage11:non_numeric_outlier_statistics:'+r['primary_gse_id'])
    else:
        problems.append('Stage11:study_outlier_table_missing:'+str(study_path))

    if 'atlas_core' in table:
        columns=table['atlas_core']['columns']
        for prefix in ('pan_','cs_'):
            for k in COMPS:
                check(None,prefix+k in columns,f'atlas_core_missing_{prefix+k}',problems)
        for col in ('cross_cancer_primary_eligible_strict','within_cancer_primary_eligible_strict',
                    'recommended_cross_cancer_score_space','recommended_within_cancer_score_space'):
            check(None,col in columns,f'atlas_core_missing_{col}',problems)

    for k in ('pan_scale_cpm','pan_scale_tmm','pan_scale_uq'):
        if k in table:
            cols=table[k]['columns']
            for col in ('model_gene','weighted_anchor_mean_log2','weighted_anchor_sd_log2'):
                check(None,col in cols,k+':missing_'+col,problems)

    out={
      'audit':'ARCHS4_09_13_FROZEN_OUTPUT_CONTRACT_READONLY_V1_1',
      'result':'PASS_REFERENCE_OUTPUT_CONTRACT' if not problems else 'HOLD_REVIEW_PROBLEMS',
      'problems':problems,
      'status_summary':{k:{a:status[k].get(a) for a in ('status','script_version','qc_failures','canonical_samples','cancer_classes','directional_features','Tier_A_samples','Tier_B_samples','Tier_C_samples','pan_Tier_A_samples','pan_Tier_B_samples','pan_Tier_C_samples')} for k in status},
      'frozen_source_sha256':sources,
      'stage13_upstream_score_replay_max_abs':score_replay,
      'stage11_study_outlier_replay':outlier_check,
      'table_summary':{k:{a:v for a,v in x.items() if a not in ('sample_ids','columns','cancer_counts')} for k,x in table.items()},
      'note':'Verifies frozen ARCHS4 metadata/score output and Stage11 outlier contracts only. Does not prove execution provenance from historical job logs, nor recount3 model transfer or cross-assay equivalence.'
    }
    print(json.dumps(out,indent=2,sort_keys=True))
    if a.output:
        p=Path(a.output).resolve()
        if frozen==p or frozen in p.parents:raise SystemExit('Refusing to write under frozen release')
        if p.exists():raise SystemExit('Refusing to overwrite existing output '+str(p))
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print('WROTE',p)
    return 0 if not problems else 2
if __name__=='__main__':raise SystemExit(main())
