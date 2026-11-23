#!/usr/bin/env python3
"""Targeted READ-ONLY model-gene and atlas identity probe for HPC.

No expression values are loaded, no datasets/models are modified, no gene remapping
is accepted automatically. Writes only candidate reports into a new directory.
"""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,csv,json,re,hashlib,gzip

def read(path):
    with open(path,newline='',encoding='utf-8') as f: return list(csv.DictReader(f,delimiter='\t'))

def base(x):return re.sub(r'\.\d+$','',str(x).strip())
def clean(x):
    if isinstance(x,bytes):return x.decode('utf-8',errors='replace').strip()
    return str(x).strip()
def write(path,rows,columns):
    with open(path,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)
def hashfile(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(4<<20),b''):h.update(chunk)
    return h.hexdigest()

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--root',type=Path,required=True)
    a.add_argument('--source-evidence',type=Path,help='optional previously extracted 3602 accession source evidence TSV')
    a.add_argument('--out',type=Path,required=True)
    v=a.parse_args();root=v.root.resolve();out=v.out.resolve()
    if out.exists():raise ValueError('Refuse to overwrite existing output: '+str(out))
    refs={r['role']:r for r in read(root/'stage4_contract_output/03_frozen_reference_file_contract.tsv')}
    h5path=Path(refs['ARCHS4_source_H5_EXTERNAL']['path'])
    atlas=Path(refs['frozen_9465_projection_atlas']['path'])
    mapfile=root/'stage5_output/00_frozen_W_to_G026_candidate_mapping.tsv'
    genefile=root/'stage5_output/00_gene_row_mapping.tsv'
    wmap=read(mapfile);gro=read(genefile)
    if len(wmap)!=1497 or len(gro)!=63856:raise ValueError('Unexpected W-map or G026 rows, investigate upstream first')
    if not {'row_zero_based','source_gene_id','gene_symbol','mapping_status'}.issubset(gro[0]):raise ValueError('G026 mapping schema drift')
    if any(int(r['row_zero_based'])!=i for i,r in enumerate(gro)):raise ValueError('G026 row-order drift')
    gby=defaultdict(list)
    for r in gro:gby[clean(r['gene_symbol'])].append(r)
    wby=defaultdict(list)
    for r in wmap:wby[r['W_model_gene']].append(r)
    if len(wby)!=499 or any(len(rows)!=3 for rows in wby.values()):raise ValueError('W 499x3 contract drift')
    held=sorted(k for k,rows in wby.items() if any(r['mapping_status']!='UNAMBIGUOUS_ONE_TO_ONE' for r in rows))
    # 11 other held are unambiguous in G026 and ambiguous only in original H5
    # determine all 34 from original frozen W diagnostic, not from wmap alone.
    bysym=defaultdict(list)
    import h5py
    with h5py.File(h5path,'r') as f:
        sym=f['meta/genes/symbol'][:]
        ens=f['meta/genes/ensembl_gene'][:]
        if len(sym)!=67186 or len(ens)!=67186:raise ValueError('Original ARCHS4 gene axis drift')
        if f['data/expression'].shape[0]!=67186:raise ValueError('Original ARCHS4 gene count drift')
    for i,(s,e) in enumerate(zip(sym,ens)):
        bysym[clean(s)].append((i,base(clean(e))))
    for s,entries in wby.items():
        if s in held:continue
        ids={(base(r['matched_source_gene_id']),r['source_row_zero_based']) for r in entries}
        if len(ids)!=1:
            held.append(s);continue
        ensg=next(iter(ids))[0]
        if len([x for x in bysym[s] if x[1]==ensg])!=1:held.append(s)
    held=sorted(set(held))
    if len(held)!=34:raise ValueError('Expected 34 held W genes, observed '+str(len(held)))
    out.mkdir(parents=True)
    report=[]
    for s in held:
        candidates=gby.get(s,[])
        h5=byssym=bysym.get(s,[])
        both=[]
        for a1 in candidates:
            for idx,eid in h5:
                if base(a1['source_gene_id'])==eid:
                    both.append((a1['row_zero_based'],base(a1['source_gene_id']),idx))
        report.append({
          'W_model_gene':s,
          'G026_W_map_statuses':';'.join(sorted(set(x['mapping_status'] for x in wby[s]))),
          'G026_symbol_row_count':len(candidates),
          'G026_symbol_gene_IDs':';'.join(sorted(set(base(x['source_gene_id']) for x in candidates))),
          'G026_candidate_row_ID_pairs':';'.join(x['row_zero_based']+':'+base(x['source_gene_id'])+':'+x['mapping_status'] for x in candidates),
          'H5_original_symbol_row_count':len(h5),
          'H5_original_symbol_row_ID_pairs':';'.join(str(idx)+':'+eid for idx,eid in h5),
          'exact_G026_H5_ID_pair_count':len(both),
          'exact_G026_H5_ID_pair_candidates':';'.join(':'.join(map(str,x)) for x in both),
          'manual_mapping_decision':'HOLD',
          'scientific_warning':'MULTIPLE_IDENTITY_POSSIBILITIES_REQUIRE_REFERENCE_RESOLUTION;NO_AUTOMATIC_MAPPING'
        })
    write(out/'01_W34_metadata_candidates.tsv',report,list(report[0]))
    opener=gzip.open if str(atlas).endswith('.gz') else open
    with opener(atlas,'rt',encoding='utf-8',newline='') as f:
        rr=csv.DictReader(f,delimiter='\t');header=rr.fieldnames or []
        atlasrows=list(rr)
    if len(atlasrows)!=9465 or 'sample_id' not in header:raise ValueError('Frozen atlas sample IDs/rows unexpectedly changed')
    sampleids=[r['sample_id'] for r in atlasrows]
    if len(set(sampleids))!=9465:raise ValueError('Frozen atlas sample IDs not unique')
    geo_atlas={x for x in sampleids if re.fullmatch(r'GSM\d+',x)}
    srcids=[];direct=[]
    if v.source_evidence:
        evidence=read(v.source_evidence)
        if len(evidence)!=3602:raise ValueError('SRA source evidence not 3602 candidate libraries')
        for a2 in evidence:
            gsms=set(re.findall(r'\bGSM\d+\b',(a2['source_experiment_attributes']+' '+a2['source_experiment_title']).upper()))
            for g in sorted(gsms):
                srcids.append(g)
                if g in geo_atlas:
                    direct.append({'library_id':a2['library_id'],'study':a2['study'],
                         'source_GSM':g,'frozen_atlas_sample_id':g,
                         'overlap_review':'CANDIDATE_DIRECT_GSM_ALIAS_MATCH_NOT_FULL_9465_CROSSWALK'})
        write(out/'02_GSM_direct_overlap_candidates.tsv',direct,
              ['library_id','study','source_GSM','frozen_atlas_sample_id','overlap_review'])
    summary={
       'status':'W34_CANDIDATES_AND_ATLAS_SCHEMA_ONLY_NOT_RELEASE_CERTIFICATION',
       'n_unresolved_W_genes':len(held),'n_with_exact_G026_H5_candidate_pairs':sum(x['exact_G026_H5_ID_pair_count']>0 for x in report),
       'n_with_single_exact_G026_H5_pair_candidate':sum(x['exact_G026_H5_ID_pair_count']==1 for x in report),
       'n_with_zero_exact_candidate_pairs':sum(x['exact_G026_H5_ID_pair_count']==0 for x in report),
       'frozen_atlas_path':str(atlas),'frozen_atlas_columns':header,
       'n_frozen_atlas_rows':len(atlasrows),
       'n_sample_ids_that_are_direct_GSM':len(geo_atlas),
       'n_recount3_candidate_GSM_accessions_seen':len(set(srcids)),
       'n_direct_GSM_overlap_candidate_libraries':len(set(x['library_id'] for x in direct)),
       'original_W_499_projection_validated':False,
       'crossatlas_independence_certified':False,
       'note':'Only identifier metadata and atlas TSV accessed; no ARCHS4 expression values read.'}
    (out/'03_readonly_probe_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (out/'SHA256SUMS.txt').open('w') as f:
        for p in sorted(out.iterdir()):
            if p.is_file() and p.name!='SHA256SUMS.txt':f.write(f'{hashfile(p)}  {p.name}\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
