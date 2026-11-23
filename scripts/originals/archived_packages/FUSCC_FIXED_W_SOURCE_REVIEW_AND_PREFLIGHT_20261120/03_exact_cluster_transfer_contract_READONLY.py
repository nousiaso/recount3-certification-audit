#!/usr/bin/env python3
"""Read-only metadata and reference gate for frozen ARCHS4/recount3 transfer.

Usage: python 03_exact_cluster_transfer_contract_READONLY.py /maps/projects/.../ARCHS4_V2_6
No output files, no modification of originals, no large matrix reads.
"""
import sys, csv, gzip, json, hashlib
from pathlib import Path
from collections import Counter

ROOT=Path(sys.argv[1]).resolve()
REL=ROOT/'RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0'
P=REL/'PROJECT_TREE'
REC=ROOT/'26_recount3_external_primary_v1'
EXPECT_W='39a96b598ec74ccbe67bea49254ac94ae568ca22ec162241ab7849976b52cdcf'
print('ROOT',ROOT,'RELEASE',REL,'RECOUNT3',REC,sep='\n')

def existing(*rel_candidates):
    for rel in rel_candidates:
        p=Path(rel) if str(rel).startswith('/') else ROOT/rel
        if p.is_file():return p
    return None

def table_header_and_count(path, keep_fields=()):
    if not path or not path.is_file():print('MISSING_TABLE',path);return None,None
    opener=gzip.open if path.name.endswith('.gz') else open
    with opener(path,'rt',newline='') as f:
        rd=csv.DictReader(f,delimiter='\t')
        cols=rd.fieldnames or []
        n=0; counts={x:Counter() for x in keep_fields if x in cols}
        examples=[]
        for row in rd:
            n+=1
            if n<=2: examples.append({k:row.get(k,'') for k in cols[:12]})
            for k in counts:counts[k][row.get(k,'')]+=1
    print('TABLE',path,'ROWS',n,'COLS',cols)
    print('EXAMPLES',examples)
    for k,v in counts.items():print('COUNTS',k,v.most_common(12))
    return cols,n

def locate_original(prefix,file):
    return existing(f'{prefix}/{file}',f'RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE/{prefix}/{file}')

model=REL/'02_MODEL/Bulk_SE_Recovery_W.rds'
if model.is_file():
    got=hashlib.sha256(model.read_bytes()).hexdigest()
    print('FROZEN_W_SHA256',got,'EXPECTED_MATCH',got==EXPECT_W)
else:print('MISSING_W',model)

for folder,filename,check in [
('10_preprojection_normalization_v1_0_1','08_model_gene_availability.tsv',()),
('11_pan_cancer_fixed_anchor_projection_v1_0_4','03_pan_cancer_anchor_scaling_CPM.tsv.gz',()),
('11_pan_cancer_fixed_anchor_projection_v1_0_4','04_pan_cancer_anchor_scaling_TMM_CPM.tsv.gz',()),
('11_pan_cancer_fixed_anchor_projection_v1_0_4','05_pan_cancer_anchor_scaling_UQ_CPM.tsv.gz',()),
('13_cancer_specific_fixed_anchor_projection_v1_0_3','03_cancer_specific_scaling_parameters_all_branches.tsv.gz',('analysis_cancer_class','branch')),
]:
    path=locate_original(folder,filename)
    table_header_and_count(path,check)

for folder,filename in [
('11_pan_cancer_fixed_anchor_projection_v1_0_4','17_primary_CPM_fraction_scores_with_metadata.tsv.gz'),
('13_cancer_specific_fixed_anchor_projection_v1_0_3','17_primary_CPM_cancer_specific_fraction_scores_with_metadata.tsv.gz')
]:
    path=locate_original(folder,filename)
    if path:
        with gzip.open(path,'rt') as f:print('SCORE_HEADER',path,f.readline().rstrip('\n')[:600])
    else:print('MISSING_SCORE',folder,filename)

for s in ('SRP074349','SRP157974'):
    path=REC/f'stage2_output/studies/{s}/selected_run_coverage_G026.h5'
    print('RECOUNT3_H5',s,'EXISTS',path.is_file(), 'SIZE_BYTES',path.stat().st_size if path.is_file() else '-')
    if not path.is_file():continue
    try:
        import h5py
        with h5py.File(path,'r') as h:
            print('H5_ROOT_KEYS',list(h.keys()))
            print('H5_ROOT_ATTRIBUTES',{k:str(v)[:100] for k,v in h.attrs.items()})
            items=[]
            def visit(name,obj):
                if len(items)>=36:return
                if isinstance(obj,h5py.Dataset):
                    items.append((name,tuple(obj.shape),str(obj.dtype),list(obj.attrs.keys())))
                elif len(items)<10:
                    items.append((name,'GROUP','',list(obj.attrs.keys())))
            h.visititems(visit)
            for item in items:print('H5_AXIS_ITEM',item)
    except ImportError:print('H5PY_MISSING_SCHEMA_NOT_CHECKED')
    except Exception as e:print('H5_SCHEMA_ERROR',type(e).__name__,str(e))

p=REC/'stage5_output/03_bridge_165_availability.tsv'
if p.is_file():
    with p.open(newline='') as f:
        rd=csv.DictReader(f,delimiter='\t');cols=rd.fieldnames or []
        rows=list(rd)
    print('BRIDGE_SCHEMA',cols,'ROWS',len(rows))
    print('BRIDGE_STUDIES',Counter(r.get('study','') for r in rows))
    print('BRIDGE_CLASSES',Counter(r.get('bridge_class','') for r in rows))
else:print('BRIDGE_LIST_MISSING',p)
print('READ_ONLY_FROZEN_TRANSFER_CONTRACT_COMPLETE')
