#!/usr/bin/env python3
"""Read-only FUSCC hold and metadata check; stdlib only, no outputs written."""
import csv, gzip, re, sys
from pathlib import Path
from collections import Counter
base=Path(sys.argv[1]); root=base/'26_recount3_external_primary_v1'
rd=root/'cohort_readiness_v1/preview/01_per_library_readiness_3602.tsv'
raw=root/'raw_studies/SRP157974/metadata/sra.sra.SRP157974.MD.gz'
def load(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8-sig',errors='replace',newline='') as f:
        r=csv.DictReader(f,delimiter='\t'); return list(r),r.fieldnames
readiness,cols=load(rd)
f=[r for r in readiness if r.get('study')=='SRP157974']
h=[r for r in f if r.get('historical_run_hold')=='YES']
non=[r for r in f if r.get('historical_run_hold')!='YES']
print('=== A. Original run/readiness evidence ===')
print('FUSCC',len(f),'held',len(h),'nonheld',len(non))
assert (len(f),len(h),len(non))==(360,26,334)
print('Readiness columns:',cols)
fields=[c for c in cols if re.search(r'hold|run|qc|qual|eligib|fail|reason|spot|base|coverage|depth|flag|status|source|map|library',c,re.I)]
print('Candidate QC/hold fields:',fields)
for c in fields:
    counts=Counter((r.get(c) or '').strip() for r in f)
    if len(counts)<12:
        print('FIELD',c,'values',dict(counts))
print('=== B. 26 held experiment records ===')
for r in h:
    vals={c:(r.get(c) or '')[:100] for c in fields if r.get(c)}
    print(r.get('library_id'),vals)
print('=== C. Source specimen labels and protocol ===')
rawrows,rawcols=load(raw)
rx={}
for r in rawrows:
    e=r.get('experiment_acc',''); rx.setdefault(e,[]).append(r)
print('Distinct source experiments',len(rx),'candidate source match',sum(r.get('library_id') in rx for r in f))
print('Source columns:',rawcols)
print('Source sample name first 5:',[(r.get('library_id'),rx[r['library_id']][0].get('sample_name')) for r in f[:5]])
print('Candidate source labels:',dict(Counter('TUMOUR' if 'tumortissue' in rx[r['library_id']][0].get('sample_name','').lower() else 'OTHER' for r in f)))
print('Original library selection:',dict(Counter(rx[r['library_id']][0].get('library_selection','') for r in f)))
print('Protocols present:',sum(bool(rx[r['library_id']][0].get('library_construction_protocol','').strip()) for r in f))
print('=== D. Original QC evidence directories ===')
for rel in ('cohort_readiness_v1','stage2_output','stage3_output','stage5_output'):
    p=root/rel
    print(rel,'exists',p.exists())
    if p.exists():
        print('Files',*[str(x.relative_to(root)) for x in sorted(p.rglob('*')) if x.is_file()][:18],sep='\n  ')
print('READ_ONLY_FUSCC_HOLD_AUDIT_COMPLETE; no original files modified')
