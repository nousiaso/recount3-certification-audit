"""V6 SRP218935: read-only reproducibility check of original source, V5 selection and run QC.
No clinical admission, patient identifiers printed, file writes or model inference.
This is a NEW consolidated diagnostic, not a claim about bytes of prior inline scripts.
"""
import csv, hashlib, io, re, sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

root = Path(sys.argv[1])
rec = root / '26_recount3_external_primary_v1'
repo = root / 'RECOUNT3_CERTIFICATION_GIT'
v5 = rec / 'clinical_specimen_patient_v5_unsigned_20261122'

sources = {
    'inventory': (v5 / '01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_CROSSWALK_UNSIGNED.tsv',
                  'ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a'),
    'screen': (v5 / '03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv',
               '06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49'),
    'v6': (repo / 'diagnostics/20261123_v6_evidence_precedence_dryrun/proposed_evidence_2055_UNSIGNED.tsv',
           '08efe2c8f469edcbdfb98231dc87f705d801434e78dcbea85f57c4b94da88eb0'),
    'selected': (repo / 'diagnostics/20261124_v6_remaining1097_source_v1/1097_original_source_evidence_LOCAL_ONLY.tsv',
                 '194a214e48b0105e9e6a40b1b72a8d6c5c7cf08abf4a728d274f80575f04bf1c'),
    'readiness': (rec / 'cohort_readiness_v1/preview/01_per_library_readiness_3602.tsv',
                  '404c66603c4d6934e42cfbd6a21eaf0233b02ead9edb568193ba26b144b33d0f'),
    'run_holds': (root / '29_recount3_adjudication_20261119/RECOUNT3_CERTIFICATION_ADJUDICATION_20261119/REVIEW_OUTPUT/05_run_hold_queue_409.tsv',
                  'd275c77591e7b93582708746e45d69ae817c8e2d13b3bb5b7f63e2e56fa56f0c'),
}

def load(name):
    path, expected = sources[name]
    assert path.is_file(), f'Missing documented input: {name}'
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    assert h == expected, f'Input SHA256 changed: {name}'
    with path.open(encoding='utf-8-sig', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        assert reader.fieldnames, f'Empty schema: {name}'
        return list(reader)

D = {k: load(k) for k in sources}
assert tuple(len(D[k]) for k in ('inventory','screen','v6','selected','readiness','run_holds')) == (5241,2055,2055,1097,3602,409)
STUDY = 'SRP218935'
by_study = lambda k, idfield: {r[idfield]: r for r in D[k] if r['study'] == STUDY}
original = by_study('inventory', 'library_id')
screen = by_study('screen', 'library_id')
v6 = by_study('v6', 'library_id')
source = by_study('selected', 'experiment_acc')
ready = by_study('readiness', 'library_id')
assert len(original) == 76 and len(screen) == len(v6) == len(source) == 75 and len(ready) == 76
selected = set(source)
held = set(original) - set(screen)
assert selected == set(screen) == set(v6) and len(held) == 1
held_exp = next(iter(held))
assert original[held_exp]['stage41_partition_NOT_ADMISSION'] == 'QC_HOLD_NOT_ADMITTED'
assert original[held_exp]['recommended_provisional_disposition'] == 'HOLD_TECHNICAL_QC'

fields = ['experiment_accession', 'run_accession', 'sample_accession', 'secondary_sample_accession',
          'sample_title', 'experiment_title', 'sample_alias', 'library_strategy']
url = 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + urlencode({
    'accession': 'PRJNA553096', 'result': 'read_run', 'fields': ','.join(fields), 'format': 'tsv'
})
with urlopen(Request(url, headers={'User-Agent':'recount3-source-certification-readonly/1.0'}), timeout=60) as fh:
    ena = list(csv.DictReader(io.StringIO(fh.read().decode('utf-8-sig')), delimiter='\t'))
assert ena and set(fields) <= set(ena[0]), 'ENA fields changed'
by_exp = defaultdict(list)
for r in ena:
    if r['experiment_accession']:
        by_exp[r['experiment_accession']].append(r)

def title_class(rows):
    titles = [((r['sample_title'] or '')+' '+(r['experiment_title'] or '')).casefold() for r in rows]
    if any('normal' in t for t in titles):
        return 'NORMAL_TITLE'
    if any('primary' in t and ('breast' in t or 'brca' in t) for t in titles):
        return 'PRIMARY_TITLE'
    return 'OTHER_OR_UNRESOLVED'

cats = {k: title_class(v) for k,v in by_exp.items()}
primary = {k for k,v in cats.items() if v == 'PRIMARY_TITLE'}
normal = {k for k,v in cats.items() if v == 'NORMAL_TITLE'}
assert len(ena) == len(by_exp) == 95 and len(primary) == 76 and len(normal) == 19
assert selected <= primary and primary - selected == held and not (normal & set(original))

# Compare sample accessions and source material IDs to the ENA sample aliases.
def norm(x):
    return ' '.join(str(x or '').strip().casefold().split())
def attrs(s):
    d = defaultdict(set)
    for item in (s or '').split('|'):
        if ';;' in item:
            k,v = item.split(';;',1)
            if k.strip() and v.strip(): d[norm(k)].add(norm(v))
    return d

def trailing_number(x):
    m = re.search(r'(\d+)\s*$', x)
    return int(m.group(1)) if m else None

sample_matches = 0
alias_matches = 0
material_ids = []
for exp in selected:
    src = source[exp]
    erows = by_exp[exp]
    possible_samples = {x for r in erows for x in
                        (r['sample_accession'],r['secondary_sample_accession']) if x}
    sample_matches += src['source_sample_acc_LOCAL_ONLY'] in possible_samples
    a = attrs(src['source_sample_attributes_LOCAL_ONLY'])
    vals = a.get('source_materal_id', set())
    assert len(vals) == 1
    material = next(iter(vals))
    material_ids.append(material)
    ali = {trailing_number(r['sample_alias']) for r in erows}
    ali.discard(None)
    alias_matches += len(ali) == 1 and trailing_number(material) in ali
assert len(set(material_ids)) == 75 and all(x.isdigit() for x in material_ids)
assert sample_matches == alias_matches == 75

status = Counter(ready[x]['historical_run_hold'] for x in original)
assert status == {'NO':75, 'YES':1}
assert all(ready[x]['historical_run_hold'] == 'NO' for x in selected)
assert ready[held_exp]['historical_run_hold'] == 'YES'
matched_holds = [r for r in D['run_holds'] if r.get('library_id') == held_exp
                 and (not r.get('study') or r['study'] == STUDY)]
assert len(matched_holds) == 1
h = matched_holds[0]
expected_reasons = {'TOP10_COVERAGE_GT_0_50', 'MAPPED_LEN_LT_40BP', 'WITHIN_STUDY_OUTLIER'}
reasons = set(h['reason_codes'].split(';'))
assert expected_reasons <= reasons
assert float(h['top10_coverage_fraction']) > .50
assert float(h['average_mapped_length_bp']) < 40

print('V6_SRP218935_REPRODUCED_DIAGNOSTIC')
print('FROZEN_INPUT_SHA256_ALL_PASS', len(sources))
print('ENA_RUN_ROWS',len(ena))
print('ENA_DISTINCT_EXPERIMENTS',len(by_exp))
print('ENA_PRIMARY_TITLE',len(primary))
print('ENA_NORMAL_TITLE',len(normal))
print('SELECTED_EXACT_ENA_MATCHES',len(selected & set(by_exp)))
print('SELECTED_SRA_SAMPLE_ALIAS_MATCHES',sample_matches)
print('MATERIAL_VS_ENA_ALIAS_MATCHES',alias_matches)
print('V5_ALL_STUDY',len(original))
print('V5_PRIMARY_SCREEN',len(screen))
print('EXCLUDED_PRIMARY_QC_HOLD',len(held))
print('NORMAL_EXPERIMENTS_OUTSIDE_V5',len(normal - set(original)))
print('READINESS_HISTORICAL_RUN_HOLD',dict(sorted(status.items())))
print('HELD_RUN_TOP10_COVERAGE_FRACTION',h['top10_coverage_fraction'])
print('HELD_RUN_AVERAGE_MAPPED_LENGTH_BP',h['average_mapped_length_bp'])
print('HELD_RUN_GENE_DETECTION_FRACTION',h['gene_detection_fraction'])
print('HELD_RUN_REASON_CODES',';'.join(sorted(reasons)))
for field in ('historical_tissue_grade','patient_grade_UNREVIEWED',
              'treatment_grade_UNREVIEWED', 'unsigned_evidence_bucket'):
    print('SELECTED_'+field.upper(), dict(sorted(Counter(v6[x][field] for x in selected).items())))
print('CLINICAL_APPROVALS',0)
print('NO_PATIENT_OR_SAMPLE_IDENTIFIERS_PRINTED')
print('NO_FILES_WRITTEN_BY_DIAGNOSTIC')
print('NO_FROZEN_DATA_CHANGED')
print('V6_SRP218935_REPRODUCED_DIAGNOSTIC_COMPLETE')
