#!/usr/bin/env python3
"""Evidence-only V5/Stage41 reconciliation. No clinical admission."""
import sys, csv, gzip, hashlib, json
from pathlib import Path
from collections import Counter

root = Path(sys.argv[1])
rec = root / "26_recount3_external_primary_v1"
v5 = rec / "clinical_specimen_patient_v5_unsigned_20261122"
out = root / "RECOUNT3_CERTIFICATION_GIT/diagnostics/20261123_v6_evidence_precedence_dryrun"

expected = {
 'all': (
  '01_ALL_5241_WITH_EXACT_SPECIMEN_PATIENT_SOURCE_CROSSWALK_UNSIGNED.tsv',
  'ea3064b884a823c854bc375e5fe2b38f8b4f187d995ffb6c30674c700e618f7a'),
 'screen': (
  '03_2055_ADMIN_PRIMARY_PROPOSED_RESTRICTIONS.tsv',
  '06b7d0443c99b365d0c76671e5bf18dbc9c2f6d5fc09ff9c167d0e2e9cf49e49'),
 'normal': (
  '04_21_SOURCE_SUPPORTED_NON_TUMOUR_IN_PRIMARY_SCREEN.tsv',
  '58d4a4dc718e1b9cc85b31f6dfe6c539536dac44d81170656652ce80b4fa78c5'),
 'sarcoma': (
  '05_52_SARCOMA_GEO_METASTATIC_DISEASE_FLAGS_IN_PRIMARY_SCREEN.tsv',
  '0cf102698e2c94faf660d647507e42ab0261eff35a64dbaadeff12e51018dcec'),
 'studies': (
  '11_176_STUDY_SOURCE_LEDGER_UNCHANGED_FROM_V4.tsv',
  'b9e1ae8d91307cc47c5411f1b8b201db2a6cfe92a8db0d0f7d79917abaef2db7'),
}

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def readtsv(p):
    with p.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))

t = {}
for name, (filename, digest) in expected.items():
    p = v5 / filename
    assert p.is_file() and sha(p) == digest, (
        f'Input absent or SHA changed: {p}'
    )
    t[name] = readtsv(p)

screen = t['screen']
by = {r['library_id']: r for r in screen}

assert len(t['all']) == len({
    r['library_id'] for r in t['all']
}) == 5241
assert len(screen) == len(by) == 2055
assert len({r['study'] for r in screen}) == 77
assert len(t['studies']) == len({
    r['study'] for r in t['studies']
}) == 176
assert all(
    r['FINAL_CLINICAL_APPROVAL'] == 'NO_UNSIGNED'
    for r in screen
)

v5_normal = {
    r['library_id'] for r in t['normal']
}
old_normal = {
    r['library_id'] for r in screen
    if r['tissue_evidence_grade'] == 'SOURCE_NON_TUMOUR'
    and r['recommended_provisional_disposition']
        == 'EXCLUDE_PRIMARY_NON_TUMOUR'
}
sarcoma = {
    r['library_id'] for r in t['sarcoma']
}

assert len(v5_normal) == 21
assert len(old_normal) == 2
assert len(sarcoma) == 52
assert v5_normal <= by.keys()
assert sarcoma <= by.keys()

normal = v5_normal | old_normal
assert len(normal) == 22
assert normal - v5_normal == {'SRX106046'}
assert 'normal' in by['SRX106046'][
    'stage41_tissue_source_evidence'
].lower()
assert not (sarcoma & normal)

# Confirm exact SRP060016 sample tissue attributes.
# Do not infer specimen class from the study-wide description.
sra = rec / (
    'raw_studies/SRP060016/metadata/'
    'sra.sra.SRP060016.MD.gz'
)

with gzip.open(sra, 'rt', newline='',
               errors='replace') as f:
    raw = list(csv.DictReader(f, delimiter='\t'))

assert len(raw) == 20

labels = {}
for r in raw:
    accession = r['experiment_acc']
    attrs = dict(
        x.split(';;', 1)
        for x in r['sample_attributes'].split('|')
        if ';;' in x
    )
    assert accession not in labels
    labels[accession] = attrs.get('tissue', '')

assert Counter(labels.values()) == Counter({
    'Normal colon': 5,
    'Primary colon cancer': 5,
    'Liver metastases': 5,
    'normal liver': 5,
})

primary = {
    k for k, value in labels.items()
    if value == 'Primary colon cancer'
}

assert primary == {
    f'SRX227306{i}' for i in range(2, 7)
}
assert all(
    k in by
    and by[k]['study'] == 'SRP060016'
    and by[k]['tissue_evidence_grade']
        == 'TISSUE_CONFLICT_UNRESOLVED'
    for k in primary
)

# GEO metastasis describes clinical disease status,
# not the anatomical site of the sequenced specimen.
assert all(
    by[k]['v5_source_metastasis_annotation']
    == 'GEO_METASTASIS_YES_LESION_SITE_NOT_PROVEN'
    for k in sarcoma
)

counts = Counter()
result = []

for r in screen:
    k = r['library_id']

    if k in normal:
        bucket = 'EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE'
        basis = (
            'V5_EXACT_SRX' if k in v5_normal
            else 'PRIOR_STAGE41_EXACT_SRX'
        )
    elif k in primary:
        bucket = 'PRIMARY_TUMOUR_SPECIMEN_ONLY_OTHER_GATES_OPEN'
        basis = 'SRA_EXPERIMENT_TISSUE_ATTRIBUTE'
    elif k in sarcoma:
        bucket = 'HOLD_SAMPLED_SITE_METASTASIS_DISEASE_FLAG'
        basis = 'V5_DISEASE_NOT_LESION_SITE'
    elif r['v5_specimen_evidence_class']:
        bucket = 'V5_TUMOUR_OR_CONTEXT_SPECIMEN_GATES_OPEN'
        basis = 'V5_SOURCE_NOT_FULL_ADMISSION'
    else:
        bucket = 'HISTORICAL_EVIDENCE_REVIEW_PENDING'
        basis = 'STAGE41_V4_NO_NEW_V5_CLASS'

    counts[bucket] += 1

    result.append({
        'study': r['study'],
        'library_id': k,
        'unsigned_evidence_bucket': bucket,
        'evidence_basis': basis,
        'historical_tissue_grade': r['tissue_evidence_grade'],
        'v5_specimen_class': r['v5_specimen_evidence_class'],
        'treatment_grade_UNREVIEWED': r['treatment_evidence_grade'],
        'patient_grade_UNREVIEWED': r['patient_evidence_grade'],
        'FINAL_APPROVAL': 'NO_UNSIGNED',
    })

assert counts == Counter({
    'EXCLUDE_PRIMARY_NON_TUMOUR_EVIDENCE': 22,
    'PRIMARY_TUMOUR_SPECIMEN_ONLY_OTHER_GATES_OPEN': 5,
    'HOLD_SAMPLED_SITE_METASTASIS_DISEASE_FLAG': 52,
    'V5_TUMOUR_OR_CONTEXT_SPECIMEN_GATES_OPEN': 280,
    'HISTORICAL_EVIDENCE_REVIEW_PENDING': 1696,
}), counts

summary = {
    'status': 'EVIDENCE_ONLY_NOT_CERTIFIED',
    'primary_screen': 2055,
    'v5_non_tumour': 21,
    'older_non_tumour': 2,
    'cumulative_union': 22,
    'additional_non_tumour_experiment': 'SRX106046',
    'remaining_admin_NOT_APPROVED': 2033,
    'exact_primary_SRP060016': 5,
    'sarcoma_sampled_site_hold': 52,
    'buckets': dict(sorted(counts.items())),
    'source_sra_sha256': sha(sra),
    'v5_input_sha256': {
        k: v[1] for k, v in expected.items()
    },
}

# Writes ONLY within the new Git diagnostic directory.
out.mkdir(parents=True, exist_ok=True)

with (out / 'proposed_evidence_2055_UNSIGNED.tsv').open(
    'w', newline=''
) as f:
    w = csv.DictWriter(
        f, fieldnames=list(result[0]), delimiter='\t'
    )
    w.writeheader()
    w.writerows(result)

(out / 'summary.json').write_text(
    json.dumps(summary, indent=2, sort_keys=True) + '\n'
)

print('PASS: V5 hashes, 20 SRA records, 2055 experiments')
for k, n in sorted(counts.items()):
    print(k, n, sep='\t')
print('CUMULATIVE_NON_TUMOUR:', len(normal))
print('ADDITIONAL_OLDER:', sorted(normal - v5_normal))
print('ADMIN_REMAINDER_NOT_APPROVED:',
      len(screen) - len(normal))
print('OUTPUT_SHA256:',
      sha(out / 'proposed_evidence_2055_UNSIGNED.tsv'))
print('ZERO CLINICAL ADMISSIONS; NO SOURCE DATA MODIFIED')
