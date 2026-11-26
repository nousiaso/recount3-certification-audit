#!/usr/bin/env python3
"""Read-only DRP005276 treatment and cross-study evidence triage."""

import csv
import gzip
import hashlib
import re
import sys
from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
repo = root / 'RECOUNT3_CERTIFICATION_GIT'

files = {
    'atlas': root / (
        '14_projection_use_analysis_ready_atlas_v1_0_3/'
        '03_projection_analysis_ready_core.tsv.gz'
    ),
    'candidate': root / (
        '41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/'
        '03_candidate_integrated_evidence_3602.tsv'
    ),
    'full': root / (
        '41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/'
        '01_full_source_experiment_identity_5221.tsv'
    ),
    'selected': repo / (
        'diagnostics/20261123_v6_primary_unattested_review_v1_1/'
        'source_evidence_599_LOCAL_ONLY.tsv'
    )
}

hashes = {
    'atlas': '7f6852e829f7aaa11a26f86b44247c9c8c0840958494b7766db1ad6abd180092',
    'candidate': 'aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81',
    'full': '95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa',
    'selected': '8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb'
}

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def table(path):
    opener = gzip.open if path.suffix == '.gz' else open

    with opener(
        str(path), 'rt', encoding='utf-8-sig', newline=''
    ) as f:
        reader = csv.DictReader(f, delimiter='\t')
        return reader.fieldnames or [], list(reader)

def safe_status(value):
    value = (value or '').strip()

    if not value:
        return '[EMPTY]'

    if re.fullmatch(r'[A-Z0-9_/-]{1,100}', value) and not re.search(
        r'(GSM|SRX|DRX|ERX|SAMN|SAMD|SAMEA)\d+', value
    ):
        return value

    return '[OTHER_SHA256_' + hashlib.sha256(
        value.encode('utf-8')
    ).hexdigest()[:12] + ']'

print('DRP005276_CLINICAL_EVIDENCE_TRIAGE_V1')

for key, path in files.items():
    assert path.is_file()
    assert sha256(path) == hashes[key], 'HASH_MISMATCH_' + key

print('FROZEN_INPUT_HASHES_VERIFIED', len(files))

fields = {}
rows = {}

for key, path in files.items():
    fields[key], rows[key] = table(path)

assert [
    len(rows[k])
    for k in ('atlas', 'candidate', 'full', 'selected')
] == [9465, 3602, 5221, 599]

can = {
    r['library_id']: r
    for r in rows['candidate']
    if r['study'] == 'DRP005276'
}

full = {
    r['experiment_acc']: r
    for r in rows['full']
    if r['study'] == 'DRP005276'
}

selected = {
    r['experiment_acc']
    for r in rows['selected']
    if r['study'] == 'DRP005276'
}

assert len(can) == len(full) == 103
assert len(selected) == 83
assert set(can) == set(full)
assert selected <= set(can)

held = set(can) - selected
assert len(held) == 20

print(
    'MEMBERSHIP_VALIDATED',
    len(can), len(selected), len(held)
)

# ---------------------------------------------------------
# A. Which treatment-related metadata fields actually exist?
# ---------------------------------------------------------

print('\nTREATMENT_RELATED_COLUMN_NAMES')

for source in ('candidate', 'full', 'selected'):
    matching = [
        name for name in fields[source]
        if re.search(
            r'pre.?treat|therap|chem|expos|collect|'
            r'time.?point|recurr',
            name, re.I
        )
    ]

    print(
        'SOURCE', source,
        'COLUMN_NAMES',
        ';'.join(matching) if matching else '[NONE]'
    )

# ---------------------------------------------------------
# B. Existing clinical and certification statuses
# ---------------------------------------------------------

print('\nCATEGORICAL_EVIDENCE')

for source, lookup in (
    ('candidate', can),
    ('full', full)
):
    if source == 'candidate':
        names = (
            'pretreatment_source_text_status',
            'pretreatment_certified',
            'patient_independence_certified',
            'primary_tumour_certified',
            'atlas_independence_certified',
            'cohort_approved',
            'tissue_review_status',
            'study_detected_atlas_overlap'
        )
    else:
        names = (
            'study_exposure_detected',
            'identity_certified'
        )

    for name in names:
        if name not in fields[source]:
            print('FIELD_ABSENT', source, name)
            continue

        for group, ids in (
            ('RETAINED83', selected),
            ('HELD20', held)
        ):
            counts = Counter(
                safe_status(lookup[x].get(name))
                for x in ids
            )

            print(
                'FIELD', source, name, group,
                'COUNTS',
                repr(dict(sorted(counts.items())))
            )

# ---------------------------------------------------------
# C. Evidence availability, without printing patient data
# ---------------------------------------------------------

print('\nEVIDENCE_POINTER_COVERAGE')

for source, lookup, pointers in (
    (
        'candidate', can,
        (
            'pretreatment_evidence_pointers',
            'tissue_evidence_pointers',
            'patient_group_evidence_pointers',
            'source_record_pointer',
            'GEO_source_series',
            'GEO_GSM_aliases'
        )
    ),
    (
        'full', full,
        (
            'source_patient_evidence_pointers',
            'biosample_identity_evidence',
            'study_exposure_detected',
            'GEO_exact_relation_GSM_aliases'
        )
    )
):
    for field in pointers:
        if field not in fields[source]:
            print('POINTER_FIELD_ABSENT', source, field)
            continue

        for group, ids in (
            ('RETAINED83', selected),
            ('HELD20', held)
        ):
            n = sum(
                bool((lookup[x].get(field) or '').strip())
                for x in ids
            )

            print(
                'POINTER', source, field, group,
                'POPULATED', n, 'OF', len(ids)
            )

# ---------------------------------------------------------
# D. Study-level GEO association only (not patient matching)
# ---------------------------------------------------------

print('\nSTUDY_LEVEL_GEO_CHECK')

gse = re.compile(
    r'(?<![A-Za-z0-9])GSE\d+(?![A-Za-z0-9])',
    re.I
)

def gses(value):
    return {
        x.upper()
        for x in gse.findall(value or '')
    }

source_gse = set()

for exp in can:
    source_gse.update(
        gses(can[exp].get('GEO_source_series'))
    )

atlas_gse = set()
bladder_gse = set()
bladder_rows = 0

for row in rows['atlas']:
    matched = (
        gses(row.get('primary_gse_id')) |
        gses(row.get('all_gse_ids'))
    )

    atlas_gse.update(matched)

    if re.search(
        r'bladd|urothel',
        row.get('analysis_cancer_class') or '',
        re.I
    ):
        bladder_rows += 1
        bladder_gse.update(matched)

print('DRP005276_GEO_SERIES_DISTINCT', len(source_gse))
print('ATLAS_GEO_SERIES_DISTINCT', len(atlas_gse))
print('ATLAS_BLADDER_LIKE_ROWS', bladder_rows)
print(
    'ATLAS_BLADDER_LIKE_GEO_SERIES_DISTINCT',
    len(bladder_gse)
)
print(
    'DRP_ATLAS_GEO_SERIES_INTERSECTION',
    len(source_gse & atlas_gse)
)

print('STUDY_LEVEL_ALIAS_CHECK_NOT_PATIENT_INDEPENDENCE')
print('PATIENT_INDEPENDENCE_CERTIFIED NO')
print('PRETREATMENT_CERTIFIED NO')
print('CLINICAL_ADMISSIONS 0')
print('NO_INDIVIDUAL_IDENTIFIERS_PRINTED')
print('NO_FROZEN_DATA_MODIFIED')
print('DRP005276_CLINICAL_EVIDENCE_TRIAGE_V1_COMPLETE')
