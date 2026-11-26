#!/usr/bin/env python3
"""Reader availability, frozen Parquet copy hashes, and Stage41 103/83/20.

No Parquet row loading. No patient/sample values printed.
Python 3.6 compatible; standard library only.
"""

import csv
import hashlib
import importlib.util
import sys

from collections import Counter
from pathlib import Path

root = Path(sys.argv[1])
repo = root / 'RECOUNT3_CERTIFICATION_GIT'

frozen = root / (
    'RELEASES/ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0/PROJECT_TREE'
)

sources = (
    '04_candidate_evidence_v1_0_4/04_candidate_sample_full_metadata.parquet',
    '05_study_sample_adjudication_v1_0_2_1/18_complete_adjudication_manifest.parquet',
)

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for buf in iter(lambda: f.read(1048576), b''):
            h.update(buf)
    return h.hexdigest()

def load(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        d = csv.DictReader(f, delimiter='\t')
        return d.fieldnames, list(d)

print('ARCHS4_PARQUET_READER_MEMBERSHIP_PREFLIGHT_V1')

# ------------------------------------------------------------
# 1. Available libraries; do not install or change environments
# ------------------------------------------------------------

print('READERS_IN_CURRENT_ENV')

for mod in ('pyarrow', 'duckdb', 'fastparquet', 'pandas'):
    print(
        'READER', mod,
        'AVAILABLE', importlib.util.find_spec(mod) is not None
    )

# Look for existing pyarrow installations in other conda environments.
# Do not import or execute them.
base = Path('/home/tpn890/orestis/conda_envs')
hits = []

if base.is_dir():
    for env in sorted(base.iterdir()):
        if not env.is_dir():
            continue
        if not (env / 'bin/python').exists():
            continue

        for path in env.glob(
            'lib/python*/site-packages/pyarrow/__init__.py'
        ):
            hits.append(env.name)
            break

print('EXISTING_CONDA_ENVIRONMENTS_WITH_PYARROW', len(hits))

for env in hits[:12]:
    print('PYARROW_ENV_NAME', env)

print('ENV_NAMES_NOT_PRINTED', max(0, len(hits) - 12))

# ------------------------------------------------------------
# 2. Stage04/05 copies: SHA256 identity, not size alone
# ------------------------------------------------------------

print('\nFROZEN_PARQUET_COPY_IDENTITY')

for name in sources:
    a, b = root / name, frozen / name

    assert a.is_file() and b.is_file(), 'PARQUET_COPY_MISSING'

    ah, bh = digest(a), digest(b)

    print(
        'SOURCE', Path(name).name,
        'BYTES_A', a.stat().st_size,
        'BYTES_B', b.stat().st_size,
        'SHA256_A', ah,
        'SHA256_B', bh,
        'BYTE_IDENTICAL', ah == bh
    )

# ------------------------------------------------------------
# 3. Reconcile the 103 Stage41 records with later frozen QC
# ------------------------------------------------------------

selected_file = repo / (
    'diagnostics/20261123_v6_primary_unattested_review_v1_1/'
    'source_evidence_599_LOCAL_ONLY.tsv'
)

assert digest(selected_file) == (
    '8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb'
)

candidates_file = root / (
    '41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/'
    '03_candidate_integrated_evidence_3602.tsv'
)

full_file = root / (
    '41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/'
    '01_full_source_experiment_identity_5221.tsv'
)

assert digest(candidates_file) == (
    'aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81'
)

assert digest(full_file) == (
    '95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa'
)

_, source = load(selected_file)
_, candidates = load(candidates_file)
_, full = load(full_file)

assert len(source) == 599
assert len(candidates) == 3602
assert len(full) == 5221

selected = {
    r['experiment_acc']
    for r in source
    if r['study'] == 'DRP005276'
}

can = {
    r['library_id']: r
    for r in candidates
    if r['study'] == 'DRP005276'
}

all_source = {
    r['experiment_acc']: r
    for r in full
    if r['study'] == 'DRP005276'
}

assert len(selected) == 83
assert len(can) == len(all_source) == 103
assert selected <= set(can)
assert set(can) == set(all_source)

held = set(can) - selected
assert len(held) == 20

print('\nDRP005276_STAGE41_RECONCILIATION')
print('STAGE41_CANDIDATE_RECORDS', len(can))
print('STAGE41_FULL_SOURCE_RECORDS', len(all_source))
print('LATER_FROZEN_TECHNICALLY_SELECTED', len(selected))
print('LATER_HISTORICAL_QC_HELD', len(held))
print('STAGE41_TO_599_SELECTED_EXACT_SET_JOIN_PASS')

# Categorical distributions only; no individual records printed.
fields = (
    'technical_hold',
    'stage41_partition_NOT_ADMISSION',
    'atlas_independence_certified',
    'cohort_approved'
)

for field in fields:
    assert field in candidates[0], 'MISSING_CANDIDATE_FIELD_' + field

    for group, ids in (
        ('SELECTED83', selected),
        ('HELD20', held)
    ):
        raw = Counter(
            (can[x].get(field) or '').strip()
            for x in ids
        )

        safe = {}

        for value, n in raw.items():
            allowed = (
                len(value) <= 70
                and all(
                    z.isupper() or z.isdigit() or z in '_-'
                    for z in value
                )
            )

            key = value if allowed else (
                'OTHER_VALUE_SHA256_' +
                hashlib.sha256(value.encode('utf-8')).hexdigest()[:12]
            )

            safe[key if key else '[EMPTY]'] = n

        print(
            'STATUS_FIELD', field,
            'GROUP', group,
            'COUNTS', repr(dict(sorted(safe.items())))
        )

print('NO_PARQUET_ROWS_READ')
print('NO_PRIVATE_ACCESSION_VALUES_PRINTED')
print('ORIGINAL_JOIN_NOT_EXECUTED')
print('PATIENT_INDEPENDENCE_UNSIGNED')
print('PRETREATMENT_UNSIGNED')
print('CLINICAL_ADMISSIONS 0')

print('ARCHS4_PARQUET_READER_MEMBERSHIP_PREFLIGHT_V1_COMPLETE')
