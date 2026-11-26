#!/usr/bin/env python3
"""Read-only targeted DRP005276 exact-accession join against frozen ARCHS4."""

import csv
import gzip
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path
import pyarrow.parquet as pq

root = Path(sys.argv[1])
repo = root / 'RECOUNT3_CERTIFICATION_GIT'

files = {
 'atlas': root/'14_projection_use_analysis_ready_atlas_v1_0_3/03_projection_analysis_ready_core.tsv.gz',
 's04': root/'04_candidate_evidence_v1_0_4/04_candidate_sample_full_metadata.parquet',
 's05': root/'05_study_sample_adjudication_v1_0_2_1/18_complete_adjudication_manifest.parquet',
 'cand': root/'41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/03_candidate_integrated_evidence_3602.tsv',
 'full': root/'41_recount3_evidence_to_adjudication_unsigned_v1_0_1_20261120/01_full_source_experiment_identity_5221.tsv',
 'selected': repo/'diagnostics/20261123_v6_primary_unattested_review_v1_1/source_evidence_599_LOCAL_ONLY.tsv',
}

expected = {
 'atlas': '7f6852e829f7aaa11a26f86b44247c9c8c0840958494b7766db1ad6abd180092',
 's04': '2e74aabd85a62c3f80835691a2ffa0c11f374bf262142349f579fd4e6f4fc36a',
 's05': '736e89b2969693f504d298c8e7274a1bcf5773a0f18667b114953c575ffd4a32',
 'cand': 'aea32ae0a63eafa8b0f6abb82c970d0fe7e2b202d39c50247f99103a9f3c2b81',
 'full': '95f450a91531e5943dd08c0a4edfda137c959b292d757a155d84e333bdc083fa',
 'selected': '8b6f017dd575ec6fc3917b784740d37c8556fe567474e6a5e69dcad066ef35eb',
}

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

for name, path in files.items():
    assert path.is_file() and digest(path) == expected[name], (
        'FROZEN_INPUT_MISMATCH_' + name
    )

print('DRP005276_ATLAS_EXACT_JOIN_V1')
print('FROZEN_INPUT_HASHES_VERIFIED', len(files))

EXP = re.compile(
    r'(?<![A-Za-z0-9])(?:SRX|ERX|DRX)\d+(?![A-Za-z0-9])', re.I
)
BIO = re.compile(
    r'(?<![A-Za-z0-9])(?:SAMN|SAMEA|SAMD)\d+(?![A-Za-z0-9])', re.I
)

def tokens(pattern, value):
    return set(
        x.upper() for x in pattern.findall(str(value or ''))
    )

def read_tsv(path, zipped=False):
    opener = gzip.open if zipped else open
    with opener(str(path), 'rt', encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))

atlas_rows = read_tsv(files['atlas'], True)
atlas = {x['sample_id'].strip().upper() for x in atlas_rows}
assert len(atlas_rows) == len(atlas) == 9465

cand = read_tsv(files['cand'])
full = read_tsv(files['full'])
allsel = read_tsv(files['selected'])
assert (len(cand), len(full), len(allsel)) == (3602, 5221, 599)

cr = [r for r in cand if r['study'] == 'DRP005276']
fr = [r for r in full if r['study'] == 'DRP005276']
sr = [r for r in allsel if r['study'] == 'DRP005276']

ci = {r['library_id']: r for r in cr}
fi = {r['experiment_acc']: r for r in fr}
selected = {r['experiment_acc'] for r in sr}

assert len(cr) == len(fr) == len(ci) == len(fi) == 103
assert len(sr) == len(selected) == 83
assert set(ci) == set(fi) and selected <= set(ci)

held = set(ci) - selected
assert len(held) == 20
assert all(ci[e]['technical_hold'] == 'NO' for e in selected)
assert all(ci[e]['technical_hold'] == 'YES' for e in held)

print('DRP005276_MEMBERSHIP SOURCE 103 SELECTED 83 HELD 20')

def graph_from_parquet(path, key):
    pf = pq.ParquetFile(str(path))
    assert key in pf.schema_arrow.names
    assert 'relation' in pf.schema_arrow.names

    seen = {}
    n = duplicate = conflict = 0

    for batch in pf.iter_batches(
        batch_size=8192, columns=[key, 'relation']
    ):
        data = batch.to_pydict()

        for acc, rel in zip(data[key], data['relation']):
            if isinstance(acc, bytes):
                acc = acc.decode('utf-8', 'replace')

            gsm = str(acc or '').strip().upper()

            if gsm not in atlas:
                continue

            n += 1
            pair = (
                frozenset(tokens(EXP, rel)),
                frozenset(tokens(BIO, rel))
            )

            if gsm in seen:
                duplicate += 1
                if pair != seen[gsm]:
                    conflict += 1
            else:
                seen[gsm] = pair

    return seen, n, duplicate, conflict

graphs = []

for label, field in [
    ('s04', 'geo_accession'),
    ('s05', 'sample_id')
]:
    g, n, d, c = graph_from_parquet(files[label], field)

    print(
        'ATLAS_SOURCE', label,
        'MATCHING_ROWS', n,
        'DISTINCT_GSMS', len(g),
        'DUPLICATE_GSM_ROWS', d,
        'CONFLICTS_WITHIN_SOURCE', c,
        'WITH_EXPERIMENT', sum(bool(x[0]) for x in g.values()),
        'WITH_BIOSAMPLE', sum(bool(x[1]) for x in g.values()),
        'WITH_BOTH', sum(bool(x[0] and x[1]) for x in g.values())
    )
    graphs.append((g, n, d, c))

g0, n0, d0, c0 = graphs[0]
g1, n1, d1, c1 = graphs[1]

disagreements = sum(
    g0[gsm] != g1[gsm]
    for gsm in atlas if gsm in g0 and gsm in g1
)

print('PARSED_STAGE04_STAGE05_DISAGREEMENTS', disagreements)

complete = (
    len(g0) == len(g1) == 9465
    and n0 == n1 == 9465
    and d0 == d1 == c0 == c1 == disagreements == 0
    and all(x[0] and x[1] for x in g0.values())
)

print('FULL_EXPERIMENT_AND_BIOSAMPLE_GRAPH_COVERAGE', complete)

if complete:
    exp_index = defaultdict(set)
    bio_index = defaultdict(set)

    for gsm, (exps, bios) in g0.items():
        for x in exps:
            exp_index[x].add(gsm)
        for x in bios:
            bio_index[x].add(gsm)

    print('ATLAS_UNIQUE_EXPERIMENT_IDS', len(exp_index))
    print('ATLAS_UNIQUE_BIOSAMPLE_IDS', len(bio_index))

    for label, ids in [
        ('SELECTED83', selected),
        ('HELD20', held),
        ('ALL103', set(fi))
    ]:
        by_exp = by_bio_full = by_bio_cand = any_hit = 0
        bios_discord = full_bio_empty = cand_bio_empty = 0

        for exp in ids:
            bfull = tokens(
                BIO, fi[exp].get('biosamples_all_sources')
            )
            bcand = tokens(
                BIO, ci[exp].get('biosample_IDs')
            )

            bios_discord += (bfull != bcand)
            full_bio_empty += (not bfull)
            cand_bio_empty += (not bcand)

            e_hit = bool(exp_index.get(exp))
            f_hit = any(bio_index.get(x) for x in bfull)
            c_hit = any(bio_index.get(x) for x in bcand)

            by_exp += e_hit
            by_bio_full += f_hit
            by_bio_cand += c_hit
            any_hit += (e_hit or f_hit or c_hit)

        print(
            'TARGET_GROUP', label,
            'N', len(ids),
            'EXPERIMENT_HITS', by_exp,
            'FULL_SOURCE_BIOSAMPLE_HITS', by_bio_full,
            'CANDIDATE_BIOSAMPLE_HITS', by_bio_cand,
            'UNION_EXACT_HITS', any_hit,
            'BIOSAMPLE_FIELD_DISCORDANCES', bios_discord,
            'FULL_SOURCE_BIOSAMPLE_EMPTY', full_bio_empty,
            'CANDIDATE_BIOSAMPLE_EMPTY', cand_bio_empty
        )

    print('EXACT_ACCESSION_SCREEN COMPLETED')
else:
    print('EXACT_ACCESSION_SCREEN NOT_RUN_INCOMPLETE_ATLAS_GRAPH')

print('PATIENT_INDEPENDENCE_CERTIFIED NO')
print('PRETREATMENT_CERTIFIED NO')
print('CLINICAL_ADMISSIONS 0')
print('NO_INDIVIDUAL_ACCESSIONS_PRINTED')
print('NO_SOURCE_DATA_MODIFIED')
print('DRP005276_ATLAS_EXACT_JOIN_V1_COMPLETE')
