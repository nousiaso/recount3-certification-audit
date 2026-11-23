#!/usr/bin/env python3
"""Recount3 accession-level SRA evidence extraction, v1.0.0.

Reads existing frozen evidence and source metadata. Creates a NEW, independent
review-only source evidence directory; never updates v1.0.1 approval templates.
No keyword match is an exclusion decision or a patient identity attestation.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)

SAMPLE_FIELDS = (
    'sample_acc', 'sample_name', 'sample_title', 'sample_description',
    'sample_attributes', 'experiment_title', 'experiment_attributes',
    'design_description', 'library_strategy', 'library_source',
    'library_selection', 'library_layout', 'library_construction_protocol',
)
# NOTE: these are high-recall keyword alerts, NOT specimen classifications.
ALERTS = {
    'normal_or_control': r'\bnormal\b|\badjacent\b|\bhealthy\b|\bturbinate\b|\bnon[- ]?tumou?r\b|\bnon[- ]?cancer\b',
    'cell_line_or_culture': r'\bcell[- ]?lines?\b|\bhepg2\b|\bpassages?\b|\bcultured?\b',
    'ascites': r'\bascites\b',
    'metastasis': r'\bmetast\w*\b',
    'treatment_or_timepoint': r'\b(?:pre|post)[-_ ]?(?:treatment|therapy|guadecitabine)\b|\btreated\b|\btherapy\b|\bchemotherap\w*\b|\bguadecitabine\b|\bneoadjuvant\b',
    'multiple_regions': r'\bmultiregion\b|\bregion\s*\d+\b|\btissue\s+\d+\b',
    'isolated_or_nonbulk': r'\bmicrodissect\w*\b|\bflow[- ]?sort\w*\b|\bFACS\b|\blaser[- ]?captur\w*\b|\bsorted cells\b',
    'explicit_primary_tumour': r'\bprimary\s+(?:breast\s+)?(?:tumou?r|carcinoma|cancer)\b|\bprimary\s+\w+\s+(?:tumou?r|carcinoma)\b',
}
PATIENT_KEYS = {
    'patient', 'patient id', 'patient_id', 'patient number', 'patient no',
    'subject', 'subject id', 'subject_id', 'donor id', 'donor_id', 'case id', 'case_id',
}

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for b in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def readtsv(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def outtsv(path, records, fields):
    with path.open('w', newline='', encoding='utf-8') as stream:
        w = csv.DictWriter(stream, delimiter='\t', fieldnames=fields, extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(records)


def parse_attributes(raw):
    attrs = collections.defaultdict(list)
    for item in (raw or '').split('|'):
        if ';;' not in item:
            continue
        name, value = item.split(';;', 1)
        name = re.sub(r'\s+', ' ', name.strip().lower())
        value = value.strip()
        if name:
            attrs[name].append(value)
    return dict(attrs)


def analyze(root, out, expect_libraries=3639, expect_candidates=3602, expect_source_rows=4166):
    root = root.resolve()
    out = out.resolve()
    if out.exists():
        raise RuntimeError(f'Output already exists, refusing overwrite: {out}')
    ledger_path = root / 'stage5_output/04_library_preanalysis_decision_3639.tsv'
    selected_runs_path = root / 'manifest/selected_runs.tsv'
    studies_path = root / 'manifest/studies.tsv'
    libs = readtsv(ledger_path)
    runs = readtsv(selected_runs_path)
    studies = readtsv(studies_path)
    if len(libs) != expect_libraries:
        raise RuntimeError(f'Expected {expect_libraries} libraries, found {len(libs)}')
    keys = [(r['study'], r['library_id']) for r in libs]
    if len(keys) != len(set(keys)):
        raise RuntimeError('Repeated library ID within study in decision ledger')
    ledger = dict(zip(keys, libs))
    candidate = {k for k, v in ledger.items() if v['preanalysis_disposition'] != 'EXCLUDE_CONVENTIONAL_DGE_TAG'}
    if len(candidate) != expect_candidates:
        raise RuntimeError(f'Expected {expect_candidates} candidates, got {len(candidate)}')
    manifest_runs = {(r['study'], r['library_id'], r['run_id']) for r in runs}
    if len(manifest_runs) != len(runs):
        raise RuntimeError('Duplicate selected run in run manifest')
    if any((s, l) not in ledger for s, l, _ in manifest_runs):
        raise RuntimeError('Selected run whose library is missing from ledger')
    study_ids = {r['study'] for r in studies}
    if study_ids != {s for s, _ in keys}:
        raise RuntimeError('Study manifest and selected libraries disagree')

    data = collections.defaultdict(list)
    source_runs = set()
    source_paths = []
    all_errors = []
    for study in sorted(study_ids):
        path = root / 'raw_studies' / study / 'metadata' / f'sra.sra.{study}.MD.gz'
        if not path.is_file():
            raise RuntimeError(f'Source SRA metadata missing: {path}')
        source_paths.append(path)
        with gzip.open(path, 'rt', encoding='utf-8', errors='replace', newline='') as stream:
            rd = csv.DictReader(stream, delimiter='\t')
            for line_num, row in enumerate(rd, 2):
                lib = row.get('experiment_acc', '').strip()
                key = (study, lib)
                if key not in ledger:
                    continue
                rid = row.get('external_id', '').strip()
                rkey = (study, lib, rid)
                if rkey in source_runs:
                    all_errors.append(f'Duplicate run source record: {rkey}')
                source_runs.add(rkey)
                if rkey not in manifest_runs:
                    all_errors.append(f'Unrecognized source run for selected library: {rkey}')
                row['_source_line'] = str(line_num)
                data[key].append(row)

    if len(source_runs) != expect_source_rows:
        all_errors.append(f'Expected {expect_source_rows} unique source runs, found {len(source_runs)}')
    if source_runs != manifest_runs:
        all_errors.append(f'Run accession disagreement; missing={len(manifest_runs-source_runs)} extra={len(source_runs-manifest_runs)}; examples missing={list(sorted(manifest_runs-source_runs))[:5]}')
    if set(data) != set(ledger):
        all_errors.append(f'Library linkage incomplete: missing={list(sorted(set(ledger)-set(data)))[:5]}')
    for k, rr in data.items():
        for f in SAMPLE_FIELDS:
            vals = {str(r.get(f, '') or '').strip() for r in rr}
            if len(vals) > 1:
                all_errors.append(f'Conflicting {f} within {k}: {len(vals)} values')
    if all_errors:
        raise RuntimeError('Source evidence validation FAILED; no output written:\n' + '\n'.join(all_errors[:35]))

    evidence = []
    alerts_by_study = collections.defaultdict(lambda: collections.Counter())
    attr_counts = collections.Counter()
    for key in sorted(candidate):
        study, lib = key
        r = data[key][0]
        l = ledger[key]
        attrs = parse_attributes(r.get('sample_attributes'))
        attr_counts.update(attrs.keys())
        # Only accession/sample-specific metadata; study title NOT used to assign class.
        text = ' | '.join(r.get(f, '') or '' for f in (
            'sample_attributes', 'sample_name', 'sample_title',
            'sample_description', 'experiment_title', 'experiment_attributes'))
        hits = [name for name, regex in ALERTS.items() if re.search(regex, text, re.I)]
        for hit in hits:
            alerts_by_study[study][hit] += 1
        ptokens = {k: attrs[k] for k in attrs if k in PATIENT_KEYS}
        row = {
            'study': study, 'library_id': lib, 'cancer_from_existing_curation': l.get('cancer', ''),
            'manifest_biosample': l.get('biosample', ''),
            'existing_patient_key': l.get('patient_key', ''),
            'existing_patient_key_basis': l.get('patient_key_basis', ''),
            'source_sample_acc': r.get('sample_acc', ''),
            'source_sample_name': r.get('sample_name', ''),
            'source_sample_title': r.get('sample_title', ''),
            'source_sample_description': r.get('sample_description', ''),
            'source_experiment_title': r.get('experiment_title', ''),
            'source_experiment_attributes': r.get('experiment_attributes', ''),
            'source_sample_attributes_raw': r.get('sample_attributes', ''),
            'parsed_sample_attributes_json': json.dumps(attrs, ensure_ascii=False, sort_keys=True),
            'source_patient_token_candidates_json': json.dumps(ptokens, ensure_ascii=False, sort_keys=True),
            'source_library_strategy': r.get('library_strategy', ''),
            'source_library_selection': r.get('library_selection', ''),
            'source_design_description': r.get('design_description', ''),
            'source_library_construction_protocol': r.get('library_construction_protocol', ''),
            'selected_run_count': str(len(data[key])),
            'selected_run_accessions': ';'.join(sorted(x['external_id'] for x in data[key])),
            'metadata_source_path': str((Path('raw_studies') / study / 'metadata' / f'sra.sra.{study}.MD.gz')),
            'metadata_source_sha256': '',  # filled below
            'source_row_record_indices': ';'.join(x['_source_line'] for x in data[key]),
            'keyword_alerts_NOT_DECISIONS': ';'.join(hits),
            'specimen_class_provisional': 'UNRESOLVED_REVIEW_REQUIRED',
            'patient_independence_provisional': 'NOT_CERTIFIED',
            'release_decision': 'NO_APPROVAL_FROM_THIS_TABLE',
        }
        evidence.append(row)

    hashes = {p.parent.parent.name: sha(p) for p in source_paths}
    for r in evidence:
        r['metadata_source_sha256'] = hashes[r['study']]
    assert len(evidence) == expect_candidates
    study_counts = collections.Counter(r['study'] for r in evidence)
    summary = {
        'tool_version': 'recount3_accession_source_evidence_v1_0_0',
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'status': 'EXACT_SOURCE_LINKAGE_VERIFIED_EVIDENCE_ONLY_NOT_CERTIFIED',
        'studies_total': len(study_ids),
        'source_run_rows': len(source_runs),
        'selected_libraries': len(ledger),
        'candidate_libraries': len(candidate),
        'protocol_excluded_libraries': len(ledger)-len(candidate),
        'candidate_metadata_with_sample_attributes': sum(bool(r['source_sample_attributes_raw']) for r in evidence),
        'candidate_with_direct_patient_attribute': sum(r['source_patient_token_candidates_json'] != '{}' for r in evidence),
        'candidate_keyword_alert_counts': dict(collections.Counter(x for r in evidence for x in r['keyword_alerts_NOT_DECISIONS'].split(';') if x)),
        'study_candidate_counts': dict(sorted(study_counts.items())),
        'per_study_alert_counts': {k: dict(v) for k, v in sorted(alerts_by_study.items())},
        'restriction': 'No automatic class, patient, run, source attestation, or model-transfer approvals.',
        'next_use': 'Reviewer/source-backed adjudication and separately controlled entry into v1.0.1 review ledgers.',
    }
    out.mkdir(parents=True, exist_ok=False)
    outtsv(out/'accession_source_evidence_3602.tsv', evidence, list(evidence[0]))
    with (out/'source_metadata_fingerprints.tsv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['relative_source_path','sha256'])
        for p in sorted(source_paths):
            writer.writerow([str(p.relative_to(root)), hashes[p.parent.parent.name]])
        for p in (ledger_path, selected_runs_path, studies_path):
            writer.writerow([str(p.relative_to(root)), sha(p)])
    (out/'source_evidence_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    with (out/'SHA256SUMS.txt').open('w') as f:
        for p in (out/'accession_source_evidence_3602.tsv', out/'source_metadata_fingerprints.tsv', out/'source_evidence_summary.json'):
            f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps({k:summary[k] for k in ('status','studies_total','source_run_rows','selected_libraries','candidate_libraries','candidate_metadata_with_sample_attributes','candidate_with_direct_patient_attribute','candidate_keyword_alert_counts')}, indent=2))
    print('NEW evidence-only directory:', out)
    print('No certification approval tables changed.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--expect-libraries', type=int, default=3639)
    p.add_argument('--expect-candidates', type=int, default=3602)
    p.add_argument('--expect-source-rows', type=int, default=4166)
    a = p.parse_args()
    analyze(a.root, a.out, a.expect_libraries, a.expect_candidates, a.expect_source_rows)

if __name__ == '__main__':
    main()

