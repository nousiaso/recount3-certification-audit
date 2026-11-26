#!/usr/bin/env python3
"""Inspect original crosswalk source and search existing logs; never run join."""
import ast
import hashlib
import os
import re
import sys
import zipfile
from pathlib import Path

root = Path(sys.argv[1])
archive = root / '26_recount3_external_primary_v1/ARCHS4_RECOUNT3_EXACT_IDENTITY_JOIN_READONLY_20261120.zip'
expected_zip = '20c22c10e07ec8d9bc5847b42c65eae3786614eabdf675afd46ad21b48502c05'
expected_code = 'ff713734d101b34269eaabbf4397740ab1972fc355e64090d7526f72476c1765'

def hash_file(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

assert archive.is_file()
assert hash_file(archive) == expected_zip

with zipfile.ZipFile(str(archive)) as z:
    code_bytes = z.read('archs4_recount3_exact_identity_join_READONLY.py')

assert hashlib.sha256(code_bytes).hexdigest() == expected_code
code = code_bytes.decode('utf-8-sig')
ast.parse(code)

print('ARCHS4_IDENTITY_LOG_INPUT_PROBE_V1')
print('VERIFIED_HISTORICAL_ZIP_SHA256', expected_zip)
print('VERIFIED_SCRIPT_SHA256', expected_code)
print('ORIGINAL_JOIN_EXECUTED NO')

print('\nSOURCE_CODE_PATH_DEFINITIONS')

# Only the first 25 lines of the original program are displayed.
# They describe the historical input paths and definitions.
for i, line in enumerate(code.splitlines()[:25], 1):
    line = re.sub(
        r'\b(?:DRX|DRR|DRS|SRX|SRR|SRS|ERX|ERR|ERS|GSM|SAMN|SAMD|SAMEA)\d+\b',
        '[ACCESSION]', line, flags=re.I
    )
    print('SOURCE_LINE', i, line.rstrip()[:240])

print('\nHISTORICAL_LOG_SEARCH')

roots = [
    (root / '26_recount3_external_primary_v1', 5),
    (root / '27_recount3_final_cert_v1_0_0', 4),
    (root / '28_recount3_final_cert_v1_0_1', 4),
    (root / '29_recount3_adjudication_20261119', 4),
    (root / '31_frozen_execution_contract_audit_v1_1_20261119', 4),
    (root, 2),
]

prune = {
    '.git', 'raw_studies', 'RECOUNT3_CERTIFICATION_GIT',
    'RELEASES', 'matrices', 'counts', 'fastq',
    'gene_expression', '__pycache__'
}

signatures = (
    b'READ_ONLY_EXACT_IDENTITY_JOIN_COMPLETE',
    b'Atlas GSMs:',
    b'Agreed source GSMs:',
    b'Full source experiments:',
)

seen = set()
scanned = 0
skipped_large = 0
candidates = 0
cap = 12000

for base, maxdepth in roots:
    if not base.is_dir():
        print('MISSING_SEARCH_ROOT', str(base.relative_to(root)))
        continue

    local_scanned = 0

    for parent, dirs, files in os.walk(str(base)):
        depth = len(Path(parent).relative_to(base).parts)

        dirs[:] = sorted(
            x for x in dirs
            if x not in prune and depth < maxdepth
        )

        for name in sorted(files):
            if not name.lower().endswith(
                ('.txt', '.log', '.out', '.md', '.stdout')
            ):
                continue

            p = Path(parent) / name
            if p in seen:
                continue

            seen.add(p)
            scanned += 1
            local_scanned += 1

            try:
                size = p.stat().st_size

                if size > 8000000:
                    skipped_large += 1
                    continue

                with p.open('rb') as f:
                    data = f.read(8000001)

                signals = sum(s in data for s in signatures)

                if signatures[0] in data:
                    candidates += 1
                    print(
                        'POSSIBLE_LOG',
                        str(p.relative_to(root)),
                        'BYTES', size,
                        'SHA256', hashlib.sha256(data).hexdigest(),
                        'MARKER_SIGNATURES', signals
                    )

            except (OSError, ValueError):
                print('UNREADABLE_LOG_FILE')

            if scanned >= cap:
                break

        if scanned >= cap:
            break

    print(
        'SEARCH_ROOT_SUMMARY',
        str(base.relative_to(root)) if base != root else 'ROOT',
        'TEXT_FILES_VISITED', local_scanned
    )

    if scanned >= cap:
        break

print('TEXT_FILES_VISITED', scanned)
print('LARGE_FILES_SKIPPED', skipped_large)
print('MARKER_CANDIDATES', candidates)
print('SEARCH_CAPPED', scanned >= cap)

print('ORIGINAL_RESULTS_AUTHENTICATED NO')
print('DRP005276_CROSS_ATLAS_INDEPENDENCE NOT_CERTIFIED')
print('DRP005276_PRETREATMENT NOT_CERTIFIED')
print('CLINICAL_ADMISSIONS 0')
print('NO_INDIVIDUAL_SAMPLE_VALUES_PRINTED')
print('NO_FROZEN_FILES_MODIFIED')
print('ARCHS4_IDENTITY_LOG_INPUT_PROBE_V1_COMPLETE')
