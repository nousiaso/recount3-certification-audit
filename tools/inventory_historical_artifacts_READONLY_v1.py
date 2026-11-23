#!/usr/bin/env python3
"""Read-only inventory: known historical ARCHS4/recount3 scripts and reports.

Only lists paths, sizes and SHA-256 for files <=8 MiB. No scientific execution,
no reading patient rows, no writes, and no opening archive contents. A missing
file here means not found in the inspected roots, not that it never existed.
"""
import argparse
import hashlib
import os
from collections import defaultdict
from pathlib import Path

TARGETS = {
 'FROZEN_EXECUTION_CONTRACT': {'script': ['audit_frozen_outputs_READONLY_v1_1.py'],
   'result': ['frozen_output_contract_v1_1.json'],
   'package': ['ARCHS4_EXECUTION_CONTRACT_AUDIT_V1_1_20261119.zip']},
 'FROZEN_ORIGINAL_SOURCE': {'script': [
   '06_audit_ARCHS4_primary_expression_compatibility_v1.py',
   '09_prepare_ARCHS4_locked_expression_and_normalization_v1_0_1.py',
   '10_run_ARCHS4_pan_cancer_fixed_anchor_projection_v1_0_4.R',
   '11_build_ARCHS4_layer1_adequacy_robustness_atlas_v1_0_2_1.R',
   '12_run_ARCHS4_cancer_specific_fixed_anchor_projection_v1_0_3.R',
   '13_build_ARCHS4_projection_use_analysis_ready_atlas_v1_0_3.R']},
 'ATLAS_IDENTITY_OVERLAP': {'script': ['probe_W34_and_ARCHS4_atlas_READONLY.py'],
   'result': ['01_W34_metadata_candidates.tsv', '02_GSM_direct_overlap_candidates.tsv', '03_readonly_probe_summary.json'],
   'package': ['RECOUNT3_W34_ATLAS_CANDIDATES.tar.gz', 'RECOUNT3_W34_OVERLAP_EVIDENCE_REVIEW_20261119.zip']},
 'RECOUNT3_SOURCE_REVIEW': {'script': ['recount3_accession_source_evidence_v1_0_0.py'],
   'package': ['RECOUNT3_SOURCE_EVIDENCE_3602.tar.gz', 'RECOUNT3_CERTIFICATION_ADJUDICATION_20261119.zip',
     'RECOUNT3_CERTIFICATION_V1_0_1_20261118.zip']},
 'STAGE41_SCREEN': {'package': ['RECOUNT3_STAGE41_V101_REVIEW_RESULTS.tar.gz']},
 'FUSCC_AND_TRANSFER_PREFLIGHT': {'package': ['FUSCC_COHORT_RECONCILIATION_UNSIGNED_20261120.zip',
   'FUSCC_MODEL_COMPATIBILITY_EVIDENCE.tar.gz', 'FUSCC_FROZEN_W_PROJECTION_METHODS.tar.gz',
   'FUSCC_ORIGINAL_PROJECTION_SOURCE.tar.gz', 'FUSCC_FIXED_W_SOURCE_REVIEW_AND_PREFLIGHT_20261120.zip']},
 'W499_MEMBERSHIP': {'result': ['16A_frozen_W_gene_symbol_harmonization.tsv'],
   'package': ['ARCHS4_FROZEN_GENE_PROVENANCE_SOURCE.tar.gz']},
 'TRANSFER_INPUT_PREFLIGHT': {'result': ['03_bridge_165_availability.tsv', '00_gene_row_mapping.tsv',
   '04_library_preanalysis_decision_3639.tsv']},
 'V5_SOURCE_CROSSWALK': {'script': ['build_source_crosswalk_v5.py'],
   'result': ['08_V5_QA.json', '09_OUTPUT_SHA256SUMS.txt'],
   'package': ['RECOUNT3_CLINICAL_SPECIMEN_PATIENT_CROSSWALK_V5_20261122.zip']},
}
SKIP_DIRS = {'.git', '__pycache__', 'raw_studies', 'stage2_output', '00_download',
 'RECOUNT3_CERTIFICATION_GIT', '.venv', 'venv', 'node_modules'}

def filehash(path, size):
    if size > 8 * 1024 * 1024:
        return 'NOT_HASHED_GT_8MiB'
    h = hashlib.sha256()
    try:
        with path.open('rb') as f:
            for part in iter(lambda: f.read(1024 * 1024), b''):
                h.update(part)
        return h.hexdigest()
    except OSError as e:
        return 'READ_ERROR_' + type(e).__name__

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('--secondary-root', type=Path)
    a = p.parse_args()
    roots = [a.root]
    if a.secondary_root and a.secondary_root != a.root:
        roots.append(a.secondary_root)
    names = {n for roles in TARGETS.values() for pats in roles.values() for n in pats}
    found = defaultdict(set)
    scopes = []
    scanned = 0
    for root in roots:
        if not root.is_dir():
            print('MISSING_ROOT\t' + str(root)); continue
        scoped = [root / '26_recount3_external_primary_v1',
                  root / '31_frozen_execution_contract_audit_v1_1_20261119',
                  root / 'RELEASES' / 'ARCHS4_SPATIALECOTYPER_FROZEN_v1_0_0' / 'PROJECT_TREE']
        scopes.append(str(root) + ' [TOP LEVEL ONLY]')
        for path in root.iterdir():
            if path.is_file():
                scanned += 1
                if path.name in names and not path.is_symlink():
                    found[path.name].add(path)
        for scope in scoped:
            if not scope.is_dir():
                continue
            scopes.append(str(scope))
            for folder, dirs, files in os.walk(scope, followlinks=False):
                dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not d.startswith('.')]
                for filename in files:
                    scanned += 1
                    if filename in names:
                        candidate = Path(folder) / filename
                        if not candidate.is_symlink():
                            found[filename].add(candidate)
    print('=== HISTORICAL AUDIT FILE INVENTORY V1 ===')
    for scope in scopes:
        print('SEARCH_SCOPE\t' + scope)
    print('SCANNED_FILENAMES\t' + str(scanned))
    print('NOTE\tFOUND does not mean scientifically verified or executed')
    print('NOTE\tNOT_FOUND only describes selected search locations')
    print('audit_id\trole\tstatus\tfilename\tbytes\tsha256_or_note\tpath')
    nfound = nall = 0
    for audit, roles in TARGETS.items():
        for role, patterns in roles.items():
            for filename in patterns:
                nall += 1
                matches = sorted(found[filename])
                if not matches:
                    print(f'{audit}\t{role}\tNOT_FOUND\t{filename}\t-\t-\t-')
                    continue
                nfound += 1
                for path in matches[:6]:
                    try:
                        size = path.stat().st_size
                        checksum = filehash(path, size)
                    except OSError as e:
                        size, checksum = '-', 'STAT_ERROR_' + type(e).__name__
                    print(f'{audit}\t{role}\tFOUND\t{filename}\t{size}\t{checksum}\t{path}')
                if len(matches) > 6:
                    print(f'{audit}\t{role}\tMORE_MATCHES\t{filename}\t{len(matches)-6}\t-\t-')
    print(f'SUMMARY\tFOUND_PATTERNS={nfound}/{nall}')
    print('READ_ONLY_COMPLETE')

if __name__ == '__main__':
    main()

