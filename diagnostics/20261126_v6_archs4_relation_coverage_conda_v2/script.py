import csv
import gzip
import sys
from pathlib import Path

import pyarrow
import pyarrow.parquet as pq

root = Path(sys.argv[1])

atlas = root / (
    '14_projection_use_analysis_ready_atlas_v1_0_3/'
    '03_projection_analysis_ready_core.tsv.gz'
)

sources = (
    (
        'STAGE04',
        root / '04_candidate_evidence_v1_0_4/'
        '04_candidate_sample_full_metadata.parquet',
        'geo_accession'
    ),
    (
        'STAGE05',
        root / '05_study_sample_adjudication_v1_0_2_1/'
        '18_complete_adjudication_manifest.parquet',
        'sample_id'
    ),
)

print('ARCHS4_RELATION_COVERAGE_CONDA_V2')
print('PYTHON', sys.version.split()[0],
      'PYARROW', pyarrow.__version__)

with gzip.open(
    str(atlas), 'rt',
    encoding='utf-8-sig', newline=''
) as fh:
    rows = list(csv.DictReader(fh, delimiter='\t'))

atlas_ids = {
    row['sample_id'].strip().upper()
    for row in rows
}

assert len(rows) == len(atlas_ids) == 9465
print('FROZEN_GSMS', len(atlas_ids))

for label, path, accession_col in sources:

    parquet = pq.ParquetFile(str(path))
    schema = parquet.schema_arrow

    assert accession_col in schema.names, (
        label, 'MISSING_ACCESSION_COLUMN'
    )
    assert 'relation' in schema.names, (
        label, 'MISSING_RELATION_COLUMN'
    )

    found = set()
    related = set()
    matched_rows = 0

    for batch in parquet.iter_batches(
        batch_size=8192,
        columns=[accession_col, 'relation']
    ):
        data = batch.to_pydict()

        for acc, relation in zip(
            data[accession_col], data['relation']
        ):
            if isinstance(acc, bytes):
                acc = acc.decode('utf-8', 'replace')

            gsm = str(acc or '').strip().upper()

            if gsm not in atlas_ids:
                continue

            matched_rows += 1
            found.add(gsm)

            if relation is not None and (
                bool(relation.strip())
                if isinstance(relation, str)
                else bool(relation)
            ):
                related.add(gsm)

    print(
        'SOURCE', label,
        'ROWS', parquet.metadata.num_rows,
        'ROW_GROUPS', parquet.metadata.num_row_groups
    )

    print(
        'SOURCE', label,
        'ACCESSION_TYPE',
        str(schema.field(accession_col).type),
        'RELATION_TYPE',
        str(schema.field('relation').type)
    )

    print(
        'SOURCE', label,
        'ATLAS_MATCHING_ROWS', matched_rows,
        'DISTINCT_GSMS', len(found),
        'GSMS_WITH_RELATION', len(related)
    )

    print(
        'SOURCE', label,
        'MISSING_GSMS', len(atlas_ids - found),
        'GSMS_WITHOUT_RELATION', len(atlas_ids - related)
    )

print('ORIGINAL_IDENTITY_JOIN_EXECUTED NO')
print('PATIENT_INDEPENDENCE_NOT_CERTIFIED')
print('PRETREATMENT_NOT_CERTIFIED')
print('CLINICAL_ADMISSIONS 0')
print('NO_INDIVIDUAL_ACCESSIONS_PRINTED')
print('NO_FROZEN_FILES_MODIFIED')
print('ARCHS4_RELATION_COVERAGE_CONDA_V2_COMPLETE')
