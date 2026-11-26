# ARCHS4 frozen identity input schema: read-only preflight

Date: 2026-11-26
Decision: SOURCE_FILES_LOCATED_PARQUET_SCHEMA_UNREADABLE_IN_CURRENT_ENV

The original historical identity-join input paths are present under the
ARCHS4 project root. Stage04 and Stage05 Parquet copies are also located
inside the frozen RELEASES tree. File sizes agree, but byte identity of
these two pairs was not established by this diagnostic.

Frozen atlas: 9,465 rows, 73 fields, known SHA256.
Stage41 candidate table: 3,602 rows, 52 fields; DRP005276 has 103 rows.
Stage41 full-source table: 5,221 rows, 21 fields; DRP005276 has 103 rows.
All 103 DRP005276 candidate records have a populated biosample_IDs field;
zero have recovered GEO GSM aliases or populated direct atlas matches.

CRITICAL: Stage41's 103 DRP005276 records are the source/candidate
universe, not the later 83 technically retained profiles. The original
20 historical technical holds are not reinstated. The atlas table has no
study column, so a zero count from searching `study=DRP005276` in that
file is not a genuine identity-overlap comparison.

The selected Python environment lacks pyarrow, so the actual Parquet
schema, relation-column evidence and 9,465 atlas-to-source mapping were
NOT inspected here. Stage04/05 GEO relation fields represent the same
underlying GEO origin and are not independent corroboration.

Next: hash the two frozen-copy Parquet pairs, inventory available
Parquet readers without installing packages, and reconcile the exact
103 Stage41 members against the 83/20 frozen technical selection.

No original identity join rerun. No model rerun. No clinical admissions.
Atlas/patient independence: UNSIGNED. Pretreatment: UNSIGNED.
