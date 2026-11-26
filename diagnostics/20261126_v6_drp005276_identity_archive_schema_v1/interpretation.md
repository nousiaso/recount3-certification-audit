# DRP005276 — archive format and provenance inspection

Date: 2026-11-26

Decision: IDENTITY_ZIP_SOURCE_PACKAGE_ONLY_RESULT_FILES_NOT_PRESENT

## Verified historical package

ARCHS4_RECOUNT3_EXACT_IDENTITY_JOIN_READONLY_20261120.zip

SHA256:
20c22c10e07ec8d9bc5847b42c65eae3786614eabdf675afd46ad21b48502c05

Its three members pass ZIP CRC and uncompressed SHA256 inspection:

1. Original read-only Python crosswalk script (8,301 bytes).
2. README for the script (2,644 bytes).
3. SHA256SUMS manifest (226 bytes).

**No generated identity-crosswalk or atlas-source records are
packaged in this ZIP.**

The interpretation that the ZIP itself holds 9,465 resolved
atlas identities is corrected.

The historically reported 9,465 identity links and overlap
counts must be traced to original execution outputs, not
inferred from this code archive.

## Earlier final-cert review files

- 27_recount3_final_cert_v1_0_0/draft/ARCHS4_overlap_status.tsv
- 28_recount3_final_cert_v1_0_1/draft/ARCHS4_overlap_status.tsv

These are byte-identical, SHA256:
15c3562b3b83415c9f7e875d5a2b7e9b3396e825931e6d27ebb9b761a95e0132

Both contain 3,602 rows and six status fields.

They lack a study column; EXPLICIT_DRP005276_STUDY_ROWS 0
does not mean that DRP005276 is absent.

They are not the complete atlas-side identity graph.

## Open gates and next steps

Inspect the packaged README and Python source as inert text
to discover historical output locations, source inputs and
expected checksums.

Recover original execution outputs before a targeted
DRP005276-to-atlas join.

Do not rerun the archived atlas-wide program.

Clinical admissions: ZERO.
Pretreatment status: UNSIGNED.
Cross-study/atlas independence: UNSIGNED.
Original 20 QC-held experiments remain excluded.
No frozen data or model files modified.
