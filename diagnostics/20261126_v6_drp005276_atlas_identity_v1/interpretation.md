# DRP005276 — frozen ARCHS4 identifier-screen audit

Date: 2026-11-26
Decision: DIRECT_IDENTIFIER_SCREEN_NO_MATCHES_INDEPENDENCE_UNSIGNED

## Audited universes

- Original DRP005276 experiments: 103; 83 selected, 20 historical QC-held.
- Frozen V5 cross-study inventory: 5,241 libraries.
- Frozen analysis-ready ARCHS4 atlas: 9,465 GEO GSM sample identifiers.
- Original SRA, selected-source and V5 inventory inputs verified by SHA256.
- Frozen atlas file SHA256:
  `7f6852e829f7aaa11a26f86b44247c9c8c0840958494b7766db1ad6abd180092`.

## Results and limits

- **Zero** GSM aliases were recovered for any of 103 DRP005276 experiments
  from the inspected original SRA/selected-source fields.
- Consequently, zero exact GSM matches are **not informative evidence**
  of independence from the frozen atlas; GSM-link coverage is zero.
- No original DRP005276 experiment accession repeats under another
  study in the 5,241-row V5 inventory.
- No direct original sample-name collisions were found using the **one**
  eligible named sample field available in the other-study V5 records.
- No source SRA sample-accession string matched values in the 6 broadly
  selected accession-like V5 fields. This screen mixes source accession
  columns with status/conflict fields; it is not a complete, schema-
  validated BioSample or sample-accession crosswalk, nor an atlas link.
- Eight V5 overlap-related fields were inventoried, but several have
  categorical values other than literal YES/NO, and some were empty.
  Their meanings were not decoded in this diagnostic.
- The screen did not directly link DRP005276 to the original ARCHS4
  SRA experiment/BioSample source graph or all independent depositions.

## Scientific disposition

No demonstrated overlap in the identifier channels examined. This does
**not** establish cross-study patient independence or atlas independence.
No biological admissions, no technical reinstatements, and no changes
to frozen V5/V6, ARCHS4, or SpatialEcoTyper assets.

Clinical admissions: ZERO.
Pretreatment: UNSIGNED.
Atlas/cross-study independence: UNSIGNED.

Next: inspect the exact frozen archive identity-crosswalk resources and
V5 overlap-field value semantics; then carry out a targeted, read-only
source accession graph check using documented paths and schemas.

The executed code, captured stdout and SHA256 manifest are archived
adjacent to this interpretation. No sample identifiers in public stdout.
