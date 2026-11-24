# V6 — Exact GEO–SRA provenance audit

Date: 2026-11-24
Status: EXACT SOURCE RELATION PASS; CLINICAL LOCK OPEN.

All 64 selected GSM records directly identify their expected SRX:

- SRP094795 / GSE91065: 23/23
- SRP072494 / GSE79671: 20/20
- SRP194626 / GSE130469: 4/4
- SRP112751 / GSE101588: 17/17

Zero source-relation mismatches. Zero selected cell lines flagged by the implemented check.

Full GSE79671 subject-token multiplicities: 16 doubletons, four singletons.

Distinct breast title prefixes do not certify distinct patients.
Pre-bevacizumab is not equivalent to treatment-naive.

Original GEO snapshot SHA256 values:

- GSE91065: `4dd7d20aabbaafb4c2d0b16d541da000787fa880d652e2db0b1aa1b52a7ac730`
- GSE79671: `25d6a138fa91e90b75bbefc375b06cdb34e01b25608e98a38fa6b89b6266eedc`
- GSE130469: `44317d378636086bfa304756520251776709efa1f1a36eacbd9d74d3effc9c58`
- GSE101588: `2974a161614d170a59d5d640192ce808ffb414c8f8ca6d2f25c577866cc1c407`

Original audit outputs:

- `diagnostics/20261124_v6_priority_64_source_linkage/64_source_linkage_LOCAL_ONLY.tsv`: `eb9ccdab319e9531f35b9c5ca2133feb9952fb8f56979b86ae426aea6a244e6b`
- `diagnostics/20261124_v6_priority_64_source_linkage/summary.json`: `d8f7bd376f4926b271114f244b5fb85de6adc84468f9135ec007d9865a8aa14c`
- `diagnostics/20261124_v6_priority_64_geo_exact/64_geo_exact_source_evidence_LOCAL_ONLY.tsv`: `71f788ee1fd34edb3adb8b7fe26dc0cd06067a3a87b31e377b9d5124748fdae5`
- `diagnostics/20261124_v6_priority_64_geo_exact/summary.json`: `46939760ddddae2e55ebd20996a14026dc506f76f2c313d3818000a03170b1ee`
- `diagnostics/20261124_v6_priority_64_geo_exact/geo_soft_LOCAL_ONLY/GSE91065_family.soft.gz`: `4dd7d20aabbaafb4c2d0b16d541da000787fa880d652e2db0b1aa1b52a7ac730`
- `diagnostics/20261124_v6_priority_64_geo_exact/geo_soft_LOCAL_ONLY/GSE79671_family.soft.gz`: `25d6a138fa91e90b75bbefc375b06cdb34e01b25608e98a38fa6b89b6266eedc`
- `diagnostics/20261124_v6_priority_64_geo_exact/geo_soft_LOCAL_ONLY/GSE130469_family.soft.gz`: `44317d378636086bfa304756520251776709efa1f1a36eacbd9d74d3effc9c58`
- `diagnostics/20261124_v6_priority_64_geo_exact/geo_soft_LOCAL_ONLY/GSE101588_family.soft.gz`: `2974a161614d170a59d5d640192ce808ffb414c8f8ca6d2f25c577866cc1c407`

Scripts:

- `tools/v6_priority_64_source_linkage_READONLY_v1.py`: `e76d01e52d109066923a0645a0dad7e4e44c3ff5205545284cca957769658975`
- `tools/v6_priority_64_geo_exact_READONLY_v1.py`: `23893925c9528a3c852cb219e0cca2ab720a35e9719a8174da5bda781726d6d5`

All 64 records remain NOT_ADJUDICATED and NO_UNSIGNED.
No frozen source files or model outputs changed.
