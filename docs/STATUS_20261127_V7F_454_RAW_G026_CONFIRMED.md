# V7F completed — the full 454-run calibration source exists locally

Date: 27 November 2026. Source is the user's executed V7F cluster stdout.

The original SRP074349 compressed recount3 G026 gene-sums file is available:

- File: `26_recount3_external_primary_v1/raw_studies/SRP074349/gene_sums/sra.gene_sums.SRP074349.G026.gz`.
- Size: 56,142,028 compressed bytes.
- SHA256: `b056331eb16d6ec612409fd37a266abfe10c4c15329b486dbb9df213612797c5`.
- The first two decompressed lines are non-matrix preamble lines. **The third line (0-based index 2) contains the tab-separated matrix run header.**
- Header width: 669 fields (gene ID plus **668 run columns**).
- All **454 distinct bridge SRRs** appear in the real matrix header.
- 17 of the first inspected postheader lines contain recognizable Ensembl gene IDs.
- Earlier bridge report recorded use of all 454 run columns and 143 ARCHS4-paired specimens for the exploratory **465-gene** comparison (median across-specimen W-gene Pearson r 0.9990371069). That diagnostic did not project the complete frozen 499-gene ecotype model.

**Scientific correction:** the V7F field `TAB_FIRST_LINE_BRIDGE_RUN_OVERLAP 0` is NOT a failed run-membership check; the first line is a preamble, not the header. The meaningful result is `DELIMITER TAB BEST_LINE_0_BASED 2 RUN_VALUES 668 BRIDGE_OVERLAP 454`. Thus the matching original gene-body coverage matrix has **already been downloaded**. There is **no reason to redownload** SRP074349 coverage or to use the unrelated 144-run selected HDF5.

Next: the new V7G workflow discovers the exact 668-run header, extracts all 454 matched run columns, verifies the full 63,856-gene row axis against the frozen Stage5 map, **sums raw G026 base-pair coverage within the 143 SRX experiments**, and extracts the matched ARCHS4 Kallisto pseudocount columns from the original authoritative HDF5.

Both assays are independently normalized with the original frozen Stage09 CPM/TMM/UQ methods over **all valid symbol groups** (TMM/UQ factors estimated on the defined 143-pair calibration subset), then projected using frozen Stage10 pan-cancer anchor means/SDs and unchanged official W, for all 499 genes. This design does not require every matched GSM to belong to the frozen 9,465-sample atlas; when frozen H overlap exists, CPM replay is checked separately.

V7G includes descriptive cross-assay H/fraction/reconstruction readouts, normalization robustness, and a **nine-annotation-difference sensitivity** that sets the corresponding frozen-Z features to zero in both assays while retaining unchanged W. This is an ablation sensitivity, NOT repair or replacement of source annotation.

Samples remain technical calibration controls with **zero clinical admissions**. Do not infer patient independence, treatment history or external biological replication from this work.
