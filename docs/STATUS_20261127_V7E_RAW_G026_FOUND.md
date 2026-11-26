# V7E completed — raw SRP074349 G026 matrix is the next asset to check

Date: 27 November 2026. Evidence is the user's executed V7E stdout.

- All **129** Stage2 selected G026 HDF5 files were opened successfully (zero HDF5 read errors).
- **0/454** SRP074349 bridge-calibration runs occur in the Stage2 run QC inventory (4,166 rows).
- **0/454** occur in any of the 129 selected Stage2 HDF5 run axes.
- Original SRA SRP074349 metadata exists and recognizes all **143** source bridge experiments (668 metadata rows).
- Exact original gene-sums source file exists: `26_recount3_external_primary_v1/raw_studies/SRP074349/gene_sums/sra.gene_sums.SRP074349.G026.gz`, **56,142,028 compressed bytes**.
- The raw gene-sums G026 gzip **has not yet been checked** for the 454 runs, and no expression has been extracted for them in this work.
- Stage5 earlier 465-gene `bridge_SRP074349` diagnostic artifacts exist, but are not equivalent to complete 499-gene model transfer.

**Decisive next step:** execute `scripts/run_v7f_raw_srp074349_gene_sums_header.sh`. This reads the first 20 lines of the original compressed G026 file, reports only aggregate accession/header overlap and layout diagnostics, and checks its SHA256. It does not modify, download or generate any expression data.

If all 454 run accessions are present in the original file, design a versioned **run-preserving, then experiment-aggregated** targeted extraction using verified file schema. Do not re-download the already available original source or substitute the unrelated 144-run selected HDF5. If not present, investigate source format or genuine missing data before considering any downloads.

Preserve unsuccessful V7C provenance. Keep clinical releases at zero; these are same-sample technical calibration controls, not independent validation specimens.
