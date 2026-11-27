# V7K source inspection PASS; next inspect lower-level NMF solver

Date: 27 November 2026. Evidence: user-supplied full original supekgate V7K stdout.

## Installed source facts

- Installed SpatialEcoTyper version 1.0.3; R 4.4.3.
- Actual namespace DeconvoluteSE is the publicly exposed function.
- The function filters input rows to W model genes and may apply an **automatic log2(dat+1)** transformation only when all input elements are nonnegative and max(dat) exceeds 80. The actual condition on the frozen 143 and 9465 columns has not yet been checked numerically.
- Original Stage10 explicitly uses scale=FALSE, W fixed, nsample.per.run=10000, sum2one=FALSE, ncores=1.
- NMFpredict also has scaling disabled and partitions the matrix only when ncol(testdat)>10000. Both 143 and 9465 fall below that threshold; **the explicit sample-chunk branch is not the cause of the observed 143-versus-9465 H difference.**
- For both cases NMFpredict directly calls its lower-level '.nmf.predict(W,testdat,scale=FALSE,normalize=FALSE)' function. V7K did not inspect the implementation of that function.

## Scientific interpretation

V7J previously passed original raw G026 preparation (454/454 SRR runs mapped to 143 SRX experiments, 63856 gene rows; 0 coverage aggregation mass loss), original frozen ARCHS4 CPM-Z replay (max abs 2.04e-14), gene harmonization 499, and archived H/fraction consistency (8.88e-16). It **correctly stopped** when a new 143-column ARCHS4 projection differed from the immutable 9465-sample H: max abs H 0.0126240960966; max fraction TVD 0.00959058827769.

The installed wrapper code rules out explicit splitting at nsample.per.run=10000 but does not establish the underlying cause. **Do not infer a specific algorithmic mechanism before examining '.nmf.predict'.** The automatic log branch also needs a bounded numerical check using actual inputs.

No recount3 H/fraction scores, cross-assay agreement metrics or clinical admissions have been produced by V7J.

## V7L: inspect true namespace-local solver

Versioned files:
- tools/v7l_inspect_installed_nmf_predict_READONLY_v1.R
- scripts/run_v7l_installed_nmf_predict_READONLY.sh

Run after pulling main:
bash scripts/run_v7l_installed_nmf_predict_READONLY.sh

Output:
diagnostics/20261127_v7l_installed_nmf_predict_source_v1/01_installed_nmf_predict_stdout.txt

The V7L script prints the actual lower-level solver implementation, its function arguments and namespace helper details. It does not read expression, run a model or change source files. Preserve stopped V7J logs and private data. Only after reviewing V7L should we design the smallest numerical experiment that distinguishes a data-dependent solver from differences introduced by input filtering or preprocessing.

All 143 SRP074349 specimens remain same-source technical calibration controls. Clinical admissions: 0.
