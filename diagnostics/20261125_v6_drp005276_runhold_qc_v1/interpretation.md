# DRP005276 — historical QC exclusion reasons (2026-11-25)

Decision: ORIGINAL_QC_HOLDS_20_OF_20_RECONCILED_CLINICAL_UNSIGNED.

## Source and run-level checks

- Frozen V5 original study inventory: 103 distinct experiments.
- V5 administrative screen: 83, all readiness historical_run_hold=NO.
- Outside screen: 20, all QC_HOLD_NOT_ADMITTED and historical_run_hold=YES.
- Original 409-row run-hold queue: 20 exact experiment matches, precisely one row for each excluded experiment, zero mismatched study labels, zero selected-experiment matches.
- Original QC reason-code experiment counts: TOP10_COVERAGE_GT_0_50 = 19; DETECTION_LT_0_10 = 1; WITHIN_STUDY_OUTLIER = 2. Some reasons overlap; 22 flags belong to 20 excluded experiments.
- Among held experiments, original gene detection fraction: minimum 0.0698133, maximum 0.54377; 1 below 0.10.
- Original top-ten coverage fraction: minimum 0.302583, maximum 0.915882; 19 above 0.50.
- Original average mapped length: minimum 70.1 bp, maximum 70.97 bp; zero below 40 bp.

The diagnostic field named MIN_MEDIAN_MAX computes the element at index floor(n/2), i.e. the 11th order statistic for 20 observations; this is not the conventional median (average of the 10th and 11th). We do not interpret its printed middle value as a conventional median. This does not affect the thresholds, reason codes or exclusion counts.

## Interpretation and open gates

Original historical technical hold membership and run-level reasons are fully reconciled. No experiment is reinstated. QC evidence supports exclusion; it does not certify that the selected 83 are independent patients, pathologically primary bladder tumours, pretreatment, atlas-independent or eligible for fixed-W SpatialEcoTyper projection.

The study remains CLINICAL_UNSIGNED; clinical admissions = 0. Frozen ARCHS4, V5, V6, W and downstream biological inference unchanged.

## Provenance

Executed original read-only diagnostic script and its captured stdout are archived adjacent to this interpretation as script.py and stdout.txt; exact byte hashes are in SHA256SUMS.txt. Input hashes validated against the frozen V5 inventory (5241), screen (2055), readiness (3602) and original run-hold queue (409). Only aggregate data are included in the public repository; no individual accession or patient information was printed.

Next: original DRP005276 study publication, sample-to-patient relationship and treatment/primary pathology evidence. Do not repeat the original QC audit on unchanged inputs.
