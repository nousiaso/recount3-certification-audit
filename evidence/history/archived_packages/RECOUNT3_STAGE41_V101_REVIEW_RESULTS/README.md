# Stage 41 v1.0.1 — FUSCC identity evidence carry-forward, unsigned adjudication integration

**Purpose:** consolidate prior Stage 37/39/40 evidence and original SRA study records
into one independently auditable human review worksheet for all 3,602 recount3
candidate experiments. The 9,465-sample frozen ARCHS4 model and all prior stage
files remain untouched. All 3,602 clinical and patient decisions remain EMPTY.

## Identity semantics
- SRA SRX/ERX/DRX = experiment, SRR/ERR/DRR = sequencing run; SRS/ERS/DRS = SRA sample.
- BioSample SAMN/SAMEA/SAMD = deposited biological material, not automatically person.
- GEO GSM = submitted GEO sample; an experiment can have multiple GSM aliases.
- Explicit SRA and GEO patient IDs are **proposals**, namespaced by study;
  shared IDs indicate possible multiple samples per person, not independent patients.
- RR-#### isolates: **suspected related specimens only**; not attested people.
- The 334 FUSCC source crosswalk profiles (279 corroborated GEO patient tokens,
  55 RNA isolate-only records) are carried forward as UNSIGNED evidence.
  'LOCKED' in an evidence filename does not certify patient or treatment eligibility.
  Source WES cross-assay mismatch warnings remain reviewer flags.
- An accession-shaped SRA run_alias is only an alias TOKEN, NOT a verified
  official run accession; it remains separate from the run-accession evidence gate.
- Shared SRS or BioSample across experiments is evidence of related deposited
  material, not proof of two independent patient observations.
- Nonmatching GEO aliases are NOT proof of complete GEO or patient independence.
- All study catalog completeness flags remain NOT CERTIFIED.

## Reviewer-facing evidence
01: Full 5,221 source experiments, observed BioSample/run, all GSM aliases and patient proposals.
02: Every retained Stage40 edge plus proposed-patient, BioSample, run and RR family edges.
03: One row per 3,602 candidates: evidence, conservative clinical flags and partition.
04: Independent **UNSIGNED** reviewer worksheet; every reviewer field blank.
05: Original GEO clinical evidence for candidates; not transformed into patient/treatment approval.
06: Explicit source-level patient-group proposals and RR suspect groups (not certified).
07: Per-candidate prioritized review exceptions, including overlap, repeated patients, tissue.
08: All 126 study summaries of identities, pretreatment leads and coverage.
09: Full source experiment-level missing run/BioSample/GSM and identity conflicts.
10/11: Summary and SHA256 ledger of all consumed original files; SHA256SUMS checks output.

## Clinical adjudication
To approve a sample for core validation *later*, a reviewer must establish all of:
(a) patient group and independence from other biological units;
(b) unequivocal primary invasive tumour and acceptable pretreatment status;
(c) technical QC, proper assay chemistry, and independent atlas identity;
(d) study-held-out (strict core), considering all overlapping experiments,
    not just the selected recount3 libraries.
No automatic decision follows from a keyword, study title, unique BioSample,
absence of a GSM alias, or a zero match in an incomplete atlas-side crosswalk.
A source 'metastasis: Yes' field may mean an outcome, not biopsied metastasis.
An oligo-dT cDNA priming step after ribosomal depletion is not poly(A) RNA enrichment.

## Graph completeness
An input containing 5,221 experiments does NOT imply 5,221 verified runs/BioSamples.
Stage40 provided BioSample mapping for only its 3,602 candidates and no run edges.
Provide separately verified TSVs with study,experiment_acc,run_acc and/or
study,experiment_acc,biosample if available. Missing relationships remain explicitly
unknown rather than fabricated. Atlas-side exact accession crosswalk may be supplied
with sample_id,experiment_acc,biosample; partial input can detect more matches but
cannot certify negative patient overlap.

## Governance
The script NEVER signs decisions, writes to an existing directory, overwrites
prior source/protocol files, or creates a certified projection manifest.
Passing hashes is integrity, not biological certification. Reviewer decisions are
subject to human attestation before a separately authorized cohort lock.
