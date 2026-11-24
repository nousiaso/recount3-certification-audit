# V6 NSCLC82 — TCIA clinical provenance

Date: 2026-11-24
Study: SRP117020
GEO series: GSE103584

## Completed source checks

- Exact selected SRX-to-GEO-to-TCIA clinical Case ID links: 82/82.
- Distinct matched TCIA subjects: 82.
- GEO/TCIA histology concordance: 82/82.
- Adenocarcinoma: 64.
- Squamous-cell carcinoma: 18.
- Pathological M0: 77.
- Pathological M1a: 1.
- Pathological M1b: 4.
- Adjuvant treatment Yes: 26.
- Chemotherapy Yes: 26.
- Radiation Yes: 10.
- No explicit preoperative treatment fields identified by the audit.

## Source provenance

Original TCIA CSV SHA256:
0a1f7cdf02cef58450db17f116a3fb0ecfd4d737ca6843bb1939dbe421071476

Private 82-row clinical crosswalk SHA256:
dbfc06c14511ee031c8cccaec8d020eeb4370fd6ad8cd1b4f7d8b89a79be3c79

Saved summary JSON SHA256:
a48c19d266f514cfb58bd4127664d62c0e2e2fef0d5f2cfd1034d172d2594e18

## Output-capture incident

The original diagnostic completed and printed its final status.
Its summary JSON and private clinical crosswalk were written.

The shell's tee process failed because its destination directory
did not exist when tee attempted to open stdout.txt.

Consequently, the original stdout file was not captured.
The saved summary and source hashes are the archived evidence.
No replacement stdout is presented as original execution output.

## Interpretation

The exact original clinical patient crosswalk is supported.

Pathological M1 does not independently establish that the sequenced
specimen came from a metastatic lesion. All five M1 cases require
sampled-lesion/site review.

Postoperative adjuvant therapy must not be interpreted as
preoperative treatment.

The published treatment-naive statement is cohort-level evidence.
Individual specimen timing remains a separate certification gate.

Clinical approval: NO_UNSIGNED.
Frozen V5: unchanged.
SpatialEcoTyper: paused.

Status:
NSCLC82_TCIA_LINK_AND_HISTOLOGY_CONCORDANT_NOT_CLINICALLY_CERTIFIED
