# Raw SRA metadata coverage (NOT clinical certification)

77 study metadata files existed. All 2,055 screen experiments were matched. 1,984 had at
least one field matching broad keyword selection. This is evidence of data availability only.

Potential weakness 1: `has_specific_attributes` is a misleading variable name: `organism`
contains `organ`, `disease`, `sample type=Frozen`, etc.; none proves specimen status.
Potential weakness 2: `previously_unattested_or_ambiguous` tallies any broad attribute among
old tissue-grade categories, even if the attribute has no bearing on lesion identity.
Potential weakness 3: the non-primary regex tests **only keys containing `tissue`**, and only
for `normal|healthy|benign|adjacent|metasta`. V5 already has source-confirmed non-tumours
that these conditions need not detect, e.g. prostate U/uninvolved, generic salivary tissue.
Zero printed flags must **not** be interpreted as absence of nonprimary samples.
Potential weakness 4: `first-line therapy` is not necessarily a sample collection timepoint.

No repeat of the same broad coverage scan is recommended. Any future review should use
actual study and sample-specific evidence and reconcile existing V5 proposed restrictions.

The script is transcribed from this conversation (not matched to a saved cluster script).
The stdout text is transcribed from the user's pasted terminal result, not a cluster file;
this difference in provenance is recorded deliberately.
