# Evidence and decision policy

**Observed**: directly present in the supplied cluster output or source file.
**Verified checksum**: file byte identity established, not biological truth.
**Source-supported**: explicit experiment/specimen attribute or publication evidence.
**Provisional**: source data suggest an interpretation, but one or more release gates remain.
**Unresolved**: insufficient or contradictory evidence; never fill with assumptions.
**Final approved**: requires explicit clinical + technical + patient + atlas + model-transfer signoff;
**none has been granted**.

The full clinical universe is 5,241 records / 176 studies. The Stage41 administrative
primary screen has 2,055 rows, of which V5 flags 21 as source-supported non-tumour.
The resulting **2,034** are not a newly certified cohort. The raw physical screen file
retains 2,055 rows; filters are proposed decisions only.

Cross-study independence is not established by accession non-overlap. Study-scoped patient
IDs must not be treated as globally unique individuals.

The source code for earlier historical audits must be captured from original cluster files
if possible. Titles and references from the README are a historic inventory, not code bytes.

Data governance: keep the hosted repository private. These example logs include public
SRA sample attributes (some contain age/sex); verify institutional rules before pushing any
sample-level output to a third party. If required, retain full logs only on the cluster and
commit only redacted summaries/hashes. Never commit private clinical records, tokens or passwords.
