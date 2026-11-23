# SRP060016 source-specific tissue identities

The provided SRA source metadata explicitly labels 20 experiments across five patients:
N1–N5 = normal colon; C1–C5 = primary colon cancer; CL1–CL5 = liver metastases;
L1–L5 = normal liver. Exactly five primary-screen records C1–C5 are supported as
primary colorectal tumour *specimens*; do not infer pretreatment, QC pass or final admission.

Five 4-tissue matched sets are present; observations should not be treated as 20 patients.
`stdout.md` is the actual uploaded cluster printout; `script.py` was transcribed from
assistant-provided inline code, not independently hashed against a saved cluster script.

Historical V5 review queue was correct to hold these until exact raw per-experiment attributes
were inspected. The data now resolve only specimen/site identity for C1–C5.
