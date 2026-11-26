# DRP005276 — identity ZIP schema checkpoint

Date: 2026-11-26

Status: SOURCE_PACKAGE_ONLY_COMPLETE_CROSSWALK_OUTPUTS_NOT_INSIDE_ZIP

- ZIP SHA256 verified; three members: Python script,
  README and SHA256 manifest.
- ZIP contains no generated 9,465-row atlas crosswalk tables.
- Older v1.0.0/v1.0.1 overlap TSVs are identical 3,602-row
  status tables.
- Complete historical identity graph and DRP005276-specific
  join: not yet recovered.
- V5 no-overlap flags remain limited technical evidence.
- Pretreatment and global patient independence: UNSIGNED.
- Clinical admissions: ZERO.
- Frozen data and model assets unchanged.

Evidence:
diagnostics/20261126_v6_drp005276_identity_archive_schema_v1/

Next: inspect original ZIP README/script and locate execution outputs.
