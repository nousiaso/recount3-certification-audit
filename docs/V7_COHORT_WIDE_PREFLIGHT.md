# V7: Cohort-wide clinical gate preflight (READ ONLY)

Purpose: combine frozen V5 (2,055 primary-screen experiments), the previously archived V6 evidence-precedence classification, and Stage41 technical/identity evidence for all 77 contributing studies.

**This is an evidence-baseline preflight, not an automatic clinical certification or model run.** A study-wide claim or a unique BioSample/patient proxy does not establish specimen-level primary-tumour status, treatment-naive status, patient independence, or external validation eligibility.

Run on supekgate after pulling the repository:

```bash
cd /maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT
git pull --ff-only origin main
bash scripts/run_v7_integrated_preflight.sh || echo "V7_FAILED_SSH_SESSION_PRESERVED"
```

The launcher runs in a child Bash session (do not source it) and validates the script SHA256 before execution. Input source files are never modified. Local-only outputs go under `private_data/v7_integrated_candidate_gates_v1/` (excluded from Git), and safe aggregate stdout/checksums go under `diagnostics/20261127_v7_integrated_candidate_gates_v1/`.

Expected baseline categories are 22 source-supported non-tumour exclusions, 5 exact primary specimens with other gates open, 52 sampled-site-uncertain sarcoma specimens, 280 source-context records, and 1,696 historical-evidence records. These categories total 2,055 **administrative** candidates; zero samples are yet released for definitive validation.

Next: review the aggregate 77-study output, prioritize specimen-specific evidence and technical assay compatibility, then explicitly lock a reproducible first eligible cohort. No SpatialEcoTyper inferential analysis should be claimed from this preflight.
