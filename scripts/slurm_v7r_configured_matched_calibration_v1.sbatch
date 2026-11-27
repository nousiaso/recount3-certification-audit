#!/usr/bin/env bash
#SBATCH --job-name=v7r_matched500
#SBATCH --partition=normal_prio
#SBATCH --account=no_limit_account
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT/private_data/v7r_slurm_%j.out
#SBATCH --error=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT/private_data/v7r_slurm_%j.err
set -euo pipefail
cd /maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT
bash -n scripts/run_v7r_configured_matched_calibration_v1.sh
bash scripts/run_v7r_configured_matched_calibration_v1.sh "${1:?Study required}"
