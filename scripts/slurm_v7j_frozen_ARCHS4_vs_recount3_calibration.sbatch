#!/usr/bin/env bash
#SBATCH --job-name=v7j_143_499_transfer
#SBATCH --partition=normal_prio
#SBATCH --account=no_limit_account
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT/private_data/v7j_slurm_%j.out
#SBATCH --error=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT/private_data/v7j_slurm_%j.err

set -euo pipefail
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
cd "$REPO"
bash scripts/run_v7j_frozen_ARCHS4_vs_recount3_calibration.sh
