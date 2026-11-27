#!/usr/bin/env bash
#SBATCH --job-name=v7o_143_fixedW_KL
#SBATCH --partition=normal_prio
#SBATCH --account=no_limit_account
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT/private_data/v7o_slurm_%j.out
#SBATCH --error=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT/private_data/v7o_slurm_%j.err

set -euo pipefail
ROOT=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6
REPO="$ROOT/RECOUNT3_CERTIFICATION_GIT"
cd "$REPO"
bash scripts/run_v7o_143_frozen_ARCHS4_500_vs_496_KL_READONLY.sh
