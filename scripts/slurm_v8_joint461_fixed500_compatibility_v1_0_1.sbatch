#!/usr/bin/env bash
#SBATCH --job-name=v8_fixed501
#SBATCH --partition=normal_prio
#SBATCH --account=no_limit_account
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT/private_data/v8_fixed500_v101_slurm_%j.out
#SBATCH --error=/maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT/private_data/v8_fixed500_v101_slurm_%j.err
set -euo pipefail
cd /maps/projects/andersen_group/people/tpn890/ST_ECOTYPES/ARCHS4_V2_6/RECOUNT3_CERTIFICATION_GIT
bash scripts/run_v8_joint461_fixed500_compatibility_v1_0_1.sh
