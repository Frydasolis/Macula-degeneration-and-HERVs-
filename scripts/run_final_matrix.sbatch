#!/bin/bash

#SBATCH --job-name=GeneTE
#SBATCH --partition=prod
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/storage/lemus_g/roldan/ARMD/log/GeneTE_%x_%j.out
#SBATCH --error=/storage/lemus_g/roldan/ARMD/log/GeneTE_%x_%j.err

set -eo pipefail

BASE="/storage/lemus_g/roldan/ARMD"

echo "=================================================="
echo "Final Gene + TE matrix"
echo "=================================================="
echo "Sample: $SAMPLE"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"
echo "=================================================="

source /software/anaconda/etc/profile.d/conda.sh
conda activate stellarscope_protocol_env

python \
"$BASE/scripts/build_final_gene_TE_matrix.py" \
--sample "$SAMPLE"

echo "Finished: $(date)"
