#!/bin/bash

#SBATCH --job-name=Seurat17
#SBATCH --partition=prod
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00

#SBATCH --output=/storage/lemus_g/roldan/ARMD/log/Seurat17_%j.out
#SBATCH --error=/storage/lemus_g/roldan/ARMD/log/Seurat17_%j.err

set -eo pipefail

source /software/anaconda/etc/profile.d/conda.sh
conda activate seurat_env

echo "=============================================="
echo "Creating combined Seurat object"
echo "Job: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start: $(date)"
echo "=============================================="

which R
R --version

Rscript \
/storage/lemus_g/roldan/ARMD/scripts/create_Seurat_SRP413248.R

echo "=============================================="
echo "Finished: $(date)"
echo "=============================================="
