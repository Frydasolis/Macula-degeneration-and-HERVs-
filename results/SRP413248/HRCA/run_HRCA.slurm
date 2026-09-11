#!/bin/bash
#SBATCH --job-name=HRCA_annotation
#SBATCH --partition=prod
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=HRCA_annotation_%j.out
#SBATCH --error=HRCA_annotation_%j.err

source ~/.bashrc
conda activate HRCA_scvi120

cd /storage/lemus_g/roldan/ARMD/results/SRP413248/HRCA

python annotate_HRCA.py
