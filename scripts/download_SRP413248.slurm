#!/bin/bash

#SBATCH --job-name=SRP413248
#SBATCH --output=/storage/lemus_g/roldan/ARMD/logs/SRP413248/SRP413248_%A_%a.out
#SBATCH --error=/storage/lemus_g/roldan/ARMD/logs/SRP413248/SRP413248_%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

BASE=/storage/lemus_g/roldan/ARMD
SRR_LIST=$BASE/resourses/SRP413248_SRR.txt
SRA_DIR=$BASE/resourses/SRP413248/sra

source ~/.bashrc
conda activate sra211

mkdir -p "$SRA_DIR"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")

echo "=========================================="
echo "Study:       SRP413248"
echo "SRR:         $SRR"
echo "Job ID:      $SLURM_JOB_ID"
echo "Array ID:    $SLURM_ARRAY_TASK_ID"
echo "Node:        $(hostname)"
echo "Start:       $(date)"
echo "=========================================="

prefetch "$SRR" \
    --output-directory "$SRA_DIR" \
    --max-size 100G

STATUS=$?

echo "=========================================="
echo "SRR:         $SRR"
echo "Exit status: $STATUS"
echo "End:         $(date)"
echo "=========================================="

exit $STATUS
