#!/bin/bash

#SBATCH --job-name=SRP413248_FQ
#SBATCH --output=/storage/lemus_g/roldan/ARMD/logs/SRP413248/fastqdump_%A_%a.out
#SBATCH --error=/storage/lemus_g/roldan/ARMD/logs/SRP413248/fastqdump_%A_%a.err

#SBATCH --cpus-per-task=8

BASE=/storage/lemus_g/roldan/ARMD

SRR_LIST=$BASE/resourses/SRP413248_SRR.txt
SRA_DIR=$BASE/resourses/SRP413248/sra
FASTQ_DIR=$BASE/results/SRP413248/fastq

source ~/.bashrc
conda activate sra211

mkdir -p "$FASTQ_DIR"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SRR_LIST")

SRA_FILE="$SRA_DIR/$SRR/$SRR.sra"

echo "=========================================="
echo "Study:       SRP413248"
echo "SRR:         $SRR"
echo "SRA:         $SRA_FILE"
echo "Job ID:      $SLURM_JOB_ID"
echo "Array ID:    $SLURM_ARRAY_TASK_ID"
echo "Node:        $(hostname)"
echo "Start:       $(date)"
echo "=========================================="

if [ ! -f "$SRA_FILE" ]; then
    echo "ERROR: SRA file not found"
    exit 1
fi

fastq-dump \
    "$SRA_FILE" \
    --split-spot \
    --split-files \
    --outdir "$FASTQ_DIR"

STATUS=$?

echo "=========================================="
echo "SRR:         $SRR"
echo "Exit status: $STATUS"
echo "End:         $(date)"
echo "=========================================="

exit $STATUS
