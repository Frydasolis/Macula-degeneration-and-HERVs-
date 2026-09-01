#!/bin/bash
#SBATCH --job-name=STAR_index_ARMD
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/lemus_g/roldan/ARMD/logs/SRP413248/STAR_index_%j.out
#SBATCH --error=/storage/lemus_g/roldan/ARMD/logs/SRP413248/STAR_index_%j.err

set -euo pipefail

echo "=========================================="
echo "STAR genome index - SRP413248"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "=========================================="

STAR="/home/roldan/.conda/envs/star_env/bin/STAR"

GENOME_DIR="/storage/lemus_g/roldan/ARMD/resourses/SRP413248/reference/STAR_index"
FASTA="/storage/lemus_g/roldan/ARMD/resourses/SRP413248/reference/GRCh38.primary_assembly.genome.fa"
GTF="/storage/lemus_g/roldan/ARMD/resourses/SRP413248/reference/gencode.v38.annotation.gtf"

mkdir -p "$GENOME_DIR"

"$STAR" \
    --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$FASTA" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 97 \
    --genomeSAindexNbases 14

echo "=========================================="
echo "INDEX COMPLETED"
echo "End: $(date)"
echo "=========================================="
