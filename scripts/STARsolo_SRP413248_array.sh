#!/bin/bash
#SBATCH --job-name=STARsolo_ARMD
#SBATCH --array=1-53%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/lemus_g/roldan/ARMD/logs/SRP413248/STARsolo_%A_%a.out
#SBATCH --error=/storage/lemus_g/roldan/ARMD/logs/SRP413248/STARsolo_%A_%a.err

set -euo pipefail

BASE="/storage/lemus_g/roldan/ARMD"

STAR="/home/roldan/.conda/envs/star_env/bin/STAR"

GENOME_DIR="$BASE/resourses/SRP413248/reference/STAR_index"
TABLE="$BASE/resourses/SRP413248_STARsolo_samples.tsv"

WHITELIST_V2="$BASE/resourses/SRP413248/whitelist/737K-august-2016.txt"
WHITELIST_V3="$BASE/resourses/SRP413248/whitelist/3M-february-2018.txt"

FASTQ_DIR="$BASE/results/SRP413248/fastq"
OUT_BASE="$BASE/results/SRP413248/STARsolo"

TASK_ID="${SLURM_ARRAY_TASK_ID}"

echo "=============================================="
echo "STARsolo SRP413248"
echo "Task: $TASK_ID"
echo "Job: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo "=============================================="

# ------------------------------------------------
# Read sample information from table
# ------------------------------------------------

LINE=$(sed -n "$((TASK_ID + 1))p" "$TABLE")

SRR=$(echo "$LINE" | awk '{print $1}')
GSM=$(echo "$LINE" | awk '{print $2}')
CHEMISTRY=$(echo "$LINE" | awk '{print $3}')
CBLEN=$(echo "$LINE" | awk '{print $4}')
UMILEN=$(echo "$LINE" | awk '{print $5}')

CB_UMI="$FASTQ_DIR/${SRR}_2.fastq"
CDNA="$FASTQ_DIR/${SRR}_3.fastq"

OUTDIR="$OUT_BASE/$SRR"

mkdir -p "$OUTDIR"

# ------------------------------------------------
# Select whitelist
# ------------------------------------------------

if [[ "$CHEMISTRY" == "V2" ]]; then
    WHITELIST="$WHITELIST_V2"
elif [[ "$CHEMISTRY" == "V3" ]]; then
    WHITELIST="$WHITELIST_V3"
else
    echo "ERROR: Unknown chemistry: $CHEMISTRY"
    exit 1
fi

echo ""
echo "Sample information:"
echo "SRR:       $SRR"
echo "GSM:       $GSM"
echo "Chemistry: $CHEMISTRY"
echo "CB length: $CBLEN"
echo "UMI length:$UMILEN"

echo ""
echo "Input files:"
echo "CB/UMI: $CB_UMI"
echo "cDNA:   $CDNA"

echo ""
echo "Whitelist:"
echo "$WHITELIST"

# ------------------------------------------------
# Checks
# ------------------------------------------------

echo ""
echo "Checking STAR..."
"$STAR" --version

echo ""
echo "Checking genome index..."
test -s "$GENOME_DIR/Genome"
test -s "$GENOME_DIR/SA"
test -s "$GENOME_DIR/SAindex"
echo "Genome index: OK"

echo ""
echo "Checking whitelist..."
test -s "$WHITELIST"
echo "Whitelist: OK"

echo ""
echo "Checking FASTQ..."
test -s "$CB_UMI"
test -s "$CDNA"
echo "FASTQ files: OK"

# ------------------------------------------------
# Run STARsolo
# ------------------------------------------------

echo ""
echo "Starting STARsolo..."
echo "Start: $(date)"

"$STAR" \
    --runThreadN 16 \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$CDNA" "$CB_UMI" \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 \
    --soloCBlen "$CBLEN" \
    --soloUMIstart 17 \
    --soloUMIlen "$UMILEN" \
    --soloCBwhitelist "$WHITELIST" \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloMultiMappers EM \
    --outFilterMultimapNmax 500 \
    --outFilterMultimapScoreRange 5 \
    --outFilterScoreMin 30 \
    --clipAdapterType CellRanger4 \
    --limitOutSJcollapsed 5000000 \
    --outSAMattributes NH HI AS NM nM MD CR CY UR UY CB UB GX GN sS sQ sM \
    --outSAMunmapped Within \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$OUTDIR/"

echo ""
echo "=============================================="
echo "STARsolo finished successfully"
echo "SRR: $SRR"
echo "GSM: $GSM"
echo "End: $(date)"
echo "Output: $OUTDIR"
echo "=============================================="
