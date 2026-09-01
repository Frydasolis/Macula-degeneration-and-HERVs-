#!/bin/bash

BASE=/storage/lemus_g/roldan/ARMD
METADATA=$BASE/resourses/SRP413248_metadata.tsv
FASTQ=$BASE/results/SRP413248/fastq
OUT=$BASE/results/SRP413248/cellranger_fastq

mkdir -p "$OUT"

tail -n +2 "$METADATA" | while IFS=$'\t' read -r SRR EXPERIMENT GSM
do

    mkdir -p "$OUT/$GSM"

    ln -sf "$FASTQ/${SRR}_2.fastq" \
        "$OUT/$GSM/${SRR}_R1.fastq"

    ln -sf "$FASTQ/${SRR}_3.fastq" \
        "$OUT/$GSM/${SRR}_R2.fastq"

done

echo "=========================================="
echo "FASTQ organization completed"
echo "=========================================="

find "$OUT" -mindepth 1 -maxdepth 1 -type d | sort
