#!/usr/bin/env python3

import os
import pysam
import pandas as pd
import argparse


BASE = "/storage/lemus_g/roldan/ARMD"

STAR_DIR = f"{BASE}/results/SRP413248/STARsolo_17samples"
TE_DIR = f"{BASE}/results/SRP413248/Stellarscope/assign"
OUT_DIR = f"{BASE}/results/SRP413248/Stellarscope/overlap_simple"


def extract_star_molecules(star_bam):

    molecules = set()

    bam = pysam.AlignmentFile(star_bam, "rb")

    for read in bam:

        if read.has_tag("CB") and read.has_tag("UB"):

            cb = read.get_tag("CB")
            ub = read.get_tag("UB")

            molecules.add(f"{cb}_{ub}")

    bam.close()

    return molecules


def find_te_overlap(te_bam, star_molecules):

    overlap = []

    bam = pysam.AlignmentFile(te_bam, "rb")

    for read in bam:

        if read.has_tag("CB") and read.has_tag("UB"):

            cb = read.get_tag("CB")
            ub = read.get_tag("UB")

            molecule_id = f"{cb}_{ub}"

            if molecule_id in star_molecules:

                overlap.append({
                    "Cell": cb,
                    "UMI": ub
                })

    bam.close()

    return pd.DataFrame(overlap)


def main():

    parser = argparse.ArgumentParser(
        description="STARsolo/Stellarscope CB+UMI overlap"
    )

    parser.add_argument(
        "--sample",
        required=True
    )

    args = parser.parse_args()

    sample = args.sample

    star_bam = (
        f"{STAR_DIR}/{sample}/"
        "Aligned.sortedByCoord.out.bam"
    )

    te_bam = (
        f"{TE_DIR}/{sample}/"
        f"{sample}-updated.bam"
    )

    os.makedirs(
        OUT_DIR,
        exist_ok=True
    )

    output = (
        f"{OUT_DIR}/"
        f"{sample}_gene_TE_overlap.tsv"
    )

    print("Extracting STARsolo molecules...")

    star_molecules = extract_star_molecules(
        star_bam
    )

    print(
        f"STARsolo molecules detected: "
        f"{len(star_molecules):,}"
    )

    print("Searching TE overlap...")

    overlap = find_te_overlap(
        te_bam,
        star_molecules
    )

    overlap = overlap.drop_duplicates()

    overlap.to_csv(
        output,
        sep="\t",
        index=False
    )

    print(
        f"Overlapping molecules: "
        f"{len(overlap):,}"
    )

    print(
        f"Saved correction table: {output}"
    )


if __name__ == "__main__":
    main()
