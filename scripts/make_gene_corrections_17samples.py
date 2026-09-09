#!/usr/bin/env python3

import os
import pandas as pd


BASE = "/storage/lemus_g/roldan/ARMD"

OVERLAP_DIR = (
    f"{BASE}/results/SRP413248/Stellarscope/overlap"
)

SAMPLES = [
    "GSM6841143",
    "GSM6841144",
    "GSM6841145",
    "GSM6841146",
    "GSM6841147",
    "GSM6841148",
    "GSM6841149",
    "GSM6841150",
    "GSM6841151",
    "GSM6841152",
    "GSM6841153",
    "GSM6841154",
    "GSM6841155",
    "GSM6841156",
    "GSM6841157",
    "GSM6841158",
    "GSM6841159",
]


def main():

    summary = []

    for sample in SAMPLES:

        print("=" * 70)
        print(sample)

        input_file = (
            f"{OVERLAP_DIR}/{sample}_overlap.tsv"
        )

        output_file = (
            f"{OVERLAP_DIR}/"
            f"{sample}_gene_correction.tsv"
        )

        summary_file = (
            f"{OVERLAP_DIR}/"
            f"{sample}_gene_correction_summary.txt"
        )

        if not os.path.exists(input_file):

            print(
                f"WARNING: missing {input_file}"
            )

            continue

        df = pd.read_csv(
            input_file,
            sep="\t",
            dtype=str
        )

        print(
            f"Input rows: {len(df):,}"
        )

        # ------------------------------------------------------
        # Normalize column names
        # ------------------------------------------------------

        if "Sample" in df.columns:

            # Standard format:
            # Sample Cell UMI Gene Gene_name TE

            rename = {
                "Sample": "GSM",
                "Gene": "Gene_ID"
            }

            df = df.rename(
                columns=rename
            )

        elif "STARsolo_GX" in df.columns:

            # GSM6841149 special format:
            # Cell UMI STARsolo_GX STARsolo_GN Stellarscope_TE

            df = df.rename(
                columns={
                    "STARsolo_GX": "Gene_ID",
                    "STARsolo_GN": "Gene_name",
                    "Stellarscope_TE": "TE"
                }
            )

            if "GSM" not in df.columns:
                df["GSM"] = sample

        else:

            raise RuntimeError(
                f"Unexpected columns for {sample}: "
                f"{list(df.columns)}"
            )

        required = [
            "Cell",
            "UMI",
            "Gene_ID",
            "Gene_name"
        ]

        missing = [
            x for x in required
            if x not in df.columns
        ]

        if missing:

            raise RuntimeError(
                f"{sample}: missing columns {missing}"
            )

        # ------------------------------------------------------
        # One correction per unique Cell + UMI + Gene
        # ------------------------------------------------------

        molecule_gene = (
            df[
                [
                    "Cell",
                    "UMI",
                    "Gene_ID",
                    "Gene_name"
                ]
            ]
            .drop_duplicates()
        )

        print(
            "Unique Cell+UMI+Gene molecules:",
            f"{len(molecule_gene):,}"
        )

        # ------------------------------------------------------
        # Count overlap molecules per Cell + Gene
        # ------------------------------------------------------

        correction = (
            molecule_gene
            .groupby(
                [
                    "Cell",
                    "Gene_ID",
                    "Gene_name"
                ],
                as_index=False
            )
            .size()
            .rename(
                columns={
                    "size": "Overlap_UMIs"
                }
            )
        )

        correction = correction.sort_values(
            [
                "Cell",
                "Gene_ID"
            ]
        )

        correction.to_csv(
            output_file,
            sep="\t",
            index=False
        )

        total = int(
            correction["Overlap_UMIs"].sum()
        )

        with open(summary_file, "w") as f:

            f.write(
                f"Sample\t{sample}\n"
            )

            f.write(
                f"Input_overlap_rows\t"
                f"{len(df)}\n"
            )

            f.write(
                f"Unique_Cell_UMI_Gene_molecules\t"
                f"{len(molecule_gene)}\n"
            )

            f.write(
                f"Cells_with_overlap\t"
                f"{molecule_gene['Cell'].nunique()}\n"
            )

            f.write(
                f"Genes_with_overlap\t"
                f"{molecule_gene['Gene_ID'].nunique()}\n"
            )

            f.write(
                f"Total_gene_UMIs_to_remove\t"
                f"{total}\n"
            )

            f.write(
                f"Correction_Cell_Gene_rows\t"
                f"{len(correction)}\n"
            )

        summary.append(
            {
                "Sample": sample,
                "Input_rows": len(df),
                "Unique_Cell_UMI_Gene": len(
                    molecule_gene
                ),
                "Cells": molecule_gene[
                    "Cell"
                ].nunique(),
                "Genes": molecule_gene[
                    "Gene_ID"
                ].nunique(),
                "Gene_UMIs_to_remove": total,
                "Correction_Cell_Gene_rows":
                    len(correction)
            }
        )

        print(
            "Gene UMIs to remove:",
            f"{total:,}"
        )

    # ----------------------------------------------------------
    # Global summary
    # ----------------------------------------------------------

    summary_df = pd.DataFrame(summary)

    global_summary = (
        f"{OVERLAP_DIR}/"
        "ALL_samples_gene_correction_summary.tsv"
    )

    summary_df.to_csv(
        global_summary,
        sep="\t",
        index=False
    )

    print("\n" + "=" * 70)
    print("ALL SAMPLES COMPLETE")
    print("=" * 70)

    print(global_summary)


if __name__ == "__main__":
    main()
