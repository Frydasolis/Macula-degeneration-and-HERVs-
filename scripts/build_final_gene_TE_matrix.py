#!/usr/bin/env python3

import os
import argparse
import pandas as pd

from scipy.io import mmread, mmwrite
from scipy.sparse import coo_matrix, vstack, lil_matrix


BASE = "/storage/lemus_g/roldan/ARMD"

STAR_DIR = f"{BASE}/results/SRP413248/STARsolo_17samples"
TE_DIR = f"{BASE}/results/SRP413248/Stellarscope/assign"
OVERLAP_DIR = f"{BASE}/results/SRP413248/Stellarscope/overlap"

OUT_DIR = f"{BASE}/results/SRP413248/final_gene_TE"


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


def read_lines(path):
    with open(path) as f:
        return [x.rstrip("\n") for x in f]


def main(sample):

    print("=" * 70)
    print("FINAL GENE + TE MATRIX")
    print("=" * 70)
    print("Sample:", sample)

    gene_dir = (
        f"{STAR_DIR}/{sample}/Solo.out/Gene/filtered"
    )

    gene_mtx = f"{gene_dir}/matrix.mtx"
    gene_features = f"{gene_dir}/features.tsv"
    gene_barcodes = f"{gene_dir}/barcodes.tsv"

    te_dir = f"{TE_DIR}/{sample}"

    te_mtx = f"{te_dir}/{sample}-TE_counts.mtx"
    te_features = f"{te_dir}/{sample}-features.tsv"
    te_barcodes = f"{te_dir}/{sample}-barcodes.tsv"

    correction_file = (
        f"{OVERLAP_DIR}/{sample}_gene_correction.tsv"
    )

    for path in [
        gene_mtx,
        gene_features,
        gene_barcodes,
        te_mtx,
        te_features,
        te_barcodes,
        correction_file,
    ]:
        if not os.path.exists(path):
            raise FileNotFoundError(path)

    out_dir = f"{OUT_DIR}/{sample}"
    os.makedirs(out_dir, exist_ok=True)

    # ---------------------------------------------------------
    # Gene matrix
    # ---------------------------------------------------------

    print("\nLoading GeneFull matrix...")

    genes = mmread(gene_mtx).tocsr()

    gf = pd.read_csv(
        gene_features,
        sep="\t",
        header=None,
        dtype=str
    )

    gene_ids = gf.iloc[:, 0].tolist()
    gene_names = gf.iloc[:, 1].tolist()

    if gf.shape[1] >= 3:
        gene_types = gf.iloc[:, 2].tolist()
    else:
        gene_types = ["Gene"] * len(gene_ids)

    barcodes = read_lines(gene_barcodes)

    print("Gene matrix:", genes.shape)

    # ---------------------------------------------------------
    # TE matrix
    # ---------------------------------------------------------

    print("\nLoading TE matrix...")

    tes = mmread(te_mtx).tocsr()

    te_names_all = read_lines(te_features)
    te_barcodes = read_lines(te_barcodes)

    print("TE matrix:", tes.shape)

    # ---------------------------------------------------------
    # Align TE cells to STARsolo cells
    # ---------------------------------------------------------

    te_cell_index = {
        bc: i
        for i, bc in enumerate(te_barcodes)
    }

    te_order = [
        te_cell_index.get(bc, None)
        for bc in barcodes
    ]

    missing_cells = sum(
        x is None
        for x in te_order
    )

    print(
        "STARsolo cells missing from Stellarscope:",
        missing_cells
    )

    # Keep all STARsolo cells.
    # Cells absent from Stellarscope receive zero TE counts.

    te_aligned = lil_matrix(
        (
            tes.shape[0],
            len(barcodes)
        ),
        dtype=tes.dtype
    )

    for new_col, old_col in enumerate(te_order):

        if old_col is not None:
            te_aligned[:, new_col] = tes[:, old_col]

    tes = te_aligned.tocsr()

    # ---------------------------------------------------------
    # Load correction
    # ---------------------------------------------------------

    print("\nLoading correction table...")

    correction = pd.read_csv(
        correction_file,
        sep="\t",
        dtype={
            "Cell": str,
            "Gene": str,
            "Gene_name": str,
            "Overlap_UMIs": int
        }
    )

    print(
        "Correction rows:",
        f"{len(correction):,}"
    )

    # ---------------------------------------------------------
    # Dictionaries
    # ---------------------------------------------------------

    gene_index = {
        gene: i
        for i, gene in enumerate(gene_ids)
    }

    cell_index = {
        cell: i
        for i, cell in enumerate(barcodes)
    }

    # ---------------------------------------------------------
    # Apply correction safely
    # ---------------------------------------------------------

    corr_rows = []
    corr_cols = []
    corr_values = []

    requested_total = 0
    applied_total = 0
    capped_entries = 0

    conflict_rows = []

    for row in correction.itertuples(index=False):

        requested = int(row.Overlap_UMIs)

        requested_total += requested

        if row.Gene_ID not in gene_index:
            continue

        if row.Cell not in cell_index:
            continue

        r = gene_index[row.Gene_ID]
        c = cell_index[row.Cell]

        original = int(genes[r, c])

        # Never subtract more than the existing GeneFull count.
        applied = min(requested, original)

        if applied < requested:

            capped_entries += 1

            conflict_rows.append({
                "Cell": row.Cell,
                "Gene": row.Gene_ID,
                "Gene_name": row.Gene_name,
                "GeneFull_count": original,
                "Requested": requested,
                "Applied": applied,
                "Not_applied": requested - applied
            })

        if applied > 0:

            corr_rows.append(r)
            corr_cols.append(c)
            corr_values.append(applied)

            applied_total += applied

    correction_matrix = coo_matrix(
        (
            corr_values,
            (corr_rows, corr_cols)
        ),
        shape=genes.shape,
        dtype="int32"
    ).tocsr()

    correction_matrix.sum_duplicates()

    genes_corrected = (
        genes - correction_matrix
    ).tocsr()

    genes_corrected.eliminate_zeros()

    # ---------------------------------------------------------
    # Remove __no_feature from TE matrix
    # ---------------------------------------------------------

    keep_te = [
        i for i, te in enumerate(te_names_all)
        if te != "__no_feature"
    ]

    tes = tes[keep_te, :].tocsr()

    te_names = [
        te_names_all[i]
        for i in keep_te
    ]

    # ---------------------------------------------------------
    # Build combined matrix
    # ---------------------------------------------------------

    combined = vstack(
        [
            genes_corrected,
            tes
        ],
        format="csr"
    )

    combined.eliminate_zeros()

    # ---------------------------------------------------------
    # Features
    # ---------------------------------------------------------

    feature_ids = (
        gene_ids +
        [f"TE_{x}" for x in te_names]
    )

    feature_names = (
        gene_names +
        te_names
    )

    feature_types = (
        gene_types +
        ["TE"] * len(te_names)
    )

    feature_df = pd.DataFrame({
        "Feature_ID": feature_ids,
        "Feature_Name": feature_names,
        "Feature_Type": feature_types
    })

    # ---------------------------------------------------------
    # Write outputs
    # ---------------------------------------------------------

    matrix_out = f"{out_dir}/matrix.mtx"
    features_out = f"{out_dir}/features.tsv"
    barcodes_out = f"{out_dir}/barcodes.tsv"

    summary_out = f"{out_dir}/summary.txt"

    conflicts_out = (
        f"{out_dir}/gene_correction_conflicts.tsv"
    )

    mmwrite(
        matrix_out,
        combined
    )

    feature_df.to_csv(
        features_out,
        sep="\t",
        header=False,
        index=False
    )

    with open(barcodes_out, "w") as f:
        for bc in barcodes:
            f.write(bc + "\n")

    if conflict_rows:

        pd.DataFrame(
            conflict_rows
        ).to_csv(
            conflicts_out,
            sep="\t",
            index=False
        )

    else:

        with open(conflicts_out, "w") as f:
            f.write(
                "Cell\tGene\tGene_name\t"
                "GeneFull_count\tRequested\t"
                "Applied\tNot_applied\n"
            )

    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------

    with open(summary_out, "w") as f:

        f.write(f"Sample\t{sample}\n")
        f.write(f"Cells\t{len(barcodes)}\n")
        f.write(f"Canonical_genes\t{len(gene_ids)}\n")
        f.write(f"TE_features\t{len(te_names)}\n")
        f.write(f"Combined_features\t{len(feature_ids)}\n")

        f.write(
            f"GeneFull_total_before\t"
            f"{int(genes.sum())}\n"
        )

        f.write(
            f"Gene_UMIs_requested\t"
            f"{requested_total}\n"
        )

        f.write(
            f"Gene_UMIs_applied\t"
            f"{applied_total}\n"
        )

        f.write(
            f"GeneFull_total_after\t"
            f"{int(genes_corrected.sum())}\n"
        )

        f.write(
            f"TE_total\t"
            f"{int(tes.sum())}\n"
        )

        f.write(
            f"Combined_total\t"
            f"{int(combined.sum())}\n"
        )

        f.write(
            f"Correction_entries_capped\t"
            f"{capped_entries}\n"
        )

    print("\nDONE")
    print("Matrix:", matrix_out)
    print("Summary:", summary_out)
    print(
        "Gene UMIs requested:",
        f"{requested_total:,}"
    )
    print(
        "Gene UMIs applied:",
        f"{applied_total:,}"
    )
    print(
        "Correction entries capped:",
        f"{capped_entries:,}"
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample",
        required=True
    )

    args = parser.parse_args()

    main(args.sample)
