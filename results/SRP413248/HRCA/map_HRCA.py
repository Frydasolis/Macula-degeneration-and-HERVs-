import scanpy as sc
import pandas as pd
import numpy as np
import os

input_file = "SRP413248_HRCA_annotated.h5ad"

print("="*70)
print("LOADING HRCA ANNOTATED DATA")
print("="*70)

adata = sc.read_h5ad(input_file)

print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)
print("Obs:")
print(adata.obs.columns.tolist())

# ------------------------------------------------------------
# Find HRCA annotation column
# ------------------------------------------------------------

possible = [
    "HRCA_celltype",
    "scANVI_predictions",
    "predicted_celltype",
    "celltype",
    "CellType"
]

annotation = None

for col in possible:
    if col in adata.obs.columns:
        annotation = col
        break

if annotation is None:
    raise ValueError(
        "No HRCA annotation column found. Available columns:\n"
        + "\n".join(adata.obs.columns)
    )

print("\nUsing annotation:", annotation)

print("\nNumber of cell types:", adata.obs[annotation].nunique())

print("\nTop cell types:")
print(adata.obs[annotation].value_counts().head(30))

# ------------------------------------------------------------
# Check whether dimensional reductions already exist
# ------------------------------------------------------------

print("\nExisting embeddings:")
print(list(adata.obsm.keys()))

# ------------------------------------------------------------
# Calculate PCA/UMAP if necessary
# ------------------------------------------------------------

if "X_umap" not in adata.obsm:

    print("\nNo UMAP found.")

    if "X_pca" not in adata.obsm:

        print("Calculating PCA...")

        # Normalize only for visualization
        adata_map = adata.copy()

        sc.pp.normalize_total(
            adata_map,
            target_sum=1e4
        )

        sc.pp.log1p(adata_map)

        sc.pp.highly_variable_genes(
            adata_map,
            n_top_genes=3000,
            flavor="seurat"
        )

        sc.tl.pca(
            adata_map,
            n_comps=50,
            use_highly_variable=True
        )

        sc.pp.neighbors(
            adata_map,
            n_neighbors=15,
            n_pcs=30
        )

        sc.tl.umap(
            adata_map,
            random_state=42
        )

        adata.obsm["X_pca"] = adata_map.obsm["X_pca"]
        adata.obsm["X_umap"] = adata_map.obsm["X_umap"]

    else:

        print("Using existing PCA.")

        sc.pp.neighbors(
            adata,
            n_neighbors=15,
            n_pcs=30
        )

        sc.tl.umap(
            adata,
            random_state=42
        )

else:

    print("Using existing UMAP.")

# ------------------------------------------------------------
# GLOBAL HRCA MAP
# ------------------------------------------------------------

print("\nGenerating global HRCA map...")

sc.pl.umap(
    adata,
    color=annotation,
    legend_loc="right margin",
    legend_fontsize=6,
    frameon=False,
    size=8,
    alpha=0.7,
    title="SRP413248 — HRCA cell-type annotation",
    show=False,
    save="_HRCA_global.png"
)

sc.pl.umap(
    adata,
    color=annotation,
    legend_loc="right margin",
    legend_fontsize=5,
    frameon=False,
    size=8,
    alpha=0.7,
    title="SRP413248 — HRCA cell-type annotation",
    show=False,
    save="_HRCA_global.pdf"
)

print("\nDONE")

print("\nFiles created:")
print("  figures/umap_HRCA_global.png")
print("  figures/umap_HRCA_global.pdf")

