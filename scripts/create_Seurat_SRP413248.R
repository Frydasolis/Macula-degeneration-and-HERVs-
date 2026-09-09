#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

BASE <- "/storage/lemus_g/roldan/ARMD"

INPUT_DIR <- file.path(
    BASE,
    "results/SRP413248/final_gene_TE"
)

OUTPUT_DIR <- file.path(
    BASE,
    "results/SRP413248/Seurat"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

SAMPLES <- c(
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
    "GSM6841159"
)

cat("====================================================\n")
cat("Creating Seurat object: SRP413248\n")
cat("====================================================\n\n")

seurat_objects <- list()

feature_reference <- NULL

for (sample in SAMPLES) {

    cat("-----------------------------------------------\n")
    cat("Sample:", sample, "\n")
    cat("-----------------------------------------------\n")

    sample_dir <- file.path(
        INPUT_DIR,
        sample
    )

    matrix_file <- file.path(
        sample_dir,
        "matrix.mtx"
    )

    features_file <- file.path(
        sample_dir,
        "features.tsv"
    )

    barcodes_file <- file.path(
        sample_dir,
        "barcodes.tsv"
    )

    required <- c(
        matrix_file,
        features_file,
        barcodes_file
    )

    missing <- required[
        !file.exists(required)
    ]

    if (length(missing) > 0) {
        stop(
            paste(
                "Missing files for",
                sample,
                ":",
                paste(missing, collapse = ", ")
            )
        )
    }

    # ------------------------------------------------------
    # Read matrix
    # ------------------------------------------------------

    counts <- readMM(
        matrix_file
    )

    counts <- as(
        counts,
        "dgCMatrix"
    )

    # ------------------------------------------------------
    # Read features
    # ------------------------------------------------------

    features <- read.delim(
        features_file,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE
    )

    if (ncol(features) < 2) {
        stop(
            paste(
                "Unexpected features.tsv format:",
                sample
            )
        )
    }

    feature_id <- features[[1]]
    feature_name <- features[[2]]

    if (ncol(features) >= 3) {
        feature_type <- features[[3]]
    } else {
        feature_type <- rep(
            "Gene",
            length(feature_id)
        )
    }

    # ------------------------------------------------------
    # Ensure feature IDs are unique
    # ------------------------------------------------------

    feature_id <- make.unique(
        feature_id
    )

    rownames(counts) <- feature_id

    # ------------------------------------------------------
    # Read barcodes
    # ------------------------------------------------------

    barcodes <- readLines(
        barcodes_file
    )

    if (nrow(counts) != length(feature_id)) {
        stop(
            paste(
                sample,
                ": number of matrix rows does not match features.tsv"
            )
        )
    }

    if (ncol(counts) != length(barcodes)) {
        stop(
            paste(
                sample,
                ": number of matrix columns does not match barcodes.tsv"
            )
        )
    }

    # ------------------------------------------------------
    # Store feature annotation from first sample
    # ------------------------------------------------------

    current_features <- data.frame(
        Feature_ID = feature_id,
        Feature_Name = feature_name,
        Feature_Type = feature_type,
        stringsAsFactors = FALSE
    )

    if (is.null(feature_reference)) {

        feature_reference <- current_features

    } else {

        if (!identical(
            feature_reference$Feature_ID,
            current_features$Feature_ID
        )) {

            stop(
                paste(
                    "Feature IDs differ between samples:",
                    sample
                )
            )
        }
    }

    # ------------------------------------------------------
    # Make cell names sample-specific
    # ------------------------------------------------------

    original_barcodes <- barcodes

    cell_names <- paste(
        sample,
        barcodes,
        sep = "_"
    )

    colnames(counts) <- cell_names

    # ------------------------------------------------------
    # Create Seurat object
    # ------------------------------------------------------

    obj <- CreateSeuratObject(
        counts = counts,
        assay = "RNA",
        project = sample,
        min.cells = 0,
        min.features = 0
    )

    # ------------------------------------------------------
    # Add basic sample metadata
    # ------------------------------------------------------

    obj$sample <- sample

    obj$GSM <- sample

    obj$original_barcode <- original_barcodes

    # ------------------------------------------------------
    # Save individual Seurat object
    # ------------------------------------------------------

    individual_file <- file.path(
        OUTPUT_DIR,
        paste0(sample, "_Seurat.rds")
    )

    saveRDS(
        obj,
        individual_file
    )

    cat(
        "Cells:",
        ncol(obj),
        "\n"
    )

    cat(
        "Features:",
        nrow(obj),
        "\n"
    )

    cat(
        "Total counts:",
        sum(GetAssayData(obj, layer = "counts")),
        "\n"
    )

    seurat_objects[[sample]] <- obj
}

cat("\n====================================================\n")
cat("Merging 17 samples\n")
cat("====================================================\n")

# ----------------------------------------------------------
# Merge
# ----------------------------------------------------------

combined <- merge(
    x = seurat_objects[[1]],
    y = seurat_objects[2:length(seurat_objects)],
    merge.data = TRUE
)

# Seurat v5 creates multiple count layers after merge.
# Join them before accessing the combined counts layer.

combined[["RNA"]] <- JoinLayers(
    combined[["RNA"]]
)

# ----------------------------------------------------------
# Add feature annotations as a separate object
# ----------------------------------------------------------

feature_annotation_file <- file.path(
    OUTPUT_DIR,
    "feature_annotation.tsv"
)

write.table(
    feature_reference,
    feature_annotation_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# ----------------------------------------------------------
# Identify TE and gene features
# ----------------------------------------------------------

gene_features <- feature_reference$Feature_ID[
    feature_reference$Feature_Type != "TE"
]

te_features <- feature_reference$Feature_ID[
    feature_reference$Feature_Type == "TE"
]

# Save lists
writeLines(
    gene_features,
    file.path(
        OUTPUT_DIR,
        "canonical_gene_features.txt"
    )
)

writeLines(
    te_features,
    file.path(
        OUTPUT_DIR,
        "TE_features.txt"
    )
)

# ----------------------------------------------------------
# Global metadata
# ----------------------------------------------------------

combined$nCount_RNA <- Matrix::colSums(
    GetAssayData(
        combined,
        assay = "RNA",
        layer = "counts"
    )
)

combined$nFeature_RNA <- Matrix::colSums(
    GetAssayData(
        combined,
        assay = "RNA",
        layer = "counts"
    ) > 0
)

# ----------------------------------------------------------
# Save combined object
# ----------------------------------------------------------

combined_file <- file.path(
    OUTPUT_DIR,
    "SRP413248_GeneTE_merged.rds"
)

saveRDS(
    combined,
    combined_file
)

# ----------------------------------------------------------
# Save metadata
# ----------------------------------------------------------

metadata_file <- file.path(
    OUTPUT_DIR,
    "SRP413248_cell_metadata.tsv"
)

meta <- combined@meta.data

meta$cell <- rownames(meta)

write.table(
    meta,
    metadata_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

# ----------------------------------------------------------
# Summary
# ----------------------------------------------------------

summary_file <- file.path(
    OUTPUT_DIR,
    "SRP413248_Seurat_summary.txt"
)

cat(
    paste0(
        "Samples\t", length(SAMPLES), "\n",
        "Cells\t", ncol(combined), "\n",
        "Features\t", nrow(combined), "\n",
        "Canonical_gene_features\t", length(gene_features), "\n",
        "TE_features\t", length(te_features), "\n",
        "Total_counts\t",
        sum(
            GetAssayData(
                combined,
                assay = "RNA",
                layer = "counts"
            )
        ),
        "\n"
    ),
    file = summary_file
)

cat("\n====================================================\n")
cat("DONE\n")
cat("====================================================\n")

cat(
    "Seurat object:\n",
    combined_file,
    "\n\n"
)

cat(
    "Metadata:\n",
    metadata_file,
    "\n\n"
)

cat(
    "Feature annotation:\n",
    feature_annotation_file,
    "\n\n"
)

cat(
    "Summary:\n",
    summary_file,
    "\n"
)
