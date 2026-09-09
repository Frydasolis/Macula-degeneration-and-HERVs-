#!/usr/bin/env Rscript

library(Seurat)

BASE <- "/storage/lemus_g/roldan/ARMD"

seurat_file <- file.path(
    BASE,
    "results/SRP413248/Seurat/SRP413248_GeneTE_merged.rds"
)

metadata_file <- file.path(
    BASE,
    "resourses/SraRunTable-7.csv"
)

output_file <- file.path(
    BASE,
    "results/SRP413248/Seurat/SRP413248_GeneTE_with_metadata.rds"
)

# ------------------------------------------------------------
# Load Seurat
# ------------------------------------------------------------

obj <- readRDS(seurat_file)

# ------------------------------------------------------------
# Read SRA metadata
# ------------------------------------------------------------

meta <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# ------------------------------------------------------------
# One row per GSM
# ------------------------------------------------------------

meta <- meta[
    !duplicated(meta$`Sample Name`),
]

# ------------------------------------------------------------
# Rename relevant columns
# ------------------------------------------------------------

meta_clean <- meta[, c(
    "Sample Name",
    "BioSample",
    "AGE",
    "disease_state",
    "cell_type",
    "tissue"
)]

colnames(meta_clean) <- c(
    "sample",
    "BioSample",
    "Age",
    "Disease",
    "CellType",
    "Tissue"
)

rownames(meta_clean) <- meta_clean$sample

# ------------------------------------------------------------
# Check GSM matching
# ------------------------------------------------------------

gsm_obj <- unique(obj$sample)
gsm_meta <- meta_clean$sample

missing_in_meta <- setdiff(
    gsm_obj,
    gsm_meta
)

missing_in_obj <- setdiff(
    gsm_meta,
    gsm_obj
)

if (length(missing_in_meta) > 0) {
    stop(
        paste(
            "GSM missing from metadata:",
            paste(missing_in_meta, collapse = ", ")
        )
    )
}

if (length(missing_in_obj) > 0) {
    warning(
        paste(
            "Metadata GSM not present in Seurat:",
            paste(missing_in_obj, collapse = ", ")
        )
    )
}

# ------------------------------------------------------------
# Match metadata to cells
# ------------------------------------------------------------

cell_metadata <- meta_clean[
    match(obj$sample, rownames(meta_clean)),
    ,
    drop = FALSE
]

rownames(cell_metadata) <- colnames(obj)

# Do not add sample twice
cell_metadata$sample <- obj$sample

# ------------------------------------------------------------
# Add metadata
# ------------------------------------------------------------

obj <- AddMetaData(
    object = obj,
    metadata = cell_metadata[
        ,
        c(
            "BioSample",
            "Age",
            "Disease",
            "CellType",
            "Tissue"
        )
    ]
)

# ------------------------------------------------------------
# Convert Disease to factor
# ------------------------------------------------------------

obj$Disease <- factor(
    obj$Disease,
    levels = c(
        "Healthy Control",
        "Dry AMD",
        "Wet AMD"
    )
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
    obj,
    output_file
)

# ------------------------------------------------------------
# Print checks
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Metadata successfully added\n")
cat("========================================\n")

cat("\nCells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")

cat("\nDisease:\n")
print(table(obj$Disease, useNA = "ifany"))

cat("\nAge by sample:\n")
print(
    unique(
        obj@meta.data[
            ,
            c("sample", "Age", "Disease")
        ]
    )
)

cat("\nSaved:\n")
cat(output_file, "\n")
