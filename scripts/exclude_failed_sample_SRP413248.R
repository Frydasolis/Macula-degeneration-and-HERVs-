#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
})

BASE <- "/storage/lemus_g/roldan/ARMD"

INPUT <- file.path(
    BASE,
    "results/SRP413248/Seurat/SRP413248_GeneTE_with_metadata.rds"
)

OUTPUT <- file.path(
    BASE,
    "results/SRP413248/Seurat/SRP413248_GeneTE_QC_noGSM6841158.rds"
)

EXCLUDED_SAMPLE <- "GSM6841158"

# ---------------------------------------------------------
# Load object
# ---------------------------------------------------------

obj <- readRDS(INPUT)

cat("============================================\n")
cat("Sample exclusion\n")
cat("============================================\n")

cat(
    "Original cells:",
    ncol(obj),
    "\n"
)

cat(
    "Original samples:",
    length(unique(obj$sample)),
    "\n\n"
)

# ---------------------------------------------------------
# Verify sample exists
# ---------------------------------------------------------

if (!EXCLUDED_SAMPLE %in% unique(obj$sample)) {
    stop(
        paste(
            "Sample not found:",
            EXCLUDED_SAMPLE
        )
    )
}

# ---------------------------------------------------------
# Exclude failed sample
# ---------------------------------------------------------

obj_qc <- subset(
    obj,
    subset = sample != EXCLUDED_SAMPLE
)

# ---------------------------------------------------------
# Add QC exclusion information
# ---------------------------------------------------------

obj_qc$QC_sample_status <- "Included"

# ---------------------------------------------------------
# Verify exclusion
# ---------------------------------------------------------

if (
    EXCLUDED_SAMPLE %in%
    unique(obj_qc$sample)
) {
    stop(
        "Excluded sample is still present."
    )
}

cat(
    "Excluded sample:",
    EXCLUDED_SAMPLE,
    "\n"
)

cat(
    "Cells removed:",
    sum(obj$sample == EXCLUDED_SAMPLE),
    "\n"
)

cat(
    "Remaining cells:",
    ncol(obj_qc),
    "\n"
)

cat(
    "Remaining samples:",
    length(unique(obj_qc$sample)),
    "\n\n"
)

cat("Cells by disease:\n")
print(
    table(
        obj_qc$Disease,
        useNA = "ifany"
    )
)

cat("\nCells by sample:\n")
print(
    table(
        obj_qc$sample
    )
)

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

saveRDS(
    obj_qc,
    OUTPUT
)

cat("\nSaved:\n")
cat(
    OUTPUT,
    "\n"
)

# ---------------------------------------------------------
# Exclusion record
# ---------------------------------------------------------

record <- file.path(
    dirname(OUTPUT),
    "SRP413248_sample_exclusion.txt"
)

cat(
    paste0(
        "Excluded sample\t", EXCLUDED_SAMPLE, "\n",
        "Reason\tSample-level technical QC failure\n",
        "Evidence\t79.36% reads unmapped as too short; ",
        "median 35 UMI/cell; median 33 genes/cell; ",
        "24.4% valid barcodes\n",
        "Original_cells\t",
        sum(obj$sample == EXCLUDED_SAMPLE),
        "\n",
        "Remaining_cells\t",
        ncol(obj_qc),
        "\n"
    ),
    file = record
)

cat(
    "\nExclusion record:\n",
    record,
    "\n"
)
