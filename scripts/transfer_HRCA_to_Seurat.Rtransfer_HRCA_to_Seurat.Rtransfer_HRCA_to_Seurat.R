
############################################################
# HRCA annotation transfer to Seurat
# ARMD - SRP413248
#
# Steps:
# 1. Transfer HRCA cell-type predictions to Seurat
# 2. Group HRCA cell types into major retinal classes
# 3. Generate UMAPs using the original Seurat UMAP
# 4. Save final annotated Seurat object
############################################################


############################
# 0. Libraries
############################

library(Seurat)
library(ggplot2)


############################
# 1. Paths
############################

seurat_file <- paste0(
    "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
    "Seurat/downstream/",
    "SRP413248_GeneTE_QC_normalized_PCA_UMAP_clusters.rds"
)

hrca_predictions_file <- paste0(
    "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
    "HRCA/SRP413248_HRCA_predictions.csv"
)

output_dir <- paste0(
    "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
    "HRCA/figures"
)

output_seurat <- paste0(
    "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
    "Seurat/downstream/",
    "SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)


############################
# 2. Load original Seurat object
############################

cat("Loading Seurat object...\n")

seu <- readRDS(seurat_file)

cat("Cells:", ncol(seu), "\n")
cat("Genes:", nrow(seu), "\n")

cat("\nAvailable reductions:\n")
print(Reductions(seu))


############################
# 3. Load HRCA predictions
############################

cat("\nLoading HRCA predictions...\n")

hrca <- read.csv(
    hrca_predictions_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# First column contains the Seurat-compatible cell IDs
hrca$cell_id <- hrca[[1]]

cat("HRCA cells:", nrow(hrca), "\n")

cat("\nHRCA prediction columns:\n")
print(colnames(hrca))


############################
# 4. Match HRCA cells to Seurat cells
############################

cat("\nMatching HRCA annotations to Seurat...\n")

seurat_cells <- colnames(seu)
hrca_cells <- hrca$cell_id

common_cells <- intersect(
    seurat_cells,
    hrca_cells
)

cat("Seurat cells:", length(seurat_cells), "\n")
cat("HRCA cells:", length(hrca_cells), "\n")
cat("Common cells:", length(common_cells), "\n")
cat(
    "Seurat cells without HRCA annotation:",
    sum(!seurat_cells %in% hrca_cells),
    "\n"
)


############################
# 5. Transfer HRCA cell types
############################

cat("\nTransferring HRCA annotations...\n")

hrca_celltype <- setNames(
    hrca$HRCA_celltype,
    hrca$cell_id
)

hrca_probability <- setNames(
    hrca$HRCA_prediction_probability,
    hrca$cell_id
)

seu$HRCA_celltype <- hrca_celltype[colnames(seu)]

seu$HRCA_prediction_probability <-
    hrca_probability[colnames(seu)]


############################
# 6. Validate transferred annotations
############################

cat("\nAnnotation validation:\n")

cat(
    "Annotated Seurat cells:",
    sum(!is.na(seu$HRCA_celltype)),
    "\n"
)

cat(
    "Unannotated Seurat cells:",
    sum(is.na(seu$HRCA_celltype)),
    "\n"
)

cat("\nHRCA cell-type distribution:\n")

print(
    sort(
        table(seu$HRCA_celltype),
        decreasing = TRUE
    )
)


############################
# 7. Define major retinal classes
############################

cat("\nAssigning major retinal classes...\n")

seu$HRCA_majorclass <- NA_character_


#-----------------------------------------------------------
# Photoreceptors
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype %in% c(
        "Rod",
        "ML_Cone",
        "S_Cone"
    )
] <- "Photoreceptor"


#-----------------------------------------------------------
# Horizontal cells
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype %in% c(
        "H1",
        "H2"
    )
] <- "Horizontal"


#-----------------------------------------------------------
# Bipolar cells
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype %in% c(
        "BB",
        "DB1",
        "DB2",
        "DB3a",
        "DB3b",
        "DB4a",
        "DB4b",
        "DB5",
        "DB6",
        "FMB",
        "GB",
        "IMB",
        "RB",
        "OFFx"
    )
] <- "Bipolar"


#-----------------------------------------------------------
# Amacrine cells
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype %in% c(
        "AII",
        "AII*",
        "CAI",
        "CAII",
        "nNOS",
        "OFF-SAC",
        "OFF-SAC*",
        "ON-SAC",
        "PENK",
        "SEG",
        "SEG*",
        "VG3",
        "VG3*",
        "VIP"
    )
] <- "Amacrine"

# HAC subtypes
seu$HRCA_majorclass[
    grepl(
        "^HAC",
        seu$HRCA_celltype
    )
] <- "Amacrine"


#-----------------------------------------------------------
# Retinal ganglion cells
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype %in% c(
        "MG_ON",
        "MG_OFF",
        "PG_ON",
        "PG_OFF",
        "ipRGC"
    )
] <- "RGC"

# HRGC subtypes
seu$HRCA_majorclass[
    grepl(
        "^HRGC",
        seu$HRCA_celltype
    )
] <- "RGC"


#-----------------------------------------------------------
# Müller glia
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype == "MG"
] <- "Muller glia"


#-----------------------------------------------------------
# Astrocytes
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype == "Astrocyte"
] <- "Astrocyte"


#-----------------------------------------------------------
# Microglia
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype == "Microglia"
] <- "Microglia"


#-----------------------------------------------------------
# Retinal pigment epithelium
#-----------------------------------------------------------

seu$HRCA_majorclass[
    seu$HRCA_celltype == "RPE"
] <- "RPE"


#-----------------------------------------------------------
# Remaining / unknown cells
#-----------------------------------------------------------

seu$HRCA_majorclass[
    is.na(seu$HRCA_majorclass)
] <- "Other"


############################
# 8. Check major classes
############################

cat("\nMajor cell-class distribution:\n")

major_counts <- sort(
    table(seu$HRCA_majorclass),
    decreasing = TRUE
)

print(major_counts)

cat(
    "\nTotal cells:",
    ncol(seu),
    "\n"
)

cat(
    "Assigned cells:",
    sum(!is.na(seu$HRCA_majorclass)),
    "\n"
)

cat(
    "Other cells:",
    sum(seu$HRCA_majorclass == "Other"),
    "\n"
)


############################
# 9. Create output directory
############################

if (!dir.exists(output_dir)) {
    dir.create(
        output_dir,
        recursive = TRUE
    )
}


############################
# 10. UMAP - HRCA cell types
############################

cat("\nGenerating HRCA cell-type UMAP...\n")

p_celltype <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "HRCA_celltype",
    label = FALSE,
    raster = TRUE
) +
    ggtitle(
        "ARMD - HRCA cell-type annotation"
    )


############################
# 11. Save cell-type UMAP
############################

ggsave(
    filename = file.path(
        output_dir,
        "UMAP_HRCA_celltypes.png"
    ),
    plot = p_celltype,
    width = 14,
    height = 10,
    dpi = 300,
    limitsize = FALSE
)

ggsave(
    filename = file.path(
        output_dir,
        "UMAP_HRCA_celltypes.pdf"
    ),
    plot = p_celltype,
    width = 14,
    height = 10,
    limitsize = FALSE
)


############################
# 12. UMAP - major classes
############################

cat("\nGenerating major-class UMAP...\n")

p_major <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "HRCA_majorclass",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
) +
    ggtitle(
        "ARMD - HRCA major cell classes"
    )


############################
# 13. Save major-class UMAP
############################

ggsave(
    filename = file.path(
        output_dir,
        "UMAP_HRCA_majorclasses.png"
    ),
    plot = p_major,
    width = 14,
    height = 10,
    dpi = 300,
    limitsize = FALSE
)

ggsave(
    filename = file.path(
        output_dir,
        "UMAP_HRCA_majorclasses.pdf"
    ),
    plot = p_major,
    width = 14,
    height = 10,
    limitsize = FALSE
)


############################
# 14. Save final Seurat object
############################

cat("\nSaving final Seurat object...\n")

saveRDS(
    seu,
    output_seurat
)


############################
# 15. Final validation
############################

cat("\n========================================\n")
cat("HRCA ANNOTATION PIPELINE COMPLETE\n")
cat("========================================\n\n")

cat(
    "Cells:",
    ncol(seu),
    "\n"
)

cat(
    "Genes:",
    nrow(seu),
    "\n"
)

cat(
    "HRCA annotations:",
    sum(!is.na(seu$HRCA_celltype)),
    "\n"
)

cat(
    "Major classes:",
    sum(!is.na(seu$HRCA_majorclass)),
    "\n"
)

cat(
    "Other:",
    sum(seu$HRCA_majorclass == "Other"),
    "\n"
)

cat(
    "\nFinal Seurat object:\n",
    output_seurat,
    "\n"
)

cat(
    "\nFigures:\n",
    output_dir,
    "\n"
)
