library(Seurat)
library(Matrix)
library(dplyr)
library(readr)

# ---------------------------------------------
# Input
# ---------------------------------------------

seurat_file <- "/storage/lemus_g/roldan/ARMD/results/SRP413248/Seurat/downstream/SRP413248_GeneTE_QC_UMAP_HRCA_annotated_geneSymbols.rds"

output_file <- "/storage/lemus_g/roldan/ARMD/results/SRP413248/CellChat/network_metrics/LINE1_expression_per_donor.csv"

# ---------------------------------------------
# Load Seurat object
# ---------------------------------------------

obj <- readRDS(seurat_file)

# ---------------------------------------------
# Identify TE-L1 features
# ---------------------------------------------

L1_features <- grep("^TE-L1", rownames(obj), value = TRUE)

cat("Total TE-L1 features:", length(L1_features), "\n")

# ---------------------------------------------
# Extract normalized expression
# ---------------------------------------------

L1_expr <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)[L1_features, , drop = FALSE]

# ---------------------------------------------
# Calculate total TE-L1 expression per cell
# ---------------------------------------------

L1_per_cell <- Matrix::colSums(L1_expr)

# ---------------------------------------------
# Create cell-level table
# ---------------------------------------------

L1_cell_df <- data.frame(
  Cell = names(L1_per_cell),
  LINE1_expression = as.numeric(L1_per_cell),
  GSM = obj$GSM[names(L1_per_cell)]
)

# ---------------------------------------------
# Calculate mean LINE-1 expression per donor
# ---------------------------------------------

L1_donor <- L1_cell_df %>%
  group_by(GSM) %>%
  summarise(
    LINE1_expression = mean(LINE1_expression, na.rm = TRUE),
    Cells = n(),
    .groups = "drop"
  )

# ---------------------------------------------
# Add disease state
# ---------------------------------------------

condition_df <- obj@meta.data %>%
  select(GSM, disease_state) %>%
  distinct()

L1_donor <- L1_donor %>%
  left_join(condition_df, by = "GSM") %>%
  mutate(
    Condition = case_when(
      disease_state == "Healthy Control" ~ "Healthy",
      disease_state == "Wet AMD" ~ "Wet",
      disease_state == "Dry AMD" ~ "Dry"
    )
  )

# ---------------------------------------------
# Save
# ---------------------------------------------

write_csv(L1_donor, output_file)

cat("\nOutput:\n")
cat(output_file, "\n")

cat("\nNumber of donors:", nrow(L1_donor), "\n")

print(L1_donor)
