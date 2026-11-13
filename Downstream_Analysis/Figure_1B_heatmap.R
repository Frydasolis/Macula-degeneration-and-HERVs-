# ============================================================
# FIGURE 1B
# TE family abundance across retinal cell populations in AMD
#
# 61 TE families
# 7 retinal cell types
# 3 disease conditions
#
# Sample-level CPM -> mean CPM by group ->
# log2(CPM + 1) -> row-wise Z-score
#
# Descriptive heatmap; no significance annotations.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
})

# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

base <- file.path(
  "results",
  "SRP413248",
  "TE_family_subclustering",
  "Figure1B_pseudobulk"
)

outdir <- file.path(base, "Figure1B_FINAL")

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load sample-level TE family abundance
# ------------------------------------------------------------

dat <- read.csv(
  file.path(base, "TE_family_sample_abundance_CPM.csv"),
  stringsAsFactors = FALSE
)

required_columns <- c(
  "TE_family",
  "sample",
  "CPM",
  "Cell_type",
  "Disease"
)

stopifnot(all(required_columns %in% colnames(dat)))

# ------------------------------------------------------------
# 3. Define cell types and disease conditions
# ------------------------------------------------------------

celltypes <- c(
  "Amacrine",
  "Astrocyte",
  "Bipolar",
  "Horizontal",
  "Muller glia",
  "Photoreceptor",
  "RGC"
)

conditions <- c(
  "Healthy Control",
  "Dry AMD",
  "Wet AMD"
)

dat <- dat %>%
  filter(
    Cell_type %in% celltypes,
    Disease %in% conditions
  ) %>%
  mutate(
    Cell_type = factor(Cell_type, levels = celltypes),
    Disease = factor(Disease, levels = conditions)
  )

# ------------------------------------------------------------
# 4. Average sample-level CPM within each group
# ------------------------------------------------------------

summary_dat <- dat %>%
  group_by(TE_family, Cell_type, Disease) %>%
  summarise(
    mean_CPM = mean(CPM, na.rm = TRUE),
    median_CPM = median(CPM, na.rm = TRUE),
    n_samples = n_distinct(sample),
    .groups = "drop"
  )

write.csv(
  summary_dat,
  file.path(outdir, "TE_family_group_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Build TE family x cell type/disease matrix
# ------------------------------------------------------------

summary_dat <- summary_dat %>%
  mutate(
    Group = paste(
      as.character(Cell_type),
      as.character(Disease),
      sep = "__"
    )
  )

mat_df <- summary_dat %>%
  select(TE_family, Group, mean_CPM) %>%
  pivot_wider(
    names_from = Group,
    values_from = mean_CPM
  )

mat <- mat_df %>%
  column_to_rownames("TE_family") %>%
  as.matrix()

expected_columns <- unlist(
  lapply(celltypes, function(ct) {
    paste(ct, conditions, sep = "__")
  })
)

missing_columns <- setdiff(
  expected_columns,
  colnames(mat)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing cell type/disease combinations: ",
    paste(missing_columns, collapse = ", ")
  )
}

mat <- mat[, expected_columns, drop = FALSE]

# A missing group is not equivalent to zero expression.
if (anyNA(mat)) {
  stop(
    "Missing TE-family/group values detected. ",
    "Check sample coverage before generating the heatmap."
  )
}

# ------------------------------------------------------------
# 6. Log transformation
# ------------------------------------------------------------

mat_log <- log2(mat + 1)

# ------------------------------------------------------------
# 7. Row-wise Z-score
# ------------------------------------------------------------

row_means <- rowMeans(mat_log)

row_sds <- apply(mat_log, 1, sd)

mat_z <- sweep(
  mat_log,
  1,
  row_means,
  "-"
)

nonzero_sd <- row_sds > 0 & is.finite(row_sds)

mat_z[nonzero_sd, ] <- sweep(
  mat_z[nonzero_sd, , drop = FALSE],
  1,
  row_sds[nonzero_sd],
  "/"
)

# Families with zero variance have no relative pattern.
mat_z[!nonzero_sd, ] <- 0

# ------------------------------------------------------------
# 8. Save matrices
# ------------------------------------------------------------

write.csv(
  mat,
  file.path(outdir, "TE_family_mean_CPM_matrix.csv")
)

write.csv(
  mat_log,
  file.path(outdir, "TE_family_log2CPM_matrix.csv")
)

write.csv(
  mat_z,
  file.path(outdir, "TE_family_Zscore_matrix.csv")
)

# ------------------------------------------------------------
# 9. Column annotations
# ------------------------------------------------------------

annotation_col <- data.frame(
  Cell_type = rep(celltypes, each = 3),
  Disease = rep(conditions, times = 7)
)

rownames(annotation_col) <- colnames(mat_z)

annotation_col$Cell_type <- factor(
  annotation_col$Cell_type,
  levels = celltypes
)

annotation_col$Disease <- factor(
  annotation_col$Disease,
  levels = conditions
)

disease_colors <- c(
  "Healthy Control" = "#607D8B",
  "Dry AMD" = "#E5A04B",
  "Wet AMD" = "#B85B69"
)

celltype_colors <- setNames(
  brewer.pal(7, "Set2"),
  celltypes
)

annotation_colors <- list(
  Cell_type = celltype_colors,
  Disease = disease_colors
)

# ------------------------------------------------------------
# 10. Heatmap colors
# ------------------------------------------------------------

heat_colors <- colorRampPalette(
  c(
    "#2166AC",
    "#67A9CF",
    "#F7F7F7",
    "#EF8A62",
    "#B2182B"
  )
)(101)

# Fixed symmetric scale
breaks <- seq(-2.5, 2.5, length.out = 102)

# ------------------------------------------------------------
# 11. Display labels
# ------------------------------------------------------------

column_labels <- rep(
  c("Healthy", "Dry", "Wet"),
  times = length(celltypes)
)

# ------------------------------------------------------------
# 12. Generate heatmap
# ------------------------------------------------------------

draw_heatmap <- function(filename, device) {

  pheatmap(
    mat_z,
    color = heat_colors,
    breaks = breaks,

    cluster_rows = TRUE,
    cluster_cols = FALSE,

    annotation_col = annotation_col,
    annotation_colors = annotation_colors,

    gaps_col = seq(3, 18, by = 3),

    labels_col = column_labels,

    show_rownames = TRUE,
    show_colnames = TRUE,

    fontsize = 10,
    fontsize_row = 7,
    fontsize_col = 9,

    angle_col = 45,

    border_color = NA,

    treeheight_row = 45,
    treeheight_col = 0,

    cellwidth = 24,
    cellheight = 10,

    main = "TE family abundance across retinal cell populations",

    legend = TRUE,
    annotation_legend = TRUE,

    filename = filename,

    width = 12,
    height = 12
  )
}

draw_heatmap(
  file.path(outdir, "Figure1B_TE_family_Zscore.pdf"),
  "pdf"
)

draw_heatmap(
  file.path(outdir, "Figure1B_TE_family_Zscore.png"),
  "png"
)

# ------------------------------------------------------------
# 13. Report
# ------------------------------------------------------------

cat("\n====================================\n")
cat("FIGURE 1B COMPLETED\n")
cat("====================================\n")

cat("TE families:", nrow(mat_z), "\n")
cat("Cell type/condition groups:", ncol(mat_z), "\n")
cat("Z-score range:",
    round(min(mat_z), 3), "to",
    round(max(mat_z), 3), "\n")

cat("Zero-variance families:",
    sum(!nonzero_sd), "\n")

cat("Output directory:", outdir, "\n")
