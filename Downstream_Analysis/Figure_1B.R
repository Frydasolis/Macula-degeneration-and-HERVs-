# ============================================================
# FIGURE 1B — TE FAMILY HEATMAP
# Age-related macular degeneration
#
# No additional filtering of TE families or nuclei.
# 61 TE families x 39,329 nuclei
# 9 retinal cell types x 3 conditions
#
# Sample-level CPM -> mean CPM -> log2(CPM+1)
# -> row-wise Z-score
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

base <- "results/SRP413248"

outdir <- file.path(
  base,
  "TE_family_subclustering",
  "Figure1B_FINAL"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load original data
# ------------------------------------------------------------

te <- readRDS(file.path(
  base,
  "TE_family_subclustering",
  "TE_family_by_cell_sparse.rds"
))

seu <- readRDS(file.path(
  base,
  "Seurat/downstream",
  "SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
))

stopifnot(
  inherits(te, "dgCMatrix"),
  nrow(te) == 61,
  ncol(te) == 39329,
  !anyDuplicated(colnames(te)),
  all(colnames(te) %in% colnames(seu)),
  all(te@x >= 0),
  all(abs(te@x - round(te@x)) < 1e-8)
)

# ------------------------------------------------------------
# 3. Metadata
# ------------------------------------------------------------

meta <- seu@meta.data[
  colnames(te),
  c(
    "sample",
    "Disease",
    "HRCA_majorclass",
    "nCount_RNA"
  ),
  drop = FALSE
]

stopifnot(
  identical(rownames(meta), colnames(te)),
  !anyNA(meta),
  all(is.finite(meta$nCount_RNA)),
  all(meta$nCount_RNA > 0)
)

meta$sample <- as.character(meta$sample)
meta$Disease <- as.character(meta$Disease)
meta$HRCA_majorclass <- as.character(meta$HRCA_majorclass)

celltypes <- c(
  "Amacrine",
  "Astrocyte",
  "Bipolar",
  "Horizontal",
  "Microglia",
  "Muller glia",
  "Photoreceptor",
  "RGC",
  "RPE"
)

conditions <- c(
  "Healthy Control",
  "Dry AMD",
  "Wet AMD"
)

stopifnot(
  all(meta$HRCA_majorclass %in% celltypes),
  all(meta$Disease %in% conditions)
)

# ------------------------------------------------------------
# 4. Aggregate all nuclei by sample and cell type
# ------------------------------------------------------------

meta$Group <- paste(
  meta$sample,
  meta$HRCA_majorclass,
  sep = "__"
)

groups <- unique(meta$Group)

membership <- sparse.model.matrix(
  ~ 0 + factor(meta$Group, levels = groups)
)

counts <- as.matrix(
  te %*% membership
)

colnames(counts) <- groups

# Total RNA counts for each sample-cell type
libsize <- as.numeric(
  rowsum(
    meta$nCount_RNA,
    group = factor(meta$Group, levels = groups),
    reorder = FALSE
  )[, 1]
)

stopifnot(
  length(libsize) == ncol(counts),
  all(libsize > 0)
)

# ------------------------------------------------------------
# 5. Calculate CPM for all TE families
# ------------------------------------------------------------

cpm <- sweep(
  counts,
  2,
  libsize,
  "/"
) * 1e6

group_meta <- meta[
  match(groups, meta$Group),
  c("sample", "Disease", "HRCA_majorclass"),
  drop = FALSE
]

sample_cpm <- data.frame(
  TE_family = rep(rownames(cpm), times = ncol(cpm)),
  Group = rep(colnames(cpm), each = nrow(cpm)),
  CPM = as.vector(cpm)
)

sample_cpm$sample <- group_meta$sample[
  match(sample_cpm$Group, groups)
]

sample_cpm$Disease <- group_meta$Disease[
  match(sample_cpm$Group, groups)
]

sample_cpm$Cell_type <- group_meta$HRCA_majorclass[
  match(sample_cpm$Group, groups)
]

write.csv(
  sample_cpm,
  file.path(outdir, "TE_family_sample_CPM.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Mean CPM per cell type and condition
# ------------------------------------------------------------

summary_dat <- sample_cpm %>%
  group_by(TE_family, Cell_type, Disease) %>%
  summarise(
    mean_CPM = mean(CPM),
    median_CPM = median(CPM),
    n_samples = n_distinct(sample),
    .groups = "drop"
  )

write.csv(
  summary_dat,
  file.path(outdir, "TE_family_group_summary.csv"),
  row.names = FALSE
)

summary_dat <- summary_dat %>%
  mutate(
    Group = paste(Cell_type, Disease, sep = "__")
  )

mat <- summary_dat %>%
  select(TE_family, Group, mean_CPM) %>%
  pivot_wider(
    names_from = Group,
    values_from = mean_CPM
  ) %>%
  column_to_rownames("TE_family") %>%
  as.matrix()

expected <- unlist(lapply(
  celltypes,
  function(ct) paste(ct, conditions, sep = "__")
))

stopifnot(all(expected %in% colnames(mat)))

mat <- mat[, expected, drop = FALSE]

# Missing groups are not assumed to be zero.
if (anyNA(mat)) {
  stop(
    "A cell type-condition group has no observed nuclei. ",
    "Cannot calculate its Z-score without imputation."
  )
}

# ------------------------------------------------------------
# 7. Log2 transformation and row Z-score
# ------------------------------------------------------------

mat_log <- log2(mat + 1)

mat_z <- t(scale(t(mat_log)))

# Constant rows have no relative variation.
constant_rows <- apply(mat_log, 1, sd) == 0

mat_z[constant_rows, ] <- 0

stopifnot(all(is.finite(mat_z)))

# ------------------------------------------------------------
# 8. Export matrices
# ------------------------------------------------------------

write.csv(
  mat,
  file.path(outdir, "TE_family_mean_CPM.csv")
)

write.csv(
  mat_log,
  file.path(outdir, "TE_family_log2CPM.csv")
)

write.csv(
  mat_z,
  file.path(outdir, "TE_family_Zscore.csv")
)

# ------------------------------------------------------------
# 9. Heatmap annotations
# ------------------------------------------------------------

annotation_col <- data.frame(
  Cell_type = factor(
    rep(celltypes, each = 3),
    levels = celltypes
  ),
  Disease = factor(
    rep(conditions, times = 9),
    levels = conditions
  ),
  row.names = colnames(mat_z)
)

celltype_colors <- setNames(
  grDevices::hcl.colors(9, "Set 2"),
  celltypes
)

disease_colors <- c(
  "Healthy Control" = "#607D8B",
  "Dry AMD" = "#E5A04B",
  "Wet AMD" = "#B85B69"
)

annotation_colors <- list(
  Cell_type = celltype_colors,
  Disease = disease_colors
)

# ------------------------------------------------------------
# 10. Heatmap
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

plot_args <- list(
  mat = mat_z,
  color = heat_colors,
  breaks = seq(-2.5, 2.5, length.out = 102),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_col = seq(3, 24, by = 3),
  labels_col = rep(
    c("Healthy", "Dry", "Wet"),
    times = 9
  ),
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 45,
  border_color = NA,
  cellwidth = 23,
  cellheight = 10,
  main = "TE family abundance across retinal cell populations",
  width = 14,
  height = 12
)

do.call(
  pheatmap,
  c(
    plot_args,
    list(
      filename = file.path(
        outdir,
        "Figure1B_TE_family_Zscore.pdf"
      )
    )
  )
)

do.call(
  pheatmap,
  c(
    plot_args,
    list(
      filename = file.path(
        outdir,
        "Figure1B_TE_family_Zscore.png"
      )
    )
  )
)

# ------------------------------------------------------------
# 11. Summary
# ------------------------------------------------------------

cat("\n=== FIGURE 1B COMPLETE ===\n")
cat("TE families:", nrow(mat_z), "\n")
cat("Nuclei:", ncol(te), "\n")
cat("Cell types:", length(celltypes), "\n")
cat("Conditions:", length(conditions), "\n")
cat("Heatmap columns:", ncol(mat_z), "\n")
cat("Sample-celltype groups:", length(groups), "\n")
cat("Constant TE families:", sum(constant_rows), "\n")
cat("Output:", outdir, "\n")
