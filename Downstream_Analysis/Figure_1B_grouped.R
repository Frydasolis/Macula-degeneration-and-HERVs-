
# ============================================================
# FIGURE 1B - TE FAMILY HEATMAP GROUPED BY LINEAGE
#
# 61 TE families
# 9 retinal cell types
# 3 disease conditions
#
# No filtering
# Uses previously calculated row-wise Z-scores
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

base <- file.path(
  "results",
  "SRP413248",
  "TE_family_subclustering",
  "Figure1B_FINAL"
)

outdir <- file.path(base, "Grouped")

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load original Z-score matrix
# ------------------------------------------------------------

mat <- as.matrix(
  read.csv(
    file.path(base, "TE_family_Zscore.csv"),
    row.names = 1,
    check.names = FALSE
  )
)

storage.mode(mat) <- "numeric"

stopifnot(
  nrow(mat) == 61,
  ncol(mat) == 27,
  !anyDuplicated(rownames(mat)),
  all(is.finite(mat))
)

# ------------------------------------------------------------
# 3. Define TE lineage groups
# ------------------------------------------------------------

classify_te <- function(x) {

  x <- toupper(x)

  # HERV-K and HML-related labels
  if (grepl("^HERVK|^HML", x)) {
    return("HERV-K / HML")
  }

  # HERV-H-related labels
  if (grepl("^HERVH", x)) {
    return("HERV-H")
  }

  # HERV-L-related labels
  if (grepl("^HERVL", x)) {
    return("HERV-L")
  }

  # Other explicitly HERV-labeled families
  if (grepl("^HERV", x)) {
    return("Other HERV")
  }

  # ERV-labeled families
  if (grepl("^ERV", x)) {
    return("Other ERV")
  }

  # LINE-associated labels
  if (grepl("^L1", x)) {
    return("LINE / L1")
  }

  # LTR-labeled families
  if (grepl("^LTR", x)) {
    return("LTR-labeled")
  }

  # MER-labeled families
  if (grepl("^MER", x)) {
    return("MER-labeled")
  }

  # Names requiring independent annotation
  return("Other / unclassified")
}

te_annotation <- data.frame(
  TE_family = rownames(mat),
  TE_group = vapply(
    rownames(mat),
    classify_te,
    character(1)
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 4. Define biological group order
# ------------------------------------------------------------

group_order <- c(
  "HERV-K / HML",
  "HERV-H",
  "HERV-L",
  "Other HERV",
  "Other ERV",
  "LINE / L1",
  "LTR-labeled",
  "MER-labeled",
  "Other / unclassified"
)

te_annotation$TE_group <- factor(
  te_annotation$TE_group,
  levels = group_order
)

# ------------------------------------------------------------
# 5. Order families within each group
# ------------------------------------------------------------

te_annotation <- te_annotation %>%
  arrange(TE_group, TE_family)

mat_ordered <- mat[
  te_annotation$TE_family,
  ,
  drop = FALSE
]

stopifnot(
  nrow(mat_ordered) == nrow(mat),
  setequal(rownames(mat_ordered), rownames(mat))
)

# ------------------------------------------------------------
# 6. Row annotations
# ------------------------------------------------------------

annotation_row <- data.frame(
  TE_group = te_annotation$TE_group,
  row.names = te_annotation$TE_family
)

group_colors <- c(
  "HERV-K / HML" = "#8E6C9E",
  "HERV-H" = "#D17A83",
  "HERV-L" = "#4C9F8A",
  "Other HERV" = "#5C87B2",
  "Other ERV" = "#D5A253",
  "LINE / L1" = "#A46F51",
  "LTR-labeled" = "#6BA3B0",
  "MER-labeled" = "#B18D57",
  "Other / unclassified" = "#999999"
)

# ------------------------------------------------------------
# 7. Column annotations
# ------------------------------------------------------------

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

expected <- unlist(lapply(
  celltypes,
  function(ct) paste(ct, conditions, sep = "__")
))

stopifnot(
  setequal(colnames(mat_ordered), expected)
)

mat_ordered <- mat_ordered[
  ,
  expected,
  drop = FALSE
]

annotation_col <- data.frame(
  Cell_type = factor(
    rep(celltypes, each = 3),
    levels = celltypes
  ),
  Disease = factor(
    rep(conditions, times = 9),
    levels = conditions
  ),
  row.names = colnames(mat_ordered)
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
  TE_group = group_colors,
  Cell_type = celltype_colors,
  Disease = disease_colors
)

# ------------------------------------------------------------
# 8. Gaps between TE groups
# ------------------------------------------------------------

group_sizes <- table(te_annotation$TE_group)

group_sizes <- group_sizes[group_sizes > 0]

gaps_row <- head(
  cumsum(as.integer(group_sizes)),
  -1
)

# ------------------------------------------------------------
# 9. Heatmap colors
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

# ------------------------------------------------------------
# 10. Generate heatmap
# ------------------------------------------------------------

plot_args <- list(
  mat = mat_ordered,

  color = heat_colors,
  breaks = seq(-2.5, 2.5, length.out = 102),

  cluster_rows = FALSE,
  cluster_cols = FALSE,

  annotation_row = annotation_row,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,

  gaps_row = gaps_row,
  gaps_col = seq(3, 24, by = 3),

  labels_col = rep(
    c("Healthy", "Dry", "Wet"),
    times = 9
  ),

  show_rownames = TRUE,
  show_colnames = TRUE,

  fontsize_row = 7,
  fontsize_col = 8,

  angle_col = 45,
  border_color = NA,

  cellwidth = 23,
  cellheight = 10,

  main = "TE family abundance across retinal cell populations",

  legend = TRUE,
  annotation_legend = TRUE,

  width = 15,
  height = 13
)

do.call(
  pheatmap,
  c(
    plot_args,
    list(
      filename = file.path(
        outdir,
        "Figure1B_TE_grouped.pdf"
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
        "Figure1B_TE_grouped.png"
      )
    )
  )
)

# ------------------------------------------------------------
# 11. Export annotation and ordered matrix
# ------------------------------------------------------------

write.csv(
  te_annotation,
  file.path(outdir, "TE_family_classification.csv"),
  row.names = FALSE
)

write.csv(
  mat_ordered,
  file.path(outdir, "TE_family_grouped_Zscore.csv")
)

# ------------------------------------------------------------
# 12. Summary
# ------------------------------------------------------------

cat("\n=== FIGURE 1B GROUPED COMPLETE ===\n")

cat("TE families:", nrow(mat_ordered), "\n")
cat("Heatmap columns:", ncol(mat_ordered), "\n")

cat("\n=== TE GROUPS ===\n")
print(table(te_annotation$TE_group))

cat("\nOutput:", outdir, "\n")
