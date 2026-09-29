############################################################
# 01_Global_LINE1_vs_CellChat.R
#
# FIGURE A-B
#
# Global L1FLnI expression vs donor-level CellChat metrics
#
# Unit of inference = DONOR
#
# A: Global L1FLnI expression vs network density
# B: Global L1FLnI expression vs total communication strength
#
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

############################################################
# PATHS
############################################################

root <- "/storage/lemus_g/roldan/ARMD"

seurat_file <- file.path(
  root,
  "results/SRP413248/Seurat/downstream",
  "SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

network_dir <- file.path(
  root,
  "results/SRP413248/CellChat/donor_level_symbol",
  "donor_figure_analysis"
)

out_dir <- file.path(
  network_dir,
  "locus_edge_correlations/paper_figures/Global_LINE1"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# LOAD SEURAT
############################################################

cat("\nLoading Seurat object...\n")

seu <- readRDS(seurat_file)

DefaultAssay(seu) <- "RNA"

expr <- GetAssayData(
  seu,
  assay = "RNA",
  layer = "data"
)

meta <- seu@meta.data
meta$cell <- rownames(meta)

############################################################
# IDENTIFY ALL L1FLnI FEATURES
############################################################

line1_features <- grep(
  "^TE-L1FLnI-",
  rownames(expr),
  value = TRUE
)

cat("\n============================================\n")
cat("GLOBAL L1FLnI FEATURES\n")
cat("============================================\n")

cat("Number of L1FLnI loci:", length(line1_features), "\n")

if (length(line1_features) == 0) {
  stop("No TE-L1FLnI features found.")
}

############################################################
# DONOR-LEVEL GLOBAL L1FLnI EXPRESSION
#
# mean_expression:
# mean normalized expression across:
#   all L1FLnI loci x all retinal cells in donor
#
# detection_breadth:
# fraction of L1FLnI loci detected at least once in donor
############################################################

donors <- unique(meta$orig.ident)

global_list <- list()

for (d in donors) {

  cells_d <- meta$cell[
    meta$orig.ident == d
  ]

  mat <- expr[
    line1_features,
    cells_d,
    drop = FALSE
  ]

  mean_expression <- mean(mat)

  locus_detected <- Matrix::rowSums(mat > 0) > 0

  detection_breadth <- mean(locus_detected)

  condition <- unique(
    meta$disease_state[
      meta$orig.ident == d
    ]
  )[1]

  global_list[[d]] <- data.frame(
    donor = d,
    condition = condition,
    n_cells = length(cells_d),
    global_L1FLnI_expression = mean_expression,
    detection_breadth = detection_breadth,
    stringsAsFactors = FALSE
  )
}

global_expr <- bind_rows(global_list)

############################################################
# LOAD DONOR CELLCHAT NETWORK METRICS
############################################################

candidate_files <- c(
  file.path(network_dir, "donor_edge_qc.csv"),
  file.path(network_dir, "condition_network_summary.csv"),
  file.path(network_dir, "donor_network_metrics.csv"),
  file.path(network_dir, "donor_metrics.csv")
)

metric_file <- candidate_files[file.exists(candidate_files)][1]

if (is.na(metric_file)) {

  cat("\nAvailable CSV files in network directory:\n")

  print(
    list.files(
      network_dir,
      pattern = "\\.csv$",
      full.names = TRUE
    )
  )

  stop(
    "Could not automatically identify donor-level network metrics file."
  )
}

cat("\nUsing network metrics file:\n")
cat(metric_file, "\n")

metrics <- read.csv(
  metric_file,
  stringsAsFactors = FALSE
)

cat("\nMetric columns:\n")
print(names(metrics))

############################################################
# STANDARDIZE DONOR COLUMN
############################################################

if (!"donor" %in% names(metrics)) {

  donor_candidates <- c(
    "orig.ident",
    "sample",
    "sample_id",
    "Donor",
    "GSM"
  )

  found <- donor_candidates[
    donor_candidates %in% names(metrics)
  ]

  if (length(found) == 0) {
    stop("Could not identify donor column.")
  }

  names(metrics)[
    names(metrics) == found[1]
  ] <- "donor"
}

############################################################
# VERIFY REQUIRED NETWORK METRICS
############################################################

required <- c(
  "donor",
  "density",
  "total_strength"
)

missing_cols <- setdiff(
  required,
  names(metrics)
)

if (length(missing_cols) > 0) {

  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

############################################################
# JOIN
############################################################

dat <- global_expr %>%
  left_join(
    metrics %>%
      select(
        donor,
        density,
        total_strength
      ),
    by = "donor"
  ) %>%
  filter(
    !is.na(density),
    !is.na(total_strength)
  )

cat("\n============================================\n")
cat("DONORS INCLUDED\n")
cat("============================================\n")

print(
  dat %>%
    select(
      donor,
      condition,
      global_L1FLnI_expression,
      detection_breadth,
      density,
      total_strength
    )
)

cat("\nn =", nrow(dat), "donors\n")

############################################################
# SAVE MASTER TABLE
############################################################

write.csv(
  dat,
  file.path(
    out_dir,
    "Global_L1FLnI_Donor_CellChat.csv"
  ),
  row.names = FALSE
)

############################################################
# PEARSON TESTS
############################################################

cor_density <- cor.test(
  dat$global_L1FLnI_expression,
  dat$density,
  method = "pearson"
)

cor_strength <- cor.test(
  dat$global_L1FLnI_expression,
  dat$total_strength,
  method = "pearson"
)

stats <- data.frame(
  outcome = c(
    "Network density",
    "Total communication strength"
  ),

  n = c(
    nrow(dat),
    nrow(dat)
  ),

  Pearson_r = c(
    unname(cor_density$estimate),
    unname(cor_strength$estimate)
  ),

  P_value = c(
    cor_density$p.value,
    cor_strength$p.value
  )
)

stats$R2 <- stats$Pearson_r^2

write.csv(
  stats,
  file.path(
    out_dir,
    "Global_L1FLnI_Pearson_Statistics.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("PEARSON RESULTS\n")
cat("============================================\n")

print(stats)

############################################################
# CONDITION COLORS / SHAPES
############################################################

condition_colors <- c(
  "Healthy Control" = "#4C78A8",
  "Dry AMD" = "#F2A51A",
  "Wet AMD" = "#E45756"
)

condition_shapes <- c(
  "Healthy Control" = 16,
  "Dry AMD" = 17,
  "Wet AMD" = 15
)

############################################################
# STAT LABEL FUNCTION
############################################################

stat_label <- function(test, n) {

  r <- unname(test$estimate)
  p <- test$p.value

  p_txt <- ifelse(
    p < 0.001,
    format(
      p,
      scientific = TRUE,
      digits = 2
    ),
    sprintf("%.3f", p)
  )

  paste0(
    "Pearson r = ",
    sprintf("%.2f", r),
    "\nP = ",
    p_txt,
    "\nn = ",
    n,
    " donors"
  )
}

############################################################
# COMMON THEME
############################################################

paper_theme <- theme_classic(
  base_size = 13
) +
  theme(
    axis.title = element_text(
      face = "bold",
      size = 13
    ),

    axis.text = element_text(
      size = 11
    ),

    plot.title = element_text(
      face = "bold",
      size = 14
    ),

    plot.subtitle = element_text(
      size = 10
    ),

    legend.title = element_blank(),

    legend.position = "bottom",

    legend.text = element_text(
      size = 10
    ),

    plot.margin = margin(
      12, 15, 12, 12
    )
  )

############################################################
# PANEL A — NETWORK DENSITY
############################################################

pA <- ggplot(
  dat,
  aes(
    x = global_L1FLnI_expression,
    y = density
  )
) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey30",
    fill = "grey80",
    linewidth = 0.8,
    alpha = 0.30
  ) +

  geom_point(
    aes(
      color = condition,
      shape = condition
    ),
    size = 4,
    stroke = 1
  ) +

  scale_color_manual(
    values = condition_colors
  ) +

  scale_shape_manual(
    values = condition_shapes
  ) +

  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = stat_label(
      cor_density,
      nrow(dat)
    ),
    hjust = -0.05,
    vjust = 1.15,
    size = 4
  ) +

  labs(
    title = "A",
    x = "Global L1FLnI mean normalized expression",
    y = "Network density"
  ) +

  paper_theme

############################################################
# PANEL B — TOTAL COMMUNICATION STRENGTH
############################################################

pB <- ggplot(
  dat,
  aes(
    x = global_L1FLnI_expression,
    y = total_strength
  )
) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "grey30",
    fill = "grey80",
    linewidth = 0.8,
    alpha = 0.30
  ) +

  geom_point(
    aes(
      color = condition,
      shape = condition
    ),
    size = 4,
    stroke = 1
  ) +

  scale_color_manual(
    values = condition_colors
  ) +

  scale_shape_manual(
    values = condition_shapes
  ) +

  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = stat_label(
      cor_strength,
      nrow(dat)
    ),
    hjust = -0.05,
    vjust = 1.15,
    size = 4
  ) +

  labs(
    title = "B",
    x = "Global L1FLnI mean normalized expression",
    y = "Total communication strength"
  ) +

  paper_theme

############################################################
# COMBINE A + B
############################################################

combined <- (
  pA +
  pB +
  plot_layout(
    guides = "collect"
  )
) &
  theme(
    legend.position = "bottom"
  )

############################################################
# SAVE INDIVIDUAL PANELS
############################################################

ggsave(
  file.path(
    out_dir,
    "Panel_A_Global_L1FLnI_vs_Density.png"
  ),
  pA,
  width = 6.5,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Panel_A_Global_L1FLnI_vs_Density.pdf"
  ),
  pA,
  width = 6.5,
  height = 5.5,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Panel_B_Global_L1FLnI_vs_Strength.png"
  ),
  pB,
  width = 6.5,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Panel_B_Global_L1FLnI_vs_Strength.pdf"
  ),
  pB,
  width = 6.5,
  height = 5.5,
  bg = "white"
)

############################################################
# SAVE COMBINED FIGURE
############################################################

ggsave(
  file.path(
    out_dir,
    "Figure_Global_L1FLnI_vs_CellChat.png"
  ),
  combined,
  width = 12,
  height = 5.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure_Global_L1FLnI_vs_CellChat.pdf"
  ),
  combined,
  width = 12,
  height = 5.5,
  bg = "white"
)

############################################################
# DONE
############################################################

cat("\n============================================\n")
cat("GLOBAL L1FLnI FIGURE COMPLETE\n")
cat("============================================\n")

cat("\nFigures saved in:\n")
cat(out_dir, "\n")

cat("\nFiles:\n")
print(
  list.files(
    out_dir,
    pattern = "\\.(png|pdf|csv)$"
  )
)
