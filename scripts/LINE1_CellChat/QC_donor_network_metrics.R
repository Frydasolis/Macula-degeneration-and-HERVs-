############################################################
# DONOR-LEVEL CELLCHAT NETWORK QC
#
# Input:
#   network_metrics/donor_network_metrics.csv
#
# Outputs:
#   network_metrics/QC/
#
# Figure:
#   A. Network metrics by condition
#   B. Correlation between network metrics
#   C. n_cells vs total_strength
#   D. n_cells vs density
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(patchwork)
})

############################################################
# PATHS
############################################################

input_file <- "network_metrics/donor_network_metrics.csv"

outdir <- "network_metrics/QC"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# LOAD DATA
############################################################

df <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)

df$Condition <- factor(
  df$Condition,
  levels = c(
    "Healthy Control",
    "Dry AMD",
    "Wet AMD"
  )
)

cat("\n============================================\n")
cat("DONOR NETWORK QC\n")
cat("============================================\n")
cat("Donors:", nrow(df), "\n\n")

print(df)

############################################################
# VARIABLES OF INTEREST
############################################################

network_vars <- c(
  "total_LR_interactions",
  "active_edges",
  "density",
  "total_strength",
  "mean_edge_strength",
  "mean_LR_per_edge",
  "n_pathways"
)

############################################################
# PANEL A
# DISTRIBUTION OF NETWORK METRICS BY CONDITION
############################################################

plot_df <- df %>%
  select(
    Donor,
    Condition,
    density,
    total_strength,
    mean_edge_strength
  ) %>%
  pivot_longer(
    cols = c(
      density,
      total_strength,
      mean_edge_strength
    ),
    names_to = "Metric",
    values_to = "Value"
  )

plot_df$Metric <- factor(
  plot_df$Metric,
  levels = c(
    "density",
    "total_strength",
    "mean_edge_strength"
  ),
  labels = c(
    "Network density",
    "Total communication strength",
    "Mean edge strength"
  )
)

pA <- ggplot(
  plot_df,
  aes(
    x = Condition,
    y = Value,
    fill = Condition
  )
) +

  geom_boxplot(
    width = 0.55,
    alpha = 0.25,
    outlier.shape = NA
  ) +

  geom_jitter(
    aes(color = Condition),
    width = 0.12,
    size = 2.8,
    alpha = 0.9
  ) +

  facet_wrap(
    ~ Metric,
    scales = "free_y",
    nrow = 1, ncol = 3,
  ) +

  labs(
    title = "A",
    x = NULL,
    y = NULL
  ) +

  theme_classic(base_size = 12) +

  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      size = 16
    )
  )

############################################################
# PANEL B
# SPEARMAN CORRELATION BETWEEN NETWORK METRICS
############################################################

cor_mat <- cor(
  df[, network_vars],
  method = "spearman",
  use = "pairwise.complete.obs"
)

write.csv(
  cor_mat,
  file.path(
    outdir,
    "network_metric_Spearman_correlations.csv"
  )
)

cor_df <- as.data.frame(
  as.table(cor_mat)
)

colnames(cor_df) <- c(
  "Metric1",
  "Metric2",
  "rho"
)

pretty_names <- c(
  total_LR_interactions = "LR interactions",
  active_edges = "Active edges",
  density = "Density",
  total_strength = "Total strength",
  mean_edge_strength = "Mean edge strength",
  mean_LR_per_edge = "LR per edge",
  n_pathways = "Pathways"
)

cor_df$Metric1 <- pretty_names[
  as.character(cor_df$Metric1)
]

cor_df$Metric2 <- pretty_names[
  as.character(cor_df$Metric2)
]

metric_order <- unname(pretty_names)

cor_df$Metric1 <- factor(
  cor_df$Metric1,
  levels = metric_order
)

cor_df$Metric2 <- factor(
  cor_df$Metric2,
  levels = rev(metric_order)
)

pB <- ggplot(
  cor_df,
  aes(
    x = Metric1,
    y = Metric2,
    fill = rho
  )
) +

  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +

  geom_text(
    aes(
      label = sprintf("%.2f", rho)
    ),
    size = 3.2
  ) +

  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-1, 1)
  ) +

  coord_fixed() +

  labs(
    title = "B",
    x = NULL,
    y = NULL,
    fill = "Spearman\nrho"
  ) +

  theme_minimal(base_size = 10) +

  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold",
      size = 16
    )
  )

############################################################
# FUNCTION:
# SPEARMAN TEST + LABEL
############################################################

make_spearman_label <- function(x, y) {

  test <- cor.test(
    x,
    y,
    method = "spearman",
    exact = FALSE
  )

  rho <- unname(test$estimate)
  p <- test$p.value

  p_text <- ifelse(
    p < 0.001,
    format(
      p,
      scientific = TRUE,
      digits = 2
    ),
    sprintf("%.3f", p)
  )

  paste0(
    "Spearman ρ = ",
    sprintf("%.2f", rho),
    "\nP = ",
    p_text
  )
}

############################################################
# PANEL C
# NUMBER OF CELLS VS TOTAL STRENGTH
############################################################

label_C <- make_spearman_label(
  df$n_cells,
  df$total_strength
)

pC <- ggplot(
  df,
  aes(
    x = n_cells,
    y = total_strength,
    color = Condition
  )
) +

  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +

  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "grey35",
    fill = "grey80",
    linewidth = 0.8
  ) +

  geom_text_repel(
    aes(label = Donor),
    size = 2.7,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    show.legend = FALSE
  ) +

  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_C,
    hjust = 1.1,
    vjust = 1.3,
    size = 3.6
  ) +

  labs(
    title = "C",
    x = "Number of cells",
    y = "Total communication strength"
  ) +

  theme_classic(base_size = 12) +

  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(
      face = "bold",
      size = 16
    )
  )

############################################################
# PANEL D
# NUMBER OF CELLS VS NETWORK DENSITY
############################################################

label_D <- make_spearman_label(
  df$n_cells,
  df$density
)

pD <- ggplot(
  df,
  aes(
    x = n_cells,
    y = density,
    color = Condition
  )
) +

  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +

  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "grey35",
    fill = "grey80",
    linewidth = 0.8
  ) +

  geom_text_repel(
    aes(label = Donor),
    size = 2.7,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    show.legend = FALSE
  ) +

  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_D,
    hjust = 1.1,
    vjust = 1.3,
    size = 3.6
  ) +

  labs(
    title = "D",
    x = "Number of cells",
    y = "Network density"
  ) +

  theme_classic(base_size = 12) +

  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 16
    )
  )

############################################################
# COMBINE FIGURE
############################################################

final_plot <- wrap_plots(
  pA,
  pB,
  pC,
  pD,
  design = '
  AAAAAA
  BBBBBB
  CCCDDD
  ',
  heights = c(1, 1.3, 1)
)

############################################################
# SAVE
############################################################

ggsave(
  file.path(
    outdir,
    "Donor_CellChat_network_QC.png"
  ),
  final_plot,
  width = 14,
  height = 14,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(
    outdir,
    "Donor_CellChat_network_QC.pdf"
  ),
  final_plot,
  width = 14,
  height = 14
)

############################################################
# ADDITIONAL QC:
# CORRELATION OF CELL NUMBER WITH ALL NETWORK METRICS
############################################################

qc_results <- lapply(
  network_vars,
  function(v) {

    test <- cor.test(
      df$n_cells,
      df[[v]],
      method = "spearman",
      exact = FALSE
    )

    data.frame(
      Metric = v,
      rho = unname(test$estimate),
      p_value = test$p.value
    )
  }
)

qc_results <- do.call(
  rbind,
  qc_results
)

qc_results$FDR <- p.adjust(
  qc_results$p_value,
  method = "BH"
)

write.csv(
  qc_results,
  file.path(
    outdir,
    "cell_number_vs_network_metrics.csv"
  ),
  row.names = FALSE
)

############################################################
# PRINT RESULTS
############################################################

cat("\n============================================\n")
cat("CELL NUMBER vs NETWORK METRICS\n")
cat("============================================\n\n")

print(
  qc_results,
  row.names = FALSE
)

cat("\n============================================\n")
cat("SPEARMAN CORRELATION MATRIX\n")
cat("============================================\n\n")

print(
  round(
    cor_mat,
    3
  )
)

cat("\n============================================\n")
cat("FILES CREATED\n")
cat("============================================\n")

cat(
  file.path(
    outdir,
    "Donor_CellChat_network_QC.png"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "Donor_CellChat_network_QC.pdf"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "network_metric_Spearman_correlations.csv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "cell_number_vs_network_metrics.csv"
  ),
  "\n"
)

cat("\nDONE\n")
