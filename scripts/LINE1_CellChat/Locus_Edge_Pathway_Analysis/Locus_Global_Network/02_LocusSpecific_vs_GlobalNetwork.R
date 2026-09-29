############################################################
# 02_LocusSpecific_vs_GlobalNetwork.R
#
# PANEL C
#
# Candidate locus-specific L1FLnI expression
# versus GLOBAL donor-level CellChat network metrics
#
# Unit of inference = DONOR
#
# For each of the 12 candidate loci:
#
#   locus mean expression vs Network Density
#   locus mean expression vs Total Communication Strength
#
# Pearson = primary analysis
#
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

############################################################
# PATHS
############################################################

root <- "/storage/lemus_g/roldan/ARMD"

locus_dir <- file.path(
  root,
  "results/SRP413248/CellChat/donor_level_symbol",
  "donor_figure_analysis/locus_edge_correlations"
)

global_dir <- file.path(
  locus_dir,
  "paper_figures/Global_LINE1"
)

out_dir <- file.path(
  locus_dir,
  "paper_figures/Locus_GlobalNetwork"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# INPUT FILES
############################################################

locus_file <- file.path(
  locus_dir,
  "Donor_Global_Locus_Expression.csv"
)

network_file <- file.path(
  global_dir,
  "Global_L1FLnI_Donor_CellChat.csv"
)

############################################################
# LOAD
############################################################

cat("\nLoading donor-level locus expression...\n")

locus <- read.csv(
  locus_file,
  stringsAsFactors = FALSE
)

cat("Loading donor-level network metrics...\n")

network <- read.csv(
  network_file,
  stringsAsFactors = FALSE
)

############################################################
# CHECK COLUMNS
############################################################

required_locus <- c(
  "donor",
  "condition",
  "locus",
  "mean_expression"
)

required_network <- c(
  "donor",
  "density",
  "total_strength"
)

missing_locus <- setdiff(
  required_locus,
  names(locus)
)

missing_network <- setdiff(
  required_network,
  names(network)
)

if (length(missing_locus) > 0) {
  stop(
    paste(
      "Missing locus columns:",
      paste(missing_locus, collapse = ", ")
    )
  )
}

if (length(missing_network) > 0) {
  stop(
    paste(
      "Missing network columns:",
      paste(missing_network, collapse = ", ")
    )
  )
}

############################################################
# JOIN
############################################################

dat <- locus %>%
  select(
    donor,
    condition,
    locus,
    mean_expression
  ) %>%
  left_join(
    network %>%
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

cat("\nDonors:", n_distinct(dat$donor), "\n")
cat("Loci:", n_distinct(dat$locus), "\n")
cat("Rows:", nrow(dat), "\n")

############################################################
# SAFE PEARSON
############################################################

safe_pearson <- function(x, y) {

  keep <- complete.cases(x, y)

  x <- x[keep]
  y <- y[keep]

  n <- length(x)

  if (
    n < 3 ||
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(
      data.frame(
        n = n,
        r = NA_real_,
        p = NA_real_
      )
    )
  }

  test <- tryCatch(
    cor.test(
      x,
      y,
      method = "pearson"
    ),
    error = function(e) NULL
  )

  if (is.null(test)) {
    return(
      data.frame(
        n = n,
        r = NA_real_,
        p = NA_real_
      )
    )
  }

  data.frame(
    n = n,
    r = unname(test$estimate),
    p = test$p.value
  )
}

############################################################
# CALCULATE CORRELATIONS
############################################################

results_list <- list()

counter <- 1

for (loc in unique(dat$locus)) {

  df <- dat %>%
    filter(
      locus == loc
    )

  ##########################################################
  # DENSITY
  ##########################################################

  stat_density <- safe_pearson(
    df$mean_expression,
    df$density
  )

  results_list[[counter]] <- data.frame(
    locus = loc,
    metric = "Network Density",
    n = stat_density$n,
    r = stat_density$r,
    p = stat_density$p,
    stringsAsFactors = FALSE
  )

  counter <- counter + 1

  ##########################################################
  # TOTAL STRENGTH
  ##########################################################

  stat_strength <- safe_pearson(
    df$mean_expression,
    df$total_strength
  )

  results_list[[counter]] <- data.frame(
    locus = loc,
    metric = "Total Communication Strength",
    n = stat_strength$n,
    r = stat_strength$r,
    p = stat_strength$p,
    stringsAsFactors = FALSE
  )

  counter <- counter + 1
}

results <- bind_rows(
  results_list
)

############################################################
# MULTIPLE-TESTING CORRECTION
#
# BH correction separately for each network metric
# across the 12 candidate loci.
############################################################

results <- results %>%
  group_by(
    metric
  ) %>%
  mutate(
    FDR = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()

############################################################
# CLEAN LOCUS LABELS
############################################################

results <- results %>%
  mutate(

    locus_clean = sub(
      "^TE-L1FLnI-",
      "",
      locus
    ),

    significance = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      p < 0.05 ~ "P < 0.05",
      TRUE ~ "Not significant"
    ),

    significance = factor(
      significance,
      levels = c(
        "FDR < 0.05",
        "P < 0.05",
        "Not significant"
      )
    )
  )

############################################################
# SAVE STATISTICS
############################################################

write.csv(
  results,
  file.path(
    out_dir,
    "Locus_vs_GlobalNetwork_Pearson.csv"
  ),
  row.names = FALSE
)

############################################################
# ORDER LOCI
#
# Order using strongest absolute correlation observed
# across the two network metrics.
############################################################

locus_order <- results %>%
  group_by(
    locus_clean
  ) %>%
  summarise(
    max_abs_r = max(
      abs(r),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    max_abs_r
  ) %>%
  pull(
    locus_clean
  )

results$locus_clean <- factor(
  results$locus_clean,
  levels = locus_order
)

############################################################
# COLORS / SIGNIFICANCE
#
# Blue = not significant
# Red  = significant
############################################################

results <- results %>%
  mutate(

    significance_plot = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      p < 0.05 ~ "P < 0.05",
      TRUE ~ "Not significant"
    ),

    significance_plot = factor(
      significance_plot,
      levels = c(
        "Not significant",
        "P < 0.05",
        "FDR < 0.05"
      )
    )
  )

sig_colors <- c(
  "Not significant" = "#4C78A8",
  "P < 0.05" = "#D73027",
  "FDR < 0.05" = "#D73027"
)

############################################################
# PANEL C
############################################################

pC <- ggplot(
  results,
  aes(
    x = r,
    y = locus_clean
  )
) +

  ##########################################################
  # ZERO REFERENCE
  ##########################################################

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey60"
  ) +

  ##########################################################
  # LINE FROM ZERO TO PEARSON r
  ##########################################################

  geom_segment(
    aes(
      x = 0,
      xend = r,
      y = locus_clean,
      yend = locus_clean,
      color = significance_plot
    ),
    linewidth = 1.2,
    alpha = 0.75
  ) +

  ##########################################################
  # POINTS
  ##########################################################

  geom_point(
    aes(
      color = significance_plot,
      size = significance_plot
    ),
    alpha = 0.95
  ) +

  ##########################################################
  # PEARSON r LABEL
  ##########################################################

  geom_text(
    aes(
      label = sprintf("%.2f", r),
      color = significance_plot
    ),
    hjust = ifelse(
      results$r >= 0,
      -0.45,
      1.45
    ),
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  ) +

  ##########################################################
  # NETWORK METRICS
  ##########################################################

  facet_wrap(
    ~ metric,
    nrow = 1
  ) +

  ##########################################################
  # COLORS
  ##########################################################

  scale_color_manual(
    values = sig_colors,
    name = NULL
  ) +

  scale_size_manual(
    values = c(
      "Not significant" = 3.5,
      "P < 0.05" = 4.5,
      "FDR < 0.05" = 5.2
    ),
    name = NULL
  ) +

  ##########################################################
  # X AXIS
  ##########################################################

  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(
      -1,
      1,
      by = 0.25
    )
  ) +

  ##########################################################
  # LABELS
  ##########################################################

  labs(
    title = "C",
    subtitle = "Locus-specific L1FLnI expression vs global retinal communication",
    x = "Pearson correlation (r)",
    y = "L1FLnI locus"
  ) +

  ##########################################################
  # THEME
  ##########################################################

  theme_classic(
    base_size = 13
  ) +

  theme(

    axis.title = element_text(
      face = "bold",
      size = 12
    ),

    axis.text.x = element_text(
      size = 10
    ),

    axis.text.y = element_text(
      face = "bold",
      size = 10
    ),

    strip.text = element_text(
      face = "bold",
      size = 13
    ),

    strip.background = element_rect(
      fill = "grey96",
      color = NA
    ),

    plot.title = element_text(
      face = "bold",
      size = 16
    ),

    plot.subtitle = element_text(
      size = 11
    ),

    legend.position = "bottom",

    legend.text = element_text(
      size = 10
    ),

    panel.spacing = unit(
      2,
      "lines"
    )
  )

############################################################
# SAVE
############################################################

ggsave(
  file.path(
    out_dir,
    "Panel_C_Locus_vs_GlobalNetwork.png"
  ),
  pC,
  width = 11,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Panel_C_Locus_vs_GlobalNetwork.pdf"
  ),
  pC,
  width = 11,
  height = 7,
  bg = "white"
)

############################################################
# CONSOLE RESULTS
############################################################

cat("\n============================================\n")
cat("PANEL C COMPLETE\n")
cat("============================================\n")

cat("\nPearson correlations:\n")

print(
  results %>%
    arrange(
      metric,
      p
    ) %>%
    select(
      locus_clean,
      metric,
      n,
      r,
      p,
      FDR,
      significance
    ),
  n = Inf
)

cat("\nNominal P < 0.05:\n")

print(
  results %>%
    filter(
      p < 0.05
    ) %>%
    arrange(
      FDR,
      p
    ),
  n = Inf
)

cat("\nFDR < 0.05:\n")

print(
  results %>%
    filter(
      FDR < 0.05
    ) %>%
    arrange(
      FDR
    ),
  n = Inf
)

cat("\nOutput:\n")
cat(out_dir, "\n")
