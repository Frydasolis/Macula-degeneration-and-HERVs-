
# ============================================================
# FIGURE 3 - GLOBAL TE BURDEN STATISTICS
# ARMD / SRP413248
#
# Calculates cell-type-specific TE burden statistics:
#   - Healthy Control vs Dry AMD
#   - Healthy Control vs Wet AMD
#   - Dry AMD vs Wet AMD
#
# Statistical test:
#   Wilcoxon rank-sum test
#
# Multiple testing:
#   Benjamini-Hochberg FDR
#
# Figures:
#   Violin + boxplot + significance brackets
#
# NOTE:
# These are cell-level exploratory statistics.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({

  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)

})


# ============================================================
# 2. PATHS
# ============================================================

project_dir <- "/storage/lemus_g/roldan/ARMD"

seurat_file <- file.path(
  project_dir,
  "results/SRP413248/Seurat/downstream/",
  "SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

figure3_dir <- file.path(
  project_dir,
  "results/SRP413248/Seurat/downstream/TE_analysis/Figure3"
)

dir.create(
  figure3_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. LOAD SEURAT OBJECT
# ============================================================

cat("\n============================================\n")
cat("Loading Seurat object\n")
cat("============================================\n")

seu <- readRDS(seurat_file)

cat("Cells:", ncol(seu), "\n")
cat("Features:", nrow(seu), "\n")


# ============================================================
# 4. FIND CONDITION COLUMN
# ============================================================

possible_condition_columns <- c(
  "condition",
  "Condition",
  "disease",
  "Disease",
  "diagnosis",
  "Diagnosis"
)

condition_col <- possible_condition_columns[
  possible_condition_columns %in% colnames(seu@meta.data)
][1]

if (is.na(condition_col)) {

  stop(
    "Could not find a disease/condition metadata column."
  )

}

cat("\nCondition column:", condition_col, "\n")


# ============================================================
# 5. CREATE WORKING DATAFRAME
# ============================================================

meta <- seu@meta.data

meta$condition_raw <- as.character(
  meta[[condition_col]]
)

meta$HRCA_majorclass <- as.character(
  meta$HRCA_majorclass
)

meta$percent_TE <- as.numeric(
  meta$percent_TE
)


# ============================================================
# 6. STANDARDIZE DISEASE LABELS
# ============================================================

meta$condition <- case_when(

  meta$condition_raw %in% c(
    "Healthy",
    "Healthy Control",
    "Control",
    "HC",
    "healthy",
    "healthy control"
  ) ~ "Healthy Control",

  meta$condition_raw %in% c(
    "Dry",
    "Dry AMD",
    "dry AMD",
    "Dry_AMD"
  ) ~ "Dry AMD",

  meta$condition_raw %in% c(
    "Wet",
    "Wet AMD",
    "wet AMD",
    "Wet_AMD"
  ) ~ "Wet AMD",

  TRUE ~ meta$condition_raw

)


# ============================================================
# 7. FACTOR ORDER
# ============================================================

condition_levels <- c(
  "Healthy Control",
  "Dry AMD",
  "Wet AMD"
)

meta$condition <- factor(
  meta$condition,
  levels = condition_levels
)

meta$HRCA_majorclass <- factor(
  meta$HRCA_majorclass
)


# ============================================================
# 8. CHECK CONDITIONS
# ============================================================

cat("\n============================================\n")
cat("Conditions detected\n")
cat("============================================\n")

print(
  table(
    meta$condition,
    useNA = "ifany"
  )
)

cat("\nUnique conditions:\n")
print(
  unique(
    as.character(meta$condition)
  )
)


# ============================================================
# 9. PLOT DATA
# ============================================================

plot_df <- meta %>%

  filter(
    !is.na(HRCA_majorclass),
    !is.na(condition),
    !is.na(percent_TE),
    is.finite(percent_TE)
  ) %>%

  filter(
    condition %in% condition_levels
  )


# ============================================================
# 10. SUMMARY STATISTICS
# ============================================================

cat("\n============================================\n")
cat("Summary statistics\n")
cat("============================================\n")

summary_stats <- plot_df %>%

  group_by(
    HRCA_majorclass,
    condition
  ) %>%

  summarise(

    n = n(),

    mean_percent_TE =
      mean(
        percent_TE,
        na.rm = TRUE
      ),

    median_percent_TE =
      median(
        percent_TE,
        na.rm = TRUE
      ),

    sd_percent_TE =
      sd(
        percent_TE,
        na.rm = TRUE
      ),

    IQR_percent_TE =
      IQR(
        percent_TE,
        na.rm = TRUE
      ),

    .groups = "drop"

  )

print(summary_stats)


# ============================================================
# 11. SAVE SUMMARY
# ============================================================

write_csv(
  summary_stats,
  file.path(
    figure3_dir,
    "TE_burden_summary_by_celltype_condition.csv"
  )
)


# ============================================================
# 12. PAIRWISE COMPARISONS
# ============================================================

comparisons <- list(

  c(
    "Healthy Control",
    "Dry AMD"
  ),

  c(
    "Healthy Control",
    "Wet AMD"
  ),

  c(
    "Dry AMD",
    "Wet AMD"
  )

)


# ============================================================
# 13. FUNCTION FOR SIGNIFICANCE LABEL
# ============================================================

get_significance <- function(fdr) {

  if (is.na(fdr)) {

    return("NA")

  }

  if (fdr < 0.0001) {

    return("****")

  }

  if (fdr < 0.001) {

    return("***")

  }

  if (fdr < 0.01) {

    return("**")

  }

  if (fdr < 0.05) {

    return("*")

  }

  return("ns")

}


# ============================================================
# 14. CALCULATE WILCOXON TESTS
# ============================================================

cat("\n============================================\n")
cat("PAIRWISE STATISTICS\n")
cat("============================================\n")


stats_list <- list()

counter <- 1


for (ct in levels(plot_df$HRCA_majorclass)) {

  dat_ct <- plot_df %>%

    filter(
      HRCA_majorclass == ct
    )


  for (cmp in comparisons) {

    group1 <- cmp[1]
    group2 <- cmp[2]


    dat1 <- dat_ct %>%

      filter(
        condition == group1
      ) %>%

      pull(percent_TE)


    dat2 <- dat_ct %>%

      filter(
        condition == group2
      ) %>%

      pull(percent_TE)


    # Skip if either group is absent
    if (
      length(dat1) == 0 ||
      length(dat2) == 0
    ) {

      next

    }


    # Wilcoxon rank-sum test
    wt <- wilcox.test(
      dat1,
      dat2,
      exact = FALSE
    )


    stats_list[[counter]] <- tibble(

      cell_type = ct,

      comparison =
        paste(
          group2,
          "vs",
          group1
        ),

      group1 = group1,

      group2 = group2,

      n_group1 = length(dat1),

      n_group2 = length(dat2),

      median_group1 =
        median(
          dat1,
          na.rm = TRUE
        ),

      median_group2 =
        median(
          dat2,
          na.rm = TRUE
        ),

      mean_group1 =
        mean(
          dat1,
          na.rm = TRUE
        ),

      mean_group2 =
        mean(
          dat2,
          na.rm = TRUE
        ),

      median_difference =
        median(
          dat2,
          na.rm = TRUE
        ) -
        median(
          dat1,
          na.rm = TRUE
        ),

      mean_difference =
        mean(
          dat2,
          na.rm = TRUE
        ) -
        mean(
          dat1,
          na.rm = TRUE
        ),

      W = as.numeric(
        wt$statistic
      ),

      p_value = wt$p.value

    )


    counter <- counter + 1

  }

}


stats <- bind_rows(
  stats_list
)


# ============================================================
# 15. FDR CORRECTION
# ============================================================

stats <- stats %>%

  mutate(

    FDR = p.adjust(
      p_value,
      method = "BH"
    ),

    significance =
      vapply(
        FDR,
        get_significance,
        character(1)
      )

  )


# ============================================================
# 16. PRINT STATISTICS
# ============================================================

print(stats)


# ============================================================
# 17. SAVE ALL PAIRWISE RESULTS
# ============================================================

stats_file <- file.path(
  figure3_dir,
  "TE_burden_pairwise_Wilcoxon_by_celltype.csv"
)

write_csv(
  stats,
  stats_file
)

cat(
  "\nSaved:\n",
  stats_file,
  "\n"
)


# ============================================================
# 18. SAVE SIGNIFICANT RESULTS
# ============================================================

stats_significant <- stats %>%

  filter(
    FDR < 0.05
  )

write_csv(

  stats_significant,

  file.path(
    figure3_dir,
    "TE_burden_pairwise_Wilcoxon_significant.csv"
  )

)


# ============================================================
# 19. COLORS
# ============================================================

condition_colors <- c(

  "Healthy Control" = "#6BAED6",

  "Dry AMD" = "#F4A261",

  "Wet AMD" = "#D95F59"

)


# ============================================================
# 20. FUNCTION TO CREATE PLOT
# ============================================================

make_plot <- function(ct) {


  # ----------------------------------------------------------
  # Data for this cell type
  # ----------------------------------------------------------

  dat <- plot_df %>%

    filter(
      HRCA_majorclass == ct
    )


  stat_ct <- stats %>%

    filter(
      cell_type == ct
    )


  # ----------------------------------------------------------
  # Maximum value
  # ----------------------------------------------------------

  ymax <- max(
    dat$percent_TE,
    na.rm = TRUE
  )


  # ----------------------------------------------------------
  # Spacing
  # ----------------------------------------------------------

  step <- max(
    ymax * 0.12,
    0.5
  )


  # ----------------------------------------------------------
  # Annotation positions
  # ----------------------------------------------------------

  annotation_df <- stat_ct %>%

    mutate(

      xmin =
        match(
          group1,
          condition_levels
        ),

      xmax =
        match(
          group2,
          condition_levels
        ),

      y.position =
        ymax +
        seq_len(n()) * step

    )


  # ----------------------------------------------------------
  # Base violin plot
  # ----------------------------------------------------------

  p <- ggplot(

    dat,

    aes(

      x = condition,

      y = percent_TE,

      fill = condition

    )

  ) +

    geom_violin(

      trim = FALSE,

      alpha = 0.65,

      color = "black"

    ) +

    geom_boxplot(

      width = 0.15,

      outlier.shape = NA,

      fill = "white",

      color = "black"

    ) +

    scale_fill_manual(

      values = condition_colors

    ) +

    labs(

      title = ct,

      x = NULL,

      y = "TE burden (%)"

    ) +

    theme_classic(

      base_size = 12

    ) +

    theme(

      plot.title =
        element_text(

          face = "bold",

          hjust = 0.5

        ),

      legend.position = "none"

    )


  # ----------------------------------------------------------
  # ADD STATISTICAL BRACKETS
  # ----------------------------------------------------------

  if (
    nrow(annotation_df) > 0
  ) {


    bracket_height <-
      step * 0.25


    text_height <-
      step * 0.10


    # Horizontal bracket
    p <- p +

      geom_segment(

        data = annotation_df,

        aes(

          x = xmin,

          xend = xmax,

          y = y.position,

          yend = y.position

        ),

        inherit.aes = FALSE,

        linewidth = 0.4

      )


    # Left vertical bracket
    p <- p +

      geom_segment(

        data = annotation_df,

        aes(

          x = xmin,

          xend = xmin,

          y = y.position,

          yend =
            y.position -
            bracket_height

        ),

        inherit.aes = FALSE,

        linewidth = 0.4

      )


    # Right vertical bracket
    p <- p +

      geom_segment(

        data = annotation_df,

        aes(

          x = xmax,

          xend = xmax,

          y = y.position,

          yend =
            y.position -
            bracket_height

        ),

        inherit.aes = FALSE,

        linewidth = 0.4

      )


    # Significance label
    p <- p +

      geom_text(

        data = annotation_df,

        aes(

          x =
            (xmin + xmax) / 2,

          y =
            y.position +
            text_height,

          label = significance

        ),

        inherit.aes = FALSE,

        size = 4

      )


    # Expand y-axis
    p <- p +

      scale_y_continuous(

        expand =
          expansion(

            mult =
              c(
                0.05,
                0.25
              )

          )

      )

  }


  return(p)

}


# ============================================================
# 21. GENERATE INDIVIDUAL PLOTS
# ============================================================

cat("\n============================================\n")
cat("Generating individual plots\n")
cat("============================================\n")


plots <- list()


for (
  ct in levels(
    plot_df$HRCA_majorclass
  )
) {


  if (
    is.na(ct)
  ) {

    next

  }


  p <- make_plot(
    ct
  )


  plots[[ct]] <- p


  safe_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    ct
  )


  # ----------------------------------------------------------
  # PDF
  # ----------------------------------------------------------

  ggsave(

    file.path(

      figure3_dir,

      paste0(

        "TE_burden_violin_",

        safe_name,

        ".pdf"

      )

    ),

    p,

    width = 5.5,

    height = 5.5

  )


  # ----------------------------------------------------------
  # PNG
  # ----------------------------------------------------------

  ggsave(

    file.path(

      figure3_dir,

      paste0(

        "TE_burden_violin_",

        safe_name,

        ".png"

      )

    ),

    p,

    width = 5.5,

    height = 5.5,

    dpi = 300

  )


  cat(
    "Saved:",
    ct,
    "\n"
  )

}


# ============================================================
# 22. COMBINED FIGURE
# ============================================================

if (
  requireNamespace(
    "patchwork",
    quietly = TRUE
  )
) {


  library(
    patchwork
  )


  combined <- wrap_plots(

    plots,

    ncol = 3

  )


  # ----------------------------------------------------------
  # Combined PDF
  # ----------------------------------------------------------

  ggsave(

    file.path(

      figure3_dir,

      "Figure3_TE_burden_by_celltype_with_statistics.pdf"

    ),

    combined,

    width = 15,

    height = 15

  )


  # ----------------------------------------------------------
  # Combined PNG
  # ----------------------------------------------------------

  ggsave(

    file.path(

      figure3_dir,

      "Figure3_TE_burden_by_celltype_with_statistics.png"

    ),

    combined,

    width = 15,

    height = 15,

    dpi = 300

  )


  cat(
    "\nCombined figure saved.\n"
  )


} else {


  cat(
    "\nNOTE: patchwork is not installed.\n"
  )

  cat(
    "Individual plots were still generated.\n"
  )

}


# ============================================================
# 23. SESSION INFO
# ============================================================

writeLines(

  capture.output(
    sessionInfo()
  ),

  file.path(

    figure3_dir,

    "sessionInfo_TE_burden_statistics.txt"

  )

)


# ============================================================
# 24. FINAL MESSAGE
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "DONE\n"
)

cat(
  "============================================\n"
)

cat(
  "\nResults directory:\n"
)

cat(
  figure3_dir,
  "\n\n"
)

