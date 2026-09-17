#!/usr/bin/env Rscript

# ============================================================
# CLEAN LINE-1 vs ALL CELLCHAT GLOBAL METRICS
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# ------------------------------------------------------------
# PATHS
# ------------------------------------------------------------

project_dir <- "/storage/lemus_g/roldan/ARMD"

seurat_file <- file.path(
  project_dir,
  "results/SRP413248/Seurat/downstream",
  "SRP413248_GeneTE_QC_UMAP_HRCA_annotated_geneSymbols.rds"
)

metrics_file <- file.path(
  project_dir,
  "results/SRP413248/CellChat/network_metrics/global_all_metrics",
  "CellChat_global_all_network_metrics.csv"
)

out_dir <- file.path(
  project_dir,
  "results/SRP413248/CellChat/network_metrics/global_all_metrics/LINE1_correlations"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# LOAD SEURAT
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("LOADING SEURAT OBJECT\n")
cat("====================================================\n")

obj <- readRDS(seurat_file)

cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")

# ------------------------------------------------------------
# IDENTIFY LINE-1 FEATURES
# ------------------------------------------------------------

features <- rownames(obj)

L1_features <- grep("^TE-L1", features, value = TRUE)

cat("\nLINE-1 features:", length(L1_features), "\n")

if (length(L1_features) == 0) {
  stop("No TE-L1 features found.")
}

# ------------------------------------------------------------
# EXTRACT NORMALIZED LINE-1 EXPRESSION
# ------------------------------------------------------------

cat("\nCalculating LINE-1 expression per cell...\n")

L1_expr <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)[L1_features, , drop = FALSE]

# Sum all LINE-1 features per cell
L1_per_cell <- Matrix::colSums(L1_expr)

L1_cell_df <- data.frame(
  Cell = names(L1_per_cell),
  LINE1_expression = as.numeric(L1_per_cell),
  GSM = obj$GSM[names(L1_per_cell)]
)

# ------------------------------------------------------------
# AGGREGATE LINE-1 PER DONOR
# ------------------------------------------------------------

L1_donor <- L1_cell_df %>%
  group_by(GSM) %>%
  summarise(
    LINE1_expression = mean(LINE1_expression, na.rm = TRUE),
    Cells_LINE1 = n(),
    .groups = "drop"
  )

cat("\nLINE-1 donor measurements:", nrow(L1_donor), "\n")

# ------------------------------------------------------------
# LOAD GLOBAL CELLCHAT METRICS
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("LOADING CELLCHAT GLOBAL METRICS\n")
cat("====================================================\n")

metrics <- read_csv(
  metrics_file,
  show_col_types = FALSE
)

cat("Rows:", nrow(metrics), "\n")
cat("Columns:", ncol(metrics), "\n")

# ------------------------------------------------------------
# CLEAN DONOR ID
# ------------------------------------------------------------

metrics <- metrics %>%
  mutate(
    GSM_clean = gsub("\\.rds$", "", donor)
  )

# ------------------------------------------------------------
# ADD CONDITION EXPLICITLY
# ------------------------------------------------------------

metrics <- metrics %>%
  mutate(
    condition = case_when(
      GSM_clean %in% c(
        "GSM6841143",
        "GSM6841144",
        "GSM6841145",
        "GSM6841146",
        "GSM6841147",
        "GSM6841148"
      ) ~ "Healthy",

      GSM_clean %in% c(
        "GSM6841149",
        "GSM6841150",
        "GSM6841151",
        "GSM6841152",
        "GSM6841153",
        "GSM6841154",
        "GSM6841155"
      ) ~ "Wet",

      GSM_clean %in% c(
        "GSM6841156",
        "GSM6841157",
        "GSM6841159"
      ) ~ "Dry",

      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# MERGE
# ------------------------------------------------------------

merged <- metrics %>%
  left_join(
    L1_donor,
    by = c("GSM_clean" = "GSM")
  )

# ------------------------------------------------------------
# VALIDATION
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("VALIDATION\n")
cat("====================================================\n")

cat("Merged donors:", nrow(merged), "\n")

cat("\nCondition counts:\n")
print(table(merged$condition, useNA = "ifany"))

if (nrow(merged) != 16) {
  stop("Expected 16 donors after merge.")
}

if (any(is.na(merged$LINE1_expression))) {
  stop("Missing LINE-1 expression for some donors.")
}

if (any(is.na(merged$condition))) {
  stop("Missing condition for some donors.")
}

# ------------------------------------------------------------
# DEFINE METRICS TO TEST
# ------------------------------------------------------------

# Remove:
# - donor identifiers
# - condition
# - duplicate cell count
# - constant signaling_genes
# - LINE-1 itself
# - Cells_LINE1

exclude_metrics <- c(
  "donor",
  "GSM_clean",
  "condition",
  "LINE1_expression",
  "Cells_LINE1",
  "signaling_genes",
  "cells",
  "Cells"
)

numeric_metrics <- names(merged)[
  sapply(merged, is.numeric)
]

metrics_to_test <- setdiff(
  numeric_metrics,
  exclude_metrics
)

cat("\n====================================================\n")
cat("METRICS INCLUDED IN FINAL ANALYSIS\n")
cat("====================================================\n")

cat("Number of metrics:", length(metrics_to_test), "\n\n")

print(metrics_to_test)

# ------------------------------------------------------------
# REMOVE CONSTANT METRICS
# ------------------------------------------------------------

metrics_to_test <- metrics_to_test[
  sapply(
    merged[metrics_to_test],
    function(x) {
      length(unique(x[!is.na(x)])) > 1
    }
  )
]

cat("\nAfter removing constant metrics:", length(metrics_to_test), "\n")

# ------------------------------------------------------------
# SPEARMAN CORRELATIONS
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("SPEARMAN CORRELATIONS\n")
cat("====================================================\n")

results <- lapply(
  metrics_to_test,
  function(metric) {

    x <- merged$LINE1_expression
    y <- merged[[metric]]

    complete <- complete.cases(x, y)

    x_complete <- x[complete]
    y_complete <- y[complete]

    if (length(unique(x_complete)) < 2 ||
        length(unique(y_complete)) < 2) {

      return(
        data.frame(
          metric = metric,
          N = sum(complete),
          rho = NA_real_,
          p = NA_real_
        )
      )
    }

    test <- suppressWarnings(
      cor.test(
        x_complete,
        y_complete,
        method = "spearman",
        exact = FALSE
      )
    )

    data.frame(
      metric = metric,
      N = sum(complete),
      rho = unname(test$estimate),
      p = test$p.value
    )
  }
)

results <- bind_rows(results)

# ------------------------------------------------------------
# FDR
# ------------------------------------------------------------

results <- results %>%
  mutate(
    FDR = p.adjust(p, method = "BH"),
    abs_rho = abs(rho),
    direction = case_when(
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ NA_character_
    ),
    strength = case_when(
      abs_rho < 0.20 ~ "Very weak",
      abs_rho < 0.40 ~ "Weak",
      abs_rho < 0.60 ~ "Moderate",
      abs_rho < 0.80 ~ "Strong",
      abs_rho >= 0.80 ~ "Very strong",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(desc(abs_rho))

# ------------------------------------------------------------
# SAVE FULL CLEAN RESULTS
# ------------------------------------------------------------

write_csv(
  results,
  file.path(
    out_dir,
    "LINE1_vs_ALL_metrics_CORRELATIONS_CLEAN.csv"
  )
)

# ------------------------------------------------------------
# NOMINAL p < 0.05
# ------------------------------------------------------------

nominal <- results %>%
  filter(!is.na(p), p < 0.05) %>%
  arrange(p)

write_csv(
  nominal,
  file.path(
    out_dir,
    "LINE1_vs_ALL_metrics_nominal_p05.csv"
  )
)

# ------------------------------------------------------------
# FDR < 0.05
# ------------------------------------------------------------

fdr_sig <- results %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  arrange(FDR)

write_csv(
  fdr_sig,
  file.path(
    out_dir,
    "LINE1_vs_ALL_metrics_FDR05.csv"
  )
)

# ------------------------------------------------------------
# PRIMARY METRICS
# ------------------------------------------------------------

primary_names <- c(
  "interactions_count",
  "interaction_weight",
  "network_density",
  "total_degree",
  "total_strength"
)

primary_results <- results %>%
  filter(metric %in% primary_names)

write_csv(
  primary_results,
  file.path(
    out_dir,
    "LINE1_vs_PRIMARY_metrics.csv"
  )
)

# ------------------------------------------------------------
# DONOR DATA
# ------------------------------------------------------------

donor_data <- merged %>%
  select(
    donor,
    GSM_clean,
    condition,
    LINE1_expression,
    Cells_LINE1,
    all_of(metrics_to_test)
  )

write_csv(
  donor_data,
  file.path(
    out_dir,
    "LINE1_vs_ALL_metrics_DONOR_DATA_CLEAN.csv"
  )
)

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

summary_table <- data.frame(
  n_donors = nrow(merged),
  n_metrics_tested = length(metrics_to_test),
  nominal_p05 = sum(results$p < 0.05, na.rm = TRUE),
  FDR_p05 = sum(results$FDR < 0.05, na.rm = TRUE)
)

write_csv(
  summary_table,
  file.path(
    out_dir,
    "LINE1_vs_ALL_metrics_ANALYSIS_SUMMARY.csv"
  )
)

# ------------------------------------------------------------
# PRINT RESULTS
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("TOP CORRELATIONS BY |rho|\n")
cat("====================================================\n\n")

print(
  results %>%
    select(
      metric,
      N,
      rho,
      p,
      FDR,
      direction,
      strength
    ) %>%
    head(15)
)

cat("\n====================================================\n")
cat("NOMINAL p < 0.05\n")
cat("====================================================\n\n")

if (nrow(nominal) == 0) {
  cat("None\n")
} else {
  print(nominal)
}

cat("\n====================================================\n")
cat("FDR < 0.05\n")
cat("====================================================\n\n")

if (nrow(fdr_sig) == 0) {
  cat("None\n")
} else {
  print(fdr_sig)
}

cat("\n====================================================\n")
cat("PRIMARY METRICS\n")
cat("====================================================\n\n")

print(
  primary_results %>%
    select(
      metric,
      N,
      rho,
      p,
      FDR,
      direction,
      strength
    )
)

cat("\n====================================================\n")
cat("FILES WRITTEN\n")
cat("====================================================\n\n")

print(
  list.files(
    out_dir,
    pattern = "CLEAN|p05|PRIMARY|SUMMARY",
    full.names = TRUE
  )
)

cat("\n====================================================\n")
cat("CLEAN ANALYSIS COMPLETED\n")
cat("====================================================\n")
