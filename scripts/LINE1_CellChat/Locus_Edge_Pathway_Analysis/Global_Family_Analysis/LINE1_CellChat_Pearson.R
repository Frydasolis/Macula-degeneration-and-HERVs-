############################################################
# LINE1_CellChat_Pearson.R
#
# Donor-level association between L1FLnI expression
# and CellChat network remodeling
#
# Statistical unit: DONOR
#
# Primary network outcomes:
#   1. density
#   2. total_strength
#
# Primary correlation:
#   Pearson
#
# Analysis levels:
#   A. Family-wide global L1FLnI
#   B. Family-wide L1FLnI by retinal cell type
#   C. Candidate-locus burden by cell type
#   D. Candidate locus-specific associations
############################################################


############################################################
# 0. Packages
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

cat("\n============================================\n")
cat("L1FLnI - CELLCHAT DONOR-LEVEL PEARSON\n")
cat("============================================\n\n")


############################################################
# 1. Paths
############################################################

seurat_file <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/Seurat/downstream/",
  "SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

network_file <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/CellChat/donor_level_symbol/",
  "network_metrics/donor_network_metrics.csv"
)

candidate_file <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/Seurat/downstream/",
  "TE_analysis/Figure4/",
  "LINE1_CellChat_candidates.csv"
)

outdir <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/CellChat/donor_level_symbol/",
  "network_metrics/LINE1_Pearson"
)

table_dir <- file.path(outdir, "tables")
figure_dir <- file.path(outdir, "figures")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


############################################################
# 2. Load data
############################################################

cat("Loading Seurat object...\n")

seu <- readRDS(seurat_file)

cat("Cells:", ncol(seu), "\n")
cat("Features:", nrow(seu), "\n\n")

cat("Loading donor CellChat metrics...\n")

net <- read.csv(
  network_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Network donors:", nrow(net), "\n\n")

cat("Loading candidate loci...\n")

cand <- read.csv(
  candidate_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
# 3. Check required metadata
############################################################

required_meta <- c(
  "orig.ident",
  "disease_state",
  "HRCA_majorclass"
)

missing_meta <- setdiff(
  required_meta,
  colnames(seu@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing Seurat metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}


############################################################
# 4. Standardize network table column names
############################################################

cat("Network columns:\n")
print(colnames(net))

if (!"Donor" %in% colnames(net)) {
  stop("Network table does not contain column 'Donor'.")
}

if (!"density" %in% colnames(net)) {
  stop("Network table does not contain 'density'.")
}

if (!"total_strength" %in% colnames(net)) {
  stop("Network table does not contain 'total_strength'.")
}

net$Donor <- as.character(net$Donor)


############################################################
# 5. Identify all L1FLnI loci
############################################################

all_features <- rownames(seu)

l1_features <- grep(
  "^TE-L1FLnI-",
  all_features,
  value = TRUE
)

cat("\n============================================\n")
cat("L1FLnI FAMILY\n")
cat("============================================\n")

cat("Total L1FLnI loci:", length(l1_features), "\n")

if (length(l1_features) == 0) {
  stop("No L1FLnI features found.")
}

write.table(
  data.frame(
    TE_feature = l1_features,
    stringsAsFactors = FALSE
  ),
  file.path(
    table_dir,
    "all_L1FLnI_features.csv"
  ),
  sep = ",",
  row.names = FALSE,
  quote = TRUE
)


############################################################
# 6. Metadata
############################################################

meta <- seu@meta.data

meta$Donor <- as.character(meta$orig.ident)
meta$Condition <- as.character(meta$disease_state)
meta$CellType <- as.character(meta$HRCA_majorclass)

cat("\nDonors:\n")
print(sort(unique(meta$Donor)))

cat("\nConditions:\n")
print(table(meta$Condition))

cat("\nCell types:\n")
print(table(meta$CellType))


############################################################
# 7. Get normalized sparse matrix
############################################################

DefaultAssay(seu) <- "RNA"

cat("\nExtracting normalized RNA layer...\n")

expr <- GetAssayData(
  seu,
  assay = "RNA",
  layer = "data"
)

l1_expr <- expr[
  l1_features,
  ,
  drop = FALSE
]

cat(
  "L1FLnI matrix:",
  nrow(l1_expr),
  "loci x",
  ncol(l1_expr),
  "cells\n"
)


############################################################
# 8. Helper: safe Pearson correlation
############################################################

safe_pearson <- function(x, y) {

  keep <- is.finite(x) & is.finite(y)

  x <- x[keep]
  y <- y[keep]

  n <- length(x)

  if (n < 3) {
    return(
      data.frame(
        n = n,
        Pearson_r = NA_real_,
        R2 = NA_real_,
        P_value = NA_real_
      )
    )
  }

  if (sd(x) == 0 || sd(y) == 0) {
    return(
      data.frame(
        n = n,
        Pearson_r = NA_real_,
        R2 = NA_real_,
        P_value = NA_real_
      )
    )
  }

  test <- cor.test(
    x,
    y,
    method = "pearson"
  )

  r <- unname(test$estimate)

  data.frame(
    n = n,
    Pearson_r = r,
    R2 = r^2,
    P_value = test$p.value
  )
}


############################################################
# 9. Helper: donor condition
############################################################

donor_condition <- unique(
  meta[
    ,
    c(
      "Donor",
      "Condition"
    )
  ]
)

if (anyDuplicated(donor_condition$Donor)) {

  check_condition <- table(
    donor_condition$Donor
  )

  if (any(check_condition > 1)) {
    stop(
      "At least one donor has more than one disease condition."
    )
  }
}


############################################################
# 10. FAMILY-WIDE GLOBAL L1FLnI BURDEN
#
# One value per donor.
#
# Expression burden:
# Mean normalized expression across ALL L1FLnI
# molecules/features and all cells from donor.
#
# Detection breadth:
# Fraction of 13,292 loci detected in >=1 cell
# from that donor.
############################################################

cat("\n============================================\n")
cat("A. FAMILY-WIDE GLOBAL L1FLnI\n")
cat("============================================\n")

global_list <- list()

donors <- sort(unique(meta$Donor))

for (d in donors) {

  cells <- rownames(meta)[
    meta$Donor == d
  ]

  mat <- l1_expr[
    ,
    cells,
    drop = FALSE
  ]

  # Mean expression of each locus across donor cells
  locus_means <- Matrix::rowMeans(mat)

  # Family-wide expression burden
  expression_burden <- mean(locus_means)

  # Number of loci detected in >=1 cell
  detected_loci <- sum(
    Matrix::rowSums(mat > 0) > 0
  )

  detection_breadth <- detected_loci /
    length(l1_features)

  global_list[[d]] <- data.frame(
    Donor = d,
    n_cells = length(cells),
    n_L1FLnI_loci = length(l1_features),
    detected_loci = detected_loci,
    detection_breadth = detection_breadth,
    L1_expression_burden = expression_burden,
    stringsAsFactors = FALSE
  )
}

global_burden <- do.call(
  rbind,
  global_list
)

rownames(global_burden) <- NULL

global_burden <- merge(
  global_burden,
  donor_condition,
  by = "Donor",
  all.x = TRUE
)

global_burden <- merge(
  global_burden,
  net[
    ,
    c(
      "Donor",
      "density",
      "total_strength"
    )
  ],
  by = "Donor",
  all.x = TRUE
)

write.csv(
  global_burden,
  file.path(
    table_dir,
    "family_global_L1FLnI_burden.csv"
  ),
  row.names = FALSE
)

cat("\nGlobal donor burden:\n")
print(global_burden)


############################################################
# 11. Global Pearson correlations
############################################################

global_results <- list()

predictors <- c(
  "L1_expression_burden",
  "detection_breadth"
)

network_metrics <- c(
  "density",
  "total_strength"
)

counter <- 1

for (pred in predictors) {

  for (metric in network_metrics) {

    res <- safe_pearson(
      global_burden[[pred]],
      global_burden[[metric]]
    )

    global_results[[counter]] <- data.frame(
      Analysis_level = "Family_global",
      CellType = "ALL",
      Feature = "ALL_L1FLnI",
      Predictor = pred,
      Network_metric = metric,
      n = res$n,
      Pearson_r = res$Pearson_r,
      R2 = res$R2,
      P_value = res$P_value,
      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

global_results <- do.call(
  rbind,
  global_results
)

global_results$FDR <- p.adjust(
  global_results$P_value,
  method = "BH"
)

write.csv(
  global_results,
  file.path(
    table_dir,
    "family_global_Pearson.csv"
  ),
  row.names = FALSE
)

cat("\nGlobal Pearson results:\n")
print(global_results)


############################################################
# 12. FAMILY-WIDE L1FLnI BY CELL TYPE
#
# For every donor x retinal major cell type:
#   - all 13,292 L1FLnI loci
#   - expression burden
#   - detection breadth
############################################################

cat("\n============================================\n")
cat("B. FAMILY-WIDE L1FLnI BY CELL TYPE\n")
cat("============================================\n")

celltypes <- sort(
  unique(meta$CellType)
)

celltype_list <- list()

counter <- 1

for (ct in celltypes) {

  cat("Processing:", ct, "\n")

  ct_donors <- sort(
    unique(
      meta$Donor[
        meta$CellType == ct
      ]
    )
  )

  for (d in ct_donors) {

    cells <- rownames(meta)[
      meta$Donor == d &
        meta$CellType == ct
    ]

    if (length(cells) == 0) {
      next
    }

    mat <- l1_expr[
      ,
      cells,
      drop = FALSE
    ]

    locus_means <- Matrix::rowMeans(mat)

    expression_burden <- mean(
      locus_means
    )

    detected_loci <- sum(
      Matrix::rowSums(mat > 0) > 0
    )

    detection_breadth <- detected_loci /
      length(l1_features)

    celltype_list[[counter]] <- data.frame(
      Donor = d,
      CellType = ct,
      n_cells = length(cells),
      n_L1FLnI_loci = length(l1_features),
      detected_loci = detected_loci,
      detection_breadth = detection_breadth,
      L1_expression_burden = expression_burden,
      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

celltype_burden <- do.call(
  rbind,
  celltype_list
)

celltype_burden <- merge(
  celltype_burden,
  donor_condition,
  by = "Donor",
  all.x = TRUE
)

celltype_burden <- merge(
  celltype_burden,
  net[
    ,
    c(
      "Donor",
      "density",
      "total_strength"
    )
  ],
  by = "Donor",
  all.x = TRUE
)

write.csv(
  celltype_burden,
  file.path(
    table_dir,
    "family_celltype_L1FLnI_burden.csv"
  ),
  row.names = FALSE
)


############################################################
# 13. Cell-type Pearson correlations
############################################################

celltype_results <- list()

counter <- 1

for (ct in sort(unique(celltype_burden$CellType))) {

  tmp <- celltype_burden[
    celltype_burden$CellType == ct,
    ,
    drop = FALSE
  ]

  for (pred in predictors) {

    for (metric in network_metrics) {

      res <- safe_pearson(
        tmp[[pred]],
        tmp[[metric]]
      )

      celltype_results[[counter]] <- data.frame(
        Analysis_level = "Family_celltype",
        CellType = ct,
        Feature = "ALL_L1FLnI",
        Predictor = pred,
        Network_metric = metric,
        n = res$n,
        Pearson_r = res$Pearson_r,
        R2 = res$R2,
        P_value = res$P_value,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }
}

celltype_results <- do.call(
  rbind,
  celltype_results
)

celltype_results$FDR <- p.adjust(
  celltype_results$P_value,
  method = "BH"
)

write.csv(
  celltype_results,
  file.path(
    table_dir,
    "family_celltype_Pearson.csv"
  ),
  row.names = FALSE
)

cat("\nFamily-wide cell-type results:\n")
print(celltype_results)


############################################################
# 14. Candidate loci
############################################################

cat("\n============================================\n")
cat("C. CANDIDATE LOCI\n")
cat("============================================\n")

cat("Candidate columns:\n")
print(colnames(cand))

# Keep L1FLnI only if TE_family exists
if ("TE_family" %in% colnames(cand)) {

  cand <- cand[
    cand$TE_family == "L1FLnI",
    ,
    drop = FALSE
  ]
}

# Determine candidate feature column
if ("TE_feature" %in% colnames(cand)) {

  feature_col <- "TE_feature"

} else if ("feature" %in% colnames(cand)) {

  feature_col <- "feature"

} else {

  stop(
    "Could not identify candidate TE feature column."
  )
}

# Determine cell-type column
if ("HRCA_majorclass" %in% colnames(cand)) {

  candidate_celltype_col <- "HRCA_majorclass"

} else if ("CellType" %in% colnames(cand)) {

  candidate_celltype_col <- "CellType"

} else {

  stop(
    "Could not identify candidate cell-type column."
  )
}

candidate_map <- unique(
  data.frame(
    CellType = as.character(
      cand[[candidate_celltype_col]]
    ),
    TE_feature = as.character(
      cand[[feature_col]]
    ),
    stringsAsFactors = FALSE
  )
)

candidate_map <- candidate_map[
  candidate_map$TE_feature %in% rownames(expr),
  ,
  drop = FALSE
]

cat(
  "Unique candidate cell-type/locus pairs:",
  nrow(candidate_map),
  "\n"
)

cat(
  "Unique candidate loci:",
  length(unique(candidate_map$TE_feature)),
  "\n"
)

print(candidate_map)

write.csv(
  candidate_map,
  file.path(
    table_dir,
    "candidate_L1FLnI_map.csv"
  ),
  row.names = FALSE
)


############################################################
# 15. Candidate locus expression per donor x cell type
############################################################

candidate_locus_list <- list()

counter <- 1

for (ct in unique(candidate_map$CellType)) {

  loci <- unique(
    candidate_map$TE_feature[
      candidate_map$CellType == ct
    ]
  )

  ct_donors <- sort(
    unique(
      meta$Donor[
        meta$CellType == ct
      ]
    )
  )

  for (d in ct_donors) {

    cells <- rownames(meta)[
      meta$Donor == d &
        meta$CellType == ct
    ]

    if (length(cells) == 0) {
      next
    }

    mat <- expr[
      loci,
      cells,
      drop = FALSE
    ]

    means <- Matrix::rowMeans(mat)

    detection <- Matrix::rowMeans(
      mat > 0
    )

    for (j in seq_along(loci)) {

      candidate_locus_list[[counter]] <- data.frame(
        Donor = d,
        CellType = ct,
        TE_feature = loci[j],
        n_cells = length(cells),
        mean_expression = means[j],
        detection_rate = detection[j],
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }
}

candidate_locus <- do.call(
  rbind,
  candidate_locus_list
)

candidate_locus <- merge(
  candidate_locus,
  donor_condition,
  by = "Donor",
  all.x = TRUE
)

candidate_locus <- merge(
  candidate_locus,
  net[
    ,
    c(
      "Donor",
      "density",
      "total_strength"
    )
  ],
  by = "Donor",
  all.x = TRUE
)

write.csv(
  candidate_locus,
  file.path(
    table_dir,
    "candidate_locus_donor_expression.csv"
  ),
  row.names = FALSE
)


############################################################
# 16. Candidate burden by donor x cell type
############################################################

candidate_groups <- unique(
  candidate_locus[
    ,
    c(
      "Donor",
      "CellType"
    )
  ]
)

candidate_burden_list <- list()

for (i in seq_len(nrow(candidate_groups))) {

  d <- candidate_groups$Donor[i]
  ct <- candidate_groups$CellType[i]

  tmp <- candidate_locus[
    candidate_locus$Donor == d &
      candidate_locus$CellType == ct,
    ,
    drop = FALSE
  ]

  candidate_burden_list[[i]] <- data.frame(
    Donor = d,
    CellType = ct,
    Condition = tmp$Condition[1],
    n_cells = tmp$n_cells[1],
    n_candidate_loci = nrow(tmp),
    candidate_expression_burden =
      mean(tmp$mean_expression),
    candidate_detection_burden =
      mean(tmp$detection_rate),
    density = tmp$density[1],
    total_strength = tmp$total_strength[1],
    stringsAsFactors = FALSE
  )
}

candidate_burden <- do.call(
  rbind,
  candidate_burden_list
)

write.csv(
  candidate_burden,
  file.path(
    table_dir,
    "candidate_burden_by_donor_celltype.csv"
  ),
  row.names = FALSE
)


############################################################
# 17. Candidate burden Pearson
############################################################

candidate_results <- list()

counter <- 1

candidate_predictors <- c(
  "candidate_expression_burden",
  "candidate_detection_burden"
)

for (ct in sort(unique(candidate_burden$CellType))) {

  tmp <- candidate_burden[
    candidate_burden$CellType == ct,
    ,
    drop = FALSE
  ]

  for (pred in candidate_predictors) {

    for (metric in network_metrics) {

      res <- safe_pearson(
        tmp[[pred]],
        tmp[[metric]]
      )

      candidate_results[[counter]] <- data.frame(
        Analysis_level = "Candidate_burden",
        CellType = ct,
        Feature = "Candidate_L1FLnI",
        Predictor = pred,
        Network_metric = metric,
        n = res$n,
        Pearson_r = res$Pearson_r,
        R2 = res$R2,
        P_value = res$P_value,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }
}

candidate_results <- do.call(
  rbind,
  candidate_results
)

candidate_results$FDR <- p.adjust(
  candidate_results$P_value,
  method = "BH"
)

write.csv(
  candidate_results,
  file.path(
    table_dir,
    "candidate_burden_Pearson.csv"
  ),
  row.names = FALSE
)


############################################################
# 18. LOCUS-SPECIFIC PEARSON
############################################################

cat("\n============================================\n")
cat("D. LOCUS-SPECIFIC PEARSON\n")
cat("============================================\n")

locus_groups <- unique(
  candidate_locus[
    ,
    c(
      "CellType",
      "TE_feature"
    )
  ]
)

locus_results <- list()

counter <- 1

for (i in seq_len(nrow(locus_groups))) {

  ct <- locus_groups$CellType[i]
  locus <- locus_groups$TE_feature[i]

  tmp <- candidate_locus[
    candidate_locus$CellType == ct &
      candidate_locus$TE_feature == locus,
    ,
    drop = FALSE
  ]

  for (pred in c(
    "mean_expression",
    "detection_rate"
  )) {

    for (metric in network_metrics) {

      res <- safe_pearson(
        tmp[[pred]],
        tmp[[metric]]
      )

      locus_results[[counter]] <- data.frame(
        Analysis_level = "Locus_specific",
        CellType = ct,
        Feature = locus,
        Predictor = pred,
        Network_metric = metric,
        n = res$n,
        Pearson_r = res$Pearson_r,
        R2 = res$R2,
        P_value = res$P_value,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }
}

locus_results <- do.call(
  rbind,
  locus_results
)

locus_results$FDR <- p.adjust(
  locus_results$P_value,
  method = "BH"
)

write.csv(
  locus_results,
  file.path(
    table_dir,
    "locus_specific_Pearson.csv"
  ),
  row.names = FALSE
)

cat("\nLocus-specific results:\n")
print(locus_results)


############################################################
# 19. MASTER RESULTS TABLE
############################################################

master <- rbind(
  global_results,
  celltype_results,
  candidate_results,
  locus_results
)

master <- master[
  order(
    master$Analysis_level,
    master$P_value
  ),
  ,
  drop = FALSE
]

write.csv(
  master,
  file.path(
    table_dir,
    "MASTER_Pearson_results.csv"
  ),
  row.names = FALSE
)


############################################################
# 20. Global scatter plot helper
############################################################

make_scatter <- function(
    dat,
    xvar,
    yvar,
    xlab,
    ylab,
    title,
    filename) {

  test <- safe_pearson(
    dat[[xvar]],
    dat[[yvar]]
  )

  subtitle_text <- paste0(
    "Pearson r = ",
    ifelse(
      is.na(test$Pearson_r),
      "NA",
      sprintf("%.3f", test$Pearson_r)
    ),
    "   P = ",
    ifelse(
      is.na(test$P_value),
      "NA",
      format.pval(
        test$P_value,
        digits = 3,
        eps = 0.001
      )
    ),
    "   n = ",
    test$n
  )

  p <- ggplot(
    dat,
    aes(
      x = .data[[xvar]],
      y = .data[[yvar]]
    )
  ) +

    geom_point(
      aes(
        shape = Condition
      ),
      size = 3,
      alpha = 0.85
    ) +

    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.8
    ) +

    labs(
      title = title,
      subtitle = subtitle_text,
      x = xlab,
      y = ylab,
      shape = "Condition"
    ) +

    theme_classic(
      base_size = 12
    ) +

    theme(
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "right"
    )

  ggsave(
    paste0(filename, ".pdf"),
    p,
    width = 6.5,
    height = 5.5
  )

  ggsave(
    paste0(filename, ".png"),
    p,
    width = 6.5,
    height = 5.5,
    dpi = 400
  )
}


############################################################
# 21. Global family-wide scatter plots
############################################################

make_scatter(
  global_burden,
  "L1_expression_burden",
  "density",
  "Global L1FLnI expression burden",
  "Network density",
  "Global L1FLnI vs network density",
  file.path(
    figure_dir,
    "Global_L1FLnI_vs_Density"
  )
)

make_scatter(
  global_burden,
  "L1_expression_burden",
  "total_strength",
  "Global L1FLnI expression burden",
  "Total communication strength",
  "Global L1FLnI vs total communication strength",
  file.path(
    figure_dir,
    "Global_L1FLnI_vs_TotalStrength"
  )
)

make_scatter(
  global_burden,
  "detection_breadth",
  "density",
  "L1FLnI detection breadth",
  "Network density",
  "L1FLnI detection breadth vs network density",
  file.path(
    figure_dir,
    "Global_L1FLnI_Breadth_vs_Density"
  )
)

make_scatter(
  global_burden,
  "detection_breadth",
  "total_strength",
  "L1FLnI detection breadth",
  "Total communication strength",
  "L1FLnI detection breadth vs total communication strength",
  file.path(
    figure_dir,
    "Global_L1FLnI_Breadth_vs_TotalStrength"
  )
)


############################################################
# 22. Cell-type heatmap
############################################################

heat <- celltype_results[
  celltype_results$Predictor ==
    "L1_expression_burden",
  ,
  drop = FALSE
]

heat$Network_metric_label <- ifelse(
  heat$Network_metric == "density",
  "Network density",
  "Total strength"
)

p_heat <- ggplot(
  heat,
  aes(
    x = Network_metric_label,
    y = CellType,
    fill = Pearson_r
  )
) +

  geom_tile(
    linewidth = 0.4
  ) +

  geom_text(
    aes(
      label = sprintf(
        "r=%.2f\nP=%.3g",
        Pearson_r,
        P_value
      )
    ),
    size = 3
  ) +

  scale_fill_gradient2(
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +

  labs(
    title = "Family-wide L1FLnI burden by retinal cell type",
    x = NULL,
    y = NULL
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    figure_dir,
    "CellType_L1FLnI_Pearson_heatmap.pdf"
  ),
  p_heat,
  width = 7,
  height = 6
)

ggsave(
  file.path(
    figure_dir,
    "CellType_L1FLnI_Pearson_heatmap.png"
  ),
  p_heat,
  width = 7,
  height = 6,
  dpi = 400
)


############################################################
# 23. Locus-specific heatmap
############################################################

locus_heat <- locus_results[
  locus_results$Predictor ==
    "mean_expression",
  ,
  drop = FALSE
]

locus_heat$Feature_clean <- sub(
  "^TE-L1FLnI-",
  "",
  locus_heat$Feature
)

locus_heat$Network_metric_label <- ifelse(
  locus_heat$Network_metric == "density",
  "Network density",
  "Total strength"
)

p_locus <- ggplot(
  locus_heat,
  aes(
    x = Network_metric_label,
    y = Feature_clean,
    fill = Pearson_r
  )
) +

  geom_tile(
    linewidth = 0.4
  ) +

  geom_text(
    aes(
      label = sprintf(
        "r=%.2f\nP=%.3g",
        Pearson_r,
        P_value
      )
    ),
    size = 2.8
  ) +

  scale_fill_gradient2(
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +

  facet_grid(
    CellType ~ .,
    scales = "free_y",
    space = "free_y"
  ) +

  labs(
    title = "Locus-specific L1FLnI associations with retinal network remodeling",
    x = NULL,
    y = "L1FLnI locus"
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    strip.text.y = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  file.path(
    figure_dir,
    "Locus_specific_L1FLnI_Pearson_heatmap.pdf"
  ),
  p_locus,
  width = 8,
  height = 8
)

ggsave(
  file.path(
    figure_dir,
    "Locus_specific_L1FLnI_Pearson_heatmap.png"
  ),
  p_locus,
  width = 8,
  height = 8,
  dpi = 400
)


############################################################
# 24. Significant / strongest results summary
############################################################

sig <- master[
  !is.na(master$P_value) &
    master$P_value < 0.05,
  ,
  drop = FALSE
]

sig <- sig[
  order(sig$P_value),
  ,
  drop = FALSE
]

write.csv(
  sig,
  file.path(
    table_dir,
    "Pearson_rawP_less_0.05.csv"
  ),
  row.names = FALSE
)

fdr_sig <- master[
  !is.na(master$FDR) &
    master$FDR < 0.05,
  ,
  drop = FALSE
]

fdr_sig <- fdr_sig[
  order(fdr_sig$FDR),
  ,
  drop = FALSE
]

write.csv(
  fdr_sig,
  file.path(
    table_dir,
    "Pearson_FDR_less_0.05.csv"
  ),
  row.names = FALSE
)


############################################################
# 25. Summary
############################################################

cat("\n\n============================================\n")
cat("ANALYSIS COMPLETED\n")
cat("============================================\n")

cat("\nDonors:", length(unique(meta$Donor)), "\n")
cat("L1FLnI loci:", length(l1_features), "\n")
cat("Cell types:", length(unique(meta$CellType)), "\n")
cat(
  "Candidate loci:",
  length(unique(candidate_map$TE_feature)),
  "\n"
)

cat("\nPrimary network outcomes:\n")
cat("  - density\n")
cat("  - total_strength\n")

cat("\nPrimary correlation:\n")
cat("  - Pearson\n")

cat("\nOutput directory:\n")
cat(outdir, "\n")

cat("\nMain tables:\n")
cat(
  file.path(
    table_dir,
    "family_global_Pearson.csv"
  ),
  "\n"
)
cat(
  file.path(
    table_dir,
    "family_celltype_Pearson.csv"
  ),
  "\n"
)
cat(
  file.path(
    table_dir,
    "candidate_burden_Pearson.csv"
  ),
  "\n"
)
cat(
  file.path(
    table_dir,
    "locus_specific_Pearson.csv"
  ),
  "\n"
)
cat(
  file.path(
    table_dir,
    "MASTER_Pearson_results.csv"
  ),
  "\n"
)

cat("\nDONE\n")
