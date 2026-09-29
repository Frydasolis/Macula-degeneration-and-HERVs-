############################################################
# LINE1_CellChat_Pearson_aligned.R
#
# Donor-level association between L1FLnI expression
# and CellChat network remodeling
#
# Statistical unit: DONOR
#
# Network outcomes:
#   1. density
#   2. total_strength
#
# Correlations:
#   PRIMARY: Pearson
#   SENSITIVITY: Spearman
#
# Conditions:
#   ALL
#   Healthy Control
#   Dry AMD
#   Wet AMD
#
# Analysis levels:
#   A. Family-wide global L1FLnI
#   B. Family-wide L1FLnI by retinal cell type
#   C. Candidate-locus burden by cell type
#   D. Candidate locus-specific associations
#
# IMPORTANT:
#   - n = number of DONORS
#   - n < 4 is marked DESCRIPTIVE_ONLY
#   - ALL should reproduce previous donor-level results
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
cat("L1FLnI - CELLCHAT ALIGNED DONOR ANALYSIS\n")
cat("Pearson primary + Spearman sensitivity\n")
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
  "network_metrics/LINE1_Pearson_aligned"
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
# 3. Check metadata
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
# 4. Check network table
############################################################

cat("Network columns:\n")
print(colnames(net))

required_net <- c(
  "Donor",
  "density",
  "total_strength"
)

missing_net <- setdiff(
  required_net,
  colnames(net)
)

if (length(missing_net) > 0) {
  stop(
    "Missing network columns: ",
    paste(missing_net, collapse = ", ")
  )
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

write.csv(
  data.frame(
    TE_feature = l1_features,
    stringsAsFactors = FALSE
  ),
  file.path(
    table_dir,
    "all_L1FLnI_features.csv"
  ),
  row.names = FALSE
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

cat("\nConditions by cells:\n")
print(table(meta$Condition))

cat("\nCell types:\n")
print(table(meta$CellType))


############################################################
# 7. Donor-condition table
############################################################

donor_condition <- unique(
  meta[, c("Donor", "Condition")]
)

check_condition <- table(
  donor_condition$Donor
)

if (any(check_condition > 1)) {
  stop(
    "At least one donor has more than one disease condition."
  )
}

cat("\nDonors by condition:\n")
print(table(donor_condition$Condition))


############################################################
# 8. Conditions to analyze
############################################################

conditions_to_test <- c(
  "ALL",
  "Healthy Control",
  "Dry AMD",
  "Wet AMD"
)

cat("\nConditions that will be analyzed:\n")
print(conditions_to_test)


############################################################
# 9. Get normalized RNA matrix
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
# 10. Safe correlation function
############################################################

safe_correlations <- function(x, y) {

  keep <- is.finite(x) & is.finite(y)

  x <- x[keep]
  y <- y[keep]

  n <- length(x)

  inference <- ifelse(
    n >= 4,
    "INFERENTIAL",
    "DESCRIPTIVE_ONLY"
  )

  empty_result <- function() {
    data.frame(
      n = n,
      Inference = inference,
      Pearson_r = NA_real_,
      Pearson_R2 = NA_real_,
      Pearson_P = NA_real_,
      Spearman_rho = NA_real_,
      Spearman_P = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (n < 3) {
    return(empty_result())
  }

  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(empty_result())
  }

  pearson_test <- tryCatch(
    cor.test(
      x,
      y,
      method = "pearson"
    ),
    error = function(e) NULL
  )

  spearman_test <- tryCatch(
    suppressWarnings(
      cor.test(
        x,
        y,
        method = "spearman",
        exact = FALSE
      )
    ),
    error = function(e) NULL
  )

  pearson_r <- if (
    is.null(pearson_test)
  ) {
    NA_real_
  } else {
    unname(pearson_test$estimate)
  }

  pearson_p <- if (
    is.null(pearson_test)
  ) {
    NA_real_
  } else {
    pearson_test$p.value
  }

  spearman_rho <- if (
    is.null(spearman_test)
  ) {
    NA_real_
  } else {
    unname(spearman_test$estimate)
  }

  spearman_p <- if (
    is.null(spearman_test)
  ) {
    NA_real_
  } else {
    spearman_test$p.value
  }

  data.frame(
    n = n,
    Inference = inference,
    Pearson_r = pearson_r,
    Pearson_R2 = pearson_r^2,
    Pearson_P = pearson_p,
    Spearman_rho = spearman_rho,
    Spearman_P = spearman_p,
    stringsAsFactors = FALSE
  )
}


############################################################
# 11. Helper: subset by condition
############################################################

subset_condition <- function(dat, condition) {

  if (condition == "ALL") {
    return(dat)
  }

  dat[
    dat$Condition == condition,
    ,
    drop = FALSE
  ]
}


############################################################
# 12. Helper: add FDR
#
# BH correction is performed separately within:
#   Analysis_level x Condition
#
# Pearson and Spearman are corrected independently.
############################################################

add_grouped_fdr <- function(dat) {

  dat$Pearson_FDR <- NA_real_
  dat$Spearman_FDR <- NA_real_

  groups <- unique(
    dat[, c(
      "Analysis_level",
      "Condition"
    )]
  )

  for (i in seq_len(nrow(groups))) {

    level_i <- groups$Analysis_level[i]
    condition_i <- groups$Condition[i]

    idx <- which(
      dat$Analysis_level == level_i &
      dat$Condition == condition_i
    )

    valid_p <- idx[
      !is.na(dat$Pearson_P[idx])
    ]

    if (length(valid_p) > 0) {
      dat$Pearson_FDR[valid_p] <- p.adjust(
        dat$Pearson_P[valid_p],
        method = "BH"
      )
    }

    valid_s <- idx[
      !is.na(dat$Spearman_P[idx])
    ]

    if (length(valid_s) > 0) {
      dat$Spearman_FDR[valid_s] <- p.adjust(
        dat$Spearman_P[valid_s],
        method = "BH"
      )
    }
  }

  dat
}


############################################################
# 13. Common network metrics
############################################################

network_metrics <- c(
  "density",
  "total_strength"
)
############################################################
# 14. A. FAMILY-WIDE GLOBAL L1FLnI
#
# One value per donor.
#
# Expression burden:
# Mean normalized expression across ALL L1FLnI loci
# and all cells from that donor.
#
# Detection breadth:
# Fraction of all L1FLnI loci detected in >=1 cell
# from that donor.
############################################################

cat("\n============================================\n")
cat("A. FAMILY-WIDE GLOBAL L1FLnI\n")
cat("============================================\n")

global_list <- list()

donors <- sort(
  unique(meta$Donor)
)

for (d in donors) {

  cells <- rownames(meta)[
    meta$Donor == d
  ]

  if (length(cells) == 0) {
    next
  }

  mat <- l1_expr[
    ,
    cells,
    drop = FALSE
  ]

  # Mean expression of every L1FLnI locus
  # across all cells from this donor
  locus_means <- Matrix::rowMeans(mat)

  # Family-wide expression burden
  expression_burden <- mean(
    locus_means
  )

  # Number of L1FLnI loci detected
  # in at least one cell
  detected_loci <- sum(
    Matrix::rowSums(mat > 0) > 0
  )

  # Fraction of all L1FLnI loci detected
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


############################################################
# Add disease condition
############################################################

global_burden <- merge(
  global_burden,
  donor_condition,
  by = "Donor",
  all.x = TRUE
)


############################################################
# Add CellChat global network metrics
############################################################

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


############################################################
# Save donor-level global burden
############################################################

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
# 15. Global correlations
#
# ALL + each disease condition
############################################################

cat("\n============================================\n")
cat("A. GLOBAL CORRELATIONS\n")
cat("============================================\n")

global_results <- list()

counter <- 1

global_predictors <- c(
  "L1_expression_burden",
  "detection_breadth"
)

for (condition_i in conditions_to_test) {

  cat("\nCondition:", condition_i, "\n")

  tmp <- subset_condition(
    global_burden,
    condition_i
  )

  cat(
    "Donors available:",
    length(unique(tmp$Donor)),
    "\n"
  )

  for (pred in global_predictors) {

    for (metric in network_metrics) {

      res <- safe_correlations(
        tmp[[pred]],
        tmp[[metric]]
      )

      global_results[[counter]] <- data.frame(
        Analysis_level = "Family_global",
        Condition = condition_i,
        CellType = "ALL",
        Feature = "ALL_L1FLnI",
        Predictor = pred,
        Network_metric = metric,
        n = res$n,
        Inference = res$Inference,
        Pearson_r = res$Pearson_r,
        Pearson_R2 = res$Pearson_R2,
        Pearson_P = res$Pearson_P,
        Spearman_rho = res$Spearman_rho,
        Spearman_P = res$Spearman_P,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }
}


global_results <- do.call(
  rbind,
  global_results
)

global_results <- add_grouped_fdr(
  global_results
)


############################################################
# Save global results
############################################################

write.csv(
  global_results,
  file.path(
    table_dir,
    "family_global_correlations.csv"
  ),
  row.names = FALSE
)

cat("\nGlobal correlation results:\n")
print(global_results)


############################################################
# 16. B. FAMILY-WIDE L1FLnI BY CELL TYPE
#
# For every donor x retinal major cell type:
#
#   - all L1FLnI loci
#   - expression burden
#   - detection breadth
#
# Statistical unit remains DONOR.
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

  cat("\nProcessing cell type:", ct, "\n")

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

    # Mean expression for each L1FLnI locus
    # within this donor/cell-type combination
    locus_means <- Matrix::rowMeans(mat)

    # Mean family-wide expression burden
    expression_burden <- mean(
      locus_means
    )

    # Number of loci detected in at least one cell
    detected_loci <- sum(
      Matrix::rowSums(mat > 0) > 0
    )

    # Fraction of all L1FLnI loci detected
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


############################################################
# Combine donor/cell-type results
############################################################

celltype_burden <- do.call(
  rbind,
  celltype_list
)

rownames(celltype_burden) <- NULL


############################################################
# Add disease condition
############################################################

celltype_burden <- merge(
  celltype_burden,
  donor_condition,
  by = "Donor",
  all.x = TRUE
)


############################################################
# Add CellChat global network metrics
############################################################

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


############################################################
# Save donor x cell-type burden table
############################################################

write.csv(
  celltype_burden,
  file.path(
    table_dir,
    "family_celltype_L1FLnI_burden.csv"
  ),
  row.names = FALSE
)


############################################################
# Basic QC
############################################################

cat("\nDonor x cell-type observations:\n")
print(
  table(
    celltype_burden$CellType,
    celltype_burden$Condition
  )
)


############################################################
# 17. Cell-type correlations
#
# ALL + Healthy + Dry + Wet
############################################################

cat("\n============================================\n")
cat("B. CELL-TYPE CORRELATIONS\n")
cat("============================================\n")

celltype_results <- list()

counter <- 1

celltype_predictors <- c(
  "L1_expression_burden",
  "detection_breadth"
)

for (ct in sort(
  unique(celltype_burden$CellType)
)) {

  cat("\nCell type:", ct, "\n")

  ct_data <- celltype_burden[
    celltype_burden$CellType == ct,
    ,
    drop = FALSE
  ]

  for (condition_i in conditions_to_test) {

    tmp <- subset_condition(
      ct_data,
      condition_i
    )

    for (pred in celltype_predictors) {

      for (metric in network_metrics) {

        res <- safe_correlations(
          tmp[[pred]],
          tmp[[metric]]
        )

        celltype_results[[counter]] <- data.frame(
          Analysis_level = "Family_celltype",
          Condition = condition_i,
          CellType = ct,
          Feature = "ALL_L1FLnI",
          Predictor = pred,
          Network_metric = metric,
          n = res$n,
          Inference = res$Inference,
          Pearson_r = res$Pearson_r,
          Pearson_R2 = res$Pearson_R2,
          Pearson_P = res$Pearson_P,
          Spearman_rho = res$Spearman_rho,
          Spearman_P = res$Spearman_P,
          stringsAsFactors = FALSE
        )

        counter <- counter + 1
      }
    }
  }
}


############################################################
# Combine results
############################################################

celltype_results <- do.call(
  rbind,
  celltype_results
)

celltype_results <- add_grouped_fdr(
  celltype_results
)


############################################################
# Save
############################################################

write.csv(
  celltype_results,
  file.path(
    table_dir,
    "family_celltype_correlations.csv"
  ),
  row.names = FALSE
)


############################################################
# Print strongest results for QC
############################################################

cat("\nStrongest family-by-cell-type Pearson results:\n")

tmp_print <- celltype_results[
  !is.na(celltype_results$Pearson_P),
  ,
  drop = FALSE
]

tmp_print <- tmp_print[
  order(tmp_print$Pearson_P),
  ,
  drop = FALSE
]

print(
  head(
    tmp_print,
    20
  )
)


############################################################
# 18. QC: ALL should reproduce original analysis
############################################################

cat("\n============================================\n")
cat("QC: FAMILY CELL-TYPE — ALL DONORS\n")
cat("============================================\n")

qc_family <- celltype_results[
  celltype_results$Condition == "ALL",
  ,
  drop = FALSE
]

qc_family <- qc_family[
  order(qc_family$Pearson_P),
  ,
  drop = FALSE
]

print(
  qc_family[
    ,
    c(
      "CellType",
      "Predictor",
      "Network_metric",
      "n",
      "Pearson_r",
      "Pearson_R2",
      "Pearson_P",
      "Pearson_FDR",
      "Spearman_rho",
      "Spearman_P"
    )
  ]
)


############################################################
# END OF BLOCK 2
############################################################
############################################################
# 19. C. CANDIDATE LOCI
############################################################

cat("\n============================================\n")
cat("C. CANDIDATE LOCI\n")
cat("============================================\n")

cat("Candidate columns:\n")
print(colnames(cand))


############################################################
# Keep L1FLnI candidates only
############################################################

if ("TE_family" %in% colnames(cand)) {

  cand <- cand[
    cand$TE_family == "L1FLnI",
    ,
    drop = FALSE
  ]
}


############################################################
# Identify candidate feature column
############################################################

if ("TE_feature" %in% colnames(cand)) {

  feature_col <- "TE_feature"

} else if ("feature" %in% colnames(cand)) {

  feature_col <- "feature"

} else {

  stop(
    "Could not identify candidate TE feature column."
  )
}


############################################################
# Identify candidate cell-type column
############################################################

if ("HRCA_majorclass" %in% colnames(cand)) {

  candidate_celltype_col <- "HRCA_majorclass"

} else if ("CellType" %in% colnames(cand)) {

  candidate_celltype_col <- "CellType"

} else {

  stop(
    "Could not identify candidate cell-type column."
  )
}


############################################################
# Build candidate map
############################################################

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


############################################################
# Keep candidates actually present in expression matrix
############################################################

candidate_map <- candidate_map[
  candidate_map$TE_feature %in% rownames(expr),
  ,
  drop = FALSE
]


cat(
  "\nUnique candidate cell-type/locus pairs:",
  nrow(candidate_map),
  "\n"
)

cat(
  "Unique candidate loci:",
  length(unique(candidate_map$TE_feature)),
  "\n"
)

cat("\nCandidate map:\n")
print(candidate_map)


############################################################
# Save candidate map
############################################################

write.csv(
  candidate_map,
  file.path(
    table_dir,
    "candidate_L1FLnI_map.csv"
  ),
  row.names = FALSE
)


############################################################
# 20. Candidate locus expression
#     donor x cell type x locus
#
# For every candidate locus:
#
#   mean_expression =
#     mean normalized expression across cells
#
#   detection_rate =
#     fraction of cells with expression > 0
############################################################

cat("\n============================================\n")
cat("C1. CANDIDATE LOCUS DONOR EXPRESSION\n")
cat("============================================\n")

candidate_locus_list <- list()

counter <- 1

for (ct in unique(candidate_map$CellType)) {

  cat("\nCell type:", ct, "\n")

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

    ########################################################
    # Mean normalized expression
    ########################################################

    means <- Matrix::rowMeans(mat)


    ########################################################
    # Detection rate
    #
    # Fraction of cells from this donor/cell type
    # expressing the locus (>0)
    ########################################################

    detection <- Matrix::rowMeans(
      mat > 0
    )


    ########################################################
    # Store one row per locus
    ########################################################

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


############################################################
# Combine candidate locus table
############################################################

candidate_locus <- do.call(
  rbind,
  candidate_locus_list
)

rownames(candidate_locus) <- NULL


############################################################
# Add disease condition
############################################################

candidate_locus <- merge(
  candidate_locus,
  donor_condition,
  by = "Donor",
  all.x = TRUE
)


############################################################
# Add global CellChat metrics
############################################################

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


############################################################
# Save
############################################################

write.csv(
  candidate_locus,
  file.path(
    table_dir,
    "candidate_locus_donor_expression.csv"
  ),
  row.names = FALSE
)


############################################################
# Candidate locus QC
############################################################

cat("\nCandidate locus observations:\n")
cat("Rows:", nrow(candidate_locus), "\n")

cat("\nObservations by cell type / condition:\n")

print(
  table(
    candidate_locus$CellType,
    candidate_locus$Condition
  )
)


############################################################
# 21. Candidate burden by donor x cell type
#
# This collapses the candidate loci belonging to each
# cell type into two donor-level measurements:
#
#   candidate_expression_burden
#   candidate_detection_burden
############################################################

cat("\n============================================\n")
cat("C2. CANDIDATE BURDEN\n")
cat("============================================\n")

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
      mean(
        tmp$mean_expression,
        na.rm = TRUE
      ),

    candidate_detection_burden =
      mean(
        tmp$detection_rate,
        na.rm = TRUE
      ),

    density = tmp$density[1],
    total_strength = tmp$total_strength[1],

    stringsAsFactors = FALSE
  )
}


candidate_burden <- do.call(
  rbind,
  candidate_burden_list
)

rownames(candidate_burden) <- NULL


############################################################
# Save candidate burden
############################################################

write.csv(
  candidate_burden,
  file.path(
    table_dir,
    "candidate_burden_by_donor_celltype.csv"
  ),
  row.names = FALSE
)


############################################################
# Candidate burden QC
############################################################

cat("\nCandidate burden observations:\n")

print(
  table(
    candidate_burden$CellType,
    candidate_burden$Condition
  )
)


############################################################
# 22. Candidate burden correlations
#
# ALL + Healthy + Dry + Wet
############################################################

cat("\n============================================\n")
cat("C3. CANDIDATE BURDEN CORRELATIONS\n")
cat("============================================\n")

candidate_results <- list()

counter <- 1

candidate_predictors <- c(
  "candidate_expression_burden",
  "candidate_detection_burden"
)


for (ct in sort(
  unique(candidate_burden$CellType)
)) {

  cat("\nCell type:", ct, "\n")

  ct_data <- candidate_burden[
    candidate_burden$CellType == ct,
    ,
    drop = FALSE
  ]


  for (condition_i in conditions_to_test) {

    tmp <- subset_condition(
      ct_data,
      condition_i
    )


    for (pred in candidate_predictors) {

      for (metric in network_metrics) {

        res <- safe_correlations(
          tmp[[pred]],
          tmp[[metric]]
        )


        candidate_results[[counter]] <- data.frame(

          Analysis_level =
            "Candidate_burden",

          Condition =
            condition_i,

          CellType =
            ct,

          Feature =
            "Candidate_L1FLnI",

          Predictor =
            pred,

          Network_metric =
            metric,

          n =
            res$n,

          Inference =
            res$Inference,

          Pearson_r =
            res$Pearson_r,

          Pearson_R2 =
            res$Pearson_R2,

          Pearson_P =
            res$Pearson_P,

          Spearman_rho =
            res$Spearman_rho,

          Spearman_P =
            res$Spearman_P,

          stringsAsFactors = FALSE
        )

        counter <- counter + 1
      }
    }
  }
}


############################################################
# Combine candidate burden results
############################################################

candidate_results <- do.call(
  rbind,
  candidate_results
)


############################################################
# Add FDR
############################################################

candidate_results <- add_grouped_fdr(
  candidate_results
)


############################################################
# Save candidate burden correlations
############################################################

write.csv(
  candidate_results,
  file.path(
    table_dir,
    "candidate_burden_correlations.csv"
  ),
  row.names = FALSE
)


############################################################
# Print strongest candidate burden associations
############################################################

candidate_print <- candidate_results[
  !is.na(candidate_results$Pearson_P),
  ,
  drop = FALSE
]

candidate_print <- candidate_print[
  order(candidate_print$Pearson_P),
  ,
  drop = FALSE
]

cat("\nStrongest candidate burden results:\n")

print(
  head(
    candidate_print,
    20
  )
)


############################################################
# 23. D. LOCUS-SPECIFIC ANALYSIS
#
# This is the key analysis for direct comparison
# with the locus-edge CellChat analysis.
#
# For each:
#
#   CellType
#   L1FLnI locus
#   Condition
#
# Test:
#
#   mean_expression -> density
#   mean_expression -> total_strength
#
#   detection_rate -> density
#   detection_rate -> total_strength
#
# Pearson = primary
# Spearman = sensitivity
############################################################

cat("\n============================================\n")
cat("D. LOCUS-SPECIFIC CORRELATIONS\n")
cat("============================================\n")


############################################################
# Unique cell type / locus combinations
############################################################

locus_groups <- unique(
  candidate_locus[
    ,
    c(
      "CellType",
      "TE_feature"
    )
  ]
)

cat(
  "\nCell type / locus combinations:",
  nrow(locus_groups),
  "\n"
)


############################################################
# Predictors
############################################################

locus_predictors <- c(
  "mean_expression",
  "detection_rate"
)


############################################################
# Run correlations
############################################################

locus_results <- list()

counter <- 1


for (i in seq_len(nrow(locus_groups))) {

  ct <- locus_groups$CellType[i]
  locus <- locus_groups$TE_feature[i]


  cat(
    "\nProcessing:",
    ct,
    "|",
    locus,
    "\n"
  )


  locus_data <- candidate_locus[
    candidate_locus$CellType == ct &
    candidate_locus$TE_feature == locus,
    ,
    drop = FALSE
  ]


  ##########################################################
  # ALL + each disease condition
  ##########################################################

  for (condition_i in conditions_to_test) {

    tmp <- subset_condition(
      locus_data,
      condition_i
    )


    ########################################################
    # mean_expression + detection_rate
    ########################################################

    for (pred in locus_predictors) {


      ######################################################
      # density + total_strength
      ######################################################

      for (metric in network_metrics) {

        res <- safe_correlations(
          tmp[[pred]],
          tmp[[metric]]
        )


        locus_results[[counter]] <- data.frame(

          Analysis_level =
            "Locus_specific",

          Condition =
            condition_i,

          CellType =
            ct,

          Feature =
            locus,

          Predictor =
            pred,

          Network_metric =
            metric,

          n =
            res$n,

          Inference =
            res$Inference,

          Pearson_r =
            res$Pearson_r,

          Pearson_R2 =
            res$Pearson_R2,

          Pearson_P =
            res$Pearson_P,

          Spearman_rho =
            res$Spearman_rho,

          Spearman_P =
            res$Spearman_P,

          stringsAsFactors = FALSE
        )


        counter <- counter + 1
      }
    }
  }
}


############################################################
# Combine locus-specific results
############################################################

locus_results <- do.call(
  rbind,
  locus_results
)


############################################################
# Add FDR
############################################################

locus_results <- add_grouped_fdr(
  locus_results
)


############################################################
# Save
############################################################

write.csv(
  locus_results,
  file.path(
    table_dir,
    "locus_specific_correlations.csv"
  ),
  row.names = FALSE
)


############################################################
# 24. LOCUS-SPECIFIC QC
############################################################

cat("\n============================================\n")
cat("LOCUS-SPECIFIC QC\n")
cat("============================================\n")

cat(
  "\nTotal locus-specific tests:",
  nrow(locus_results),
  "\n"
)

cat(
  "Pearson nominal P < 0.05:",
  sum(
    locus_results$Pearson_P < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Pearson FDR < 0.05:",
  sum(
    locus_results$Pearson_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Spearman nominal P < 0.05:",
  sum(
    locus_results$Spearman_P < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Spearman FDR < 0.05:",
  sum(
    locus_results$Spearman_FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


############################################################
# Results by condition
############################################################

cat("\nTests by condition:\n")

print(
  table(
    locus_results$Condition,
    locus_results$Inference
  )
)


############################################################
# 25. Strongest locus-specific results
############################################################

locus_print <- locus_results[
  !is.na(locus_results$Pearson_P),
  ,
  drop = FALSE
]

locus_print <- locus_print[
  order(locus_print$Pearson_P),
  ,
  drop = FALSE
]


cat("\n============================================\n")
cat("TOP LOCUS-SPECIFIC PEARSON RESULTS\n")
cat("============================================\n")

print(
  head(
    locus_print[
      ,
      c(
        "Condition",
        "CellType",
        "Feature",
        "Predictor",
        "Network_metric",
        "n",
        "Inference",
        "Pearson_r",
        "Pearson_R2",
        "Pearson_P",
        "Pearson_FDR",
        "Spearman_rho",
        "Spearman_P",
        "Spearman_FDR"
      )
    ],
    30
  ),
  row.names = FALSE
)


############################################################
# 26. SPECIAL QC
#
# These are the two loci identified in the original
# ALL-donor analysis.
#
# We specifically check that Condition = ALL reproduces
# the previous results.
############################################################

cat("\n============================================\n")
cat("QC: Xq23ab\n")
cat("============================================\n")

qc_xq23 <- locus_results[
  locus_results$Condition == "ALL" &
  locus_results$Feature ==
    "TE-L1FLnI-Xq23ab",
  ,
  drop = FALSE
]

print(
  qc_xq23[
    ,
    c(
      "Condition",
      "CellType",
      "Feature",
      "Predictor",
      "Network_metric",
      "n",
      "Pearson_r",
      "Pearson_R2",
      "Pearson_P",
      "Pearson_FDR",
      "Spearman_rho",
      "Spearman_P"
    )
  ],
  row.names = FALSE
)


cat("\n============================================\n")
cat("QC: 21q22.13b\n")
cat("============================================\n")

qc_21q <- locus_results[
  locus_results$Condition == "ALL" &
  locus_results$Feature ==
    "TE-L1FLnI-21q22.13b",
  ,
  drop = FALSE
]

print(
  qc_21q[
    ,
    c(
      "Condition",
      "CellType",
      "Feature",
      "Predictor",
      "Network_metric",
      "n",
      "Pearson_r",
      "Pearson_R2",
      "Pearson_P",
      "Pearson_FDR",
      "Spearman_rho",
      "Spearman_P"
    )
  ],
  row.names = FALSE
)


############################################################
# 27. QC EXPECTED VALUES FROM ORIGINAL ANALYSIS
############################################################

cat("\n============================================\n")
cat("EXPECTED ORIGINAL ALL-DONOR VALUES\n")
cat("============================================\n")

cat(
  "\nAstrocyte | Xq23ab | detection_rate",
  "| total_strength\n"
)

cat(
  "Expected approximately:\n",
  "n = 16\n",
  "Pearson r = 0.6024035\n",
  "P = 0.0135293\n"
)

cat(
  "\nRPE | 21q22.13b | detection_rate",
  "| density\n"
)

cat(
  "Expected approximately:\n",
  "n = 13\n",
  "Pearson r = 0.6288240\n",
  "P = 0.0213225\n"
)


############################################################
# END OF BLOCK 3
############################################################
############################################################
# 28. MASTER RESULTS TABLE
#
# Combine:
#   A. Family_global
#   B. Family_celltype
#   C. Candidate_burden
#   D. Locus_specific
############################################################

cat("\n============================================\n")
cat("BUILDING MASTER RESULTS TABLE\n")
cat("============================================\n")


master <- rbind(
  global_results,
  celltype_results,
  candidate_results,
  locus_results
)


############################################################
# Order master table
############################################################

master <- master[
  order(
    master$Analysis_level,
    master$Condition,
    master$Pearson_P
  ),
  ,
  drop = FALSE
]

rownames(master) <- NULL


############################################################
# Save master table
############################################################

write.csv(
  master,
  file.path(
    table_dir,
    "MASTER_aligned_correlations.csv"
  ),
  row.names = FALSE
)


cat(
  "\nTotal tests in MASTER:",
  nrow(master),
  "\n"
)


############################################################
# 29. SUMMARY TABLES
############################################################

cat("\n============================================\n")
cat("CREATING SUMMARY TABLES\n")
cat("============================================\n")


############################################################
# Pearson nominal P < 0.05
#
# Keep inferential results only.
############################################################

pearson_nominal <- master[
  master$Inference == "INFERENTIAL" &
  !is.na(master$Pearson_P) &
  master$Pearson_P < 0.05,
  ,
  drop = FALSE
]

pearson_nominal <- pearson_nominal[
  order(
    pearson_nominal$Pearson_P
  ),
  ,
  drop = FALSE
]

write.csv(
  pearson_nominal,
  file.path(
    table_dir,
    "Pearson_nominal_P005.csv"
  ),
  row.names = FALSE
)


############################################################
# Pearson FDR < 0.05
############################################################

pearson_fdr <- master[
  master$Inference == "INFERENTIAL" &
  !is.na(master$Pearson_FDR) &
  master$Pearson_FDR < 0.05,
  ,
  drop = FALSE
]

pearson_fdr <- pearson_fdr[
  order(
    pearson_fdr$Pearson_FDR
  ),
  ,
  drop = FALSE
]

write.csv(
  pearson_fdr,
  file.path(
    table_dir,
    "Pearson_FDR005.csv"
  ),
  row.names = FALSE
)


############################################################
# Spearman nominal P < 0.05
############################################################

spearman_nominal <- master[
  master$Inference == "INFERENTIAL" &
  !is.na(master$Spearman_P) &
  master$Spearman_P < 0.05,
  ,
  drop = FALSE
]

spearman_nominal <- spearman_nominal[
  order(
    spearman_nominal$Spearman_P
  ),
  ,
  drop = FALSE
]

write.csv(
  spearman_nominal,
  file.path(
    table_dir,
    "Spearman_nominal_P005.csv"
  ),
  row.names = FALSE
)


############################################################
# Spearman FDR < 0.05
############################################################

spearman_fdr <- master[
  master$Inference == "INFERENTIAL" &
  !is.na(master$Spearman_FDR) &
  master$Spearman_FDR < 0.05,
  ,
  drop = FALSE
]

spearman_fdr <- spearman_fdr[
  order(
    spearman_fdr$Spearman_FDR
  ),
  ,
  drop = FALSE
]

write.csv(
  spearman_fdr,
  file.path(
    table_dir,
    "Spearman_FDR005.csv"
  ),
  row.names = FALSE
)


############################################################
# Descriptive-only results
#
# Important especially for Dry AMD when n = 3.
############################################################

descriptive_only <- master[
  master$Inference == "DESCRIPTIVE_ONLY",
  ,
  drop = FALSE
]

write.csv(
  descriptive_only,
  file.path(
    table_dir,
    "DESCRIPTIVE_ONLY_results.csv"
  ),
  row.names = FALSE
)


############################################################
# 30. SUMMARY COUNTS
############################################################

cat("\n============================================\n")
cat("RESULT COUNTS\n")
cat("============================================\n")

cat(
  "\nPearson nominal P < 0.05:",
  nrow(pearson_nominal),
  "\n"
)

cat(
  "Pearson FDR < 0.05:",
  nrow(pearson_fdr),
  "\n"
)

cat(
  "Spearman nominal P < 0.05:",
  nrow(spearman_nominal),
  "\n"
)

cat(
  "Spearman FDR < 0.05:",
  nrow(spearman_fdr),
  "\n"
)

cat(
  "Descriptive-only tests:",
  nrow(descriptive_only),
  "\n"
)


############################################################
# Results by analysis level
############################################################

cat("\nPearson nominal results by analysis level:\n")

print(
  table(
    pearson_nominal$Analysis_level
  )
)


############################################################
# Results by condition
############################################################

cat("\nPearson nominal results by condition:\n")

print(
  table(
    pearson_nominal$Condition
  )
)


############################################################
# Pearson FDR results by condition
############################################################

cat("\nPearson FDR results by condition:\n")

print(
  table(
    pearson_fdr$Condition
  )
)


############################################################
# 31. TOP PEARSON RESULTS
############################################################

cat("\n============================================\n")
cat("TOP PEARSON RESULTS\n")
cat("============================================\n")


top_pearson <- master[
  master$Inference == "INFERENTIAL" &
  !is.na(master$Pearson_P),
  ,
  drop = FALSE
]

top_pearson <- top_pearson[
  order(top_pearson$Pearson_P),
  ,
  drop = FALSE
]


print(
  head(
    top_pearson[
      ,
      c(
        "Analysis_level",
        "Condition",
        "CellType",
        "Feature",
        "Predictor",
        "Network_metric",
        "n",
        "Pearson_r",
        "Pearson_R2",
        "Pearson_P",
        "Pearson_FDR",
        "Spearman_rho",
        "Spearman_P",
        "Spearman_FDR"
      )
    ],
    30
  ),
  row.names = FALSE
)


############################################################
# 32. TOP SPEARMAN RESULTS
############################################################

cat("\n============================================\n")
cat("TOP SPEARMAN RESULTS\n")
cat("============================================\n")


top_spearman <- master[
  master$Inference == "INFERENTIAL" &
  !is.na(master$Spearman_P),
  ,
  drop = FALSE
]

top_spearman <- top_spearman[
  order(top_spearman$Spearman_P),
  ,
  drop = FALSE
]


print(
  head(
    top_spearman[
      ,
      c(
        "Analysis_level",
        "Condition",
        "CellType",
        "Feature",
        "Predictor",
        "Network_metric",
        "n",
        "Pearson_r",
        "Pearson_P",
        "Spearman_rho",
        "Spearman_P",
        "Spearman_FDR"
      )
    ],
    30
  ),
  row.names = FALSE
)


############################################################
# 33. LOCUS-SPECIFIC HEATMAP FUNCTION
#
# Separate figures for:
#
#   mean_expression
#   detection_rate
#
# and for:
#
#   ALL
#   Healthy Control
#   Dry AMD
#   Wet AMD
############################################################

make_locus_heatmap <- function(
    dat,
    predictor_i,
    condition_i,
    filename_prefix) {


  ##########################################################
  # Subset
  ##########################################################

  plot_dat <- dat[
    dat$Predictor == predictor_i &
    dat$Condition == condition_i,
    ,
    drop = FALSE
  ]


  if (nrow(plot_dat) == 0) {

    cat(
      "\nNo data for heatmap:",
      predictor_i,
      condition_i,
      "\n"
    )

    return(NULL)
  }


  ##########################################################
  # Clean locus labels
  ##########################################################

  plot_dat$Feature_clean <- sub(
    "^TE-L1FLnI-",
    "",
    plot_dat$Feature
  )


  ##########################################################
  # Network metric labels
  ##########################################################

  plot_dat$Network_metric_label <- ifelse(
    plot_dat$Network_metric == "density",
    "Network density",
    "Total strength"
  )


  ##########################################################
  # Label
  ##########################################################

  plot_dat$plot_label <- ifelse(
    is.na(plot_dat$Pearson_r),
    "NA",
    paste0(
      "r=",
      sprintf(
        "%.2f",
        plot_dat$Pearson_r
      ),
      "\nP=",
      format.pval(
        plot_dat$Pearson_P,
        digits = 2,
        eps = 0.001
      )
    )
  )


  ##########################################################
  # Plot
  ##########################################################

  p <- ggplot(
    plot_dat,
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
        label = plot_label
      ),
      size = 2.7
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
      title = paste0(
        "L1FLnI locus-specific associations — ",
        condition_i
      ),
      subtitle = paste0(
        "Predictor: ",
        predictor_i
      ),
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


  ##########################################################
  # Safe condition name for filename
  ##########################################################

  condition_file <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    condition_i
  )


  ##########################################################
  # Save PDF
  ##########################################################

  ggsave(
    file.path(
      figure_dir,
      paste0(
        filename_prefix,
        "_",
        condition_file,
        ".pdf"
      )
    ),
    p,
    width = 8,
    height = 8
  )


  ##########################################################
  # Save PNG
  ##########################################################

  ggsave(
    file.path(
      figure_dir,
      paste0(
        filename_prefix,
        "_",
        condition_file,
        ".png"
      )
    ),
    p,
    width = 8,
    height = 8,
    dpi = 400
  )


  return(p)
}


############################################################
# 34. CREATE LOCUS-SPECIFIC HEATMAPS
############################################################

cat("\n============================================\n")
cat("CREATING LOCUS-SPECIFIC HEATMAPS\n")
cat("============================================\n")


for (condition_i in conditions_to_test) {

  ##########################################################
  # Mean expression
  ##########################################################

  make_locus_heatmap(
    locus_results,
    predictor_i = "mean_expression",
    condition_i = condition_i,
    filename_prefix =
      "LocusSpecific_meanExpression"
  )


  ##########################################################
  # Detection rate
  ##########################################################

  make_locus_heatmap(
    locus_results,
    predictor_i = "detection_rate",
    condition_i = condition_i,
    filename_prefix =
      "LocusSpecific_detectionRate"
  )
}


############################################################
# 35. FAMILY CELL-TYPE HEATMAP FUNCTION
############################################################

make_family_heatmap <- function(
    dat,
    predictor_i,
    condition_i,
    filename_prefix) {


  plot_dat <- dat[
    dat$Predictor == predictor_i &
    dat$Condition == condition_i,
    ,
    drop = FALSE
  ]


  if (nrow(plot_dat) == 0) {
    return(NULL)
  }


  plot_dat$Network_metric_label <- ifelse(
    plot_dat$Network_metric == "density",
    "Network density",
    "Total strength"
  )


  plot_dat$plot_label <- ifelse(
    is.na(plot_dat$Pearson_r),
    "NA",
    paste0(
      "r=",
      sprintf(
        "%.2f",
        plot_dat$Pearson_r
      ),
      "\nP=",
      format.pval(
        plot_dat$Pearson_P,
        digits = 2,
        eps = 0.001
      )
    )
  )


  p <- ggplot(
    plot_dat,
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
        label = plot_label
      ),
      size = 3
    ) +

    scale_fill_gradient2(
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson r"
    ) +

    labs(
      title = paste0(
        "Family-wide L1FLnI associations — ",
        condition_i
      ),
      subtitle = paste0(
        "Predictor: ",
        predictor_i
      ),
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


  condition_file <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    condition_i
  )


  ggsave(
    file.path(
      figure_dir,
      paste0(
        filename_prefix,
        "_",
        condition_file,
        ".pdf"
      )
    ),
    p,
    width = 7,
    height = 6
  )


  ggsave(
    file.path(
      figure_dir,
      paste0(
        filename_prefix,
        "_",
        condition_file,
        ".png"
      )
    ),
    p,
    width = 7,
    height = 6,
    dpi = 400
  )


  return(p)
}


############################################################
# 36. CREATE FAMILY-WIDE HEATMAPS
############################################################

cat("\n============================================\n")
cat("CREATING FAMILY-WIDE HEATMAPS\n")
cat("============================================\n")


for (condition_i in conditions_to_test) {


  make_family_heatmap(
    celltype_results,
    predictor_i = "L1_expression_burden",
    condition_i = condition_i,
    filename_prefix =
      "FamilyCellType_expressionBurden"
  )


  make_family_heatmap(
    celltype_results,
    predictor_i = "detection_breadth",
    condition_i = condition_i,
    filename_prefix =
      "FamilyCellType_detectionBreadth"
  )
}


############################################################
# 37. SPECIAL CANDIDATE LOCI TABLE
#
# Pull Xq23ab and 21q22.13b across ALL conditions.
#
# This will make comparison with edge-level analysis easier.
############################################################

special_loci <- c(
  "TE-L1FLnI-Xq23ab",
  "TE-L1FLnI-21q22.13b"
)


special_results <- locus_results[
  locus_results$Feature %in% special_loci,
  ,
  drop = FALSE
]


special_results <- special_results[
  order(
    special_results$Feature,
    special_results$Condition,
    special_results$Predictor,
    special_results$Network_metric
  ),
  ,
  drop = FALSE
]


write.csv(
  special_results,
  file.path(
    table_dir,
    "SPECIAL_Xq23ab_21q22.13b_all_conditions.csv"
  ),
  row.names = FALSE
)


cat("\n============================================\n")
cat("SPECIAL LOCI RESULTS\n")
cat("============================================\n")

print(
  special_results[
    ,
    c(
      "Condition",
      "CellType",
      "Feature",
      "Predictor",
      "Network_metric",
      "n",
      "Inference",
      "Pearson_r",
      "Pearson_P",
      "Pearson_FDR",
      "Spearman_rho",
      "Spearman_P"
    )
  ],
  row.names = FALSE
)


############################################################
# 38. FINAL QC:
#     Reproduce original ALL-donor results
############################################################

cat("\n\n============================================\n")
cat("FINAL QC — ORIGINAL RESULTS REPRODUCTION\n")
cat("============================================\n")


############################################################
# Xq23ab:
# Astrocyte
# detection_rate
# total_strength
############################################################

qc1 <- locus_results[
  locus_results$Condition == "ALL" &
  locus_results$CellType == "Astrocyte" &
  locus_results$Feature ==
    "TE-L1FLnI-Xq23ab" &
  locus_results$Predictor ==
    "detection_rate" &
  locus_results$Network_metric ==
    "total_strength",
  ,
  drop = FALSE
]


cat("\nXq23ab expected:\n")
cat(
  "n = 16 | r ~ 0.6024035 | P ~ 0.0135293\n"
)

cat("Observed:\n")

print(
  qc1[
    ,
    c(
      "n",
      "Pearson_r",
      "Pearson_R2",
      "Pearson_P",
      "Pearson_FDR",
      "Spearman_rho",
      "Spearman_P"
    )
  ],
  row.names = FALSE
)


############################################################
# 21q22.13b:
# RPE
# detection_rate
# density
############################################################

qc2 <- locus_results[
  locus_results$Condition == "ALL" &
  locus_results$CellType == "RPE" &
  locus_results$Feature ==
    "TE-L1FLnI-21q22.13b" &
  locus_results$Predictor ==
    "detection_rate" &
  locus_results$Network_metric ==
    "density",
  ,
  drop = FALSE
]


cat("\n21q22.13b expected:\n")
cat(
  "n = 13 | r ~ 0.6288240 | P ~ 0.0213225\n"
)

cat("Observed:\n")

print(
  qc2[
    ,
    c(
      "n",
      "Pearson_r",
      "Pearson_R2",
      "Pearson_P",
      "Pearson_FDR",
      "Spearman_rho",
      "Spearman_P"
    )
  ],
  row.names = FALSE
)


############################################################
# 39. Automatic reproduction check
############################################################

qc1_pass <- FALSE
qc2_pass <- FALSE


if (nrow(qc1) == 1) {

  qc1_pass <- (
    qc1$n == 16 &&
    abs(
      qc1$Pearson_r -
      0.602403466084913
    ) < 1e-6 &&
    abs(
      qc1$Pearson_P -
      0.0135292648243747
    ) < 1e-6
  )
}


if (nrow(qc2) == 1) {

  qc2_pass <- (
    qc2$n == 13 &&
    abs(
      qc2$Pearson_r -
      0.628823965166148
    ) < 1e-6 &&
    abs(
      qc2$Pearson_P -
      0.0213225255223915
    ) < 1e-6
  )
}


cat("\nQC Xq23ab:", qc1_pass, "\n")
cat("QC 21q22.13b:", qc2_pass, "\n")


if (qc1_pass && qc2_pass) {

  cat(
    "\nSUCCESS: ALL-donor analysis reproduces",
    "the original results.\n"
  )

} else {

  cat(
    "\nWARNING: At least one original result",
    "was not reproduced exactly.\n"
  )

  cat(
    "Do NOT interpret condition-specific",
    "results until this is checked.\n"
  )
}


############################################################
# 40. FINAL SUMMARY
############################################################

cat("\n\n============================================\n")
cat("ANALYSIS COMPLETED\n")
cat("============================================\n")


cat(
  "\nDonors:",
  length(unique(meta$Donor)),
  "\n"
)

cat(
  "L1FLnI loci:",
  length(l1_features),
  "\n"
)

cat(
  "Cell types:",
  length(unique(meta$CellType)),
  "\n"
)

cat(
  "Candidate loci:",
  length(unique(candidate_map$TE_feature)),
  "\n"
)


cat("\nConditions analyzed:\n")

for (condition_i in conditions_to_test) {

  if (condition_i == "ALL") {

    n_donors_condition <- length(
      unique(meta$Donor)
    )

  } else {

    n_donors_condition <- length(
      unique(
        donor_condition$Donor[
          donor_condition$Condition ==
            condition_i
        ]
      )
    )
  }

  cat(
    "  -",
    condition_i,
    ":",
    n_donors_condition,
    "donors\n"
  )
}


cat("\nNetwork outcomes:\n")
cat("  - density\n")
cat("  - total_strength\n")


cat("\nPredictors include:\n")
cat("  - L1_expression_burden\n")
cat("  - detection_breadth\n")
cat("  - candidate_expression_burden\n")
cat("  - candidate_detection_burden\n")
cat("  - mean_expression\n")
cat("  - detection_rate\n")


cat("\nStatistical framework:\n")
cat("  - Pearson: primary\n")
cat("  - Spearman: sensitivity\n")
cat("  - BH-FDR correction\n")
cat("  - n = DONORS\n")
cat("  - n < 4 = DESCRIPTIVE_ONLY\n")


cat("\nOutput directory:\n")
cat(outdir, "\n")


cat("\nMain tables:\n")

cat(
  "  ",
  file.path(
    table_dir,
    "MASTER_aligned_correlations.csv"
  ),
  "\n"
)

cat(
  "  ",
  file.path(
    table_dir,
    "Pearson_nominal_P005.csv"
  ),
  "\n"
)

cat(
  "  ",
  file.path(
    table_dir,
    "Pearson_FDR005.csv"
  ),
  "\n"
)

cat(
  "  ",
  file.path(
    table_dir,
    "Spearman_nominal_P005.csv"
  ),
  "\n"
)

cat(
  "  ",
  file.path(
    table_dir,
    "Spearman_FDR005.csv"
  ),
  "\n"
)

cat(
  "  ",
  file.path(
    table_dir,
    "DESCRIPTIVE_ONLY_results.csv"
  ),
  "\n"
)

cat(
  "  ",
  file.path(
    table_dir,
    "SPECIAL_Xq23ab_21q22.13b_all_conditions.csv"
  ),
  "\n"
)


cat("\nFigures directory:\n")
cat(figure_dir, "\n")


cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n")
