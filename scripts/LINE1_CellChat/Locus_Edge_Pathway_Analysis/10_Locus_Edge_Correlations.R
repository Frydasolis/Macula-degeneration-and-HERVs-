############################################################
# 10_Locus_Edge_Correlations.R
#
# LOCUS-SPECIFIC association with individual CellChat edges
#
# Unit of inference = DONOR
#
# For each:
#   locus x condition x sender x receiver
#
# Primary analysis:
#   Pearson correlation between:
#     donor-level global locus expression
#     donor-level CellChat edge weight
#
# Sensitivity analysis:
#   Spearman correlation
#
# IMPORTANT:
# - Cells are NOT treated as independent observations.
# - Dry AMD (n = 3 donors) is considered descriptive.
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
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
  "results/SRP413248/CellChat/donor_level_symbol/donor_figure_analysis"
)

candidate_file <- file.path(
  root,
  "results/SRP413248/Seurat/downstream/TE_analysis/Figure4",
  "LINE1_CellChat_candidates.csv"
)

out_dir <- file.path(
  network_dir,
  "locus_edge_correlations"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# LOAD SEURAT DATA
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

cat("Cells:", ncol(seu), "\n")
cat("Features:", nrow(expr), "\n")
cat("Donors:", length(unique(meta$orig.ident)), "\n")

############################################################
# CANDIDATE LOCI
############################################################

cand <- read.csv(
  candidate_file,
  stringsAsFactors = FALSE
)

loci <- unique(cand$TE_feature)

loci <- loci[
  loci %in% rownames(expr)
]

cat("\nCandidate loci found:", length(loci), "\n")

if (length(loci) != 12) {
  warning(
    "Expected 12 unique loci; found ",
    length(loci)
  )
}

print(loci)

############################################################
# DONOR-LEVEL GLOBAL LOCUS EXPRESSION
#
# For each donor and each candidate locus:
#
# mean_expression =
#   mean normalized expression across ALL retinal cells
#   belonging to that donor
#
# detection_rate =
#   proportion of donor cells with expression > 0
#
# IMPORTANT:
# This is GLOBAL donor-level locus expression.
# It is NOT restricted to RPE, astrocytes, Müller glia, etc.
############################################################

donors <- unique(meta$orig.ident)

expression_list <- list()

for (d in donors) {

  cells_d <- meta$cell[
    meta$orig.ident == d
  ]

  mat_d <- expr[
    loci,
    cells_d,
    drop = FALSE
  ]

  mean_expr <- Matrix::rowMeans(
    mat_d
  )

  detection <- Matrix::rowMeans(
    mat_d > 0
  )

  condition <- unique(
    meta$disease_state[
      meta$orig.ident == d
    ]
  )

  if (length(condition) != 1) {
    warning(
      "Donor ",
      d,
      " has ",
      length(condition),
      " disease-state labels."
    )
  }

  expression_list[[d]] <- data.frame(
    donor = d,
    condition = condition[1],
    locus = loci,
    mean_expression = as.numeric(mean_expr),
    detection_rate = as.numeric(detection),
    n_cells = length(cells_d),
    stringsAsFactors = FALSE
  )
}

locus_expression <- bind_rows(
  expression_list
)

############################################################
# SAVE DONOR-LEVEL LOCUS EXPRESSION
############################################################

write.csv(
  locus_expression,
  file.path(
    out_dir,
    "Donor_Global_Locus_Expression.csv"
  ),
  row.names = FALSE
)

cat(
  "\nDonor-locus expression rows:",
  nrow(locus_expression),
  "\n"
)

############################################################
# LOAD DONOR CELLCHAT EDGES
############################################################

cat("\nLoading donor CellChat edges...\n")

edges <- read.csv(
  file.path(
    network_dir,
    "donor_edges_complete.csv"
  ),
  stringsAsFactors = FALSE
)

cat("\nEdge columns:\n")
print(names(edges))

############################################################
# KEEP TESTABLE EDGES
#
# pair_testable = TRUE:
# both sender and receiver cell populations were available
# in that donor.
#
# If a cell type is absent, the edge is NOT converted to zero.
############################################################

edges2 <- edges %>%
  filter(
    pair_testable,
    !is.na(weight)
  ) %>%
  select(
    donor,
    condition,
    sender,
    receiver,
    weight
  )

cat(
  "\nTestable donor-edge observations:",
  nrow(edges2),
  "\n"
)

############################################################
# JOIN LOCUS EXPRESSION WITH CELLCHAT EDGES
#
# Every donor-locus observation is associated with every
# testable sender -> receiver CellChat edge in that donor.
############################################################

master <- locus_expression %>%
  inner_join(
    edges2,
    by = c(
      "donor",
      "condition"
    )
  )

write.csv(
  master,
  file.path(
    out_dir,
    "MASTER_Locus_Donor_CellChatEdge.csv"
  ),
  row.names = FALSE
)

cat(
  "Master donor-locus-edge rows:",
  nrow(master),
  "\n"
)

############################################################
# SAFE CORRELATION FUNCTION
#
# Pearson:
#   standard Pearson cor.test()
#
# Spearman:
#   exact = FALSE because ties are expected.
############################################################

safe_cor <- function(
  x,
  y,
  method = "pearson"
) {

  keep <- complete.cases(
    x,
    y
  )

  x <- x[keep]
  y <- y[keep]

  n <- length(x)

  ##########################################################
  # Minimum requirements
  ##########################################################

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

  ##########################################################
  # PEARSON
  ##########################################################

  if (method == "pearson") {

    test <- tryCatch(
      cor.test(
        x,
        y,
        method = "pearson"
      ),
      error = function(e) NULL
    )

  ##########################################################
  # SPEARMAN
  ##########################################################

  } else if (method == "spearman") {

    test <- tryCatch(
      cor.test(
        x,
        y,
        method = "spearman",
        exact = FALSE
      ),
      error = function(e) NULL
    )

  } else {

    stop(
      "Unsupported correlation method: ",
      method
    )
  }

  ##########################################################
  # Failed test
  ##########################################################

  if (is.null(test)) {

    return(
      data.frame(
        n = n,
        r = NA_real_,
        p = NA_real_
      )
    )
  }

  ##########################################################
  # RETURN
  ##########################################################

  data.frame(
    n = n,
    r = as.numeric(test$estimate),
    p = test$p.value
  )
}

############################################################
# ANALYSIS FUNCTION
############################################################

run_correlations <- function(
  predictor,
  method
) {

  cat(
    "\nRunning ",
    method,
    " correlations using ",
    predictor,
    "...\n",
    sep = ""
  )

  split_data <- split(
    master,
    interaction(
      master$locus,
      master$condition,
      master$sender,
      master$receiver,
      drop = TRUE
    )
  )

  res <- lapply(
    split_data,
    function(df) {

      stat <- safe_cor(
        df[[predictor]],
        df$weight,
        method = method
      )

      data.frame(
        locus = df$locus[1],
        condition = df$condition[1],
        sender = df$sender[1],
        receiver = df$receiver[1],
        predictor = predictor,
        method = method,
        n = stat$n,
        r = stat$r,
        p = stat$p,
        stringsAsFactors = FALSE
      )
    }
  )

  bind_rows(res)
}

############################################################
# RUN ALL CORRELATIONS
############################################################

pearson_mean <- run_correlations(
  predictor = "mean_expression",
  method = "pearson"
)

pearson_detect <- run_correlations(
  predictor = "detection_rate",
  method = "pearson"
)

spearman_mean <- run_correlations(
  predictor = "mean_expression",
  method = "spearman"
)

spearman_detect <- run_correlations(
  predictor = "detection_rate",
  method = "spearman"
)

results <- bind_rows(
  pearson_mean,
  pearson_detect,
  spearman_mean,
  spearman_detect
)

############################################################
# BH FDR
#
# Multiple-testing correction is performed separately
# within:
#
#   condition x predictor x method
#
# The correction therefore includes all:
#
#   loci x sender-receiver edges
#
# belonging to that analysis family.
############################################################

results <- results %>%
  group_by(
    condition,
    predictor,
    method
  ) %>%
  mutate(
    FDR = p.adjust(
      p,
      method = "BH"
    )
  ) %>%
  ungroup()

############################################################
# INFERENCE FLAG
#
# Dry AMD contains only n = 3 donors.
# Correlations are therefore descriptive/exploratory.
############################################################

results <- results %>%
  mutate(
    inference = case_when(
      condition == "Dry AMD" ~ "DESCRIPTIVE_n3",
      TRUE ~ "INFERENTIAL"
    )
  )

############################################################
# SAVE ALL RESULTS
############################################################

write.csv(
  results,
  file.path(
    out_dir,
    "Locus_Edge_Correlations_ALL.csv"
  ),
  row.names = FALSE
)

############################################################
# PRIMARY ANALYSIS
#
# Pearson correlation using mean normalized expression.
############################################################

primary <- results %>%
  filter(
    method == "pearson",
    predictor == "mean_expression"
  ) %>%
  arrange(
    condition,
    locus,
    desc(abs(r))
  )

write.csv(
  primary,
  file.path(
    out_dir,
    "Locus_Edge_Correlations_PEARSON_PRIMARY.csv"
  ),
  row.names = FALSE
)

############################################################
# NOMINAL ASSOCIATIONS
#
# Dry AMD excluded from inferential result table.
############################################################

nominal <- primary %>%
  filter(
    condition != "Dry AMD",
    !is.na(p),
    p < 0.05
  ) %>%
  arrange(
    FDR,
    p
  )

write.csv(
  nominal,
  file.path(
    out_dir,
    "Locus_Edge_Pearson_Nominal_P005.csv"
  ),
  row.names = FALSE
)

############################################################
# FDR-SIGNIFICANT ASSOCIATIONS
############################################################

fdr_sig <- primary %>%
  filter(
    condition != "Dry AMD",
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  arrange(
    FDR,
    p
  )

write.csv(
  fdr_sig,
  file.path(
    out_dir,
    "Locus_Edge_Pearson_FDR005.csv"
  ),
  row.names = FALSE
)

############################################################
# ADD EDGE LABEL
#
# Useful for downstream heatmaps / network comparison.
############################################################

fdr_sig_edges <- fdr_sig %>%
  mutate(
    edge = paste(
      sender,
      receiver,
      sep = " -> "
    ),
    direction = case_when(
      r > 0 ~ "POSITIVE",
      r < 0 ~ "NEGATIVE",
      TRUE ~ "ZERO"
    )
  ) %>%
  select(
    condition,
    locus,
    sender,
    receiver,
    edge,
    n,
    r,
    p,
    FDR,
    direction
  )

write.csv(
  fdr_sig_edges,
  file.path(
    out_dir,
    "Locus_Edge_Pearson_FDR005_ANNOTATED.csv"
  ),
  row.names = FALSE
)

############################################################
# SUMMARY BY CONDITION
############################################################

summary_condition <- primary %>%
  group_by(
    condition
  ) %>%
  summarise(
    tests = n(),
    valid = sum(!is.na(r)),
    positive = sum(
      r > 0,
      na.rm = TRUE
    ),
    negative = sum(
      r < 0,
      na.rm = TRUE
    ),
    nominal_P005 = sum(
      p < 0.05,
      na.rm = TRUE
    ),
    FDR005 = sum(
      FDR < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

############################################################
# SUMMARY OF FDR EDGES BY LOCUS
############################################################

summary_fdr <- fdr_sig %>%
  group_by(
    condition,
    locus
  ) %>%
  summarise(
    n_FDR = n(),
    positive = sum(r > 0),
    negative = sum(r < 0),
    .groups = "drop"
  ) %>%
  arrange(
    condition,
    desc(n_FDR)
  )

write.csv(
  summary_condition,
  file.path(
    out_dir,
    "Locus_Edge_Primary_Summary_ByCondition.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_fdr,
  file.path(
    out_dir,
    "Locus_Edge_FDR_Summary_ByLocus.csv"
  ),
  row.names = FALSE
)

############################################################
# FINAL CONSOLE OUTPUT
############################################################

cat("\n============================================\n")
cat("LOCUS x CELLCHAT EDGE ANALYSIS COMPLETE\n")
cat("============================================\n")

cat(
  "\nMaster rows:",
  nrow(master),
  "\n"
)

cat(
  "All correlation rows:",
  nrow(results),
  "\n"
)

cat(
  "\nPRIMARY PEARSON mean-expression results:\n"
)

print(
  summary_condition,
  row.names = FALSE
)

cat(
  "\nFDR-significant edges by locus:\n"
)

print(
  summary_fdr,
  row.names = FALSE
)

cat(
  "\nTop 30 non-Dry associations:\n"
)

print(
  primary %>%
    filter(
      condition != "Dry AMD",
      !is.na(r)
    ) %>%
    arrange(p) %>%
    select(
      locus,
      condition,
      sender,
      receiver,
      n,
      r,
      p,
      FDR
    ) %>%
    head(30),
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(out_dir, "\n")
