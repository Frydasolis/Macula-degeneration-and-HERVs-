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
# Pearson:
#   donor locus expression ~ donor CellChat edge weight
#
# Spearman also calculated as sensitivity analysis.
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
# LOAD DATA
############################################################

cat("\nLoading Seurat...\n")

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
# IMPORTANT:
# Expression is calculated across ALL retinal cells
# belonging to each donor.
#
# This matches the global locus-specific analysis.
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

  mean_expr <- Matrix::rowMeans(mat_d)

  detection <- Matrix::rowMeans(
    mat_d > 0
  )

  condition <- unique(
    meta$disease_state[
      meta$orig.ident == d
    ]
  )

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

write.csv(
  locus_expression,
  file.path(
    out_dir,
    "Donor_Global_Locus_Expression.csv"
  ),
  row.names = FALSE
)

############################################################
# LOAD DONOR CELLCHAT EDGES
############################################################

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

############################################################
# JOIN EXPRESSION WITH EVERY EDGE
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

############################################################
# SAFE CORRELATION FUNCTION
############################################################

safe_cor <- function(x, y, method = "pearson") {

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
      method = method,
      exact = FALSE
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
# ANALYSIS FUNCTION
############################################################

run_correlations <- function(
  predictor,
  method
) {

  cat(
    "\nRunning ",
    method,
    ": ",
    predictor,
    "\n",
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
# RUN
############################################################

pearson_mean <- run_correlations(
  "mean_expression",
  "pearson"
)

pearson_detect <- run_correlations(
  "detection_rate",
  "pearson"
)

spearman_mean <- run_correlations(
  "mean_expression",
  "spearman"
)

spearman_detect <- run_correlations(
  "detection_rate",
  "spearman"
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
# Correct within:
# condition x predictor x method
#
# across all loci and all sender-receiver edges
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
# DESCRIPTIVE FLAG FOR DRY
############################################################

results <- results %>%
  mutate(
    inference =
      ifelse(
        condition == "Dry AMD",
        "DESCRIPTIVE_n3",
        "INFERENTIAL"
      )
  )

############################################################
# SAVE
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
# PRIMARY RESULTS
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
# SIGNIFICANT / NOMINAL TABLES
############################################################

nominal <- primary %>%
  filter(
    condition != "Dry AMD",
    !is.na(p),
    p < 0.05
  )

fdr_sig <- primary %>%
  filter(
    condition != "Dry AMD",
    !is.na(FDR),
    FDR < 0.05
  )

write.csv(
  nominal,
  file.path(
    out_dir,
    "Locus_Edge_Pearson_Nominal_P005.csv"
  ),
  row.names = FALSE
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
# SUMMARY
############################################################

cat("\n============================================\n")
cat("LOCUS x CELLCHAT EDGE ANALYSIS COMPLETE\n")
cat("============================================\n")

cat("\nMaster rows:", nrow(master), "\n")
cat("All correlation rows:", nrow(results), "\n")

cat("\nPRIMARY PEARSON mean-expression results:\n")

print(
  primary %>%
    group_by(condition) %>%
    summarise(
      tests = n(),
      valid = sum(!is.na(r)),
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
)

cat("\nTop 30 non-Dry associations:\n")

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

cat("\nOutput:\n")
cat(out_dir, "\n")

