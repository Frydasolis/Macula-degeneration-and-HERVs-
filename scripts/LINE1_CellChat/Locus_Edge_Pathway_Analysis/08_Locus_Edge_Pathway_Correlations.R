############################################################
# 08_Locus_Edge_Pathway_Correlations.R
#
# WET AMD:
# LOCUS -> CELL-CELL EDGE -> SIGNALING PATHWAY
#
# Starting point:
#   Only locus-edge associations already significant
#   in the primary donor-level analysis (FDR < 0.05).
#
# For each significant locus-edge:
#   donor locus expression
#       vs
#   donor CellChat pathway probability
#
# Primary: Pearson
# Sensitivity: Spearman
#
# IMPORTANT:
# Unit of inference = DONOR
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

############################################################
# PATHS
############################################################

root <- "/storage/lemus_g/roldan/ARMD"

analysis_dir <- file.path(
  root,
  "results/SRP413248/CellChat/donor_level_symbol",
  "donor_figure_analysis/locus_edge_correlations"
)

pathway_dir <- file.path(
  root,
  "results/SRP413248/CellChat/donor_level_symbol",
  "donor_figure_analysis/pathway_analysis"
)

out_dir <- file.path(
  analysis_dir,
  "paper_figures/Locus_Edge_Pathway"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# INPUT 1:
# PRIMARY LOCUS-EDGE RESULTS
############################################################

stats <- read.csv(
  file.path(
    analysis_dir,
    "Locus_Edge_Correlations_PEARSON_PRIMARY.csv"
  ),
  stringsAsFactors = FALSE
)

############################################################
# INPUT 2:
# DONOR GLOBAL LOCUS EXPRESSION
############################################################

expr <- read.csv(
  file.path(
    analysis_dir,
    "Donor_Global_Locus_Expression.csv"
  ),
  stringsAsFactors = FALSE
)

############################################################
# INPUT 3:
# DONOR PATHWAY EDGES
############################################################

pathway_file <- file.path(
  pathway_dir,
  "donor_pathway_edges.csv"
)

if (!file.exists(pathway_file)) {

  stop(
    paste0(
      "Cannot find donor pathway file:\n",
      pathway_file,
      "\nRun 05_extract_donor_pathways.R first."
    )
  )
}

pathways <- read.csv(
  pathway_file,
  stringsAsFactors = FALSE
)

############################################################
# SHOW COLUMN NAMES
############################################################

cat("\n============================================\n")
cat("INPUT COLUMNS\n")
cat("============================================\n")

cat("\nExpression:\n")
print(names(expr))

cat("\nPathway data:\n")
print(names(pathways))

############################################################
# STANDARDIZE POSSIBLE PATHWAY COLUMN NAMES
############################################################

if ("source" %in% names(pathways) &&
    !"sender" %in% names(pathways)) {

  pathways$sender <- pathways$source
}

if ("target" %in% names(pathways) &&
    !"receiver" %in% names(pathways)) {

  pathways$receiver <- pathways$target
}

if ("prob" %in% names(pathways) &&
    !"pathway_probability" %in% names(pathways)) {

  pathways$pathway_probability <- pathways$prob
}

if ("probability" %in% names(pathways) &&
    !"pathway_probability" %in% names(pathways)) {

  pathways$pathway_probability <- pathways$probability
}

############################################################
# CHECK REQUIRED COLUMNS
############################################################

required_pathway <- c(
  "donor",
  "condition",
  "sender",
  "receiver",
  "pathway",
  "pathway_probability"
)

missing_pathway <- setdiff(
  required_pathway,
  names(pathways)
)

if (length(missing_pathway) > 0) {

  stop(
    paste0(
      "Missing pathway columns: ",
      paste(
        missing_pathway,
        collapse = ", "
      )
    )
  )
}

required_expr <- c(
  "donor",
  "condition",
  "locus",
  "mean_expression"
)

missing_expr <- setdiff(
  required_expr,
  names(expr)
)

if (length(missing_expr) > 0) {

  stop(
    paste0(
      "Missing expression columns: ",
      paste(
        missing_expr,
        collapse = ", "
      )
    )
  )
}

############################################################
# SELECT THE 24 PRIMARY WET AMD LOCUS-EDGE ASSOCIATIONS
############################################################

primary_edges <- stats %>%

  filter(
    condition == "Wet AMD",
    predictor == "mean_expression",
    method == "pearson",
    FDR < 0.05
  ) %>%

  mutate(
    locus_clean = sub(
      "^TE-L1FLnI-",
      "",
      locus
    )
  ) %>%

  select(
    locus,
    locus_clean,
    condition,
    sender,
    receiver,
    edge_r = r,
    edge_p = p,
    edge_FDR = FDR
  ) %>%

  distinct()

cat("\n============================================\n")
cat("PRIMARY WET LOCUS-EDGE ASSOCIATIONS\n")
cat("============================================\n")

cat("N =", nrow(primary_edges), "\n")

print(
  as.data.frame(
    primary_edges %>%
      count(
        locus_clean,
        name = "n_edges"
      )
  )
)

if (nrow(primary_edges) != 24) {

  warning(
    paste0(
      "Expected 24 Wet AMD primary associations, found ",
      nrow(primary_edges)
    )
  )
}

############################################################
# PREPARE EXPRESSION
############################################################

expr_wet <- expr %>%

  filter(
    condition == "Wet AMD"
  ) %>%

  select(
    donor,
    condition,
    locus,
    mean_expression,
    detection_rate
  )

############################################################
# PREPARE PATHWAY DATA
############################################################

pathways_wet <- pathways %>%

  filter(
    condition == "Wet AMD"
  ) %>%

  select(
    donor,
    condition,
    sender,
    receiver,
    pathway,
    pathway_probability
  )

############################################################
# RESTRICT PATHWAYS TO THE 24 PRIMARY LOCUS-EDGE PAIRS
############################################################

candidate_data <- primary_edges %>%

  inner_join(
    pathways_wet,
    by = c(
      "condition",
      "sender",
      "receiver"
    )
  ) %>%

  inner_join(
    expr_wet,
    by = c(
      "donor",
      "condition",
      "locus"
    )
  )

cat("\nCandidate donor-pathway observations:",
    nrow(candidate_data), "\n")

cat(
  "Unique locus-edge-pathway combinations:",
  n_distinct(
    paste(
      candidate_data$locus,
      candidate_data$sender,
      candidate_data$receiver,
      candidate_data$pathway
    )
  ),
  "\n"
)

write.csv(
  candidate_data,
  file.path(
    out_dir,
    "MASTER_Wet_Locus_Edge_Pathway_DonorData.csv"
  ),
  row.names = FALSE
)

############################################################
# IMPORTANT:
# COMPLETE DONOR GRID
#
# A pathway absent for a testable sender-receiver pair
# should contribute probability 0 rather than causing the
# donor to disappear from the correlation.
############################################################

wet_donors <- sort(
  unique(
    expr_wet$donor
  )
)

candidate_combinations <- candidate_data %>%

  distinct(
    locus,
    locus_clean,
    condition,
    sender,
    receiver,
    pathway,
    edge_r,
    edge_p,
    edge_FDR
  )

complete_grid <- candidate_combinations %>%

  crossing(
    donor = wet_donors
  ) %>%

  left_join(
    expr_wet,
    by = c(
      "donor",
      "condition",
      "locus"
    )
  ) %>%

  left_join(
    pathways_wet,
    by = c(
      "donor",
      "condition",
      "sender",
      "receiver",
      "pathway"
    )
  )

############################################################
# Missing pathway probability = 0
#
# NOTE:
# This assumes the donor had a testable sender-receiver pair.
# We will additionally require non-missing locus expression.
############################################################

complete_grid <- complete_grid %>%

  mutate(
    pathway_probability = ifelse(
      is.na(pathway_probability),
      0,
      pathway_probability
    )
  )

write.csv(
  complete_grid,
  file.path(
    out_dir,
    "COMPLETE_Wet_Locus_Edge_Pathway_DonorData.csv"
  ),
  row.names = FALSE
)

############################################################
# SAFE CORRELATION
############################################################

safe_cor <- function(
  x,
  y,
  method = "pearson"
) {

  ok <- is.finite(x) &
        is.finite(y)

  x <- x[ok]
  y <- y[ok]

  n <- length(x)

  if (n < 4) {

    return(
      data.frame(
        n = n,
        r = NA_real_,
        p = NA_real_
      )
    )
  }

  if (
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

  if (method == "pearson") {

    z <- suppressWarnings(
      cor.test(
        x,
        y,
        method = "pearson"
      )
    )

  } else {

    z <- suppressWarnings(
      cor.test(
        x,
        y,
        method = "spearman",
        exact = FALSE
      )
    )
  }

  data.frame(
    n = n,
    r = unname(z$estimate),
    p = z$p.value
  )
}

############################################################
# RUN CORRELATIONS
############################################################

run_correlations <- function(method_name) {

  split_data <- split(
    complete_grid,
    interaction(
      complete_grid$locus,
      complete_grid$sender,
      complete_grid$receiver,
      complete_grid$pathway,
      drop = TRUE
    )
  )

  results <- lapply(
    split_data,
    function(df) {

      z <- safe_cor(
        df$mean_expression,
        df$pathway_probability,
        method = method_name
      )

      data.frame(
        locus = df$locus[1],
        locus_clean = df$locus_clean[1],
        condition = df$condition[1],
        sender = df$sender[1],
        receiver = df$receiver[1],
        pathway = df$pathway[1],

        edge_r = df$edge_r[1],
        edge_p = df$edge_p[1],
        edge_FDR = df$edge_FDR[1],

        method = method_name,

        n = z$n,
        r = z$r,
        p = z$p,

        stringsAsFactors = FALSE
      )
    }
  )

  bind_rows(results)
}

cat("\nRunning Pearson...\n")

pearson <- run_correlations(
  "pearson"
)

cat("Running Spearman sensitivity...\n")

spearman <- run_correlations(
  "spearman"
)

results <- bind_rows(
  pearson,
  spearman
)

############################################################
# MULTIPLE TESTING
#
# Hierarchical follow-up:
# adjust pathways WITHIN each already significant
# locus-edge combination.
############################################################

results <- results %>%

  group_by(
    method,
    locus,
    condition,
    sender,
    receiver
  ) %>%

  mutate(
    pathway_FDR = p.adjust(
      p,
      method = "BH"
    )
  ) %>%

  ungroup() %>%

  mutate(
    direction = case_when(
      is.na(r) ~ NA_character_,
      r > 0 ~ "POSITIVE",
      r < 0 ~ "NEGATIVE",
      TRUE ~ "ZERO"
    ),

    pathway_significance = case_when(
      !is.na(pathway_FDR) &
        pathway_FDR < 0.05 ~ "FDR<0.05",

      !is.na(p) &
        p < 0.05 ~ "P<0.05",

      TRUE ~ "NS"
    )
  )

############################################################
# SAVE ALL
############################################################

write.csv(
  results,
  file.path(
    out_dir,
    "Locus_Edge_Pathway_Correlations_ALL.csv"
  ),
  row.names = FALSE
)

############################################################
# PRIMARY PEARSON
############################################################

primary <- results %>%

  filter(
    method == "pearson"
  ) %>%

  arrange(
    pathway_FDR,
    p
  )

write.csv(
  primary,
  file.path(
    out_dir,
    "Locus_Edge_Pathway_PEARSON.csv"
  ),
  row.names = FALSE
)

############################################################
# FDR-SIGNIFICANT PATHWAYS
############################################################

fdr_sig <- primary %>%

  filter(
    !is.na(pathway_FDR),
    pathway_FDR < 0.05
  ) %>%

  arrange(
    locus_clean,
    sender,
    receiver,
    pathway_FDR
  )

write.csv(
  fdr_sig,
  file.path(
    out_dir,
    "Locus_Edge_Pathway_PEARSON_FDR005.csv"
  ),
  row.names = FALSE
)

############################################################
# NOMINAL PATHWAYS
############################################################

nominal <- primary %>%

  filter(
    !is.na(p),
    p < 0.05
  )

write.csv(
  nominal,
  file.path(
    out_dir,
    "Locus_Edge_Pathway_PEARSON_P005.csv"
  ),
  row.names = FALSE
)

############################################################
# TOP PATHWAY PER LOCUS-EDGE
#
# Useful for network labels.
############################################################

top_pathway <- primary %>%

  filter(
    !is.na(r),
    !is.na(p)
  ) %>%

  group_by(
    locus,
    locus_clean,
    condition,
    sender,
    receiver
  ) %>%

  arrange(
    pathway_FDR,
    p,
    desc(abs(r))
  ) %>%

  slice(1) %>%

  ungroup()

write.csv(
  top_pathway,
  file.path(
    out_dir,
    "TOP_Pathway_Per_LocusEdge.csv"
  ),
  row.names = FALSE
)

############################################################
# SIGNIFICANT PATHWAY LABELS PER EDGE
############################################################

pathway_labels <- fdr_sig %>%

  group_by(
    locus,
    locus_clean,
    condition,
    sender,
    receiver
  ) %>%

  summarise(
    n_significant_pathways = n_distinct(pathway),

    significant_pathways = paste(
      sort(
        unique(pathway)
      ),
      collapse = " | "
    ),

    strongest_pathway_r = max(
      abs(r),
      na.rm = TRUE
    ),

    min_pathway_FDR = min(
      pathway_FDR,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

write.csv(
  pathway_labels,
  file.path(
    out_dir,
    "FDR_Pathway_Labels_Per_LocusEdge.csv"
  ),
  row.names = FALSE
)

############################################################
# SUMMARY BY LOCUS
############################################################

summary_locus <- fdr_sig %>%

  group_by(
    locus_clean
  ) %>%

  summarise(
    n_locus_edge_pathway_associations = n(),

    n_unique_edges = n_distinct(
      paste(
        sender,
        receiver,
        sep = " -> "
      )
    ),

    n_unique_pathways = n_distinct(
      pathway
    ),

    pathways = paste(
      sort(
        unique(pathway)
      ),
      collapse = ", "
    ),

    .groups = "drop"
  )

write.csv(
  summary_locus,
  file.path(
    out_dir,
    "Pathway_Summary_By_Locus.csv"
  ),
  row.names = FALSE
)

############################################################
# PRINT SUMMARY
############################################################

cat("\n============================================\n")
cat("LOCUS -> EDGE -> PATHWAY ANALYSIS COMPLETE\n")
cat("============================================\n\n")

cat(
  "Primary Wet locus-edge associations:",
  nrow(primary_edges),
  "\n"
)

cat(
  "Pearson locus-edge-pathway tests:",
  nrow(primary),
  "\n"
)

cat(
  "Nominal P < 0.05:",
  nrow(nominal),
  "\n"
)

cat(
  "Pathway FDR < 0.05:",
  nrow(fdr_sig),
  "\n\n"
)

cat("FDR-significant pathways by locus:\n")

if (nrow(summary_locus) > 0) {

  print(
    as.data.frame(
      summary_locus
    )
  )

} else {

  cat("None survived pathway-level FDR.\n")
}

cat("\nTop pathway per locus-edge:\n")

print(
  head(
    as.data.frame(
      top_pathway %>%
        select(
          locus_clean,
          sender,
          receiver,
          pathway,
          n,
          r,
          p,
          pathway_FDR
        )
    ),
    30
  )
)

cat("\nOutput:\n")
cat(out_dir, "\n")

