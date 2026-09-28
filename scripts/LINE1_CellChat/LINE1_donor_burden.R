############################################################
# DONOR-LEVEL L1FLnI BURDEN
#
# Statistical unit:
#   DONOR
#
# L1 loci are evaluated only in the retinal cell type
# in which they were identified.
#
# Outputs:
#   1. donor_celltype_LINE1_burden.csv
#   2. donor_locus_LINE1_expression.csv
#   3. donor_LINE1_burden_wide.csv
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

############################################################
# PATHS
############################################################

seurat_file <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "Seurat/downstream/",
  "SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

candidate_file <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "Seurat/downstream/TE_analysis/Figure4/",
  "LINE1_CellChat_candidates.csv"
)

outdir <- "network_metrics/LINE1_burden"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# LOAD DATA
############################################################

cat("\n============================================\n")
cat("LOADING DATA\n")
cat("============================================\n")

seu <- readRDS(seurat_file)

candidates <- read.csv(
  candidate_file,
  stringsAsFactors = FALSE
)

cat("Seurat cells:", ncol(seu), "\n")
cat("Candidate rows:", nrow(candidates), "\n")

############################################################
# INSPECT CANDIDATE FILE
############################################################

cat("\nCandidate columns:\n")
print(colnames(candidates))

cat("\nCandidate cell types:\n")
print(table(candidates$HRCA_majorclass))

############################################################
# KEEP L1FLnI ONLY
############################################################

candidates <- candidates %>%
  filter(TE_family == "L1FLnI")

candidate_map <- candidates %>%
  select(
    HRCA_majorclass,
    TE_feature
  ) %>%
  distinct()

cat("\n============================================\n")
cat("L1FLnI LOCI\n")
cat("============================================\n")

print(candidate_map)

cat(
  "\nUnique loci:",
  length(unique(candidate_map$TE_feature)),
  "\n"
)

############################################################
# CHECK FEATURES EXIST IN SEURAT
############################################################

features_available <- rownames(seu)

candidate_map$present <- (
  candidate_map$TE_feature %in% features_available
)

cat("\n============================================\n")
cat("FEATURE CHECK\n")
cat("============================================\n")

print(candidate_map)

if (!all(candidate_map$present)) {

  cat("\nWARNING: Missing features:\n")

  print(
    candidate_map$TE_feature[
      !candidate_map$present
    ]
  )
}

candidate_map <- candidate_map %>%
  filter(present)

############################################################
# DONOR COLUMN
############################################################

donor_col <- "orig.ident"

celltype_col <- "HRCA_majorclass"

condition_col <- "disease_state"

cat("\n============================================\n")
cat("METADATA\n")
cat("============================================\n")

cat("Donor column:", donor_col, "\n")
cat("Cell type column:", celltype_col, "\n")
cat("Condition column:", condition_col, "\n")

############################################################
# GET NORMALIZED EXPRESSION
############################################################

DefaultAssay(seu) <- "RNA"

expr <- GetAssayData(
  seu,
  assay = "RNA",
  layer = "data"
)

############################################################
# STORAGE
############################################################

locus_results <- list()

counter <- 1

############################################################
# LOOP:
# CELL TYPE -> LOCUS -> DONOR
############################################################

for (ct in unique(candidate_map$HRCA_majorclass)) {

  cat("\n============================================\n")
  cat("CELL TYPE:", ct, "\n")
  cat("============================================\n")

  loci <- candidate_map %>%
    filter(HRCA_majorclass == ct) %>%
    pull(TE_feature)

  cat("Loci:", length(loci), "\n")

  cells_ct <- colnames(seu)[
    seu@meta.data[[celltype_col]] == ct
  ]

  cat("Cells:", length(cells_ct), "\n")

  donors <- unique(
    seu@meta.data[
      cells_ct,
      donor_col
    ]
  )

  donors <- sort(donors)

  for (donor in donors) {

    cells_donor <- cells_ct[
      seu@meta.data[
        cells_ct,
        donor_col
      ] == donor
    ]

    n_cells <- length(cells_donor)

    if (n_cells == 0) {
      next
    }

    condition <- unique(
      seu@meta.data[
        cells_donor,
        condition_col
      ]
    )

    condition <- condition[1]

    for (locus in loci) {

      values <- as.numeric(
        expr[
          locus,
          cells_donor
        ]
      )

      mean_expression <- mean(
        values,
        na.rm = TRUE
      )

      median_expression <- median(
        values,
        na.rm = TRUE
      )

      detection_rate <- mean(
        values > 0,
        na.rm = TRUE
      )

      sum_expression <- sum(
        values,
        na.rm = TRUE
      )

      locus_results[[counter]] <- data.frame(

        Donor = donor,
        Condition = condition,
        CellType = ct,
        Locus = locus,

        n_cells = n_cells,

        mean_expression = mean_expression,
        median_expression = median_expression,
        detection_rate = detection_rate,
        sum_expression = sum_expression,

        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }
  }
}

############################################################
# COMBINE LOCUS RESULTS
############################################################

locus_df <- bind_rows(
  locus_results
)

############################################################
# SAVE LOCUS-LEVEL TABLE
############################################################

write.csv(
  locus_df,
  file.path(
    outdir,
    "donor_locus_LINE1_expression.csv"
  ),
  row.names = FALSE
)

############################################################
# DONOR x CELL-TYPE BURDEN
#
# Primary burden:
# mean of mean normalized expression across candidate loci
#
# Also retain:
# detection burden
############################################################

burden_df <- locus_df %>%

  group_by(
    Donor,
    Condition,
    CellType
  ) %>%

  summarise(

    n_cells = first(n_cells),

    n_loci = n_distinct(Locus),

    L1_expression_burden = mean(
      mean_expression,
      na.rm = TRUE
    ),

    L1_detection_burden = mean(
      detection_rate,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

############################################################
# QC FLAG BASED ON CELL NUMBER
#
# Do NOT remove donors yet.
# Merely flag them.
############################################################

burden_df <- burden_df %>%
  mutate(

    cell_number_flag = case_when(

      n_cells >= 50 ~ ">=50 cells",

      n_cells >= 20 ~ "20-49 cells",

      n_cells >= 10 ~ "10-19 cells",

      n_cells > 0 ~ "<10 cells",

      TRUE ~ "No cells"
    )
  )

############################################################
# SAVE LONG TABLE
############################################################

write.csv(
  burden_df,
  file.path(
    outdir,
    "donor_celltype_LINE1_burden.csv"
  ),
  row.names = FALSE
)

############################################################
# WIDE TABLE
############################################################

expression_wide <- burden_df %>%

  select(
    Donor,
    Condition,
    CellType,
    L1_expression_burden
  ) %>%

  tidyr::pivot_wider(
    names_from = CellType,
    values_from = L1_expression_burden,
    names_prefix = "L1_"
  )

detection_wide <- burden_df %>%

  select(
    Donor,
    CellType,
    L1_detection_burden
  ) %>%

  tidyr::pivot_wider(
    names_from = CellType,
    values_from = L1_detection_burden,
    names_prefix = "Detection_"
  )

wide_df <- left_join(
  expression_wide,
  detection_wide,
  by = "Donor"
)

############################################################
# SAVE WIDE TABLE
############################################################

write.csv(
  wide_df,
  file.path(
    outdir,
    "donor_LINE1_burden_wide.csv"
  ),
  row.names = FALSE
)

############################################################
# PRINT RESULTS
############################################################

cat("\n============================================\n")
cat("DONOR x CELL TYPE L1FLnI BURDEN\n")
cat("============================================\n\n")

print(
  burden_df,
  n = Inf
)

cat("\n============================================\n")
cat("CELL NUMBER QC\n")
cat("============================================\n\n")

print(
  burden_df %>%
    select(
      Donor,
      Condition,
      CellType,
      n_cells,
      cell_number_flag
    ),
  n = Inf
)

cat("\n============================================\n")
cat("WIDE DONOR TABLE\n")
cat("============================================\n\n")

print(
  wide_df,
  n = Inf
)

cat("\n============================================\n")
cat("FILES CREATED\n")
cat("============================================\n")

cat(
  file.path(
    outdir,
    "donor_locus_LINE1_expression.csv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "donor_celltype_LINE1_burden.csv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "donor_LINE1_burden_wide.csv"
  ),
  "\n"
)

cat("\nDONE\n")
