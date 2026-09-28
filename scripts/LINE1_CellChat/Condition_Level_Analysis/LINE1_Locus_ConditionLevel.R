############################################################
# LINE1_Locus_ConditionLevel.R
#
# Global locus-specific L1FLnI expression at CONDITION level
# matched to condition-level CellChat network metrics
#
# Unit:
#   Healthy Control
#   Dry AMD
#   Wet AMD
#
# IMPORTANT:
# This is a descriptive condition-level analysis (n = 3).
# It is NOT a donor-level inferential analysis.
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

############################################################
# PATHS
############################################################

seurat_file <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "Seurat/downstream/",
  "SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

cellchat_metrics_file <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "CellChat/by_condition/",
  "CellChat_condition_network_metrics.csv"
)

previous_stats_file <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "CellChat/donor_level_symbol/network_metrics/",
  "LINE1_GlobalLocus_ByDisease/",
  "Global_Locus_Pearson_ByDisease.csv"
)

outdir <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "CellChat/by_condition/",
  "LINE1_GlobalLocus_condition_level"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# GET THE 12 CANDIDATE LOCI
############################################################

previous_stats <- read.csv(
  previous_stats_file,
  stringsAsFactors = FALSE
)

candidate_loci <- sort(
  unique(previous_stats$Feature)
)

candidate_loci <- candidate_loci[
  grepl("^TE-L1FLnI-", candidate_loci)
]

cat("\n============================================\n")
cat("CANDIDATE L1FLnI LOCI\n")
cat("============================================\n")

cat("Number of loci:", length(candidate_loci), "\n")
print(candidate_loci)

############################################################
# LOAD SEURAT
############################################################

cat("\nLoading Seurat object...\n")

seu <- readRDS(seurat_file)

DefaultAssay(seu) <- "RNA"

meta <- seu@meta.data

if(!"disease_state" %in% colnames(meta)){
  stop("disease_state not found in Seurat metadata.")
}

############################################################
# CHECK LOCI
############################################################

present_loci <- intersect(
  candidate_loci,
  rownames(seu)
)

missing_loci <- setdiff(
  candidate_loci,
  rownames(seu)
)

cat("\nLoci present:", length(present_loci), "\n")

if(length(missing_loci) > 0){
  cat("\nWARNING - missing loci:\n")
  print(missing_loci)
}

############################################################
# NORMALIZED EXPRESSION
############################################################

expr <- GetAssayData(
  seu,
  assay = "RNA",
  layer = "data"
)

expr <- expr[present_loci, , drop = FALSE]

############################################################
# CONDITIONS
############################################################

conditions <- c(
  "Healthy Control",
  "Dry AMD",
  "Wet AMD"
)

results <- list()
k <- 1

############################################################
# CONDITION × LOCUS
############################################################

for(cond in conditions){

  cells <- rownames(meta)[
    meta$disease_state == cond
  ]

  cat(
    "\n",
    cond,
    ": ",
    length(cells),
    " cells\n",
    sep = ""
  )

  mat <- expr[
    ,
    cells,
    drop = FALSE
  ]

  mean_expression <- Matrix::rowMeans(mat)

  detection_rate <- Matrix::rowMeans(
    mat > 0
  )

  results[[k]] <- data.frame(
    Condition = cond,
    Feature = present_loci,
    n_cells = length(cells),
    mean_expression = as.numeric(mean_expression),
    detection_rate = as.numeric(detection_rate),
    stringsAsFactors = FALSE
  )

  k <- k + 1
}

locus_condition <- do.call(
  rbind,
  results
)

############################################################
# CELLCHAT METRICS
############################################################

network <- read.csv(
  cellchat_metrics_file,
  stringsAsFactors = FALSE
)

cat("\n============================================\n")
cat("CELLCHAT METRICS\n")
cat("============================================\n")

print(network)

############################################################
# MERGE
############################################################

final <- merge(
  locus_condition,
  network[
    ,
    c(
      "Condition",
      "density",
      "total_strength"
    )
  ],
  by = "Condition",
  all.x = TRUE
)

final <- final[
  order(
    final$Feature,
    match(final$Condition, conditions)
  ),
]

############################################################
# SAVE
############################################################

outfile <- file.path(
  outdir,
  "Global_Locus_ConditionLevel.csv"
)

write.csv(
  final,
  outfile,
  row.names = FALSE
)

cat("\n============================================\n")
cat("CONDITION-LEVEL LINE1 + CELLCHAT TABLE\n")
cat("============================================\n\n")

print(final)

cat("\nSaved:\n")
cat(outfile, "\n")

cat("\nExpected rows:", length(present_loci) * 3, "\n")
cat("Observed rows:", nrow(final), "\n")

cat("\nDONE\n")

