############################################################
# GLOBAL LOCUS-SPECIFIC L1FLnI x CELLCHAT
#
# Expression of each candidate locus across ALL retinal cells
# for each donor.
#
# Statistical unit = DONOR
# Pearson correlation within disease condition
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

base <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "CellChat/donor_level_symbol/network_metrics"
)

network_file <- file.path(
  base,
  "donor_network_metrics.csv"
)

master_file <- file.path(
  base,
  "LINE1_ByDisease/tables/MASTER_byDisease.csv"
)

outdir <- file.path(
  base,
  "LINE1_GlobalLocus_ByDisease"
)

dir.create(
  outdir,
  recursive=TRUE,
  showWarnings=FALSE
)

############################################################
# READ
############################################################

cat("Reading Seurat object...\n")

seu <- readRDS(seurat_file)

net <- read.csv(
  network_file,
  stringsAsFactors=FALSE
)

master <- read.csv(
  master_file,
  stringsAsFactors=FALSE
)

############################################################
# CANDIDATE LOCI
############################################################

# Extract the exact candidate loci used in the previous
# locus-specific analysis. Cell type is deliberately ignored
# here because this analysis measures each locus globally
# across ALL retinal cells from each donor.

candidate_loci <- unique(
  as.character(
    master$Feature[
      master$Analysis == "Locus-specific"
    ]
  )
)

candidate_loci <- candidate_loci[
  !is.na(candidate_loci) &
  candidate_loci != ""
]

cat(
  "Candidate loci from MASTER:",
  length(candidate_loci),
  "\n"
)

print(candidate_loci)

missing_loci <- setdiff(
  candidate_loci,
  rownames(seu)
)

if(length(missing_loci) > 0){

  cat(
    "\nWARNING — loci not found in Seurat:\n"
  )

  print(missing_loci)
}

candidate_loci <- intersect(
  candidate_loci,
  rownames(seu)
)

cat(
  "\nCandidate loci found in Seurat:",
  length(candidate_loci),
  "\n"
)

############################################################
# EXPRESSION MATRIX
############################################################

DefaultAssay(seu) <- "RNA"

mat <- GetAssayData(
  seu,
  assay="RNA",
  layer="data"
)

meta <- seu@meta.data

donors <- unique(
  as.character(
    meta$orig.ident
  )
)

############################################################
# GLOBAL EXPRESSION PER DONOR x LOCUS
############################################################

results <- list()

k <- 1

for(donor in donors){

  cells <- rownames(meta)[
    meta$orig.ident == donor
  ]

  condition <- unique(
    as.character(
      meta$disease_state[
        meta$orig.ident == donor
      ]
    )
  )

  if(length(condition) != 1){
    stop(
      paste(
        "Multiple conditions for donor",
        donor
      )
    )
  }

  donor_mat <- mat[
    candidate_loci,
    cells,
    drop=FALSE
  ]

  for(feature in candidate_loci){

    x <- donor_mat[
      feature,
      ,
      drop=TRUE
    ]

    results[[k]] <- data.frame(
      Donor=donor,
      Condition=condition,
      Feature=feature,
      n_cells=length(cells),

      mean_expression=
        mean(x),

      detection_rate=
        mean(x > 0),

      stringsAsFactors=FALSE
    )

    k <- k + 1
  }
}

expr <- do.call(
  rbind,
  results
)

############################################################
# NETWORK METRICS
############################################################

# Find donor column if necessary

donor_col <- c(
  "Donor",
  "donor",
  "GSM"
)

donor_col <- donor_col[
  donor_col %in% colnames(net)
][1]

if(is.na(donor_col)){
  stop("Cannot identify donor column in network metrics.")
}

net$Donor <- as.character(
  net[[donor_col]]
)

keep <- c(
  "Donor",
  "density",
  "total_strength"
)

expr <- merge(
  expr,
  net[,keep],
  by="Donor",
  all.x=TRUE
)

############################################################
# NORMALIZE CONDITION
############################################################

expr$Condition[
  expr$Condition %in%
    c("Healthy Control","HC")
] <- "Healthy"

expr$Condition[
  expr$Condition == "Dry"
] <- "Dry AMD"

expr$Condition[
  expr$Condition == "Wet"
] <- "Wet AMD"

############################################################
# SAVE DONOR TABLE
############################################################

write.csv(
  expr,
  file.path(
    outdir,
    "Global_Locus_Donor_Expression.csv"
  ),
  row.names=FALSE
)

############################################################
# PEARSON BY CONDITION
############################################################

conditions <- c(
  "Healthy",
  "Dry AMD",
  "Wet AMD"
)

predictors <- c(
  "mean_expression",
  "detection_rate"
)

outcomes <- c(
  "density",
  "total_strength"
)

stats <- list()

k <- 1

for(feature in candidate_loci){

  for(cond in conditions){

    d <- expr[
      expr$Feature == feature &
      expr$Condition == cond,
      ,
      drop=FALSE
    ]

    for(pred in predictors){

      for(outcome in outcomes){

        z <- d[
          is.finite(d[[pred]]) &
          is.finite(d[[outcome]]),
          ,
          drop=FALSE
        ]

        if(nrow(z) < 3)
          next

        # Pearson requires variation in both variables
        if(
          sd(z[[pred]]) == 0 ||
          sd(z[[outcome]]) == 0
        )
          next

        test <- cor.test(
          z[[pred]],
          z[[outcome]],
          method="pearson"
        )

        stats[[k]] <- data.frame(
          Analysis="Global locus-specific",
          Feature=feature,
          Predictor=pred,
          Outcome=outcome,
          Condition=cond,
          n=nrow(z),
          Pearson_r=unname(test$estimate),
          R2=unname(test$estimate)^2,
          P_value=test$p.value,
          Interpretation=ifelse(
            cond == "Dry AMD",
            "DESCRIPTIVE_n3",
            "INFERENTIAL"
          ),
          stringsAsFactors=FALSE
        )

        k <- k + 1
      }
    }
  }
}

stats <- do.call(
  rbind,
  stats
)

############################################################
# BH FDR
############################################################

stats$FDR <- p.adjust(
  stats$P_value,
  method="BH"
)

stats <- stats[
  order(
    stats$P_value
  ),
]

write.csv(
  stats,
  file.path(
    outdir,
    "Global_Locus_Pearson_ByDisease.csv"
  ),
  row.names=FALSE
)

############################################################
# SUMMARY
############################################################

cat("\n========================================\n")
cat("GLOBAL LOCUS ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "Loci:",
  length(unique(expr$Feature)),
  "\n"
)

cat(
  "Donors:",
  length(unique(expr$Donor)),
  "\n"
)

cat(
  "Statistical tests:",
  nrow(stats),
  "\n\n"
)

print(
  table(
    stats$Condition
  )
)

cat("\nTop results:\n")

print(
  head(
    stats,
    20
  )
)

cat(
  "\nOutput:\n",
  outdir,
  "\n",
  sep=""
)

