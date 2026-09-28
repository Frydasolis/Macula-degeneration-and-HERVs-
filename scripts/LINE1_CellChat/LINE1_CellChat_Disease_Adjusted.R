############################################################
# LINE1_CellChat_Disease_Adjusted.R
#
# Disease-adjusted donor-level analysis:
#
# Network metric ~ L1FLnI + Disease_state
#
# Statistical unit = DONOR
#
# Primary outcomes:
#   Network Density
#   Total Communication Strength
#
# Produces:
#   - disease-adjusted statistics
#   - publication-quality scatter plots
#   - summary forest plot
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

cat("\n============================================\n")
cat("DISEASE-ADJUSTED L1FLnI / CELLCHAT ANALYSIS\n")
cat("============================================\n\n")


############################################################
# 1. Paths
############################################################

base <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/CellChat/donor_level_symbol/",
  "network_metrics/LINE1_Pearson"
)

table_dir <- file.path(base, "tables")

outdir <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/CellChat/donor_level_symbol/",
  "network_metrics/LINE1_DiseaseAdjusted"
)

figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)


############################################################
# 2. Read previously generated donor tables
############################################################

global <- read.csv(
  file.path(
    table_dir,
    "family_global_L1FLnI_burden.csv"
  ),
  stringsAsFactors = FALSE
)

celltype <- read.csv(
  file.path(
    table_dir,
    "family_celltype_L1FLnI_burden.csv"
  ),
  stringsAsFactors = FALSE
)

candidate <- read.csv(
  file.path(
    table_dir,
    "candidate_locus_donor_expression.csv"
  ),
  stringsAsFactors = FALSE
)


############################################################
# 3. Standardize disease state
############################################################

fix_condition <- function(x) {

  x <- as.character(x)

  x[x %in% c(
    "Healthy",
    "Healthy Control",
    "HC"
  )] <- "Healthy"

  x[x %in% c(
    "Dry",
    "Dry AMD"
  )] <- "Dry AMD"

  x[x %in% c(
    "Wet",
    "Wet AMD"
  )] <- "Wet AMD"

  factor(
    x,
    levels = c(
      "Healthy",
      "Dry AMD",
      "Wet AMD"
    )
  )
}

global$Condition <- fix_condition(
  global$Condition
)

celltype$Condition <- fix_condition(
  celltype$Condition
)

candidate$Condition <- fix_condition(
  candidate$Condition
)


cat("Global disease distribution:\n")
print(table(global$Condition))


############################################################
# 4. Disease-adjusted model function
############################################################

adjusted_model <- function(
    dat,
    predictor,
    outcome) {

  keep <- complete.cases(
    dat[
      ,
      c(
        predictor,
        outcome,
        "Condition"
      )
    ]
  )

  d <- dat[
    keep,
    ,
    drop = FALSE
  ]

  if (nrow(d) < 6) {
    return(NULL)
  }

  if (sd(d[[predictor]]) == 0) {
    return(NULL)
  }

  form <- as.formula(
    paste(
      outcome,
      "~",
      predictor,
      "+ Condition"
    )
  )

  fit <- lm(
    form,
    data = d
  )

  co <- summary(fit)$coefficients

  if (!predictor %in% rownames(co)) {
    return(NULL)
  }

  beta <- co[
    predictor,
    "Estimate"
  ]

  se <- co[
    predictor,
    "Std. Error"
  ]

  p <- co[
    predictor,
    "Pr(>|t|)"
  ]

  ci <- confint(
    fit,
    predictor,
    level = 0.95
  )

  # standardized beta
  d_std <- d

  d_std[[predictor]] <- as.numeric(
    scale(d_std[[predictor]])
  )

  d_std[[outcome]] <- as.numeric(
    scale(d_std[[outcome]])
  )

  fit_std <- lm(
    form,
    data = d_std
  )

  beta_std <- summary(
    fit_std
  )$coefficients[
    predictor,
    "Estimate"
  ]

  # crude Pearson for reference
  pear <- cor.test(
    d[[predictor]],
    d[[outcome]],
    method = "pearson"
  )

  data.frame(
    n = nrow(d),

    beta = beta,
    SE = se,

    CI_low = ci[1],
    CI_high = ci[2],

    standardized_beta = beta_std,

    adjusted_P = p,

    crude_Pearson_r =
      unname(pear$estimate),

    crude_Pearson_P =
      pear$p.value,

    model_R2 =
      summary(fit)$r.squared,

    adjusted_R2 =
      summary(fit)$adj.r.squared,

    stringsAsFactors = FALSE
  )
}


############################################################
# 5. GLOBAL FAMILY-WIDE ANALYSIS
############################################################

cat("\n============================================\n")
cat("GLOBAL FAMILY-WIDE\n")
cat("============================================\n")

global_results <- list()

k <- 1

for (
  predictor in c(
    "L1_expression_burden",
    "detection_breadth"
  )
) {

  for (
    outcome in c(
      "density",
      "total_strength"
    )
  ) {

    res <- adjusted_model(
      global,
      predictor,
      outcome
    )

    if (!is.null(res)) {

      global_results[[k]] <- cbind(
        data.frame(
          Analysis = "Global",
          CellType = "ALL",
          Feature = "ALL_L1FLnI",
          Predictor = predictor,
          Outcome = outcome,
          stringsAsFactors = FALSE
        ),
        res
      )

      k <- k + 1
    }
  }
}

global_results <- do.call(
  rbind,
  global_results
)


############################################################
# 6. CELL-TYPE FAMILY-WIDE ANALYSIS
############################################################

cat("\n============================================\n")
cat("CELL-TYPE FAMILY-WIDE\n")
cat("============================================\n")

cell_results <- list()

k <- 1

for (
  ct in sort(
    unique(celltype$CellType)
  )
) {

  d <- celltype[
    celltype$CellType == ct,
    ,
    drop = FALSE
  ]

  for (
    predictor in c(
      "L1_expression_burden",
      "detection_breadth"
    )
  ) {

    for (
      outcome in c(
        "density",
        "total_strength"
      )
    ) {

      res <- adjusted_model(
        d,
        predictor,
        outcome
      )

      if (!is.null(res)) {

        cell_results[[k]] <- cbind(
          data.frame(
            Analysis = "Cell-type family",
            CellType = ct,
            Feature = "ALL_L1FLnI",
            Predictor = predictor,
            Outcome = outcome,
            stringsAsFactors = FALSE
          ),
          res
        )

        k <- k + 1
      }
    }
  }
}

cell_results <- do.call(
  rbind,
  cell_results
)


############################################################
# 7. LOCUS-SPECIFIC ANALYSIS
############################################################

cat("\n============================================\n")
cat("LOCUS-SPECIFIC\n")
cat("============================================\n")

groups <- unique(
  candidate[
    ,
    c(
      "CellType",
      "TE_feature"
    )
  ]
)

locus_results <- list()

k <- 1

for (
  i in seq_len(
    nrow(groups)
  )
) {

  ct <- groups$CellType[i]
  locus <- groups$TE_feature[i]

  d <- candidate[
    candidate$CellType == ct &
      candidate$TE_feature == locus,
    ,
    drop = FALSE
  ]

  for (
    predictor in c(
      "mean_expression",
      "detection_rate"
    )
  ) {

    for (
      outcome in c(
        "density",
        "total_strength"
      )
    ) {

      res <- adjusted_model(
        d,
        predictor,
        outcome
      )

      if (!is.null(res)) {

        locus_results[[k]] <- cbind(
          data.frame(
            Analysis = "Locus-specific",
            CellType = ct,
            Feature = locus,
            Predictor = predictor,
            Outcome = outcome,
            stringsAsFactors = FALSE
          ),
          res
        )

        k <- k + 1
      }
    }
  }
}

locus_results <- do.call(
  rbind,
  locus_results
)


############################################################
# 8. Combine + FDR
############################################################

master <- rbind(
  global_results,
  cell_results,
  locus_results
)

# FDR within each analysis family
master$FDR <- NA_real_

families <- unique(
  master$Analysis
)

for (a in families) {

  idx <- which(
    master$Analysis == a
  )

  master$FDR[idx] <- p.adjust(
    master$adjusted_P[idx],
    method = "BH"
  )
}

master <- master[
  order(
    master$adjusted_P
  ),
  ,
  drop = FALSE
]


############################################################
# 9. Save tables
############################################################

write.csv(
  global_results,
  file.path(
    tabdir,
    "Global_disease_adjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  cell_results,
  file.path(
    tabdir,
    "CellType_disease_adjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  locus_results,
  file.path(
    tabdir,
    "Locus_disease_adjusted.csv"
  ),
  row.names = FALSE
)

write.csv(
  master,
  file.path(
    tabdir,
    "MASTER_disease_adjusted.csv"
  ),
  row.names = FALSE
)


############################################################
# 10. Print strongest results
############################################################

cat("\n============================================\n")
cat("TOP DISEASE-ADJUSTED ASSOCIATIONS\n")
cat("============================================\n\n")

print(
  head(
    master[
      ,
      c(
        "Analysis",
        "CellType",
        "Feature",
        "Predictor",
        "Outcome",
        "n",
        "standardized_beta",
        "adjusted_P",
        "FDR",
        "crude_Pearson_r"
      )
    ],
    20
  ),
  row.names = FALSE
)


############################################################
# 11. Paper scatter function
############################################################

paper_scatter <- function(
    dat,
    predictor,
    outcome,
    title,
    xlab,
    ylab,
    filename) {

  keep <- complete.cases(
    dat[
      ,
      c(
        predictor,
        outcome,
        "Condition",
        "Donor"
      )
    ]
  )

  d <- dat[
    keep,
    ,
    drop = FALSE
  ]

  adj <- adjusted_model(
    d,
    predictor,
    outcome
  )

  pear <- cor.test(
    d[[predictor]],
    d[[outcome]],
    method = "pearson"
  )

  subtitle <- paste0(
    "Pearson r = ",
    sprintf(
      "%.2f",
      pear$estimate
    ),
    ", P = ",
    format.pval(
      pear$p.value,
      digits = 2,
      eps = 0.001
    ),
    " | disease-adjusted P = ",
    format.pval(
      adj$adjusted_P,
      digits = 2,
      eps = 0.001
    )
  )

  p <- ggplot(
    d,
    aes(
      x = .data[[predictor]],
      y = .data[[outcome]],
      color = Condition
    )
  ) +

    geom_point(
      size = 3.3,
      alpha = 0.9
    ) +

    # Overall relationship
    geom_smooth(
      aes(
        group = 1
      ),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "black",
      linewidth = 0.8,
      inherit.aes = TRUE
    ) +

    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      color = NULL
    ) +

    theme_classic(
      base_size = 13
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),

      plot.subtitle = element_text(
        size = 10
      ),

      axis.title = element_text(
        face = "bold"
      ),

      legend.position = "top"
    )

  ggsave(
    paste0(
      filename,
      ".pdf"
    ),
    p,
    width = 6.2,
    height = 5.3
  )

  ggsave(
    paste0(
      filename,
      ".png"
    ),
    p,
    width = 6.2,
    height = 5.3,
    dpi = 600
  )
}


############################################################
# 12. FIGURE A
# Global L1FLnI -> Density
############################################################

paper_scatter(
  global,

  "L1_expression_burden",

  "density",

  "Retina-wide L1FLnI expression",

  "L1FLnI expression burden",

  "Network density",

  file.path(
    figdir,
    "FigA_Global_L1FLnI_Density"
  )
)


############################################################
# 13. FIGURE B
# Global L1FLnI -> Total strength
############################################################

paper_scatter(
  global,

  "L1_expression_burden",

  "total_strength",

  "Retina-wide L1FLnI expression",

  "L1FLnI expression burden",

  "Total communication strength",

  file.path(
    figdir,
    "FigB_Global_L1FLnI_TotalStrength"
  )
)


############################################################
# 14. FIGURE C
# Photoreceptor L1FLnI -> Density
############################################################

photo <- celltype[
  celltype$CellType ==
    "Photoreceptor",
  ,
  drop = FALSE
]

paper_scatter(
  photo,

  "L1_expression_burden",

  "density",

  "Photoreceptor L1FLnI expression",

  "L1FLnI expression burden",

  "Network density",

  file.path(
    figdir,
    "FigC_Photoreceptor_L1FLnI_Density"
  )
)


############################################################
# 15. FIGURE D
# RGC breadth -> Density
############################################################

rgc <- celltype[
  celltype$CellType ==
    "RGC",
  ,
  drop = FALSE
]

paper_scatter(
  rgc,

  "detection_breadth",

  "density",

  "RGC L1FLnI detection breadth",

  "Fraction of L1FLnI loci detected",

  "Network density",

  file.path(
    figdir,
    "FigD_RGC_L1FLnI_Density"
  )
)


############################################################
# 16. FIGURE E
# Astrocyte Xq23ab
############################################################

astro_x <- candidate[
  candidate$CellType ==
    "Astrocyte" &
  candidate$TE_feature ==
    "TE-L1FLnI-Xq23ab",
  ,
  drop = FALSE
]

paper_scatter(
  astro_x,

  "detection_rate",

  "total_strength",

  "Astrocyte L1FLnI-Xq23ab",

  "Detection rate",

  "Total communication strength",

  file.path(
    figdir,
    "FigE_Astrocyte_Xq23ab_TotalStrength"
  )
)


############################################################
# 17. FIGURE F
# RPE 21q22.13b
############################################################

rpe_locus <- candidate[
  candidate$CellType ==
    "RPE" &
  candidate$TE_feature ==
    "TE-L1FLnI-21q22.13b",
  ,
  drop = FALSE
]

paper_scatter(
  rpe_locus,

  "detection_rate",

  "density",

  "RPE L1FLnI-21q22.13b",

  "Detection rate",

  "Network density",

  file.path(
    figdir,
    "FigF_RPE_21q22.13b_Density"
  )
)


############################################################
# 18. Forest plot
############################################################

forest <- master[
  master$Analysis %in%
    c(
      "Global",
      "Cell-type family"
    ),
  ,
  drop = FALSE
]

forest$Label <- ifelse(
  forest$Analysis == "Global",

  paste0(
    "Retina-wide | ",
    forest$Predictor,
    " | ",
    forest$Outcome
  ),

  paste0(
    forest$CellType,
    " | ",
    forest$Predictor,
    " | ",
    forest$Outcome
  )
)

forest$Label <- factor(
  forest$Label,
  levels = rev(
    forest$Label[
      order(
        forest$standardized_beta
      )
    ]
  )
)

p_forest <- ggplot(
  forest,
  aes(
    x = standardized_beta,
    y = Label
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_point(
    size = 2.8
  ) +

  labs(
    title =
      "Disease-adjusted associations between L1FLnI and retinal communication",

    x =
      "Standardized disease-adjusted coefficient",

    y = NULL
  ) +

  theme_classic(
    base_size = 11
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold"
      )
  )

ggsave(
  file.path(
    figdir,
    "FigG_DiseaseAdjusted_Forest.pdf"
  ),
  p_forest,
  width = 9,
  height = 9
)

ggsave(
  file.path(
    figdir,
    "FigG_DiseaseAdjusted_Forest.png"
  ),
  p_forest,
  width = 9,
  height = 9,
  dpi = 600
)


############################################################
# 19. DONE
############################################################

cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n")

cat("\nTables:\n")
cat(tabdir, "\n")

cat("\nFigures:\n")
cat(figdir, "\n")

cat("\n")
