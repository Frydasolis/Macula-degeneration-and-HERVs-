############################################################
# LINE1_CellChat_ByDisease.R
#
# Donor-level Pearson correlations stratified by disease
# Healthy / Dry AMD / Wet AMD
#
# IMPORTANT:
# Statistical unit = donor
# Dry AMD n=3 -> descriptive only
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

base <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/CellChat/donor_level_symbol/",
  "network_metrics/LINE1_Pearson"
)

table_dir <- file.path(base, "tables")

outdir <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/CellChat/donor_level_symbol/",
  "network_metrics/LINE1_ByDisease"
)

tabdir <- file.path(outdir, "tables")
figdir <- file.path(outdir, "figures")

dir.create(tabdir, recursive=TRUE, showWarnings=FALSE)
dir.create(figdir, recursive=TRUE, showWarnings=FALSE)


############################################################
# READ DATA
############################################################

global <- read.csv(
  file.path(table_dir, "family_global_L1FLnI_burden.csv"),
  stringsAsFactors=FALSE
)

celltype <- read.csv(
  file.path(table_dir, "family_celltype_L1FLnI_burden.csv"),
  stringsAsFactors=FALSE
)

candidate <- read.csv(
  file.path(table_dir, "candidate_locus_donor_expression.csv"),
  stringsAsFactors=FALSE
)


############################################################
# STANDARDIZE CONDITION
############################################################

fix_condition <- function(x){

  x <- as.character(x)

  x[x %in% c("Healthy","Healthy Control","HC")] <- "Healthy"
  x[x %in% c("Dry","Dry AMD")] <- "Dry AMD"
  x[x %in% c("Wet","Wet AMD")] <- "Wet AMD"

  factor(
    x,
    levels=c("Healthy","Dry AMD","Wet AMD")
  )
}

global$Condition    <- fix_condition(global$Condition)
celltype$Condition  <- fix_condition(celltype$Condition)
candidate$Condition <- fix_condition(candidate$Condition)


############################################################
# SAFE PEARSON
############################################################

pearson_one <- function(dat, predictor, outcome){

  keep <- complete.cases(
    dat[,c(predictor,outcome)]
  )

  d <- dat[keep,,drop=FALSE]

  n <- nrow(d)

  if(n < 3)
    return(NULL)

  if(
    sd(d[[predictor]]) == 0 ||
    sd(d[[outcome]]) == 0
  )
    return(NULL)

  ct <- cor.test(
    d[[predictor]],
    d[[outcome]],
    method="pearson"
  )

  data.frame(
    n=n,
    Pearson_r=unname(ct$estimate),
    R2=unname(ct$estimate)^2,
    P_value=ct$p.value
  )
}


############################################################
# GENERAL FUNCTION BY CONDITION
############################################################

run_by_condition <- function(
  dat,
  analysis,
  celltype_name,
  feature,
  predictor,
  outcome
){

  ans <- list()
  k <- 1

  for(cond in levels(global$Condition)){

    d <- dat[
      dat$Condition == cond,
      ,
      drop=FALSE
    ]

    res <- pearson_one(
      d,
      predictor,
      outcome
    )

    if(!is.null(res)){

      ans[[k]] <- cbind(

        data.frame(
          Analysis=analysis,
          CellType=celltype_name,
          Feature=feature,
          Predictor=predictor,
          Outcome=outcome,
          Condition=cond,
          stringsAsFactors=FALSE
        ),

        res
      )

      k <- k + 1
    }
  }

  if(length(ans)==0)
    return(NULL)

  do.call(rbind,ans)
}


############################################################
# 1. GLOBAL
############################################################

global_results <- list()
k <- 1

for(pred in c(
  "L1_expression_burden",
  "detection_breadth"
)){

  for(outcome in c(
    "density",
    "total_strength"
  )){

    tmp <- run_by_condition(
      global,
      "Global",
      "ALL",
      "ALL_L1FLnI",
      pred,
      outcome
    )

    if(!is.null(tmp)){
      global_results[[k]] <- tmp
      k <- k + 1
    }
  }
}

global_results <- do.call(
  rbind,
  global_results
)


############################################################
# 2. CELL-TYPE FAMILY-WIDE
############################################################

cell_results <- list()
k <- 1

for(ct in sort(unique(celltype$CellType))){

  d <- celltype[
    celltype$CellType == ct,
    ,
    drop=FALSE
  ]

  for(pred in c(
    "L1_expression_burden",
    "detection_breadth"
  )){

    for(outcome in c(
      "density",
      "total_strength"
    )){

      tmp <- run_by_condition(
        d,
        "Cell-type family",
        ct,
        "ALL_L1FLnI",
        pred,
        outcome
      )

      if(!is.null(tmp)){
        cell_results[[k]] <- tmp
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
# 3. LOCUS-SPECIFIC
############################################################

locus_results <- list()
k <- 1

groups <- unique(
  candidate[,c(
    "CellType",
    "TE_feature"
  )]
)

for(i in seq_len(nrow(groups))){

  ct <- groups$CellType[i]
  locus <- groups$TE_feature[i]

  d <- candidate[
    candidate$CellType == ct &
    candidate$TE_feature == locus,
    ,
    drop=FALSE
  ]

  for(pred in c(
    "mean_expression",
    "detection_rate"
  )){

    for(outcome in c(
      "density",
      "total_strength"
    )){

      tmp <- run_by_condition(
        d,
        "Locus-specific",
        ct,
        locus,
        pred,
        outcome
      )

      if(!is.null(tmp)){
        locus_results[[k]] <- tmp
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
# MASTER TABLE
############################################################

master <- rbind(
  global_results,
  cell_results,
  locus_results
)

# Dry n=3 explicitly flagged
master$Interpretation <- ifelse(
  master$n <= 3,
  "DESCRIPTIVE_n3",
  "Inferential"
)

# FDR separately within disease
master$FDR <- NA_real_

for(cond in unique(master$Condition)){

  idx <- which(
    master$Condition == cond
  )

  master$FDR[idx] <- p.adjust(
    master$P_value[idx],
    method="BH"
  )
}

master <- master[
  order(master$Condition,master$P_value),
]

write.csv(
  global_results,
  file.path(tabdir,"Global_byDisease.csv"),
  row.names=FALSE
)

write.csv(
  cell_results,
  file.path(tabdir,"CellType_byDisease.csv"),
  row.names=FALSE
)

write.csv(
  locus_results,
  file.path(tabdir,"Locus_byDisease.csv"),
  row.names=FALSE
)

write.csv(
  master,
  file.path(tabdir,"MASTER_byDisease.csv"),
  row.names=FALSE
)


############################################################
# PAPER FIGURE FUNCTION
############################################################

plot_by_disease <- function(
  dat,
  predictor,
  outcome,
  title,
  xlab,
  ylab,
  filename
){

  d <- dat[
    complete.cases(
      dat[,c(
        predictor,
        outcome,
        "Condition",
        "Donor"
      )]
    ),
  ]

  p <- ggplot(
    d,
    aes(
      x=.data[[predictor]],
      y=.data[[outcome]],
      color=Condition
    )
  ) +

    geom_point(
      size=3.4,
      alpha=0.9
    ) +

    # separate regression per condition
    geom_smooth(
      aes(group=Condition),
      method="lm",
      formula=y~x,
      se=FALSE,
      linewidth=0.9
    ) +

    labs(
      title=title,
      x=xlab,
      y=ylab,
      color=NULL
    ) +

    theme_classic(
      base_size=13
    ) +

    theme(
      plot.title=element_text(
        face="bold",
        size=14
      ),
      axis.title=element_text(
        face="bold"
      ),
      legend.position="top"
    )

  ggsave(
    paste0(filename,".pdf"),
    p,
    width=6.2,
    height=5.2
  )

  ggsave(
    paste0(filename,".png"),
    p,
    width=6.2,
    height=5.2,
    dpi=600
  )
}


############################################################
# KEY PAPER PANELS
############################################################

# Global
plot_by_disease(
  global,
  "L1_expression_burden",
  "density",
  "Retina-wide L1FLnI expression",
  "L1FLnI expression burden",
  "Network density",
  file.path(figdir,"A_Global_Density_byDisease")
)


# Photoreceptor
photo <- celltype[
  celltype$CellType=="Photoreceptor",
]

plot_by_disease(
  photo,
  "L1_expression_burden",
  "density",
  "Photoreceptor L1FLnI expression",
  "L1FLnI expression burden",
  "Network density",
  file.path(figdir,"B_Photoreceptor_Density_byDisease")
)


# RGC
rgc <- celltype[
  celltype$CellType=="RGC",
]

plot_by_disease(
  rgc,
  "detection_breadth",
  "density",
  "RGC L1FLnI detection breadth",
  "Fraction of L1FLnI loci detected",
  "Network density",
  file.path(figdir,"C_RGC_Density_byDisease")
)


# Astrocyte Xq23ab
astro <- candidate[
  candidate$CellType=="Astrocyte" &
  candidate$TE_feature=="TE-L1FLnI-Xq23ab",
]

plot_by_disease(
  astro,
  "detection_rate",
  "total_strength",
  "Astrocyte L1FLnI-Xq23ab",
  "Detection rate",
  "Total communication strength",
  file.path(figdir,"D_Astro_Xq23ab_byDisease")
)


# RPE 21q22.13b
rpe <- candidate[
  candidate$CellType=="RPE" &
  candidate$TE_feature=="TE-L1FLnI-21q22.13b",
]

plot_by_disease(
  rpe,
  "detection_rate",
  "density",
  "RPE L1FLnI-21q22.13b",
  "Detection rate",
  "Network density",
  file.path(figdir,"E_RPE_21q22_byDisease")
)


############################################################
# PRINT RESULTS
############################################################

cat("\n=========================================\n")
cat("GLOBAL RESULTS BY DISEASE\n")
cat("=========================================\n")

print(
  global_results,
  row.names=FALSE
)

cat("\n=========================================\n")
cat("TOP RESULTS BY DISEASE\n")
cat("=========================================\n")

print(
  head(
    master[
      order(master$P_value),
      c(
        "Analysis",
        "CellType",
        "Feature",
        "Predictor",
        "Outcome",
        "Condition",
        "n",
        "Pearson_r",
        "P_value",
        "FDR",
        "Interpretation"
      )
    ],
    30
  ),
  row.names=FALSE
)

cat("\nDONE\n")
cat("Tables:",tabdir,"\n")
cat("Figures:",figdir,"\n")
