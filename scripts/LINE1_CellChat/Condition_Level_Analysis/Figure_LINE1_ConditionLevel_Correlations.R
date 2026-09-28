############################################################
# Figure_LINE1_ConditionLevel_Correlations.R
#
# Condition-level association between global locus-specific
# L1FLnI expression and aggregated CellChat network metrics.
#
# THREE observations per analysis:
#   Healthy Control
#   Dry AMD
#   Wet AMD
#
# IMPORTANT:
# These correlations are DESCRIPTIVE (n = 3 conditions).
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

############################################################
# PATHS
############################################################

infile <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "CellChat/by_condition/",
  "LINE1_GlobalLocus_condition_level/",
  "Global_Locus_ConditionLevel.csv"
)

outdir <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/SRP413248/",
  "CellChat/by_condition/",
  "LINE1_GlobalLocus_condition_level/",
  "figures"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# LOAD DATA
############################################################

dat <- read.csv(
  infile,
  stringsAsFactors = FALSE
)

dat$Condition <- factor(
  dat$Condition,
  levels = c(
    "Healthy Control",
    "Dry AMD",
    "Wet AMD"
  )
)

############################################################
# COLORS / SHAPES
############################################################

condition_colors <- c(
  "Healthy Control" = "#4C78A8",
  "Dry AMD"         = "#E6A532",
  "Wet AMD"         = "#D9534F"
)

condition_shapes <- c(
  "Healthy Control" = 16,
  "Dry AMD"         = 17,
  "Wet AMD"         = 15
)

############################################################
# PRETTY LABELS
############################################################

clean_locus <- function(x){
  sub("^TE-", "", x)
}

############################################################
# FUNCTION: CALCULATE CORRELATION
############################################################

get_stats <- function(df, predictor, outcome){

  x <- df[[predictor]]
  y <- df[[outcome]]

  keep <- is.finite(x) & is.finite(y)

  x <- x[keep]
  y <- y[keep]

  if(length(x) < 3 ||
     length(unique(x)) < 2 ||
     length(unique(y)) < 2){

    return(
      data.frame(
        n = length(x),
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
    n = length(x),
    Pearson_r = r,
    R2 = r^2,
    P_value = test$p.value
  )
}

############################################################
# FUNCTION: MAKE ONE PANEL
############################################################

make_panel <- function(df,
                       predictor,
                       outcome,
                       xlab,
                       ylab){

  st <- get_stats(
    df,
    predictor,
    outcome
  )

  stat_text <- if(is.na(st$Pearson_r)){

    paste0(
      "n = ", st$n,
      " conditions\nNot estimated"
    )

  } else {

    paste0(
      "Pearson r = ",
      sprintf("%.3f", st$Pearson_r),
      "\nR² = ",
      sprintf("%.3f", st$R2),
      "\nP = ",
      format.pval(
        st$P_value,
        digits = 2,
        eps = 0.001
      ),
      "\nn = ",
      st$n,
      " conditions"
    )
  }

  ggplot(
    df,
    aes(
      x = .data[[predictor]],
      y = .data[[outcome]]
    )
  ) +

    ########################################################
    # DESCRIPTIVE REGRESSION LINE
    ########################################################

    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "grey35",
      linewidth = 0.8
    ) +

    ########################################################
    # CONDITION POINTS
    ########################################################

    geom_point(
      aes(
        color = Condition,
        shape = Condition
      ),
      size = 4
    ) +

    ########################################################
    # CONDITION LABELS
    ########################################################

    geom_text_repel(
      aes(
        label = Condition,
        color = Condition
      ),
      size = 3.4,
      show.legend = FALSE,
      box.padding = 0.5,
      point.padding = 0.4,
      min.segment.length = 0
    ) +

    scale_color_manual(
      values = condition_colors
    ) +

    scale_shape_manual(
      values = condition_shapes
    ) +

    labs(
      x = xlab,
      y = ylab,
      subtitle = stat_text
    ) +

    theme_classic(
      base_size = 12
    ) +

    theme(
      legend.position = "none",
      plot.subtitle = element_text(
        size = 10
      ),
      axis.title = element_text(
        face = "bold"
      )
    )
}

############################################################
# RUN ALL 12 LOCI
############################################################

all_stats <- list()
stat_counter <- 1

loci <- unique(dat$Feature)

cat("\n============================================\n")
cat("CONDITION-LEVEL FIGURES\n")
cat("============================================\n")
cat("Number of loci:", length(loci), "\n\n")

for(feature in loci){

  df <- dat[
    dat$Feature == feature,
    ,
    drop = FALSE
  ]

  pretty_feature <- clean_locus(feature)

  ##########################################################
  # PANEL A
  # mean expression -> density
  ##########################################################

  p1 <- make_panel(
    df,
    predictor = "mean_expression",
    outcome = "density",
    xlab = "Mean normalized expression",
    ylab = "Network density"
  )

  ##########################################################
  # PANEL B
  # detection -> density
  ##########################################################

  p2 <- make_panel(
    df,
    predictor = "detection_rate",
    outcome = "density",
    xlab = "Detection rate",
    ylab = "Network density"
  )

  ##########################################################
  # PANEL C
  # mean expression -> strength
  ##########################################################

  p3 <- make_panel(
    df,
    predictor = "mean_expression",
    outcome = "total_strength",
    xlab = "Mean normalized expression",
    ylab = "Total communication strength"
  )

  ##########################################################
  # PANEL D
  # detection -> strength
  ##########################################################

  p4 <- make_panel(
    df,
    predictor = "detection_rate",
    outcome = "total_strength",
    xlab = "Detection rate",
    ylab = "Total communication strength"
  )

  ##########################################################
  # COMBINE
  ##########################################################

  combined <- (
    p1 + p2
  ) / (
    p3 + p4
  ) +

    plot_annotation(
      title = pretty_feature,
      subtitle = paste0(
        "Aggregated condition-level analysis | ",
        "Healthy Control, Dry AMD, Wet AMD | ",
        "n = 3 conditions; descriptive"
      ),
      tag_levels = "A"
    )

  ##########################################################
  # FILE-SAFE LOCUS NAME
  ##########################################################

  safe_name <- gsub(
    "[^A-Za-z0-9._-]",
    "_",
    pretty_feature
  )

  ##########################################################
  # SAVE PNG
  ##########################################################

  ggsave(
    filename = file.path(
      outdir,
      paste0(
        safe_name,
        "_ConditionLevel_Correlation.png"
      )
    ),
    plot = combined,
    width = 11,
    height = 9,
    dpi = 400
  )

  ##########################################################
  # SAVE PDF
  ##########################################################

  ggsave(
    filename = file.path(
      outdir,
      paste0(
        safe_name,
        "_ConditionLevel_Correlation.pdf"
      )
    ),
    plot = combined,
    width = 11,
    height = 9
  )

  ##########################################################
  # COLLECT STATISTICS
  ##########################################################

  combinations <- data.frame(
    Predictor = c(
      "mean_expression",
      "detection_rate",
      "mean_expression",
      "detection_rate"
    ),
    Outcome = c(
      "density",
      "density",
      "total_strength",
      "total_strength"
    ),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(combinations))){

    st <- get_stats(
      df,
      combinations$Predictor[i],
      combinations$Outcome[i]
    )

    all_stats[[stat_counter]] <- data.frame(
      Feature = feature,
      Predictor = combinations$Predictor[i],
      Outcome = combinations$Outcome[i],
      n = st$n,
      Pearson_r = st$Pearson_r,
      R2 = st$R2,
      P_value = st$P_value,
      Analysis = "Condition-level descriptive",
      stringsAsFactors = FALSE
    )

    stat_counter <- stat_counter + 1
  }

  cat("Saved:", pretty_feature, "\n")
}

############################################################
# SAVE STATISTICS
############################################################

stats <- do.call(
  rbind,
  all_stats
)

stats$FDR <- p.adjust(
  stats$P_value,
  method = "BH"
)

stats <- stats[
  order(stats$P_value),
]

write.csv(
  stats,
  file.path(
    outdir,
    "ConditionLevel_Correlation_Statistics.csv"
  ),
  row.names = FALSE
)

############################################################
# PRINT TOP RESULTS
############################################################

cat("\n============================================\n")
cat("TOP DESCRIPTIVE CORRELATIONS\n")
cat("============================================\n\n")

print(
  head(
    stats,
    20
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("OUTPUT\n")
cat("============================================\n")

cat("Figures:", outdir, "\n")
cat("Expected PNG:", length(loci), "\n")
cat("Expected PDF:", length(loci), "\n")

cat("\nDONE\n")
