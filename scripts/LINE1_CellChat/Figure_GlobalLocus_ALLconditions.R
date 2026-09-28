############################################################
# GLOBAL LOCUS-SPECIFIC L1FLnI x CELLCHAT
# ALL CONDITIONS IN THE SAME GRAPH
#
# Each candidate locus is quantified across ALL retinal
# cells from each donor.
#
# Statistical unit = DONOR
# Pearson correlations calculated independently:
#   Healthy
#   Dry AMD
#   Wet AMD
#
# Dry AMD n=3 = DESCRIPTIVE
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

############################################################
# PATHS
############################################################

base <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/CellChat/donor_level_symbol/network_metrics/",
  "LINE1_GlobalLocus_ByDisease"
)

expr_file <- file.path(
  base,
  "Global_Locus_Donor_Expression.csv"
)

stats_file <- file.path(
  base,
  "Global_Locus_Pearson_ByDisease.csv"
)

outdir <- file.path(
  base,
  "figures_ALL_conditions"
)

dir.create(
  outdir,
  recursive=TRUE,
  showWarnings=FALSE
)

############################################################
# READ DATA
############################################################

expr <- read.csv(
  expr_file,
  stringsAsFactors=FALSE
)

stats <- read.csv(
  stats_file,
  stringsAsFactors=FALSE
)

############################################################
# STANDARDIZE CONDITIONS
############################################################

normalize_condition <- function(x){

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

  x
}

expr$Condition <- normalize_condition(
  expr$Condition
)

stats$Condition <- normalize_condition(
  stats$Condition
)

condition_order <- c(
  "Healthy",
  "Dry AMD",
  "Wet AMD"
)

expr$Condition <- factor(
  expr$Condition,
  levels=condition_order
)

############################################################
# COLORS
############################################################

condition_colors <- c(
  "Healthy" = "#4C78A8",
  "Dry AMD" = "#E6A532",
  "Wet AMD" = "#D9534F"
)

############################################################
# LINE TYPES
#
# Dry dashed to emphasize descriptive n=3 analysis
############################################################

condition_lines <- c(
  "Healthy" = "solid",
  "Dry AMD" = "dashed",
  "Wet AMD" = "solid"
)

condition_shapes <- c(
  "Healthy" = 16,
  "Dry AMD" = 17,
  "Wet AMD" = 15
)

############################################################
# HELPERS
############################################################

safe_name <- function(x){

  x <- sub(
    "^TE-",
    "",
    x
  )

  gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    x
  )
}

pretty_feature <- function(x){

  sub(
    "^TE-",
    "",
    x
  )
}

pretty_predictor <- function(x){

  if(x == "mean_expression"){
    return(
      "Mean normalized locus expression"
    )
  }

  if(x == "detection_rate"){
    return(
      "Locus detection rate"
    )
  }

  x
}

pretty_outcome <- function(x){

  if(x == "density"){
    return(
      "Network density"
    )
  }

  if(x == "total_strength"){
    return(
      "Total communication strength"
    )
  }

  x
}

############################################################
# STATISTICS LABEL
############################################################

make_stats_label <- function(
  feature,
  predictor,
  outcome
){

  lines <- character()

  for(cond in condition_order){

    z <- stats[
      stats$Feature == feature &
      stats$Predictor == predictor &
      stats$Outcome == outcome &
      stats$Condition == cond,
      ,
      drop=FALSE
    ]

    if(nrow(z) == 0){
      next
    }

    z <- z[1,,drop=FALSE]

    if(cond == "Dry AMD"){

      line <- sprintf(
        "%s: n=%d | r=%.3f | P=%.3g | FDR=%.3g [descriptive]",
        cond,
        z$n,
        z$Pearson_r,
        z$P_value,
        z$FDR
      )

    } else {

      line <- sprintf(
        "%s: n=%d | r=%.3f | P=%.3g | FDR=%.3g",
        cond,
        z$n,
        z$Pearson_r,
        z$P_value,
        z$FDR
      )
    }

    lines <- c(
      lines,
      line
    )
  }

  paste(
    lines,
    collapse="\n"
  )
}

############################################################
# LOCI
############################################################

loci <- unique(
  as.character(
    expr$Feature
  )
)

cat(
  "\nNumber of loci: ",
  length(loci),
  "\n",
  sep=""
)

print(loci)

############################################################
# GENERATE FIGURES
############################################################

plot_count <- 0

for(feature in loci){

  cat(
    "\n========================================\n",
    pretty_feature(feature),
    "\n",
    "========================================\n",
    sep=""
  )

  d_feature <- expr[
    expr$Feature == feature,
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

      ######################################################
      # DATA
      ######################################################

      d <- d_feature[
        is.finite(d_feature[[pred]]) &
        is.finite(d_feature[[outcome]]),
        ,
        drop=FALSE
      ]

      ######################################################
      # STATISTICS
      ######################################################

      stats_label <- make_stats_label(
        feature=feature,
        predictor=pred,
        outcome=outcome
      )

      ######################################################
      # BASE FIGURE
      ######################################################

      p <- ggplot(
        d,
        aes(
          x=.data[[pred]],
          y=.data[[outcome]],
          color=Condition,
          shape=Condition
        )
      ) +

        geom_point(
          size=3.8,
          alpha=0.90
        )

      ######################################################
      # HEALTHY + WET
      # regression + 95% CI
      ######################################################

      d_hw <- d[
        d$Condition %in%
          c(
            "Healthy",
            "Wet AMD"
          ),
        ,
        drop=FALSE
      ]

      if(nrow(d_hw) > 0){

        p <- p +

          geom_smooth(
            data=d_hw,
            aes(
              color=Condition,
              linetype=Condition,
              group=Condition
            ),
            method="lm",
            formula=y ~ x,
            se=TRUE,
            linewidth=0.9,
            alpha=0.15
          )
      }

      ######################################################
      # DRY AMD
      # regression only — NO confidence interval
      ######################################################

      d_dry <- d[
        d$Condition == "Dry AMD",
        ,
        drop=FALSE
      ]

      if(nrow(d_dry) >= 3){

        p <- p +

          geom_smooth(
            data=d_dry,
            aes(
              color=Condition,
              linetype=Condition,
              group=Condition
            ),
            method="lm",
            formula=y ~ x,
            se=FALSE,
            linewidth=0.9
          )
      }

      ######################################################
      # STYLE
      ######################################################

      p <- p +

        scale_color_manual(
          values=condition_colors,
          drop=FALSE
        ) +

        scale_shape_manual(
          values=condition_shapes,
          drop=FALSE
        ) +

        scale_linetype_manual(
          values=condition_lines,
          drop=FALSE
        ) +

        labs(
          title=pretty_feature(feature),

          subtitle=stats_label,

          x=pretty_predictor(pred),

          y=pretty_outcome(outcome),

          color="Condition",
          shape="Condition",
          linetype="Condition"
        ) +

        theme_classic(
          base_size=13
        ) +

        theme(
          plot.title=element_text(
            face="bold",
            size=15
          ),

          plot.subtitle=element_text(
            size=9.5,
            lineheight=1.20,
            margin=margin(
              b=8
            )
          ),

          axis.title=element_text(
            face="bold"
          ),

          legend.position="right",

          legend.title=element_text(
            face="bold"
          )
        )

      ######################################################
      # SAVE
      ######################################################

      locus_dir <- file.path(
        outdir,
        safe_name(feature)
      )

      dir.create(
        locus_dir,
        recursive=TRUE,
        showWarnings=FALSE
      )

      filename <- paste0(
        safe_name(feature),
        "_ALL_conditions_",
        pred,
        "_vs_",
        outcome
      )

      ggsave(
        file.path(
          locus_dir,
          paste0(
            filename,
            ".png"
          )
        ),
        p,
        width=8.0,
        height=6.3,
        dpi=400
      )

      ggsave(
        file.path(
          locus_dir,
          paste0(
            filename,
            ".pdf"
          )
        ),
        p,
        width=8.0,
        height=6.3
      )

      plot_count <- plot_count + 1

      cat(
        "Saved: ",
        pred,
        " -> ",
        outcome,
        "\n",
        sep=""
      )
    }
  }
}

############################################################
# COPY COMPLETE STATISTICAL TABLE
############################################################

write.csv(
  stats,
  file.path(
    outdir,
    "Global_Locus_ALL_conditions_statistics.csv"
  ),
  row.names=FALSE
)

############################################################
# DONE
############################################################

cat("\n========================================\n")
cat("ALL-CONDITION LOCUS FIGURES COMPLETE\n")
cat("========================================\n")

cat(
  "Loci: ",
  length(loci),
  "\n",
  sep=""
)

cat(
  "Figures generated: ",
  plot_count,
  "\n",
  sep=""
)

cat(
  "Expected if 12 loci: 48\n"
)

cat(
  "\nOutput:\n",
  outdir,
  "\n",
  sep=""
)

