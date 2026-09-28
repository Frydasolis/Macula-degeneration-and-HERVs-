############################################################
# Figure_LINE1_ByDisease_Paper.R
# Publication figure:
# L1FLnI vs CellChat metrics stratified by disease
#
# Statistical unit = DONOR
# Healthy n=6
# Wet AMD n=7
# Dry AMD n=3 -> descriptive only
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

############################################################
# PATHS
############################################################

base <- "network_metrics/LINE1_Pearson/tables"

outdir <- "network_metrics/LINE1_ByDisease/figures/Paper"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# READ DATA
############################################################

global <- read.csv(
  file.path(base, "family_global_L1FLnI_burden.csv"),
  stringsAsFactors = FALSE
)

celltype <- read.csv(
  file.path(base, "family_celltype_L1FLnI_burden.csv"),
  stringsAsFactors = FALSE
)

candidate <- read.csv(
  file.path(base, "candidate_locus_donor_expression.csv"),
  stringsAsFactors = FALSE
)

############################################################
# STANDARDIZE CONDITION
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

global$Condition    <- fix_condition(global$Condition)
celltype$Condition  <- fix_condition(celltype$Condition)
candidate$Condition <- fix_condition(candidate$Condition)

############################################################
# COLORS
############################################################

condition_colors <- c(
  "Healthy" = "#3B82C4",
  "Dry AMD" = "#E69F00",
  "Wet AMD" = "#C94C4C"
)

############################################################
# PEARSON FUNCTION
############################################################

get_cor <- function(dat, x, y) {

  d <- dat[
    complete.cases(dat[, c(x, y)]),
    ,
    drop = FALSE
  ]

  n <- nrow(d)

  if (n < 3) {
    return(
      data.frame(
        n = n,
        r = NA,
        p = NA
      )
    )
  }

  if (
    sd(d[[x]]) == 0 ||
    sd(d[[y]]) == 0
  ) {
    return(
      data.frame(
        n = n,
        r = NA,
        p = NA
      )
    )
  }

  z <- cor.test(
    d[[x]],
    d[[y]],
    method = "pearson"
  )

  data.frame(
    n = n,
    r = unname(z$estimate),
    p = z$p.value
  )
}

############################################################
# FORMAT P
############################################################

format_p <- function(p) {

  if (is.na(p))
    return("NA")

  if (p < 0.001)
    return(formatC(
      p,
      format = "e",
      digits = 1
    ))

  sprintf("%.3f", p)
}

############################################################
# CREATE STATS FOR ANNOTATION
############################################################

make_stats <- function(dat, x, y) {

  conditions <- c(
    "Healthy",
    "Dry AMD",
    "Wet AMD"
  )

  ans <- list()

  for (cond in conditions) {

    d <- dat[
      dat$Condition == cond,
      ,
      drop = FALSE
    ]

    z <- get_cor(d, x, y)

    ans[[cond]] <- data.frame(
      Condition = cond,
      n = z$n,
      r = z$r,
      p = z$p
    )
  }

  do.call(rbind, ans)
}

# PANEL FUNCTION
############################################################

make_panel <- function(
  dat,
  x,
  y,
  title,
  xlab,
  ylab,
  panel_letter
) {

  dat <- dat[
    complete.cases(
      dat[, c(x, y, "Condition", "Donor")]
    ),
    ,
    drop = FALSE
  ]

  stats <- make_stats(dat, x, y)

  stat_text <- function(cond, descriptive = FALSE) {

    z <- stats[stats$Condition == cond, , drop = FALSE]

    if(nrow(z) == 0 || z$n < 3 || is.na(z$r)) {

      nvalue <- ifelse(nrow(z) == 0, 0, z$n)

      out <- paste0(
        cond,
        ": correlation not estimated, n = ",
        nvalue
      )

    } else {

      out <- paste0(
        cond,
        ": r = ",
        sprintf("%.2f", z$r),
        ", P = ",
        format_p(z$p),
        ", n = ",
        z$n
      )
    }

    if(descriptive)
      out <- paste0(out, " (descriptive)")

    out
  }

  subtitle <- paste0(
    stat_text("Healthy"),
    "\n",
    stat_text("Dry AMD", TRUE),
    "\n",
    stat_text("Wet AMD")
  )

  regression_main <- dat[
    dat$Condition %in% c("Healthy", "Wet AMD"),
    ,
    drop = FALSE
  ]

  regression_dry <- dat[
    dat$Condition == "Dry AMD",
    ,
    drop = FALSE
  ]

  p <- ggplot(
    dat,
    aes(
      x = .data[[x]],
      y = .data[[y]],
      color = Condition
    )
  ) +

    geom_point(
      size = 3.4,
      alpha = 0.95
    ) +

    geom_smooth(
      data = regression_main,
      aes(group = Condition),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.9,
      alpha = 0.12
    ) +

    geom_smooth(
      data = regression_dry,
      aes(group = Condition),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.9,
      linetype = "dashed"
    ) +

    scale_color_manual(
      values = condition_colors,
      drop = FALSE
    ) +

    labs(
      title = paste0(panel_letter, "   ", title),
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      color = NULL
    ) +

    theme_classic(base_size = 11) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      plot.subtitle = element_text(
        size = 8.5,
        lineheight = 1.15
      ),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = 9),
      legend.position = "top",
      legend.text = element_text(size = 9),
      plot.margin = margin(8, 10, 8, 8)
    )

  return(p)
}

# A — GLOBAL L1FLnI
############################################################

A <- make_panel(
  global,
  "L1_expression_burden",
  "density",

  "Retina-wide L1FLnI",

  "L1FLnI expression burden",
  "Network density",

  "A"
)

############################################################
# B — PHOTORECEPTOR
############################################################

photo <- celltype[
  celltype$CellType == "Photoreceptor",
]

B <- make_panel(
  photo,
  "L1_expression_burden",
  "density",

  "Photoreceptor L1FLnI",

  "L1FLnI expression burden",
  "Network density",

  "B"
)

############################################################
# C — BIPOLAR
############################################################

bipolar <- celltype[
  celltype$CellType == "Bipolar",
]

C <- make_panel(
  bipolar,
  "detection_breadth",
  "total_strength",

  "Bipolar L1FLnI",

  "Fraction of L1FLnI loci detected",
  "Total communication strength",

  "C"
)

############################################################
# D — ASTROCYTE Xq23ab
############################################################

astro <- candidate[
  candidate$CellType == "Astrocyte" &
  candidate$TE_feature == "TE-L1FLnI-Xq23ab",
]

D <- make_panel(
  astro,
  "detection_rate",
  "total_strength",

  "Astrocyte L1FLnI-Xq23ab",

  "Detection rate",
  "Total communication strength",

  "D"
)

############################################################
# E — RPE 21q22.13b
############################################################

rpe21 <- candidate[
  candidate$CellType == "RPE" &
  candidate$TE_feature ==
    "TE-L1FLnI-21q22.13b",
]

E <- make_panel(
  rpe21,
  "detection_rate",
  "density",

  "RPE L1FLnI-21q22.13b",

  "Detection rate",
  "Network density",

  "E"
)

############################################################
# F — RPE 3q13.13ma
############################################################

rpe3q <- candidate[
  candidate$CellType == "RPE" &
  candidate$TE_feature ==
    "TE-L1FLnI-3q13.13ma",
]

F <- make_panel(
  rpe3q,
  "detection_rate",
  "total_strength",

  "RPE L1FLnI-3q13.13ma",

  "Detection rate",
  "Total communication strength",

  "F"
)

############################################################
# SAVE INDIVIDUAL PANELS
############################################################

panels <- list(
  A_Global = A,
  B_Photoreceptor = B,
  C_Bipolar = C,
  D_Astro_Xq23ab = D,
  E_RPE_21q22 = E,
  F_RPE_3q13ma = F
)

for (nm in names(panels)) {

  ggsave(
    file.path(
      outdir,
      paste0(nm, ".pdf")
    ),
    panels[[nm]],
    width = 6.5,
    height = 5
  )

  ggsave(
    file.path(
      outdir,
      paste0(nm, ".png")
    ),
    panels[[nm]],
    width = 6.5,
    height = 5,
    dpi = 600
  )
}

############################################################
# COMBINE 6 PANELS
############################################################

# Try patchwork if installed

if (
  requireNamespace(
    "patchwork",
    quietly = TRUE
  )
) {

  suppressPackageStartupMessages(
    library(patchwork)
  )

  combined <-
    (A | B) /
    (C | D) /
    (E | F) +

    plot_layout(
      guides = "collect"
    ) &

    theme(
      legend.position = "top"
    )

  ggsave(
    file.path(
      outdir,
      "Figure_LINE1_CellChat_ByDisease.pdf"
    ),
    combined,
    width = 13,
    height = 14
  )

  ggsave(
    file.path(
      outdir,
      "Figure_LINE1_CellChat_ByDisease.png"
    ),
    combined,
    width = 13,
    height = 14,
    dpi = 600
  )
}

############################################################
# SAVE PANEL STATISTICS
############################################################

extract_panel_stats <- function(
  dat,
  x,
  y,
  panel
) {

  z <- make_stats(
    dat,
    x,
    y
  )

  z$Panel <- panel
  z$Predictor <- x
  z$Outcome <- y

  z
}

stats_all <- rbind(

  extract_panel_stats(
    global,
    "L1_expression_burden",
    "density",
    "A_Global"
  ),

  extract_panel_stats(
    photo,
    "L1_expression_burden",
    "density",
    "B_Photoreceptor"
  ),

  extract_panel_stats(
    bipolar,
    "detection_breadth",
    "total_strength",
    "C_Bipolar"
  ),

  extract_panel_stats(
    astro,
    "detection_rate",
    "total_strength",
    "D_Astro_Xq23ab"
  ),

  extract_panel_stats(
    rpe21,
    "detection_rate",
    "density",
    "E_RPE_21q22"
  ),

  extract_panel_stats(
    rpe3q,
    "detection_rate",
    "total_strength",
    "F_RPE_3q13ma"
  )
)

write.csv(
  stats_all,
  file.path(
    outdir,
    "Figure_panel_statistics.csv"
  ),
  row.names = FALSE
)

############################################################
# DONE
############################################################

cat("\n============================================\n")
cat("PAPER FIGURE COMPLETE\n")
cat("============================================\n")

cat(
  "\nSaved to:\n",
  normalizePath(outdir),
  "\n"
)

cat("\nPanel statistics:\n")

print(
  stats_all,
  row.names = FALSE
)
