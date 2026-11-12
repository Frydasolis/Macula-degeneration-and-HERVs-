############################################################
# FIGURE 1A — SAMPLE-LEVEL TE BURDEN
# Retinal cell populations in AMD
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(lmerTest)
})

base <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/TE_burden_donor_analysis"
)

outdir <- file.path(base, "Figure1A_FINAL")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1. Load sample-level data
dat <- read.csv(
  file.path(base, "TE_burden_by_donor_celltype.csv"),
  stringsAsFactors = FALSE
)

celltypes <- c(
  "Amacrine",
  "Astrocyte",
  "Bipolar",
  "Horizontal",
  "Muller glia",
  "Photoreceptor",
  "RGC"
)

conditions <- c(
  "Healthy Control",
  "Dry AMD",
  "Wet AMD"
)

dat <- dat %>%
  filter(
    HRCA_majorclass %in% celltypes,
    Disease %in% conditions,
    n_cells >= 20
  ) %>%
  mutate(
    HRCA_majorclass = factor(
      HRCA_majorclass, levels = celltypes
    ),
    Disease = factor(
      Disease, levels = conditions
    )
  )

# 2. Load previously fitted mixed-effects model
model <- readRDS(
  file.path(
    base,
    "Mixed_model",
    "TE_burden_mixed_model.rds"
  )
)

anova_tab <- anova(
  model,
  type = 3,
  ddf = "Satterthwaite"
)

p_cell <- anova_tab[
  "HRCA_majorclass", "Pr(>F)"
]

p_disease <- anova_tab[
  "Disease", "Pr(>F)"
]

p_interaction <- anova_tab[
  "HRCA_majorclass:Disease", "Pr(>F)"
]

format_p <- function(x) {
  if (x < 0.001) {
    format(x, scientific = TRUE, digits = 3)
  } else {
    sprintf("%.3f", x)
  }
}

# 3. Colors
colors <- c(
  "Healthy Control" = "#607D8B",
  "Dry AMD" = "#E5A04B",
  "Wet AMD" = "#B85B69"
)

# 4. Global statistical annotation
subtitle_text <- paste0(
  "Cell type: p = ", format_p(p_cell),
  "     |     Disease: p = ", format_p(p_disease),
  "     |     Cell type × Disease: p = ",
  format_p(p_interaction)
)

# 5. Publication-quality figure
set.seed(2026)

p <- ggplot(
  dat,
  aes(
    x = Disease,
    y = mean_TE,
    fill = Disease
  )
) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.28,
    color = "#333333",
    outlier.shape = NA,
    linewidth = 0.45
  ) +
  geom_jitter(
    aes(color = Disease),
    width = 0.12,
    height = 0,
    size = 2.5,
    alpha = 0.95
  ) +
  facet_wrap(
    ~ HRCA_majorclass,
    ncol = 4,
    scales = "free_y",
    labeller = labeller(
      HRCA_majorclass = as_labeller(
        c(
          "Amacrine" = "Amacrine",
          "Astrocyte" = "Astrocyte",
          "Bipolar" = "Bipolar",
          "Horizontal" = "Horizontal",
          "Muller glia" = "Müller glia",
          "Photoreceptor" = "Photoreceptor",
          "RGC" = "RGC"
        )
      )
    )
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_x_discrete(
    labels = c(
      "Healthy Control" = "Healthy",
      "Dry AMD" = "Dry",
      "Wet AMD" = "Wet"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.12, 0.15))
  ) +
  labs(
    title = "TE-associated transcript burden across retinal cell populations",
    subtitle = subtitle_text,
    x = NULL,
    y = "Mean TE-associated transcript burden (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      size = 10,
      margin = margin(b = 15)
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.spacing = grid::unit(1.4, "lines"),
    plot.margin = margin(15, 15, 15, 15)
  )

# 6. Export
ggsave(
  file.path(outdir, "Figure1A_FINAL.pdf"),
  plot = p,
  width = 14,
  height = 8,
  device = "pdf"
)

ggsave(
  file.path(outdir, "Figure1A_FINAL.png"),
  plot = p,
  width = 14,
  height = 8,
  dpi = 600,
  bg = "white"
)

# 7. Export statistics
write.csv(
  as.data.frame(anova_tab),
  file.path(outdir, "Figure1A_mixed_model_ANOVA.csv")
)

cat("\nFIGURE 1A COMPLETED\n")
cat("Cell type p:", p_cell, "\n")
cat("Disease p:", p_disease, "\n")
cat("Interaction p:", p_interaction, "\n")
cat("Output:", outdir, "\n")
