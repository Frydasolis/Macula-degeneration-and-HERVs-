# ============================================================
# Figure 1A: TE-associated transcript burden in AMD
#
# Complete sample-level analysis:
# - Mixed-effects model
# - Type III ANOVA with Satterthwaite approximation
# - Tukey-adjusted post hoc comparisons
# - Model diagnostics
# - Publication-quality Figure 1A
#
# Input: TE_burden_by_donor_celltype.csv
# ============================================================


############################################################
# TE BURDEN MIXED-EFFECTS ANALYSIS
# Figure 1A - Retinal cell populations
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggplot2)
})

# ----------------------------------------------------------
# 1. DIRECTORIES
# ----------------------------------------------------------

base <- Sys.getenv(
  "TE_BURDEN_DATA_DIR",
  unset = file.path(
    "results",
    "SRP413248",
    "TE_burden_donor_analysis"
  )
)

outdir <- file.path(base, "Mixed_model")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# 2. LOAD DATA
# ----------------------------------------------------------

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

# Minimum nuclei per sample and cell type
min_cells <- 20

dat <- dat %>%
  filter(
    HRCA_majorclass %in% celltypes,
    Disease %in% conditions,
    n_cells >= min_cells
  ) %>%
  mutate(
    sample = factor(sample),
    HRCA_majorclass = factor(
      HRCA_majorclass,
      levels = celltypes
    ),
    Disease = factor(
      Disease,
      levels = conditions
    )
  )

cat("\n=== SAMPLE COVERAGE ===\n")

sample_coverage <- dat %>%
  count(Disease, HRCA_majorclass)

print(as.data.frame(sample_coverage), row.names = FALSE)

# ----------------------------------------------------------
# 3. MIXED-EFFECTS MODEL
# ----------------------------------------------------------

model <- lmer(
  mean_TE ~ HRCA_majorclass * Disease + (1 | sample),
  data = dat,
  REML = TRUE
)

cat("\n=== MODEL SUMMARY ===\n")
print(summary(model))

cat("\n=== TYPE III ANOVA ===\n")

anova_results <- anova(
  model,
  type = 3,
  ddf = "Satterthwaite"
)

print(anova_results)

write.csv(
  as.data.frame(anova_results),
  file.path(outdir, "Mixed_model_ANOVA.csv")
)

# ----------------------------------------------------------
# 4. POST HOC: CELL TYPES WITHIN EACH CONDITION
# ----------------------------------------------------------

emm_cell <- emmeans(
  model,
  ~ HRCA_majorclass | Disease
)

posthoc_cell <- as.data.frame(
  pairs(
    emm_cell,
    adjust = "tukey"
  )
)

cat("\n=== CELL TYPE COMPARISONS ===\n")

print(posthoc_cell, row.names = FALSE)

write.csv(
  posthoc_cell,
  file.path(outdir, "Posthoc_celltypes_by_disease.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# 5. POST HOC: DISEASE WITHIN EACH CELL TYPE
# ----------------------------------------------------------

emm_disease <- emmeans(
  model,
  ~ Disease | HRCA_majorclass
)

posthoc_disease <- as.data.frame(
  pairs(
    emm_disease,
    adjust = "tukey"
  )
)

cat("\n=== DISEASE COMPARISONS ===\n")

print(posthoc_disease, row.names = FALSE)

write.csv(
  posthoc_disease,
  file.path(outdir, "Posthoc_disease_by_celltype.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# 6. ESTIMATED MARGINAL MEANS
# ----------------------------------------------------------

emm_table <- as.data.frame(
  emmeans(
    model,
    ~ HRCA_majorclass * Disease
  )
)

write.csv(
  emm_table,
  file.path(outdir, "Estimated_marginal_means.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# 7. DIAGNOSTICS
# ----------------------------------------------------------

pdf(
  file.path(outdir, "Model_diagnostics.pdf"),
  width = 10,
  height = 5
)

par(mfrow = c(1, 2))

plot(
  fitted(model),
  resid(model),
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs fitted"
)

abline(h = 0, lty = 2)

qqnorm(
  resid(model),
  main = "Normal Q-Q plot"
)

qqline(resid(model))

dev.off()

# ----------------------------------------------------------
# 8. DESCRIPTIVE FIGURE
# ----------------------------------------------------------

p <- ggplot(
  dat,
  aes(
    x = HRCA_majorclass,
    y = mean_TE,
    color = Disease
  )
) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.7
    ),
    size = 2,
    alpha = 0.85
  ) +
  stat_summary(
    aes(group = Disease),
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.7),
    size = 3,
    shape = 18
  ) +
  labs(
    x = "Retinal cell population",
    y = "Mean TE-associated transcript burden (%)",
    color = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  file.path(outdir, "TE_burden_mixed_model.png"),
  p,
  width = 12,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(outdir, "TE_burden_mixed_model.pdf"),
  p,
  width = 12,
  height = 6
)

# ----------------------------------------------------------
# 9. SAVE MODEL
# ----------------------------------------------------------

saveRDS(
  model,
  file.path(outdir, "TE_burden_mixed_model.rds")
)

cat("\n====================================\n")
cat("MIXED MODEL ANALYSIS COMPLETED\n")
cat("====================================\n")

cat("\nOutput directory:\n", outdir, "\n")


# ============================================================
# FINAL FIGURE 1A
# Reuse the model fitted above; no external model RDS required.
# ============================================================

outdir <- file.path(base, "Figure1A_FINAL")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

anova_tab <- anova_results

p_cell <- anova_tab["HRCA_majorclass", "Pr(>F)"]
p_disease <- anova_tab["Disease", "Pr(>F)"]
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
