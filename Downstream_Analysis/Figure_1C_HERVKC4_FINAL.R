suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

root <- "results/SRP413248/TE_family_subclustering"
cpm_file <- file.path(
  root, "Figure1B_pseudobulk",
  "TE_family_sample_abundance_CPM.csv"
)
de_file <- file.path(
  root, "Figure1B_DEFINITIVE",
  "TE_family_edgeR_FINAL.csv"
)
outdir <- file.path(root, "Figure1B_DEFINITIVE")

# Load donor-level CPM and final edgeR results
cpm <- read.csv(cpm_file, check.names=FALSE)
de <- read.csv(de_file, check.names=FALSE)

df <- cpm %>%
  filter(
    TE_family == "HERVKC4",
    Cell_type == "Photoreceptor"
  ) %>%
  mutate(
    Disease = factor(
      Disease,
      levels=c("Healthy Control","Dry AMD","Wet AMD")
    )
  ) %>%
  arrange(Disease, sample)

# Verify donor-level observations
stopifnot(
  nrow(df) == 16,
  !anyNA(df$Disease),
  !anyNA(df$CPM),
  !anyDuplicated(df$sample),
  identical(
    as.integer(table(df$Disease)),
    c(6L,3L,7L)
  )
)

# Verify edgeR results
herv <- de %>%
  filter(
    TE_family == "HERVKC4",
    Cell_type == "Photoreceptor"
  )

stopifnot(nrow(herv) == 3)

# Summary for plotting
summary_df <- df %>%
  group_by(Disease) %>%
  summarise(
    n=n(),
    median_CPM=median(CPM),
    mean_CPM=mean(CPM),
    Q1=quantile(CPM,0.25),
    Q3=quantile(CPM,0.75),
    .groups="drop"
  )

# Reproducible horizontal offsets for individual donors
df <- df %>%
  group_by(Disease) %>%
  mutate(
    x=as.numeric(Disease) +
      seq(-0.12,0.12,length.out=n())
  ) %>%
  ungroup()

cols <- c(
  "Healthy Control"="#64748B",
  "Dry AMD"="#D97706",
  "Wet AMD"="#2563EB"
)



set.seed(123)

p <- ggplot(
  df,
  aes(x=Disease, y=CPM, fill=Disease)
) +
  geom_boxplot(
    width=0.48,
    alpha=0.30,
    color="black",
    linewidth=0.6,
    outlier.shape=NA
  ) +
  geom_jitter(
    aes(color=Disease),
    width=0.13,
    height=0,
    size=2.8,
    alpha=0.95
  ) +
  scale_fill_manual(values=cols) +
  scale_color_manual(values=cols) +
  scale_x_discrete(
    labels=c(
      "Healthy control\n(n = 6)",
      "Dry AMD\n(n = 3)",
      "Wet AMD\n(n = 7)"
    )
  ) +
  scale_y_continuous(
    limits=c(0, NA),
    expand=expansion(mult=c(0,0.08))
  ) +
  labs(
    title="HERVKC4",
    subtitle="Photoreceptors",
    x=NULL,
    y="HERVKC4 abundance (CPM)"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title=element_text(face="bold",size=14),
    plot.subtitle=element_text(size=11),
    axis.text.x=element_text(
      size=10,color="black",
      lineheight=1.1,
      margin=margin(t=8)
    ),
    axis.text.y=element_text(color="black"),
    axis.title.y=element_text(margin=margin(r=8)),
    legend.position="none",
    plot.margin=margin(10,12,10,10)
  )

# Export publication-ready formats
ggsave(
  file.path(outdir,"Figure1C_HERVKC4.png"),
  plot=p,width=5.8,height=4.8,
  units="in",dpi=600,bg="white"
)

ggsave(
  file.path(outdir,"Figure1C_HERVKC4.pdf"),
  plot=p,width=5.8,height=4.8,
  units="in",device=cairo_pdf
)

write.csv(
  df,
  file.path(outdir,"Figure1C_HERVKC4_donor_data.csv"),
  row.names=FALSE
)

write.csv(
  summary_df,
  file.path(outdir,"Figure1C_HERVKC4_summary.csv"),
  row.names=FALSE
)

cat("\n=== HERVKC4 DONOR SUMMARY ===\n")
print(summary_df)

cat("\n=== FINAL edgeR RESULTS ===\n")
print(
  herv %>%
    select(
      Comparison,logFC,PValue,
      FDR_contrast,FDR_celltype,FDR_global
    )
)

cat("\nFiles generated in:\n",outdir,"\n")
