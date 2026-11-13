#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

root <- "/storage/lemus_g/roldan/ARMD"

input <- file.path(
  root,
  "results/SRP413248/Seurat/downstream/",
  "SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

resdir <- file.path(
  root,
  "results/SRP413248/DE_HRCA_majorclass/",
  "TE_results/regional_enrichment/",
  "Figure5_pseudobulk_validation"
)

outdir <- file.path(resdir, "HARLEQUIN_figure")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

locus <- "TE-HARLEQUIN-Xp11.23"
celltypes <- c("Horizontal", "RGC", "Photoreceptor")

disease_levels <- c(
  "Healthy Control", "Dry AMD", "Wet AMD"
)

colors <- c(
  "Healthy Control" = "#4A4A4A",
  "Dry AMD" = "#C65A48",
  "Wet AMD" = "#437C9D"
)

theme_pub <- theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.margin = margin(10, 12, 10, 10)
  )

# ==================================================
# PANEL A: DONOR-LEVEL EXPRESSION
# ==================================================

cat("Loading Seurat object...\n")

obj <- readRDS(input)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
meta <- obj@meta.data

stopifnot(locus %in% rownames(counts))
stopifnot(identical(colnames(counts), rownames(meta)))

expression_data <- list()

for (ct in celltypes) {

  cat("Aggregating:", ct, "\n")

  cells <- rownames(meta)[
    meta$HRCA_majorclass == ct &
    meta$Disease %in% disease_levels
  ]

  m <- meta[cells, , drop = FALSE]

  donor_n <- table(m$sample)
  eligible <- names(donor_n[donor_n >= 20])

  m <- m[m$sample %in% eligible, , drop = FALSE]
  cells <- rownames(m)

  donors <- sort(unique(as.character(m$sample)))

  indicator <- sparseMatrix(
    i = seq_along(cells),
    j = match(as.character(m$sample), donors),
    x = 1,
    dims = c(length(cells), length(donors)),
    dimnames = list(cells, donors)
  )

  pb <- counts[, cells, drop = FALSE] %*% indicator

  info <- unique(m[, c("sample", "Disease")])
  rownames(info) <- info$sample
  info <- info[colnames(pb), , drop = FALSE]

  y <- DGEList(counts = pb)
  y <- calcNormFactors(y, method = "TMM")

  normalized <- cpm(y, normalized.lib.sizes = TRUE)

  expression_data[[ct]] <- data.frame(
    celltype = ct,
    sample = colnames(pb),
    Disease = as.character(info$Disease),
    n_nuclei = as.integer(donor_n[colnames(pb)]),
    raw_counts = as.numeric(pb[locus, ]),
    TMM_CPM = as.numeric(normalized[locus, ])
  )
}

expression_data <- do.call(rbind, expression_data)

expression_data$celltype <- factor(
  expression_data$celltype,
  levels = celltypes
)

expression_data$Disease <- factor(
  expression_data$Disease,
  levels = disease_levels
)

write.csv(
  expression_data,
  file.path(outdir, "HARLEQUIN_donor_TMM_CPM.csv"),
  row.names = FALSE
)

# ==================================================
# PANEL A: ACTUAL edgeR P-VALUES
# ==================================================

annotations <- list()

for (ct in celltypes) {

  for (comparison in c(
    "Dry_vs_Healthy",
    "Wet_vs_Healthy"
  )) {

    f <- file.path(
      resdir,
      paste0("DE_", ct, "_", comparison, ".csv")
    )

    d <- read.csv(f, check.names = FALSE)

    hit <- d[
      d$feature == locus,
      ,
      drop = FALSE
    ]

    if (nrow(hit) == 0) {
      cat("Locus not tested in:", ct, comparison, "\n")
      next
    }

    if (nrow(hit) != 1) {
      stop("Duplicated locus in: ", f)
    }

    annotations[[length(annotations) + 1]] <- data.frame(
      celltype = ct,
      comparison = comparison,
      p = hit$PValue,
      FDR = hit$p_val_adj
    )
  }
}

if (length(annotations) == 0) {
  stop("No statistical results found for locus")
}

annotations <- do.call(rbind, annotations)

write.csv(
  annotations,
  file.path(outdir, "HARLEQUIN_panel_A_statistics.csv"),
  row.names = FALSE
)

annotations$celltype <- factor(
  annotations$celltype,
  levels = celltypes
)

annotations$label <- paste0(
  "p = ",
  format.pval(
    annotations$p,
    digits = 2,
    eps = 1e-5
  )
)

max_y <- aggregate(
  TMM_CPM ~ celltype,
  data = expression_data,
  FUN = max
)

annotations$max_y <- max_y$TMM_CPM[
  match(annotations$celltype, max_y$celltype)
]

annotations$y <- pmax(
  annotations$max_y,
  0.01
) * ifelse(
  annotations$comparison == "Dry_vs_Healthy",
  1.25,
  1.65
)

annotations$xmin <- 1

annotations$xmax <- ifelse(
  annotations$comparison == "Dry_vs_Healthy",
  2,
  3
)


pA <- ggplot(
  expression_data,
  aes(
    x = Disease,
    y = TMM_CPM,
    fill = Disease
  )
) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.25,
    color = "black",
    linewidth = 0.45,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = Disease),
    position = position_jitter(
      width = 0.10,
      height = 0
    ),
    size = 2.7,
    alpha = 0.95
  ) +
  geom_segment(
    data = annotations,
    aes(
      x = xmin,
      xend = xmax,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.4
  ) +
  geom_text(
    data = annotations,
    aes(
      x = (xmin + xmax) / 2,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3.2
  ) +
  facet_wrap(
    ~celltype,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.05, 0.85)
    )
  ) +
  labs(
    title = "A  HARLEQUIN-Xp11.23 expression",
    subtitle = "Donor-level pseudobulk expression; p-values from edgeR",
    x = NULL,
    y = "TMM-normalized CPM"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      size = 9
    )
  )

# ==================================================
# PANEL B: edgeR LOG2 FOLD CHANGES
# ==================================================

effects <- subset(
  annotations,
  comparison == "Dry_vs_Healthy"
)

effects$logFC <- NA_real_

for (i in seq_len(nrow(effects))) {

  ct <- as.character(effects$celltype[i])

  f <- file.path(
    resdir,
    paste0("DE_", ct, "_Dry_vs_Healthy.csv")
  )

  d <- read.csv(f, check.names = FALSE)

  effects$logFC[i] <- d$logFC[
    match(locus, d$feature)
  ]
}

pB <- ggplot(
  effects,
  aes(x = celltype, y = logFC)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey65"
  ) +
  geom_point(
    size = 3,
    color = "#C65A48"
  ) +
  geom_text(
    aes(label = sprintf("%.2f", logFC)),
    hjust = -0.3,
    size = 3.5
  ) +
  coord_flip(ylim = c(-10, 1)) +
  labs(
    title = "B  Differential expression",
    x = NULL,
    y = "log2FC (Dry AMD vs Healthy)"
  ) +
  theme_pub

# ==================================================
# PANEL C: LEAVE-ONE-DONOR-OUT
# ==================================================

lodo <- read.csv(
  file.path(
    resdir,
    "HARLEQUIN_LODO_results.csv"
  )
)

lodo$celltype <- factor(
  lodo$celltype,
  levels = celltypes
)

full <- subset(
  lodo,
  excluded_donor == "NONE"
)

loo <- subset(
  lodo,
  excluded_donor != "NONE"
)

pC <- ggplot(
  loo,
  aes(
    x = celltype,
    y = logFC
  )
) +
  geom_point(
    position = position_jitter(
      width = 0.12,
      height = 0
    ),
    size = 2.5,
    alpha = 0.8,
    color = "#437C9D"
  ) +
  geom_point(
    data = full,
    aes(x = celltype, y = logFC),
    inherit.aes = FALSE,
    shape = 18,
    size = 4,
    color = "#C65A48"
  ) +
  coord_flip() +
  labs(
    title = "C  Leave-one-donor-out",
    x = NULL,
    y = "Estimated log2FC"
  ) +
  theme_pub

# ==================================================
# EXPORT
# ==================================================

ggsave(
  file.path(outdir, "Panel_A_expression_pvalues.pdf"),
  pA,
  width = 11,
  height = 4.5
)

ggsave(
  file.path(outdir, "Panel_B_edgeR_effects.pdf"),
  pB,
  width = 5.5,
  height = 4.2
)

ggsave(
  file.path(outdir, "Panel_C_LODO.pdf"),
  pC,
  width = 5.5,
  height = 4.2
)

if (requireNamespace("patchwork", quietly = TRUE)) {

  combined <- pA / (pB | pC) +
    patchwork::plot_layout(
      heights = c(1.15, 1)
    )

  ggsave(
    file.path(
      outdir,
      "Supplementary_Figure_HARLEQUIN.pdf"
    ),
    combined,
    width = 12,
    height = 9
  )

  ggsave(
    file.path(
      outdir,
      "Supplementary_Figure_HARLEQUIN.png"
    ),
    combined,
    width = 12,
    height = 9,
    dpi = 300
  )

} else {
  cat("patchwork not installed; individual panels saved\n")
}

cat("\nFIGURE COMPLETE\n")
cat("Output:", outdir, "\n")
