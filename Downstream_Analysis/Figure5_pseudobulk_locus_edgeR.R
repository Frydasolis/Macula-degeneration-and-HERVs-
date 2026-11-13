#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
})

root <- "/storage/lemus_g/roldan/ARMD"

input <- file.path(
  root,
  "results/SRP413248/Seurat/downstream/SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

outdir <- file.path(
  root,
  "results/SRP413248/DE_HRCA_majorclass/TE_results/regional_enrichment/Figure5_pseudobulk_validation"
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object...\n")
obj <- readRDS(input)
meta <- obj@meta.data

stopifnot(all(c("sample", "Disease", "HRCA_majorclass") %in% colnames(meta)))

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
stopifnot(identical(colnames(counts), rownames(meta)))

te_features <- grep("^TE-", rownames(counts), value = TRUE)
cat("Total TE loci:", length(te_features), "\n")

meta$sample <- as.character(meta$sample)
meta$Disease <- as.character(meta$Disease)
meta$HRCA_majorclass <- as.character(meta$HRCA_majorclass)

valid <- c("Healthy Control", "Dry AMD", "Wet AMD")

meta <- meta[
  !is.na(meta$sample) &
  !is.na(meta$HRCA_majorclass) &
  meta$Disease %in% valid,
  ,
  drop = FALSE
]

# Minimum 20 nuclei per donor and cell type
donor_counts <- as.data.frame(
  table(meta$HRCA_majorclass, meta$sample),
  stringsAsFactors = FALSE
)

colnames(donor_counts) <- c("celltype", "sample", "n_nuclei")
donor_counts <- donor_counts[donor_counts$n_nuclei >= 20, ]

write.csv(
  donor_counts,
  file.path(outdir, "eligible_donor_celltypes.csv"),
  row.names = FALSE
)

comparisons <- list(
  Dry_vs_Healthy = c("Healthy Control", "Dry AMD"),
  Wet_vs_Healthy = c("Healthy Control", "Wet AMD")
)

summary_rows <- list()

for (ct in sort(unique(meta$HRCA_majorclass))) {

  cat("\nCELL TYPE:", ct, "\n")

  eligible <- donor_counts$sample[donor_counts$celltype == ct]

  cells <- rownames(meta)[
    meta$HRCA_majorclass == ct &
    meta$sample %in% eligible
  ]

  if (length(cells) == 0) next

  m <- meta[cells, , drop = FALSE]
  donors <- sort(unique(m$sample))

  indicator <- Matrix::sparseMatrix(
    i = seq_along(cells),
    j = match(m$sample, donors),
    x = 1,
    dims = c(length(cells), length(donors)),
    dimnames = list(cells, donors)
  )

  # Aggregate original counts per donor
  pb <- counts[, cells, drop = FALSE] %*% indicator

  info <- unique(m[, c("sample", "Disease"), drop = FALSE])
  rownames(info) <- info$sample
  info <- info[colnames(pb), , drop = FALSE]

  for (comparison in names(comparisons)) {

    groups <- comparisons[[comparison]]
    selected <- info$Disease %in% groups
    si <- info[selected, , drop = FALSE]

    group <- factor(si$Disease, levels = groups)
    n_donors <- table(group)

    if (any(n_donors < 3)) {
      cat("SKIP", comparison, ": insufficient donors\n")
      next
    }

    # TMM normalization using all RNA features
    y_all <- DGEList(counts = pb[, selected, drop = FALSE])
    y_all <- calcNormFactors(y_all, method = "TMM")

    # Retain TE loci and original library sizes
    y <- y_all[te_features, , keep.lib.sizes = TRUE]

    design <- model.matrix(~ group)
    keep <- filterByExpr(y, design = design)

    if (sum(keep) < 10) {
      cat("SKIP", comparison, ": insufficient expressed loci\n")
      next
    }

    y <- y[keep, , keep.lib.sizes = TRUE]

    y <- estimateDisp(y, design, robust = FALSE)
    fit <- glmQLFit(y, design, robust = FALSE)
    test <- glmQLFTest(fit, coef = 2)

    res <- topTags(test, n = Inf, sort.by = "none")$table

    res$feature <- rownames(res)
    res$p_val_adj <- p.adjust(res$PValue, method = "BH")
    res$celltype <- ct
    res$comparison <- comparison

    res <- res[, c(
      "feature", "celltype", "comparison",
      "logFC", "logCPM", "F", "PValue", "p_val_adj"
    )]

    safe_ct <- gsub("[^A-Za-z0-9_-]", "_", ct)

    write.csv(
      res,
      file.path(outdir, paste0("DE_", safe_ct, "_", comparison, ".csv")),
      row.names = FALSE
    )

    n_up <- sum(
      res$p_val_adj < 0.05 & res$logFC > 0.25,
      na.rm = TRUE
    )

    n_down <- sum(
      res$p_val_adj < 0.05 & res$logFC < -0.25,
      na.rm = TRUE
    )

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      celltype = ct,
      comparison = comparison,
      n_control = as.integer(n_donors[1]),
      n_case = as.integer(n_donors[2]),
      tested_loci = nrow(res),
      upregulated = n_up,
      downregulated = n_down
    )

    cat(
      comparison,
      "tested:", nrow(res),
      "up:", n_up,
      "down:", n_down, "\n"
    )
  }
}

if (length(summary_rows) > 0) {
  summary <- do.call(rbind, summary_rows)

  write.csv(
    summary,
    file.path(outdir, "pseudobulk_DE_summary.csv"),
    row.names = FALSE
  )

  print(summary, row.names = FALSE)
} else {
  warning("No eligible comparisons")
}

writeLines(
  capture.output(sessionInfo()),
  file.path(outdir, "sessionInfo.txt")
)

cat("\nANALYSIS COMPLETE\n")
