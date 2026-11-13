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

locus <- "TE-HARLEQUIN-Xp11.23"
celltypes <- c("Horizontal", "RGC", "Photoreceptor")

cat("Loading Seurat object...\n")
obj <- readRDS(input)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
meta <- obj@meta.data

stopifnot(locus %in% rownames(counts))
stopifnot(identical(colnames(counts), rownames(meta)))

results <- list()

for (ct in celltypes) {

  cat("\nProcessing:", ct, "\n")

  cells <- rownames(meta)[
    meta$HRCA_majorclass == ct &
    meta$Disease %in% c("Healthy Control", "Dry AMD")
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

  for (excluded in c("NONE", donors)) {

    selected <- if (excluded == "NONE") {
      rep(TRUE, length(donors))
    } else {
      donors != excluded
    }

    group <- factor(
      info$Disease[selected],
      levels = c("Healthy Control", "Dry AMD")
    )

    if (any(table(group) < 2)) next

    y <- DGEList(counts = pb[, selected, drop = FALSE])
    y <- calcNormFactors(y, method = "TMM")

    design <- model.matrix(~ group)

    # Fit all expressed TE loci using the same
    # normalization strategy as the original analysis.
    te_features <- grep("^TE-", rownames(y), value = TRUE)
    y <- y[te_features, , keep.lib.sizes = TRUE]

    keep <- filterByExpr(y, design = design)
    keep[locus] <- TRUE

    y <- y[keep, , keep.lib.sizes = TRUE]

    y <- estimateDisp(y, design, robust = FALSE)
    fit <- glmQLFit(y, design, robust = FALSE)
    test <- glmQLFTest(fit, coef = 2)

    tab <- topTags(test, n = Inf, sort.by = "none")$table

    if (!locus %in% rownames(tab)) next

    results[[length(results) + 1]] <- data.frame(
      celltype = ct,
      excluded_donor = excluded,
      excluded_group = if (excluded == "NONE") {
        "NONE"
      } else {
        as.character(info[excluded, "Disease"])
      },
      n_healthy = sum(group == "Healthy Control"),
      n_dry = sum(group == "Dry AMD"),
      logFC = tab[locus, "logFC"],
      PValue = tab[locus, "PValue"],
      FDR = p.adjust(tab$PValue, method = "BH")[
        match(locus, rownames(tab))
      ]
    )

    cat(
      "Excluded:", excluded,
      "logFC:", round(tab[locus, "logFC"], 3),
      "p:", signif(tab[locus, "PValue"], 3),
      "\n"
    )
  }
}

out <- do.call(rbind, results)

write.csv(
  out,
  file.path(outdir, "HARLEQUIN_LODO_results.csv"),
  row.names = FALSE
)

print(out, row.names = FALSE)
cat("\nLODO COMPLETE\n")
