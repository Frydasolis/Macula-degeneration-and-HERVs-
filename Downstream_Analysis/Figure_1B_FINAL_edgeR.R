suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
})

root <- "results/SRP413248/TE_family_subclustering"
indir <- file.path(root, "Figure1B_reproducible")
outdir <- file.path(root, "Figure1B_DEFINITIVE")

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Original reconstructed pseudobulk counts
pb <- readRDS(file.path(
  indir, "TE_family_pseudobulk_raw_counts.rds"
))

meta <- read.csv(file.path(
  indir, "Pseudobulk_metadata.csv"
))


# Match metadata to pseudobulk columns
stopifnot(all(colnames(pb) %in% meta$group_id))
meta <- meta[match(colnames(pb), meta$group_id), , drop=FALSE]
stopifnot(identical(colnames(pb), meta$group_id))

# Obtain library sizes directly from the original Seurat object
library(Seurat)

seu <- readRDS(
  "results/SRP413248/Seurat/downstream/SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

md <- seu@meta.data
md$group_id <- paste(
  md$sample,
  md$HRCA_majorclass,
  sep="__"
)

rna_lib <- md %>%
  group_by(group_id) %>%
  summarise(
    total_nCount_RNA = sum(nCount_RNA),
    total_nCount_TE = sum(nCount_TE),
    .groups="drop"
  )

meta <- meta %>%
  left_join(rna_lib, by="group_id")

stopifnot(
  !anyNA(meta$total_nCount_RNA),
  all(meta$total_nCount_RNA > 0),
  all(colSums(pb) <= meta$total_nCount_RNA)
)

cat("\n=== LIBRARY SIZE VALIDATION ===\n")
cat("Pseudobulks:", nrow(meta), "\n")
cat("Missing RNA libraries:",
    sum(is.na(meta$total_nCount_RNA)), "\n")
cat("RNA library size range:",
    range(meta$total_nCount_RNA), "\n")

rm(seu, md)
gc()

celltypes <- c(
  "Amacrine", "Astrocyte", "Bipolar",
  "Horizontal", "Muller glia",
  "Photoreceptor", "RGC"
)

all_results <- list()
filter_log <- list()

for (ct in celltypes) {

  message("\nProcessing ", ct)

  idx <- which(meta$Cell_type == ct)

  counts <- pb[, idx, drop=FALSE]
  info <- meta[idx, , drop=FALSE]

  info$Disease <- factor(
    info$Disease,
    levels=c("Healthy Control","Dry AMD","Wet AMD")
  )

  design <- model.matrix(~0+Disease, data=info)
  colnames(design) <- c("Healthy","Dry","Wet")

  if (qr(design)$rank < ncol(design)) {
    warning("Skipping ",ct,": incomplete design")
    next
  }

  # Counts remain TE-family counts
  y <- DGEList(counts=counts)

  # Use RNA-total exposure, not TE-only library sizes
  y$samples$lib.size <- info$total_nCount_RNA
  y$samples$norm.factors <- rep(1, ncol(y))

  # Expression filtering based on RNA-total library sizes
  keep <- filterByExpr(
    y,
    design=design
  )

  filter_log[[ct]] <- data.frame(
    Cell_type=ct,
    TE_family=rownames(y),
    retained=keep
  )

  y <- y[keep, , keep.lib.sizes=TRUE]

  if (nrow(y) < 2) {
    warning("Skipping ",ct,": insufficient families")
    next
  }

  # Explicit offset based on total RNA counts
  y$offset <- matrix(
    log(info$total_nCount_RNA),
    nrow=nrow(y),
    ncol=ncol(y),
    byrow=TRUE
  )

  y <- estimateDisp(
    y, design,
    robust=TRUE
  )

  fit <- glmQLFit(
    y, design,
    robust=TRUE
  )

  contrast_matrix <- makeContrasts(
    Dry_vs_Healthy=Dry-Healthy,
    Wet_vs_Healthy=Wet-Healthy,
    Wet_vs_Dry=Wet-Dry,
    levels=design
  )

  for (contrast_name in colnames(contrast_matrix)) {

    test <- glmQLFTest(
      fit,
      contrast=contrast_matrix[,contrast_name]
    )

    tab <- topTags(
      test, n=Inf, sort.by="none"
    )$table

    tab$TE_family <- rownames(tab)
    tab$Cell_type <- ct
    tab$Comparison <- contrast_name

    all_results[[paste(ct,contrast_name)]] <- tab
  }
}

results <- bind_rows(all_results)

results <- results %>%
  group_by(Cell_type,Comparison) %>%
  mutate(FDR_contrast=p.adjust(PValue,"BH")) %>%
  ungroup() %>%
  group_by(Cell_type) %>%
  mutate(FDR_celltype=p.adjust(PValue,"BH")) %>%
  ungroup() %>%
  mutate(FDR_global=p.adjust(PValue,"BH"))

filter_table <- bind_rows(filter_log)

write.csv(
  results,
  file.path(outdir,"TE_family_edgeR_FINAL.csv"),
  row.names=FALSE
)

write.csv(
  filter_table,
  file.path(outdir,"TE_family_filtering.csv"),
  row.names=FALSE
)

write.csv(
  meta,
  file.path(outdir,"Pseudobulk_metadata_FINAL.csv"),
  row.names=FALSE
)

cat("\n=== FINAL ANALYSIS ===\n")
cat("Tests:", nrow(results), "\n")
cat("FDR contrast < 0.05:",
    sum(results$FDR_contrast < 0.05), "\n")
cat("FDR cell type < 0.05:",
    sum(results$FDR_celltype < 0.05), "\n")
cat("FDR global < 0.05:",
    sum(results$FDR_global < 0.05), "\n")

cat("\n=== TOP 15 RESULTS ===\n")
print(
  results %>%
    arrange(PValue) %>%
    select(
      TE_family,Cell_type,Comparison,
      logFC,PValue,FDR_contrast,
      FDR_celltype,FDR_global
    ) %>%
    head(15),
  n=15,
  width=120
)

cat("\n=== HERVKC4 PHOTORECEPTORS ===\n")
print(
  results %>%
    filter(
      TE_family=="HERVKC4",
      Cell_type=="Photoreceptor"
    ) %>%
    select(
      TE_family,Comparison,logFC,
      PValue,FDR_contrast,
      FDR_celltype,FDR_global
    ),
  width=120
)

cat("\nSaved to:",outdir,"\n")
