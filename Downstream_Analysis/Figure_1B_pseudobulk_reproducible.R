suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
})

root <- "results/SRP413248/TE_family_subclustering"
old_dir <- file.path(root, "Figure1B_pseudobulk")
out <- file.path(root, "Figure1B_reproducible")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

te <- readRDS(file.path(root, "TE_family_by_cell_sparse.rds"))

seu <- readRDS(
  "results/SRP413248/Seurat/downstream/SRP413248_GeneTE_QC_UMAP_HRCA_annotated.rds"
)

stopifnot(identical(colnames(te), colnames(seu)))

md <- seu@meta.data
md$sample <- as.character(md$sample)
md$Disease <- as.character(md$Disease)
md$Cell_type <- as.character(md$HRCA_majorclass)

# Align metadata to matrix
md <- md[colnames(te), , drop=FALSE]

# 1. Construct donor x cell-type pseudobulks
md$group_id <- paste(md$sample, md$Cell_type, sep="__")
group_levels <- unique(md$group_id)

design_sparse <- sparse.model.matrix(
  ~ 0 + factor(md$group_id, levels=group_levels)
)
colnames(design_sparse) <- group_levels

pb <- te %*% design_sparse
pb <- as.matrix(pb)
storage.mode(pb) <- "integer"

pb_meta <- md %>%
  distinct(group_id, sample, Disease, Cell_type) %>%
  as.data.frame()

pb_meta <- pb_meta[
  match(colnames(pb), pb_meta$group_id), , drop=FALSE
]

pb_meta$n_cells <- as.integer(table(
  factor(md$group_id, levels=group_levels)
))

# 2. Audit original CPM denominators
old <- read.csv(file.path(
  old_dir, "TE_family_sample_abundance_CPM.csv"
))

pb_long <- as.data.frame(as.table(pb))
names(pb_long) <- c("TE_family","group_id","raw_count")

pb_long <- left_join(
  pb_long, pb_meta,
  by="group_id"
)

audit <- left_join(
  pb_long,
  old,
  by=c("TE_family","sample","Cell_type","Disease")
)

audit <- audit %>%
  filter(!is.na(CPM), raw_count > 0, CPM > 0) %>%
  mutate(
    implied_library_size = raw_count * 1e6 / CPM
  )

audit_summary <- audit %>%
  group_by(group_id, sample, Cell_type) %>%
  summarise(
    median_implied_library_size =
      median(implied_library_size),
    min_implied_library_size =
      min(implied_library_size),
    max_implied_library_size =
      max(implied_library_size),
    .groups="drop"
  )

rna_totals <- md %>%
  group_by(group_id) %>%
  summarise(
    total_nCount_RNA=sum(nCount_RNA),
    total_nCount_TE=sum(nCount_TE),
    .groups="drop"
  )

audit_summary <- left_join(
  audit_summary, rna_totals, by="group_id"
)

write.csv(audit_summary,
  file.path(out,"CPM_denominator_audit.csv"),
  row.names=FALSE)

cat("\n=== CPM DENOMINATOR AUDIT ===\n")
print(head(audit_summary, 12), row.names=FALSE)

# 3. edgeR by retinal cell type
celltypes <- c(
  "Amacrine", "Astrocyte", "Bipolar",
  "Horizontal", "Muller glia",
  "Photoreceptor", "RGC"
)

contrasts <- list(
  Dry_vs_Healthy=c(-1,1,0),
  Wet_vs_Healthy=c(-1,0,1),
  Wet_vs_Dry=c(0,-1,1)
)

all_results <- list()
k <- 1

for (ct in celltypes) {
  message("Processing: ", ct)

  idx <- which(pb_meta$Cell_type == ct)
  counts <- pb[,idx,drop=FALSE]
  info <- pb_meta[idx,,drop=FALSE]

  info$Disease <- factor(
    info$Disease,
    levels=c("Healthy Control","Dry AMD","Wet AMD")
  )

  design <- model.matrix(~0+Disease, data=info)
  colnames(design) <- c("Healthy","Dry","Wet")

  if (qr(design)$rank < ncol(design)) {
    warning("Skipping ", ct, ": incomplete design")
    next
  }

  y <- DGEList(counts=counts)

  # Filter using the experimental design
  keep <- filterByExpr(y, design=design)
  y <- y[keep,,keep.lib.sizes=FALSE]

  if (nrow(y) < 2) {
    warning("Skipping ", ct, ": insufficient TE families")
    next
  }

  y <- calcNormFactors(y, method="TMM")

  y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)

  contrast_matrix <- makeContrasts(
    Dry_vs_Healthy=Dry-Healthy,
    Wet_vs_Healthy=Wet-Healthy,
    Wet_vs_Dry=Wet-Dry,
    levels=design
  )

  for (contrast_name in colnames(contrast_matrix)) {
    qlf <- glmQLFTest(
      fit,
      contrast=contrast_matrix[,contrast_name]
    )

    tab <- topTags(
      qlf,
      n=Inf,
      sort.by="none"
    )$table

    tab$TE_family <- rownames(tab)
    tab$Cell_type <- ct
    tab$Comparison <- contrast_name

    all_results[[k]] <- tab
    k <- k+1
  }
}

results <- bind_rows(all_results)

results <- results %>%
  group_by(Cell_type, Comparison) %>%
  mutate(FDR_contrast=p.adjust(PValue,"BH")) %>%
  ungroup() %>%
  group_by(Cell_type) %>%
  mutate(FDR_celltype=p.adjust(PValue,"BH")) %>%
  ungroup() %>%
  mutate(FDR_global=p.adjust(PValue,"BH"))

write.csv(results,
  file.path(out,"TE_family_edgeR_reproducible.csv"),
  row.names=FALSE)

write.csv(pb_meta,
  file.path(out,"Pseudobulk_metadata.csv"),
  row.names=FALSE)

saveRDS(pb,
  file.path(out,"TE_family_pseudobulk_raw_counts.rds"))

cat("\n=== RESULTS ===\n")
cat("Tests:",nrow(results),"\n")
cat("Contrast FDR < 0.05:",
    sum(results$FDR_contrast<0.05),"\n")
cat("Cell-type FDR < 0.05:",
    sum(results$FDR_celltype<0.05),"\n")
cat("Global FDR < 0.05:",
    sum(results$FDR_global<0.05),"\n")

cat("\n=== TOP RESULTS ===\n")
print(
  results %>%
    arrange(PValue) %>%
    select(TE_family,Cell_type,Comparison,
           logFC,PValue,FDR_contrast,FDR_global) %>%
    head(15),
  row.names=FALSE
)

cat("\n=== HERVKC4 ===\n")
print(
  results %>%
    filter(TE_family=="HERVKC4",
           Cell_type=="Photoreceptor"),
  row.names=FALSE
)

cat("\nSaved to:",out,"\n")
