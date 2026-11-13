library(dplyr)

root <- "results/SRP413248/TE_family_subclustering"
indir <- file.path(root, "Figure1B_pseudobulk")
outdir <- file.path(indir, "Validation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ab <- read.csv(file.path(indir,
  "TE_family_sample_abundance_CPM.csv"))
res <- read.csv(file.path(indir,
  "TE_family_edgeR_all_results.csv"))

# 1. Integrity of abundance table
stopifnot(all(c("TE_family","sample","CPM",
                "Cell_type","Disease") %in% names(ab)))
stopifnot(all(c("TE_family","Cell_type","Comparison",
                "PValue","FDR","FDR_celltype",
                "FDR_global") %in% names(res)))

cat("\n=== INTEGRITY ===\n")
cat("Families:", n_distinct(ab$TE_family), "\n")
cat("Donors:", n_distinct(ab$sample), "\n")
cat("Cell types:", n_distinct(ab$Cell_type), "\n")
cat("Missing CPM:", sum(is.na(ab$CPM)), "\n")
cat("Negative CPM:", sum(ab$CPM < 0, na.rm=TRUE), "\n")
cat("Duplicated family/sample/cell-type:",
    sum(duplicated(ab[c("TE_family","sample",
                        "Cell_type")])), "\n")

# 2. Sum of family CPM per pseudobulk
lib <- ab %>%
  group_by(sample, Disease, Cell_type) %>%
  summarise(
    total_family_CPM = sum(CPM, na.rm=TRUE),
    n_families = n_distinct(TE_family),
    .groups="drop"
  )

write.csv(lib, file.path(outdir,
  "CPM_sums_by_pseudobulk.csv"), row.names=FALSE)

cat("\n=== CPM SUMS ===\n")
print(summary(lib$total_family_CPM))
print(head(lib, 10), row.names=FALSE)

# 3. Validate BH corrections
res$BH_global_check <- p.adjust(res$PValue, "BH")

res <- res %>%
  group_by(Cell_type, Comparison) %>%
  mutate(BH_contrast_check = p.adjust(PValue, "BH")) %>%
  ungroup() %>%
  group_by(Cell_type) %>%
  mutate(BH_celltype_check = p.adjust(PValue, "BH")) %>%
  ungroup()

cat("\n=== FDR VALIDATION ===\n")
for (pair in list(
  c("FDR","BH_contrast_check"),
  c("FDR_celltype","BH_celltype_check"),
  c("FDR_global","BH_global_check")
)) {
  delta <- abs(res[[pair[1]]] - res[[pair[2]]])
  cat(pair[1], ": max difference =",
      max(delta, na.rm=TRUE), "\n")
}

write.csv(res, file.path(outdir,
  "edgeR_FDR_validation.csv"), row.names=FALSE)

# 4. Descriptive donor-level consistency
hc <- ab %>%
  filter(Disease == "Healthy Control") %>%
  group_by(TE_family, Cell_type) %>%
  summarise(
    HC_max = max(CPM, na.rm=TRUE),
    HC_median = median(CPM, na.rm=TRUE),
    n_HC = n_distinct(sample),
    .groups="drop"
  )

consistency <- ab %>%
  filter(Disease != "Healthy Control") %>%
  inner_join(hc, by=c("TE_family","Cell_type")) %>%
  group_by(TE_family, Cell_type, Disease) %>%
  summarise(
    n_donors = n_distinct(sample),
    mean_CPM = mean(CPM, na.rm=TRUE),
    median_CPM = median(CPM, na.rm=TRUE),
    HC_max = first(HC_max),
    HC_median = first(HC_median),
    n_HC = first(n_HC),
    n_above_HC_max = sum(CPM > HC_max, na.rm=TRUE),
    fraction_above_HC_max =
      mean(CPM > HC_max, na.rm=TRUE),
    .groups="drop"
  )

write.csv(consistency, file.path(outdir,
  "TE_family_donor_consistency.csv"),
  row.names=FALSE)

cat("\n=== HERVKC4 PHOTORECEPTOR ===\n")
print(consistency %>%
  filter(TE_family=="HERVKC4",
         Cell_type=="Photoreceptor"),
  row.names=FALSE)

cat("\n=== TOP EXPLORATORY EDGE R RESULTS ===\n")
print(res %>%
  arrange(PValue) %>%
  select(TE_family, Cell_type, Comparison,
         logFC, PValue, FDR, FDR_celltype,
         FDR_global) %>%
  head(15), row.names=FALSE)

cat("\n=== SIGNIFICANCE COUNTS ===\n")
cat("Contrast FDR < 0.05:",
    sum(res$FDR < 0.05, na.rm=TRUE), "\n")
cat("Cell-type FDR < 0.05:",
    sum(res$FDR_celltype < 0.05, na.rm=TRUE), "\n")
cat("Global FDR < 0.05:",
    sum(res$FDR_global < 0.05, na.rm=TRUE), "\n")

cat("\nValidation files:", outdir, "\n")
