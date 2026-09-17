library(readr)
library(dplyr)

input_file <- "/storage/lemus_g/roldan/ARMD/results/SRP413248/CellChat/network_metrics/CellChat_ALL_celltypes_ALL_metrics_per_donor.csv"

output_file <- "/storage/lemus_g/roldan/ARMD/results/SRP413248/CellChat/network_metrics/CellChat_metrics_ONE_value_per_donor.csv"

df <- read_csv(input_file, show_col_types = FALSE)

centrality_cols <- c(
  "outdeg_unweighted_mean",
  "indeg_unweighted_mean",
  "outdeg_mean",
  "indeg_mean",
  "hub_mean",
  "authority_mean",
  "eigen_mean",
  "page_rank_mean",
  "betweenness_mean"
)

donor_metrics <- df %>%
  group_by(Donor, Condition) %>%
  summarise(
    Cell_types_available = n_distinct(Cell_type),

    across(
      all_of(centrality_cols),
      ~ mean(.x, na.rm = TRUE)
    ),

    .groups = "drop"
  )

write_csv(donor_metrics, output_file)

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")
cat("Output:\n")
cat(output_file, "\n\n")

print(donor_metrics)
