
library(CellChat)
library(Seurat)
library(dplyr)
library(readr)

input_dir <- "/storage/lemus_g/roldan/ARMD/results/SRP413248/CellChat/donor_level_symbol"

output_dir <- "/storage/lemus_g/roldan/ARMD/results/SRP413248/CellChat/network_metrics"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

files <- list.files(
  input_dir,
  pattern = "^CellChat_GSM.*\\.rds$",
  full.names = TRUE
)

# ============================================================
# Conditions
# ============================================================

get_condition <- function(donor) {

  if (donor %in% c(
    "GSM6841143",
    "GSM6841144",
    "GSM6841145",
    "GSM6841146",
    "GSM6841147",
    "GSM6841148"
  )) return("Healthy")

  if (donor %in% c(
    "GSM6841149",
    "GSM6841150",
    "GSM6841151",
    "GSM6841152",
    "GSM6841153",
    "GSM6841154",
    "GSM6841155"
  )) return("Wet")

  if (donor %in% c(
    "GSM6841156",
    "GSM6841157",
    "GSM6841159"
  )) return("Dry")

  return(NA_character_)
}

# ============================================================
# CellChat centrality metrics
# ============================================================

metric_names <- c(
  "outdeg_unweighted",
  "indeg_unweighted",
  "outdeg",
  "indeg",
  "hub",
  "authority",
  "eigen",
  "page_rank",
  "betweenness"
)

results <- list()

# ============================================================
# LOOP OVER DONORS
# ============================================================

for (f in files) {

  donor <- sub(
    "^CellChat_(GSM[0-9]+)\\.rds$",
    "\\1",
    basename(f)
  )

  condition <- get_condition(donor)

  cat("\n")
  cat("============================================\n")
  cat("Processing:", donor, "\n")
  cat("Condition:", condition, "\n")
  cat("============================================\n")

  cellchat <- readRDS(f)

  # ----------------------------------------------------------
  # Cell types
  # ----------------------------------------------------------

  cell_types <- levels(cellchat@idents)

  cell_counts <- table(
    as.character(cellchat@idents)
  )

  cat(
    "Cell types:",
    length(cell_types),
    "\n"
  )

  # ----------------------------------------------------------
  # Network matrices
  # ----------------------------------------------------------

  count_matrix <- cellchat@net$count
  weight_matrix <- cellchat@net$weight

  # ----------------------------------------------------------
  # Compute centrality
  # ----------------------------------------------------------

  cellchat <- netAnalysis_computeCentrality(
    cellchat,
    slot.name = "netP"
  )

  centrality <- cellchat@netP$centr

  cat(
    "Pathways:",
    length(centrality),
    "\n"
  )

  donor_results <- list()

  # ==========================================================
  # LOOP OVER CELL TYPES
  # ==========================================================

  for (cell_type in cell_types) {

    # --------------------------------------------------------
    # Number of cells
    # --------------------------------------------------------

    n_cells <- as.numeric(
      cell_counts[cell_type]
    )

    # --------------------------------------------------------
    # Count
    # --------------------------------------------------------

    count_out <- sum(
      count_matrix[cell_type, ],
      na.rm = TRUE
    )

    count_in <- sum(
      count_matrix[, cell_type],
      na.rm = TRUE
    )

    # --------------------------------------------------------
    # Weight
    # --------------------------------------------------------

    weight_out <- sum(
      weight_matrix[cell_type, ],
      na.rm = TRUE
    )

    weight_in <- sum(
      weight_matrix[, cell_type],
      na.rm = TRUE
    )

    # --------------------------------------------------------
    # Store pathway-level centrality
    # --------------------------------------------------------

    metric_values <- list()

    for (metric in metric_names) {

      metric_values[[metric]] <- numeric(0)

    }

    # ========================================================
    # LOOP OVER PATHWAYS
    # ========================================================

    for (pathway in names(centrality)) {

      pathway_centr <- centrality[[pathway]]

      if (is.null(pathway_centr)) {
        next
      }

      # ------------------------------------------------------
      # Check that this metric exists
      # ------------------------------------------------------

      for (metric in metric_names) {

        if (!metric %in% names(pathway_centr)) {
          next
        }

        metric_vector <- pathway_centr[[metric]]

        # ----------------------------------------------------
        # Check cell type
        # ----------------------------------------------------

        if (!cell_type %in% names(metric_vector)) {
          next
        }

        value <- metric_vector[cell_type]

        if (
          length(value) == 1 &&
          is.finite(value)
        ) {

          metric_values[[metric]] <-
            c(
              metric_values[[metric]],
              as.numeric(value)
            )
        }
      }
    }

    # ========================================================
    # MEAN ACROSS PATHWAYS
    # ========================================================

    metric_means <- lapply(
      metric_values,
      function(x) {

        if (length(x) == 0) {

          return(NA_real_)

        }

        mean(
          x,
          na.rm = TRUE
        )
      }
    )

    # ========================================================
    # ONE ROW
    # ========================================================

    donor_results[[cell_type]] <- data.frame(

      Donor = donor,

      Condition = condition,

      Cell_type = cell_type,

      Cell_cells = n_cells,

      count_out = count_out,

      count_in = count_in,

      weight_out = weight_out,

      weight_in = weight_in,

      outdeg_unweighted_mean =
        metric_means$outdeg_unweighted,

      indeg_unweighted_mean =
        metric_means$indeg_unweighted,

      outdeg_mean =
        metric_means$outdeg,

      indeg_mean =
        metric_means$indeg,

      hub_mean =
        metric_means$hub,

      authority_mean =
        metric_means$authority,

      eigen_mean =
        metric_means$eigen,

      page_rank_mean =
        metric_means$page_rank,

      betweenness_mean =
        metric_means$betweenness,

      stringsAsFactors = FALSE
    )
  }

  results[[donor]] <-
    bind_rows(donor_results)

  cat(
    "Completed:",
    donor,
    "\n"
  )
}

# ============================================================
# COMBINE
# ============================================================

all_metrics <- bind_rows(results)

# ============================================================
# ORDER
# ============================================================

all_metrics <- all_metrics %>%
  mutate(
    Condition = factor(
      Condition,
      levels = c(
        "Healthy",
        "Dry",
        "Wet"
      )
    )
  ) %>%
  arrange(
    Condition,
    Donor,
    Cell_type
  )

all_metrics$Condition <-
  as.character(all_metrics$Condition)

# ============================================================
# SAVE
# ============================================================

output_file <- file.path(
  output_dir,
  "CellChat_ALL_celltypes_ALL_metrics_per_donor.csv"
)

write_csv(
  all_metrics,
  output_file
)

# ============================================================
# SUMMARY
# ============================================================

cat("\n")
cat("====================================================\n")
cat("CELLCHAT METRICS COMPLETED\n")
cat("====================================================\n")

cat(
  "Donors:",
  length(unique(all_metrics$Donor)),
  "\n"
)

cat(
  "Cell types:",
  length(unique(all_metrics$Cell_type)),
  "\n"
)

cat(
  "Rows:",
  nrow(all_metrics),
  "\n"
)

cat(
  "Columns:",
  ncol(all_metrics),
  "\n"
)

cat("\nOutput:\n")
cat(output_file, "\n")

cat("====================================================\n")

print(all_metrics)

