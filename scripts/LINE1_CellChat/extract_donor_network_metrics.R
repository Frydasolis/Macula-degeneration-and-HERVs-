suppressPackageStartupMessages({
  library(CellChat)
})

############################################################
# OUTPUT DIRECTORY
############################################################

outdir <- "network_metrics"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# FIND ALL DONOR CELLCHAT OBJECTS
############################################################

files <- list.files(
  pattern = "^CellChat_GSM[0-9]+\\.rds$",
  full.names = TRUE
)

cat("\n============================================\n")
cat("DONOR-LEVEL CELLCHAT NETWORK METRICS\n")
cat("============================================\n")
cat("Objects found:", length(files), "\n\n")

############################################################
# CONDITION MAP
############################################################

get_condition <- function(donor) {

  healthy <- paste0("GSM68411", 43:48)
  wet     <- paste0("GSM68411", 49:55)
  dry     <- c(
    "GSM6841156",
    "GSM6841157",
    "GSM6841159"
  )

  if (donor %in% healthy) {
    return("Healthy Control")
  }

  if (donor %in% wet) {
    return("Wet AMD")
  }

  if (donor %in% dry) {
    return("Dry AMD")
  }

  return(NA_character_)
}

############################################################
# STORAGE
############################################################

network_results <- list()
celltype_results <- list()
pathway_results <- list()

############################################################
# LOOP THROUGH DONORS
############################################################

for (f in files) {

  donor <- sub(
    "^CellChat_(GSM[0-9]+)\\.rds$",
    "\\1",
    basename(f)
  )

  condition <- get_condition(donor)

  cat("--------------------------------------------\n")
  cat("Processing:", donor, "\n")
  cat("Condition :", condition, "\n")

  cc <- readRDS(f)

  ##########################################################
  # MATRICES
  ##########################################################

  count_mat  <- cc@net$count
  weight_mat <- cc@net$weight

  if (is.null(count_mat) || is.null(weight_mat)) {

    warning(
      donor,
      " does not contain aggregated network matrices."
    )

    next
  }

  ##########################################################
  # BASIC NETWORK CHARACTERISTICS
  ##########################################################

  n_celltypes <- nrow(weight_mat)

  total_interactions <- sum(
    count_mat,
    na.rm = TRUE
  )

  total_strength <- sum(
    weight_mat,
    na.rm = TRUE
  )

  active_edges <- sum(
    weight_mat > 0,
    na.rm = TRUE
  )

  ##########################################################
  # DENSITY
  #
  # Directed network including autocrine interactions.
  ##########################################################

  possible_edges <- n_celltypes^2

  density <- active_edges / possible_edges

  ##########################################################
  # MEAN STRENGTH AMONG ACTIVE EDGES
  ##########################################################

  mean_edge_strength <- ifelse(
    active_edges > 0,
    total_strength / active_edges,
    NA
  )

  ##########################################################
  # LR INTERACTIONS PER ACTIVE CELL-CELL EDGE
  ##########################################################

  mean_LR_per_edge <- ifelse(
    active_edges > 0,
    total_interactions / active_edges,
    NA
  )

  ##########################################################
  # PATHWAYS
  ##########################################################

  n_pathways <- if (!is.null(cc@netP$pathways)) {
    length(cc@netP$pathways)
  } else {
    NA_integer_
  }

  ##########################################################
  # NUMBER OF CELLS
  ##########################################################

  n_cells <- length(cc@idents)

  ##########################################################
  # MASTER NETWORK TABLE
  ##########################################################

  network_results[[donor]] <- data.frame(

    Donor = donor,
    Condition = condition,

    n_cells = n_cells,
    n_celltypes = n_celltypes,

    total_LR_interactions = total_interactions,

    active_edges = active_edges,
    possible_edges = possible_edges,

    density = density,

    total_strength = total_strength,
    mean_edge_strength = mean_edge_strength,

    mean_LR_per_edge = mean_LR_per_edge,

    n_pathways = n_pathways,

    stringsAsFactors = FALSE
  )

  ##########################################################
  # CELL-TYPE LEVEL METRICS
  #
  # rows    = sender
  # columns = receiver
  ##########################################################

  celltypes <- rownames(weight_mat)

  outgoing_strength <- rowSums(
    weight_mat,
    na.rm = TRUE
  )

  incoming_strength <- colSums(
    weight_mat,
    na.rm = TRUE
  )

  outgoing_edges <- rowSums(
    weight_mat > 0,
    na.rm = TRUE
  )

  incoming_edges <- colSums(
    weight_mat > 0,
    na.rm = TRUE
  )

  outgoing_LR <- rowSums(
    count_mat,
    na.rm = TRUE
  )

  incoming_LR <- colSums(
    count_mat,
    na.rm = TRUE
  )

  celltype_results[[donor]] <- data.frame(

    Donor = donor,
    Condition = condition,
    CellType = celltypes,

    outgoing_strength = outgoing_strength,
    incoming_strength = incoming_strength,

    outgoing_edges = outgoing_edges,
    incoming_edges = incoming_edges,

    outgoing_LR = outgoing_LR,
    incoming_LR = incoming_LR,

    stringsAsFactors = FALSE
  )

  ##########################################################
  # PATHWAY-LEVEL METRICS
  ##########################################################

  if (
    !is.null(cc@netP$prob) &&
    !is.null(cc@netP$pathways)
  ) {

    prob_array <- cc@netP$prob
    pathways <- cc@netP$pathways

    for (k in seq_along(pathways)) {

      pmat <- prob_array[, , k]

      total_probability <- sum(
        pmat,
        na.rm = TRUE
      )

      pathway_edges <- sum(
        pmat > 0,
        na.rm = TRUE
      )

      mean_probability <- ifelse(
        pathway_edges > 0,
        total_probability / pathway_edges,
        NA
      )

      pathway_results[[paste(donor, k, sep = "_")]] <-
        data.frame(

          Donor = donor,
          Condition = condition,
          Pathway = pathways[k],

          total_probability = total_probability,
          active_edges = pathway_edges,
          mean_probability = mean_probability,

          stringsAsFactors = FALSE
        )
    }
  }

  cat(
    "  Cells:", n_cells,
    "| Cell types:", n_celltypes,
    "| LR:", total_interactions,
    "| Edges:", active_edges,
    "| Density:", round(density, 3),
    "| Strength:", round(total_strength, 3),
    "| Pathways:", n_pathways,
    "\n\n"
  )
}

############################################################
# COMBINE RESULTS
############################################################

network_df <- do.call(
  rbind,
  network_results
)

celltype_df <- do.call(
  rbind,
  celltype_results
)

pathway_df <- do.call(
  rbind,
  pathway_results
)

rownames(network_df) <- NULL
rownames(celltype_df) <- NULL
rownames(pathway_df) <- NULL

############################################################
# ORDER DONORS
############################################################

network_df <- network_df[
  order(network_df$Donor),
]

celltype_df <- celltype_df[
  order(
    celltype_df$Donor,
    celltype_df$CellType
  ),
]

pathway_df <- pathway_df[
  order(
    pathway_df$Donor,
    pathway_df$Pathway
  ),
]

############################################################
# SAVE
############################################################

write.csv(
  network_df,
  file.path(
    outdir,
    "donor_network_metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  celltype_df,
  file.path(
    outdir,
    "donor_celltype_metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  pathway_df,
  file.path(
    outdir,
    "donor_pathway_metrics.csv"
  ),
  row.names = FALSE
)

############################################################
# PRINT MASTER TABLE
############################################################

cat("\n============================================\n")
cat("MASTER DONOR NETWORK TABLE\n")
cat("============================================\n\n")

print(
  network_df,
  row.names = FALSE
)

cat("\n============================================\n")
cat("FILES CREATED\n")
cat("============================================\n")

cat(
  file.path(
    outdir,
    "donor_network_metrics.csv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "donor_celltype_metrics.csv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "donor_pathway_metrics.csv"
  ),
  "\n"
)

cat("\nDONE\n")
