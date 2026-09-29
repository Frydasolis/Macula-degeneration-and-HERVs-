############################################################
# 11_Figure_LocusSpecific_Pearson_Networks.R
#
# Locus-specific donor-level Pearson networks
#
# Each edge represents:
#   Pearson correlation between donor-level global locus
#   expression and donor-level CellChat communication weight
#
# Red  = positive association
# Blue = negative association
#
# Thick/opaque = FDR < 0.05
# Thin/transparent = nominal P < 0.05, FDR >= 0.05
#
# P >= 0.05 not displayed.
#
# Dry AMD n=3 = DESCRIPTIVE ONLY
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(scales)
})

############################################################
# PATHS
############################################################

root <- "/storage/lemus_g/roldan/ARMD"

analysis_dir <- file.path(
  root,
  "results/SRP413248/CellChat/donor_level_symbol/",
  "donor_figure_analysis/locus_edge_correlations"
)

candidate_file <- file.path(
  root,
  "results/SRP413248/Seurat/downstream/TE_analysis/Figure4",
  "LINE1_CellChat_candidates.csv"
)

out_dir <- file.path(
  analysis_dir,
  "figures_Pearson_networks"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# LOAD PEARSON RESULTS
############################################################

x <- read.csv(
  file.path(
    analysis_dir,
    "Locus_Edge_Correlations_PEARSON_PRIMARY.csv"
  ),
  stringsAsFactors = FALSE
)

cand <- read.csv(
  candidate_file,
  stringsAsFactors = FALSE
)

############################################################
# CELL TYPES
############################################################

celltypes <- c(
  "Bipolar",
  "Astrocyte",
  "Amacrine",
  "Microglia",
  "RPE",
  "RGC",
  "Photoreceptor",
  "Muller glia",
  "Horizontal"
)

############################################################
# FIXED CIRCULAR LAYOUT
############################################################

theta <- seq(
  0,
  2*pi,
  length.out = length(celltypes) + 1
)[-(length(celltypes) + 1)]

coords <- cbind(
  cos(theta),
  sin(theta)
)

rownames(coords) <- celltypes

############################################################
# COLORS
############################################################

positive_col <- "#D73027"
negative_col <- "#4575B4"
altered_col  <- "#FFD21F"
normal_col   <- "white"

############################################################
# LOCUS METADATA FROM FIGURE 4
############################################################

candidate_info <- cand %>%
  transmute(
    locus = TE_feature,
    celltype = HRCA_majorclass,
    comparison = comparison,
    direction = toupper(direction),
    log2FC = log2FC,
    FDR_DE = FDR
  )

############################################################
# FUNCTION: WHICH NODE SHOULD BE HIGHLIGHTED?
############################################################

get_node_info <- function(
  locus_name,
  condition_name
) {

  z <- candidate_info %>%
    filter(
      locus == locus_name
    )

  # For Healthy:
  # do not call the node "altered" because Healthy is reference.
  if (condition_name == "Healthy Control") {
    return(NULL)
  }

  # Wet:
  # accept Figure-4 comparisons where Wet is the TEST group.
  if (condition_name == "Wet AMD") {

    z <- z %>%
      filter(
        grepl(
          "^Wet AMD vs",
          comparison
        )
      )
  }

  # Dry:
  # accept Figure-4 comparisons where Dry is TEST group.
  if (condition_name == "Dry AMD") {

    z <- z %>%
      filter(
        grepl(
          "^Dry AMD vs",
          comparison
        )
      )
  }

  if (nrow(z) == 0) {
    return(NULL)
  }

  z
}

############################################################
# DRAW ONE LOCUS x CONDITION NETWORK
############################################################

draw_locus_network <- function(
  locus_name,
  condition_name
) {

  dat <- x %>%
    filter(
      locus == locus_name,
      condition == condition_name,
      !is.na(r),
      !is.na(p),
      p < 0.05
    )

  ##########################################################
  # VERTICES
  ##########################################################

  vertices <- data.frame(
    name = celltypes,
    stringsAsFactors = FALSE
  )

  node_info <- get_node_info(
    locus_name,
    condition_name
  )

  vertices$altered <- FALSE
  vertices$direction <- ""

  if (!is.null(node_info)) {

    for (j in seq_len(nrow(node_info))) {

      ct <- node_info$celltype[j]

      if (ct %in% vertices$name) {

        vertices$altered[
          vertices$name == ct
        ] <- TRUE

        # If multiple comparisons point to same cell type,
        # combine labels only if needed.
        old <- vertices$direction[
          vertices$name == ct
        ]

        new <- node_info$direction[j]

        if (old == "") {
          vertices$direction[
            vertices$name == ct
          ] <- new
        } else if (!grepl(new, old, fixed = TRUE)) {
          vertices$direction[
            vertices$name == ct
          ] <- paste(
            old,
            new,
            sep = "/"
          )
        }
      }
    }
  }

  ##########################################################
  # EMPTY NETWORK HANDLING
  ##########################################################

  if (nrow(dat) == 0) {

    plot(
      NA,
      xlim = c(-1.5, 1.5),
      ylim = c(-1.5, 1.5),
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = paste0(
        sub("^TE-", "", locus_name),
        "\n",
        condition_name
      )
    )

    text(
      0,
      0,
      "No Pearson edges\nwith P < 0.05",
      cex = 1.2
    )

    return(invisible(NULL))
  }

  ##########################################################
  # GRAPH
  ##########################################################

  g <- graph_from_data_frame(
    dat %>%
      select(
        sender,
        receiver,
        r,
        p,
        FDR,
        n
      ),
    directed = TRUE,
    vertices = vertices
  )

  xy <- coords[
    V(g)$name,
    ,
    drop = FALSE
  ]

  ##########################################################
  # EDGE COLORS
  ##########################################################

  raw_cols <- ifelse(
    E(g)$r > 0,
    positive_col,
    negative_col
  )

  sig_fdr <- E(g)$FDR < 0.05

  edge_cols <- ifelse(
    sig_fdr,
    adjustcolor(
      raw_cols,
      alpha.f = 0.90
    ),
    adjustcolor(
      raw_cols,
      alpha.f = 0.25
    )
  )

  ##########################################################
  # EDGE WIDTH
  ##########################################################

  # FDR edges are emphasized.
  # Nominal edges remain visible but thin.

  base_width <- rescale(
    abs(E(g)$r),
    to = c(1.0, 5.5)
  )

  edge_width <- ifelse(
    sig_fdr,
    base_width,
    pmax(
      0.6,
      base_width * 0.35
    )
  )

  ##########################################################
  # NODE COLORS / LABELS
  ##########################################################

  idx <- match(
    V(g)$name,
    vertices$name
  )

  node_cols <- ifelse(
    vertices$altered[idx],
    altered_col,
    normal_col
  )

  labels <- V(g)$name

  has_direction <-
    vertices$altered[idx] &
    vertices$direction[idx] != ""

  labels[has_direction] <- paste0(
    labels[has_direction],
    "\n",
    vertices$direction[idx][has_direction]
  )

  ##########################################################
  # TITLE
  ##########################################################

  n_donors <- max(
    dat$n,
    na.rm = TRUE
  )

  n_nominal <- sum(
    dat$p < 0.05,
    na.rm = TRUE
  )

  n_fdr <- sum(
    dat$FDR < 0.05,
    na.rm = TRUE
  )

  title1 <- paste0(
    sub("^TE-", "", locus_name),
    " — ",
    condition_name
  )

  title2 <- paste0(
    "Pearson | n=",
    n_donors,
    " donors | P<0.05: ",
    n_nominal,
    " | FDR<0.05: ",
    n_fdr
  )

  if (condition_name == "Dry AMD") {
    title2 <- paste0(
      title2,
      " | DESCRIPTIVE"
    )
  }

  ##########################################################
  # PLOT
  ##########################################################

  plot(
    g,

    layout = xy,

    vertex.color = node_cols,
    vertex.frame.color = "black",
    vertex.frame.width = 2,

    vertex.size = 27,

    vertex.label = labels,
    vertex.label.cex = 0.72,
    vertex.label.color = "black",
    vertex.label.dist = 1.30,

    edge.color = edge_cols,
    edge.width = edge_width,
    edge.arrow.size = 0.35,
    edge.curved = 0.14,

    main = paste0(
      title1,
      "\n",
      title2
    )
  )
}

############################################################
# GENERATE ALL 12 LOCI
#
# One figure = Healthy | Dry | Wet
############################################################

loci <- sort(
  unique(x$locus)
)

manifest <- list()

for (locus_name in loci) {

  safe <- gsub(
    "[^A-Za-z0-9._-]",
    "_",
    sub("^TE-", "", locus_name)
  )

  cat(
    "\nGenerating ",
    locus_name,
    "\n",
    sep = ""
  )

  ##########################################################
  # PNG
  ##########################################################

  png(
    file.path(
      out_dir,
      paste0(
        safe,
        "_Pearson_Network.png"
      )
    ),
    width = 4500,
    height = 1700,
    res = 300
  )

  par(
    mfrow = c(1,3),
    mar = c(2,2,5,2),
    oma = c(3,1,3,1)
  )

  draw_locus_network(
    locus_name,
    "Healthy Control"
  )

  draw_locus_network(
    locus_name,
    "Dry AMD"
  )

  draw_locus_network(
    locus_name,
    "Wet AMD"
  )

  mtext(
    "Locus-specific association with donor-level retinal communication",
    outer = TRUE,
    side = 3,
    line = 1,
    cex = 1.35,
    font = 2
  )

  mtext(
    "Red = positive Pearson r | Blue = negative Pearson r | Thick/opaque = FDR < 0.05 | Thin/transparent = nominal P < 0.05",
    outer = TRUE,
    side = 1,
    line = 1,
    cex = 0.78
  )

  dev.off()

  ##########################################################
  # PDF
  ##########################################################

  pdf(
    file.path(
      out_dir,
      paste0(
        safe,
        "_Pearson_Network.pdf"
      )
    ),
    width = 15,
    height = 5.7
  )

  par(
    mfrow = c(1,3),
    mar = c(2,2,5,2),
    oma = c(3,1,3,1)
  )

  draw_locus_network(
    locus_name,
    "Healthy Control"
  )

  draw_locus_network(
    locus_name,
    "Dry AMD"
  )

  draw_locus_network(
    locus_name,
    "Wet AMD"
  )

  mtext(
    "Locus-specific association with donor-level retinal communication",
    outer = TRUE,
    side = 3,
    line = 1,
    cex = 1.25,
    font = 2
  )

  mtext(
    "Red = positive Pearson r | Blue = negative Pearson r | Thick/opaque = FDR < 0.05 | Thin/transparent = nominal P < 0.05",
    outer = TRUE,
    side = 1,
    line = 1,
    cex = 0.75
  )

  dev.off()

  ##########################################################
  # MANIFEST
  ##########################################################

  tmp <- x %>%
    filter(
      locus == locus_name
    ) %>%
    group_by(condition) %>%
    summarise(
      n_valid_edges = sum(
        !is.na(r)
      ),
      n_nominal = sum(
        p < 0.05,
        na.rm = TRUE
      ),
      n_FDR = sum(
        FDR < 0.05,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      locus = locus_name
    )

  manifest[[locus_name]] <- tmp
}

############################################################
# SAVE MANIFEST
############################################################

manifest <- bind_rows(
  manifest
) %>%
  select(
    locus,
    condition,
    everything()
  )

write.csv(
  manifest,
  file.path(
    out_dir,
    "Pearson_Network_Figure_Manifest.csv"
  ),
  row.names = FALSE
)

############################################################
# SAVE EXACT EDGES DISPLAYED
############################################################

displayed <- x %>%
  filter(
    !is.na(r),
    !is.na(p),
    p < 0.05
  ) %>%
  mutate(
    edge_class = ifelse(
      FDR < 0.05,
      "FDR_significant",
      "Nominal_only"
    ),
    association = ifelse(
      r > 0,
      "Positive",
      "Negative"
    )
  )

write.csv(
  displayed,
  file.path(
    out_dir,
    "Pearson_Edges_DISPLAYED.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("PEARSON LOCUS-SPECIFIC NETWORKS COMPLETE\n")
cat("============================================\n")

cat("\nLoci:", length(loci), "\n")
cat("Expected:", length(loci), "PNG +", length(loci), "PDF\n")

cat("\nOutput:\n")
cat(out_dir, "\n")

cat("\nFDR edges by condition:\n")

print(
  x %>%
    filter(
      !is.na(FDR),
      FDR < 0.05
    ) %>%
    count(
      condition,
      locus,
      name = "n_FDR"
    ) %>%
    arrange(
      condition,
      desc(n_FDR)
    )
)

