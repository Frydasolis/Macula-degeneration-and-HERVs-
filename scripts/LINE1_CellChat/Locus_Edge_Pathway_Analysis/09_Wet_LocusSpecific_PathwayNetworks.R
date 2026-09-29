############################################################
# 09_Wet_LocusSpecific_PathwayNetworks.R
#
# WET AMD
# L1FLnI locus -> CellChat edge -> top pathway
#
# Primary locus-edge association:
#   donor-level Pearson
#   mean locus expression vs CellChat edge weight
#   FDR < 0.05
#
# Pathway annotation:
#   top associated pathway from hierarchical pathway analysis
#
# Edge width = |Pearson r| for locus-edge association
# Edge label = top associated CellChat pathway
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
})

############################################################
# PATHS
############################################################

root <- "/storage/lemus_g/roldan/ARMD"

analysis_dir <- file.path(
  root,
  "results/SRP413248/CellChat/donor_level_symbol",
  "donor_figure_analysis/locus_edge_correlations"
)

pathway_dir <- file.path(
  analysis_dir,
  "paper_figures/Locus_Edge_Pathway"
)

out_dir <- file.path(
  analysis_dir,
  "paper_figures/Wet_LocusSpecific_PathwayNetworks"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# INPUT
############################################################

stats <- read.csv(
  file.path(
    analysis_dir,
    "Locus_Edge_Correlations_PEARSON_PRIMARY.csv"
  ),
  stringsAsFactors = FALSE
)

top_pathways <- read.csv(
  file.path(
    pathway_dir,
    "TOP_Pathway_Per_LocusEdge.csv"
  ),
  stringsAsFactors = FALSE
)

############################################################
# PRIMARY WET LOCUS-EDGE ASSOCIATIONS
############################################################

wet <- stats %>%

  filter(
    condition == "Wet AMD",
    predictor == "mean_expression",
    method == "pearson",
    FDR < 0.05
  ) %>%

  mutate(
    locus_clean = sub(
      "^TE-L1FLnI-",
      "",
      locus
    )
  )

############################################################
# PREPARE PATHWAY ANNOTATIONS
############################################################

pathway_annotation <- top_pathways %>%

  filter(
    condition == "Wet AMD",
    method == "pearson"
  ) %>%

  select(
    locus,
    sender,
    receiver,
    pathway,
    pathway_r = r,
    pathway_p = p,
    pathway_FDR
  )

############################################################
# JOIN
############################################################

wet <- wet %>%

  left_join(
    pathway_annotation,
    by = c(
      "locus",
      "sender",
      "receiver"
    )
  )

cat("\n============================================\n")
cat("WET LOCUS + PATHWAY NETWORK DATA\n")
cat("============================================\n\n")

cat("Locus-edge associations:", nrow(wet), "\n")
cat("With pathway annotation:", sum(!is.na(wet$pathway)), "\n\n")

print(
  as.data.frame(
    wet %>%
      count(
        locus_clean,
        name = "n_edges"
      )
  )
)

############################################################
# SAVE EXACT DATA USED IN FIGURE
############################################################

write.csv(
  wet,
  file.path(
    out_dir,
    "Wet_Locus_Pathway_Network_DATA.csv"
  ),
  row.names = FALSE
)

############################################################
# FIXED NODE POSITIONS
############################################################

nodes <- data.frame(

  celltype = c(
    "RPE",
    "Astrocyte",
    "Muller glia",
    "RGC",
    "Horizontal",
    "Amacrine",
    "Bipolar",
    "Photoreceptor",
    "Microglia"
  ),

  x = c(
     0.00,
    -1.75,
    -0.70,
     0.75,
     1.75,
     1.55,
     0.60,
    -0.80,
    -1.70
  ),

  y = c(
     0.00,
     1.05,
     1.70,
     1.70,
     0.75,
    -0.90,
    -1.70,
    -1.70,
    -0.75
  ),

  stringsAsFactors = FALSE
)

nodes <- nodes %>%

  mutate(
    node_class = ifelse(
      celltype == "RPE",
      "RPE",
      "Other"
    )
  )

############################################################
# COLORS
############################################################

node_colors <- c(
  "RPE"   = "#F2B66D",
  "Other" = "#F1F4F6"
)

locus_colors <- c(
  "3q13.13ka" = "#A84A5B",
  "3q13.13ja" = "#D45D79",
  "3q13.13ea" = "#E58C62",
  "2p24.3g"   = "#7467A8"
)

############################################################
# ADD COORDINATES
############################################################

prepare_edges <- function(df) {

  df %>%

    left_join(
      nodes %>%
        select(
          sender = celltype,
          x_sender = x,
          y_sender = y
        ),
      by = "sender"
    ) %>%

    left_join(
      nodes %>%
        select(
          receiver = celltype,
          x_receiver = x,
          y_receiver = y
        ),
      by = "receiver"
    )
}

############################################################
# NETWORK FUNCTION
############################################################

make_network <- function(locus_name) {

  df <- wet %>%
    filter(
      locus_clean == locus_name
    ) %>%
    prepare_edges()

  col <- unname(
    locus_colors[locus_name]
  )

  normal_edges <- df %>%
    filter(
      sender != receiver
    )

  self_edges <- df %>%
    filter(
      sender == receiver
    )

  ##########################################################
  # OFFSET PATHWAY LABELS
  #
  # Move labels slightly perpendicular to each edge.
  ##########################################################

  normal_edges <- normal_edges %>%

    mutate(

      mid_x = (
        x_sender +
        x_receiver
      ) / 2,

      mid_y = (
        y_sender +
        y_receiver
      ) / 2,

      dx = x_receiver - x_sender,
      dy = y_receiver - y_sender,

      edge_length = sqrt(
        dx^2 + dy^2
      ),

      label_x =
        mid_x -
        0.13 * dy / edge_length,

      label_y =
        mid_y +
        0.13 * dx / edge_length
    )

  ##########################################################
  # PLOT
  ##########################################################

  p <- ggplot() +

    ########################################################
    # EDGES
    ########################################################

    geom_curve(
      data = normal_edges,

      aes(
        x = x_sender,
        y = y_sender,
        xend = x_receiver,
        yend = y_receiver,
        linewidth = abs(r)
      ),

      curvature = 0.13,

      color = col,

      alpha = 0.78,

      lineend = "round",

      arrow = arrow(
        type = "closed",
        length = unit(
          0.115,
          "inches"
        )
      )
    ) +

    ########################################################
    # PATHWAY LABELS
    ########################################################

    geom_label(
      data = normal_edges,

      aes(
        x = label_x,
        y = label_y,
        label = pathway
      ),

      size = 3.15,

      fontface = "bold",

      color = "#263238",

      fill = "white",

      label.size = 0.15,

      label.padding = unit(
        0.10,
        "lines"
      ),

      alpha = 0.94
    ) +

    ########################################################
    # NODES
    ########################################################

    geom_point(
      data = nodes,

      aes(
        x = x,
        y = y,
        fill = node_class
      ),

      shape = 21,

      size = 12,

      stroke = 1.15,

      color = "#45525C"
    ) +

    ########################################################
    # NODE LABELS
    ########################################################

    geom_text(
      data = nodes,

      aes(
        x = x,
        y = y,
        label = celltype
      ),

      size = 3.4,

      fontface = "bold",

      color = "#263238"
    ) +

    ########################################################
    # SCALES
    ########################################################

    scale_fill_manual(
      values = node_colors,
      guide = "none"
    ) +

    scale_linewidth_continuous(
      range = c(
        1.5,
        5.2
      ),

      limits = c(
        0.90,
        1.00
      ),

      breaks = c(
        0.90,
        0.95,
        1.00
      ),

      name = "|Pearson r|"
    ) +

    ########################################################
    # TITLES
    ########################################################

    labs(

      title = paste0(
        "L1FLnI-",
        locus_name
      ),

      subtitle = paste0(
        nrow(df),
        " FDR-significant communication edges"
      ),

      caption =
        "Labels indicate the top associated CellChat pathway"
    ) +

    coord_equal(
      xlim = c(
        -2.40,
        2.40
      ),

      ylim = c(
        -2.30,
        2.30
      ),

      clip = "off"
    ) +

    theme_void(
      base_size = 12
    ) +

    theme(

      plot.title = element_text(
        size = 17,
        face = "bold",
        hjust = 0.5,
        color = col
      ),

      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5,
        color = "#56636B"
      ),

      plot.caption = element_text(
        size = 8.5,
        hjust = 0.5,
        color = "#68757D",
        margin = margin(
          t = 8
        )
      ),

      legend.position = "bottom",

      legend.title = element_text(
        size = 9,
        face = "bold"
      ),

      legend.text = element_text(
        size = 8
      ),

      plot.margin = margin(
        12,
        18,
        12,
        18
      )
    )

  ##########################################################
  # SELF-LOOPS
  ##########################################################

  if (nrow(self_edges) > 0) {

    for (j in seq_len(nrow(self_edges))) {

      loop <- self_edges[j, ]

      theta <- seq(
        0.15 * pi,
        1.85 * pi,
        length.out = 120
      )

      radius <- 0.31

      loop_df <- data.frame(

        x =
          loop$x_sender +
          radius * cos(theta),

        y =
          loop$y_sender +
          0.33 +
          radius * sin(theta)
      )

      loop_width <- 1.5 +
        3.7 *
        (
          abs(loop$r) - 0.90
        ) / 0.10

      loop_width <- max(
        1.5,
        min(
          5.2,
          loop_width
        )
      )

      p <- p +

        geom_path(
          data = loop_df,

          aes(
            x = x,
            y = y
          ),

          inherit.aes = FALSE,

          color = col,

          linewidth = loop_width,

          alpha = 0.78,

          lineend = "round",

          arrow = arrow(
            type = "closed",
            length = unit(
              0.11,
              "inches"
            )
          )
        ) +

        annotate(
          "label",

          x = loop$x_sender,

          y = loop$y_sender + 0.72,

          label = loop$pathway,

          size = 3.15,

          fontface = "bold",

          color = "#263238",

          fill = "white"
        )
    }
  }

  return(p)
}

############################################################
# LOCUS ORDER
############################################################

loci_order <- c(
  "3q13.13ka",
  "3q13.13ja",
  "3q13.13ea",
  "2p24.3g"
)

############################################################
# CREATE INDIVIDUAL FIGURES
############################################################

plots <- list()

for (loc in loci_order) {

  if (!loc %in% wet$locus_clean) {
    next
  }

  cat(
    "Creating network:",
    loc,
    "\n"
  )

  p <- make_network(
    loc
  )

  plots[[loc]] <- p

  safe <- gsub(
    "[^A-Za-z0-9._-]",
    "_",
    loc
  )

  ggsave(
    file.path(
      out_dir,
      paste0(
        "Wet_",
        safe,
        "_PathwayNetwork.png"
      )
    ),
    p,
    width = 7,
    height = 7,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      paste0(
        "Wet_",
        safe,
        "_PathwayNetwork.pdf"
      )
    ),
    p,
    width = 7,
    height = 7,
    bg = "white"
  )
}

############################################################
# COMBINED 2 x 2
############################################################

if (
  requireNamespace(
    "patchwork",
    quietly = TRUE
  )
) {

  suppressPackageStartupMessages(
    library(patchwork)
  )

  combined <- (
    plots[["3q13.13ka"]] |
    plots[["3q13.13ja"]]
  ) /
  (
    plots[["3q13.13ea"]] |
    plots[["2p24.3g"]]
  ) +

    plot_annotation(

      title =
        "Locus-specific L1FLnI associations converge on RPE-centered signaling in wet AMD",

      subtitle =
        "Donor-level locus–CellChat associations surviving BH-FDR < 0.05",

      theme = theme(

        plot.title = element_text(
          size = 20,
          face = "bold",
          hjust = 0.5,
          color = "#263238"
        ),

        plot.subtitle = element_text(
          size = 11,
          hjust = 0.5,
          color = "#5F6B73"
        )
      )
    )

  ggsave(
    file.path(
      out_dir,
      "Wet_Four_Loci_PathwayNetworks.png"
    ),
    combined,
    width = 14,
    height = 14,
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(
      out_dir,
      "Wet_Four_Loci_PathwayNetworks.pdf"
    ),
    combined,
    width = 14,
    height = 14,
    bg = "white"
  )
}

############################################################
# PRINT EXACT PATHWAY TABLE
############################################################

figure_table <- wet %>%

  select(
    locus_clean,
    sender,
    receiver,
    r,
    FDR,
    pathway,
    pathway_r,
    pathway_FDR
  ) %>%

  arrange(
    locus_clean,
    sender,
    receiver
  )

write.csv(
  figure_table,
  file.path(
    out_dir,
    "Pathways_Shown_In_Figure.csv"
  ),
  row.names = FALSE
)

cat("\n============================================\n")
cat("PATHWAY-ANNOTATED NETWORKS COMPLETE\n")
cat("============================================\n\n")

cat(
  "Networks:",
  length(plots),
  "\n"
)

cat(
  "Total locus-edge associations:",
  nrow(wet),
  "\n"
)

cat(
  "Pathway annotations:",
  sum(!is.na(wet$pathway)),
  "\n\n"
)

print(
  as.data.frame(
    figure_table
  )
)

cat("\nOutput:\n")
cat(out_dir, "\n")

