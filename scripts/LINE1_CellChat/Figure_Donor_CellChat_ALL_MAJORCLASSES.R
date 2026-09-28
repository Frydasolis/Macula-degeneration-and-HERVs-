############################################################
# Figure_Donor_CellChat_ALL_MAJORCLASSES.R
#
# Donor-level CellChat networks
#
# 9 HRCA major classes
# ALL signaling pathways
# GLOBAL probability scale across donors
#
# Produces:
#   1. FULL interactome
#   2. HIGHLIGHTED interactome
#
# Highlighted pathways:
# APP, CADM, NCAM, NRXN, PSAP, VEGF
############################################################

suppressPackageStartupMessages({
  library(CellChat)
  library(igraph)
})

############################################################
# PATHS
############################################################

base <- paste0(
  "/storage/lemus_g/roldan/ARMD/results/",
  "SRP413248/CellChat/donor_level_symbol"
)

outdir <- file.path(
  base,
  "network_metrics/donor_networks_ALL_majorclasses"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

############################################################
# DONORS
############################################################

condition_map <- c(

  GSM6841143="Healthy",
  GSM6841144="Healthy",
  GSM6841145="Healthy",
  GSM6841146="Healthy",
  GSM6841147="Healthy",
  GSM6841148="Healthy",

  GSM6841149="Wet AMD",
  GSM6841150="Wet AMD",
  GSM6841151="Wet AMD",
  GSM6841152="Wet AMD",
  GSM6841153="Wet AMD",
  GSM6841154="Wet AMD",
  GSM6841155="Wet AMD",

  GSM6841156="Dry AMD",
  GSM6841157="Dry AMD",
  GSM6841159="Dry AMD"
)

donor_order <- names(condition_map)

############################################################
# ALL 9 MAJOR CLASSES
############################################################

cell_order <- c(
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
# FIXED MANUSCRIPT-STYLE COORDINATES
############################################################

coords <- rbind(

  "Bipolar" =
    c( 0.00,  1.05),

  "Astrocyte" =
    c( 0.72,  0.76),

  "Amacrine" =
    c( 1.08,  0.25),

  "Microglia" =
    c( 1.08, -0.35),

  "RPE" =
    c( 0.62, -0.82),

  "RGC" =
    c( 0.00, -1.05),

  "Photoreceptor" =
    c(-0.72, -0.78),

  "Muller glia" =
    c(-1.08, -0.10),

  "Horizontal" =
    c(-0.72,  0.72)
)

############################################################
# PATHWAYS TO HIGHLIGHT
############################################################

highlight_pathways <- c(
  "APP",
  "CADM",
  "NCAM",
  "NRXN",
  "PSAP",
  "VEGF"
)

highlight_colors <- c(
  APP  = "#D73027",
  CADM = "#2A9D8F",
  NCAM = "#7B2CBF",
  NRXN = "#E76F51",
  PSAP = "#F26B21",
  VEGF = "#E69F00"
)

############################################################
# LOCATE DONOR OBJECTS
############################################################

all_files <- list.files(
  base,
  pattern="\\.rds$",
  full.names=TRUE
)

donor_files <- setNames(
  rep(NA_character_, length(donor_order)),
  donor_order
)

for(donor in donor_order){

  hits <- all_files[
    grepl(donor, basename(all_files))
  ]

  if(length(hits) == 0){

    warning("No file for ", donor)
    next
  }

  if(length(hits) > 1){

    cat(
      "\nMultiple files found for ",
      donor,
      ":\n",
      paste(hits, collapse="\n"),
      "\n",
      sep=""
    )

    stop("Ambiguous donor file.")
  }

  donor_files[donor] <- hits
}

############################################################
# EXTRACT ALL PATHWAY EDGES
############################################################

extract_all_edges <- function(cc, donor){

  pathways <- cc@netP$pathways

  if(length(pathways) == 0)
    return(NULL)

  out <- list()
  k <- 1

  for(pathway in pathways){

    idx <- which(
      cc@netP$pathways == pathway
    )

    if(length(idx) == 0)
      next

    p <- cc@netP$prob[,,idx,drop=FALSE]

    if(length(dim(p)) == 3)
      p <- p[,,1]

    rn <- rownames(p)
    cn <- colnames(p)

    if(is.null(rn))
      rn <- dimnames(p)[[1]]

    if(is.null(cn))
      cn <- dimnames(p)[[2]]

    for(i in seq_len(nrow(p))){

      for(j in seq_len(ncol(p))){

        value <- p[i,j]

        if(
          is.finite(value) &&
          value > 0
        ){

          out[[k]] <- data.frame(

            Donor=donor,

            Condition=
              unname(condition_map[donor]),

            source=rn[i],

            target=cn[j],

            pathway=pathway,

            probability=
              as.numeric(value),

            stringsAsFactors=FALSE
          )

          k <- k + 1
        }
      }
    }
  }

  if(length(out) == 0)
    return(NULL)

  do.call(rbind,out)
}

############################################################
# EXTRACT FROM ALL DONORS
############################################################

edge_list <- list()

for(donor in donor_order){

  f <- donor_files[donor]

  if(is.na(f))
    next

  cat(
    "\nExtracting ",
    donor,
    " — ",
    condition_map[[donor]],
    "\n",
    sep=""
  )

  cc <- readRDS(f)

  cat(
    "Cell types in object: ",
    paste(
      dimnames(cc@net$weight)[[1]],
      collapse=", "
    ),
    "\n",
    sep=""
  )

  cat(
    "Number of pathways: ",
    length(cc@netP$pathways),
    "\n",
    sep=""
  )

  cat(
    "APP present: ",
    "APP" %in% cc@netP$pathways,
    "\n",
    sep=""
  )

  edge_list[[donor]] <-
    extract_all_edges(
      cc,
      donor
    )
}

edge_list <- edge_list[
  !vapply(
    edge_list,
    is.null,
    logical(1)
  )
]

edges <- do.call(
  rbind,
  edge_list
)

############################################################
# STANDARDIZE POSSIBLE MULLER SPELLINGS
############################################################

edges$source[
  edges$source %in%
  c(
    "Müller glia",
    "Muller_glia",
    "Muller Glia"
  )
] <- "Muller glia"

edges$target[
  edges$target %in%
  c(
    "Müller glia",
    "Muller_glia",
    "Muller Glia"
  )
] <- "Muller glia"

############################################################
# REPORT ANY UNEXPECTED CELL TYPES
############################################################

observed_cells <- sort(
  unique(
    c(
      edges$source,
      edges$target
    )
  )
)

unexpected <- setdiff(
  observed_cells,
  cell_order
)

cat("\n========================================\n")
cat("CELL TYPE CHECK\n")
cat("========================================\n")

cat(
  "Observed:\n",
  paste(observed_cells, collapse=", "),
  "\n"
)

if(length(unexpected) > 0){

  cat(
    "\nWARNING — unexpected cell types:\n",
    paste(unexpected, collapse=", "),
    "\n"
  )
}

############################################################
# KEEP THE 9 MAJOR CLASSES
############################################################

edges <- edges[
  edges$source %in% cell_order &
  edges$target %in% cell_order,
  ,
  drop=FALSE
]

############################################################
# SAVE COMPLETE EDGE TABLE
############################################################

write.csv(
  edges,
  file.path(
    outdir,
    "ALL_DONORS_ALL_PATHWAYS_EDGES.csv"
  ),
  row.names=FALSE
)

############################################################
# PATHWAY SUMMARY
############################################################

pathway_summary <- aggregate(
  probability ~
    Condition +
    Donor +
    pathway,
  data=edges,
  FUN=sum
)

write.csv(
  pathway_summary,
  file.path(
    outdir,
    "DONOR_PATHWAY_SUMMARY.csv"
  ),
  row.names=FALSE
)

############################################################
# MICROGLIA-SPECIFIC TABLE
############################################################

micro_edges <- edges[
  edges$source == "Microglia" |
  edges$target == "Microglia",
  ,
  drop=FALSE
]

write.csv(
  micro_edges,
  file.path(
    outdir,
    "MICROGLIA_ALL_DONOR_EDGES.csv"
  ),
  row.names=FALSE
)

############################################################
# MICROGLIA + APP TABLE
############################################################

micro_app <- micro_edges[
  micro_edges$pathway == "APP",
  ,
  drop=FALSE
]

write.csv(
  micro_app,
  file.path(
    outdir,
    "MICROGLIA_APP_EDGES.csv"
  ),
  row.names=FALSE
)

############################################################
# GLOBAL PROBABILITY SCALE
############################################################

q <- quantile(
  edges$probability,
  probs=c(
    0,
    .25,
    .50,
    .75,
    .90,
    .95,
    .99,
    1
  ),
  na.rm=TRUE
)

print(q)

visual_max <- as.numeric(
  q["99%"]
)

cat(
  "\nP99 visual ceiling:",
  visual_max,
  "\n"
)

############################################################
# WIDTH FUNCTIONS
############################################################

prob_width <- function(p){

  p2 <- pmin(
    p,
    visual_max
  )

  0.6 +
    7.5 *
    p2 /
    visual_max
}

prob_arrow <- function(p){

  p2 <- pmin(
    p,
    visual_max
  )

  0.20 +
    0.40 *
    p2 /
    visual_max
}

############################################################
# DRAW FUNCTION
############################################################

draw_network <- function(
  donor,
  mode=c(
    "FULL",
    "HIGHLIGHT"
  ),
  device=c(
    "png",
    "pdf"
  )
){

  mode <- match.arg(mode)
  device <- match.arg(device)

  d <- edges[
    edges$Donor == donor,
    ,
    drop=FALSE
  ]

  condition <- condition_map[[donor]]

  if(nrow(d) == 0)
    return(NULL)

  ##########################################################
  # EDGE TABLE FOR IGRAPH
  ##########################################################

  edge_df <- data.frame(
    from=d$source,
    to=d$target,
    pathway=d$pathway,
    probability=d$probability,
    stringsAsFactors=FALSE
  )

  vertices <- data.frame(
    name=cell_order,
    stringsAsFactors=FALSE
  )

  g <- graph_from_data_frame(
    edge_df,
    directed=TRUE,
    vertices=vertices
  )

  lay <- coords[
    V(g)$name,
    ,
    drop=FALSE
  ]

  ##########################################################
  # COLORS
  ##########################################################

  if(mode == "FULL"){

    # All pathways represented neutrally.
    # This prevents a huge unreadable color legend.
    edge_col <- rep(
      "#8A8A8A55",
      nrow(edge_df)
    )

  } else {

    edge_col <- rep(
      "#BDBDBD45",
      nrow(edge_df)
    )

    selected <-
      edge_df$pathway %in%
      highlight_pathways

    edge_col[selected] <-
      unname(
        highlight_colors[
          edge_df$pathway[selected]
        ]
      )
  }

  E(g)$color <- edge_col

  E(g)$width <-
    prob_width(
      edge_df$probability
    )

  E(g)$arrow.size <-
    prob_arrow(
      edge_df$probability
    )

  ##########################################################
  # OUTPUT
  ##########################################################

  mode_dir <- file.path(
    outdir,
    mode,
    gsub(" ","_",condition)
  )

  dir.create(
    mode_dir,
    recursive=TRUE,
    showWarnings=FALSE
  )

  stem <- file.path(
    mode_dir,
    paste0(
      donor,
      "_",
      mode,
      "_all_majorclasses"
    )
  )

  if(device=="png"){

    png(
      paste0(stem,".png"),
      width=4800,
      height=3200,
      res=400
    )

  } else {

    pdf(
      paste0(stem,".pdf"),
      width=12,
      height=8
    )
  }

  par(
    mar=c(1,1,2,9),
    xpd=NA,
    family="sans"
  )

  ##########################################################
  # DRAW GRAPH
  ##########################################################

  plot(
    g,

    layout=lay,

    rescale=FALSE,

    xlim=c(-1.55,2.25),
    ylim=c(-1.35,1.35),

    asp=1,

    vertex.size=29,
    vertex.color="white",
    vertex.frame.color="black",
    vertex.frame.width=4,

    vertex.label=NA,

    edge.color=E(g)$color,
    edge.width=E(g)$width,
    edge.arrow.size=E(g)$arrow.size,

    edge.curved=0.04,

    main=""
  )

  ##########################################################
  # MANUAL LABELS
  ##########################################################

  text(
    0.00,1.27,
    "Bipolar",
    font=2,cex=1.20
  )

  text(
    0.92,0.84,
    "Astrocyte",
    font=2,cex=1.20
  )

  text(
    1.28,0.28,
    "Amacrine",
    font=2,cex=1.20
  )

  text(
    1.31,-0.39,
    "Microglia",
    font=2,cex=1.20
  )

  text(
    0.77,-0.94,
    "RPE",
    font=2,cex=1.20
  )

  text(
    0.00,-1.27,
    "RGC",
    font=2,cex=1.20
  )

  text(
    -0.88,-0.91,
    "Photoreceptor",
    font=2,cex=1.20
  )

  text(
    -1.30,-0.10,
    "Muller glia",
    font=2,cex=1.20
  )

  text(
    -0.91,0.82,
    "Horizontal",
    font=2,cex=1.20
  )

  ##########################################################
  # DONOR
  ##########################################################

  text(
    -1.48,
    1.30,

    paste0(
      donor,
      " — ",
      condition
    ),

    adj=c(0,0.5),
    font=2,
    cex=1.0
  )

  ##########################################################
  # HIGHLIGHT LEGEND
  ##########################################################

  if(mode=="HIGHLIGHT"){

    text(
      1.48,
      1.28,

      "Signaling pathway",

      adj=c(0,0.5),

      font=2,

      cex=1.40
    )

    yy <- seq(
      1.03,
      -0.02,
      length.out=
        length(highlight_pathways)
    )

    donor_pathways <-
      unique(d$pathway)

    for(i in seq_along(
      highlight_pathways
    )){

      pw <-
        highlight_pathways[i]

      col_i <- if(
        pw %in% donor_pathways
      ){
        highlight_colors[pw]
      } else {
        "grey80"
      }

      arrows(
        1.49,
        yy[i],

        1.70,
        yy[i],

        length=0.07,

        lwd=3,

        col=col_i
      )

      text(
        1.79,
        yy[i],

        pw,

        adj=c(0,0.5),

        cex=1.10
      )
    }

  } else {

    text(
      1.48,
      1.25,

      "Full CellChat\ninteractome",

      adj=c(0,1),

      font=2,

      cex=1.35
    )

    text(
      1.48,
      0.72,

      paste0(
        length(unique(d$pathway)),
        " signaling pathways"
      ),

      adj=c(0,0.5),

      cex=1.05
    )
  }

  ##########################################################
  # PROBABILITY LEGEND
  ##########################################################

  legend_prob <- c(
    0.10,
    0.15,
    0.20,
    0.25
  )

  text(
    1.48,
    -0.40,

    "Communication\nprobability",

    adj=c(0,1),

    font=2,

    cex=1.30
  )

  yy2 <- c(
    -0.78,
    -0.96,
    -1.14,
    -1.32
  )

  for(i in seq_along(
    legend_prob
  )){

    p <- legend_prob[i]

    arrows(
      1.49,
      yy2[i],

      1.70,
      yy2[i],

      length=0.07,

      lwd=prob_width(p),

      col="black"
    )

    text(
      1.79,
      yy2[i],

      sprintf(
        "%.2f",
        p
      ),

      adj=c(0,0.5),

      cex=1.05
    )
  }

  dev.off()
}

############################################################
# GENERATE ALL NETWORKS
############################################################

for(donor in donor_order){

  if(
    !donor %in%
    unique(edges$Donor)
  )
    next

  cat(
    "\nDrawing ",
    donor,
    " — ",
    condition_map[[donor]],
    "\n",
    sep=""
  )

  for(mode in c(
    "FULL",
    "HIGHLIGHT"
  )){

    draw_network(
      donor,
      mode,
      "png"
    )

    draw_network(
      donor,
      mode,
      "pdf"
    )
  }
}

############################################################
# FINAL CHECK
############################################################

cat("\n========================================\n")
cat("ALL MAJOR-CLASS NETWORKS COMPLETE\n")
cat("========================================\n")

cat(
  "Donors:",
  length(unique(edges$Donor)),
  "\n"
)

cat(
  "Cell types:\n",
  paste(
    sort(
      unique(
        c(
          edges$source,
          edges$target
        )
      )
    ),
    collapse=", "
  ),
  "\n"
)

cat(
  "Total pathways:",
  length(unique(edges$pathway)),
  "\n"
)

cat(
  "Microglia edges:",
  nrow(micro_edges),
  "\n"
)

cat(
  "Microglia APP edges:",
  nrow(micro_app),
  "\n"
)

cat(
  "\nOutput:\n",
  outdir,
  "\n"
)
