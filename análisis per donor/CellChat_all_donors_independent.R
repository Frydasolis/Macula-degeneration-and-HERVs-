library(Seurat)
library(CellChat)

# ============================================================
# INPUT / OUTPUT
# ============================================================

input <- "/storage/lemus_g/roldan/ARMD/results/SRP413248/Seurat/downstream/SRP413248_GeneTE_QC_UMAP_HRCA_annotated_geneSymbols.rds"

outdir <- "/storage/lemus_g/roldan/ARMD/results/SRP413248/CellChat/donor_level_symbol"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ============================================================
# LOAD SEURAT
# ============================================================

cat("Loading Seurat object...\n")

seurat <- readRDS(input)

donors <- sort(unique(seurat$GSM))

cat("\nDonors detected:\n")
print(donors)

cat("\nNumber of donors:", length(donors), "\n")

# ============================================================
# CELLCHAT DATABASE
# ============================================================

data(CellChatDB.human)

# ============================================================
# LOOP
# ============================================================

for (donor_id in donors) {

    cat("\n\n")
    cat("############################################################\n")
    cat("PROCESSING:", donor_id, "\n")
    cat("############################################################\n")

    outfile <- file.path(
        outdir,
        paste0("CellChat_", donor_id, ".rds")
    )

    # Skip already completed
    if (file.exists(outfile)) {

        cat("Already exists. Skipping:", outfile, "\n")
        next
    }

    # ----------------------------------------------------------
    # SUBSET DONOR
    # ----------------------------------------------------------

    donor <- subset(
        seurat,
        subset = GSM == donor_id
    )

    cat("Cells:", ncol(donor), "\n")

    cat("Cell types:\n")
    print(table(donor$HRCA_majorclass))

    # ----------------------------------------------------------
    # CREATE CELLCHAT
    # ----------------------------------------------------------

    cellchat <- createCellChat(
        object = donor,
        group.by = "HRCA_majorclass"
    )

    cellchat@DB <- CellChatDB.human

    # ----------------------------------------------------------
    # SIGNALING GENES
    # ----------------------------------------------------------

    lr_genes <- unique(c(
        cellchat@DB$interaction$ligand,
        cellchat@DB$interaction$receptor
    ))

    lr_genes <- lr_genes[!is.na(lr_genes)]

    signaling_genes <- intersect(
        lr_genes,
        rownames(donor)
    )

    cat("Signaling genes:", length(signaling_genes), "\n")

    # ----------------------------------------------------------
    # EXPRESSION DATA
    # ----------------------------------------------------------

    expr <- GetAssayData(
        donor,
        assay = "RNA",
        layer = "data"
    )

    cellchat@data.signaling <- expr[
        signaling_genes,
        ,
        drop = FALSE
    ]

    cat("data.signaling:\n")
    print(dim(cellchat@data.signaling))

    # ----------------------------------------------------------
    # OVEREXPRESSED GENES
    # ----------------------------------------------------------

    cellchat <- identifyOverExpressedGenes(
        cellchat,
        thresh.pc = 0,
        thresh.fc = 0,
        thresh.p = 0.05,
        min.cells = 10
    )

    cat(
        "Overexpressed genes:",
        length(cellchat@var.features$features),
        "\n"
    )

    # ----------------------------------------------------------
    # OVEREXPRESSED INTERACTIONS
    # ----------------------------------------------------------

    cellchat <- identifyOverExpressedInteractions(cellchat)

    cat(
        "LR interactions:",
        nrow(cellchat@LR$LRsig),
        "\n"
    )

    # ----------------------------------------------------------
    # COMMUNICATION PROBABILITY
    # ----------------------------------------------------------

    set.seed(1)

    cellchat <- computeCommunProb(
        cellchat,
        type = "triMean",
        trim = 0.1,
        raw.use = TRUE,
        population.size = FALSE,
        nboot = 100,
        seed.use = 1
    )

    # ----------------------------------------------------------
    # FILTER
    # ----------------------------------------------------------

    cellchat <- filterCommunication(
        cellchat,
        min.cells = 10
    )

    # ----------------------------------------------------------
    # PATHWAYS
    # ----------------------------------------------------------

    cellchat <- computeCommunProbPathway(cellchat)

    cat(
        "Pathways:",
        length(cellchat@netP$pathways),
        "\n"
    )

    # ----------------------------------------------------------
    # AGGREGATE
    # ----------------------------------------------------------

    cellchat <- aggregateNet(cellchat)

    cat("Network dimensions:\n")
    print(dim(cellchat@net$prob))

    cat(
        "Significant edges:",
        sum(cellchat@net$pval < 0.05, na.rm = TRUE),
        "\n"
    )

    # ----------------------------------------------------------
    # SAVE
    # ----------------------------------------------------------

    saveRDS(
        cellchat,
        outfile,
        compress = TRUE
    )

    cat("\nSAVED:\n")
    cat(outfile, "\n")

    rm(cellchat)
    rm(donor)
    rm(expr)

    gc()

}

cat("\n\n")
cat("============================================================\n")
cat("ALL DONORS COMPLETED\n")
cat("============================================================\n")

print(
    list.files(
        outdir,
        pattern = "^CellChat_GSM.*\\.rds$",
        full.names = FALSE
    )
)
