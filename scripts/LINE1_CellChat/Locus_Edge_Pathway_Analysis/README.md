# L1FLnI locus–CellChat pathway analysis in wet AMD

## Overview

This analysis investigates whether expression of specific L1FLnI loci is associated with donor-level remodeling of retinal intercellular communication in wet age-related macular degeneration (AMD).

The analysis integrates donor-level LINE-1 expression with donor-level CellChat communication networks.

The overall analytical framework is:

**L1FLnI locus expression → cell-cell communication edge → signaling pathway**

The unit of statistical inference is the **donor**, not individual cells.

---

## Biological question

Are specific L1FLnI loci associated with remodeling of retinal cell-cell communication in wet AMD, and which CellChat signaling pathways characterize these locus-associated communication patterns?

---

## Primary locus-edge analysis

For each donor, global expression of individual candidate L1FLnI loci was calculated across retinal cells.

Locus expression was correlated with donor-level CellChat edge weights using Pearson correlation.

The primary analysis identified **24 FDR-significant locus-edge associations in wet AMD** involving four L1FLnI loci:

| Locus | FDR-significant edges |
|---|---:|
| L1FLnI-3q13.13ja | 7 |
| L1FLnI-3q13.13ea | 6 |
| L1FLnI-2p24.3g | 6 |
| L1FLnI-3q13.13ka | 5 |

All 24 associations were positive.

These associations showed strong convergence on communication edges involving the retinal pigment epithelium (RPE).

---

## Pathway-level decomposition

The 24 FDR-significant locus-edge associations were subsequently decomposed at the CellChat signaling-pathway level.

For each significant locus-edge pair, donor-level locus expression was correlated with donor-level CellChat pathway probability.

Pearson correlation was used as the primary analysis, with Spearman correlation as a sensitivity analysis.

Multiple-testing correction was performed using the Benjamini-Hochberg procedure within each previously selected locus-edge combination.

### Results

A total of:

- 1,440 Pearson locus-edge-pathway tests were evaluated
- 220 associations had nominal P < 0.05
- 218 associations had pathway-level FDR < 0.05

The significant pathway associations included signaling programs such as:

- APP
- NCAM
- CypA
- PTN
- FLRT
- CNTN
- ADGRL
- PTPR
- SLITRK
- COLLAGEN
- NRXN
- NOTCH
- RELN
- SLIT
- and others

Because multiple pathways were associated with individual cell-cell edges, the main network visualization displays the **top statistically supported pathway for each locus-edge association**.

All pathway-level associations are retained in the corresponding output tables.

---

## Wet AMD locus-specific networks

Four locus-specific networks were generated:

- L1FLnI-3q13.13ka
- L1FLnI-3q13.13ja
- L1FLnI-3q13.13ea
- L1FLnI-2p24.3g

In these figures:

- Nodes represent retinal cell classes.
- Directed arrows represent FDR-significant donor-level locus–CellChat edge associations.
- Arrow direction represents CellChat sender → receiver direction.
- Edge width represents the absolute Pearson correlation coefficient between locus expression and CellChat edge weight.
- Edge labels indicate the top associated CellChat signaling pathway.
- RPE is highlighted because significant associations showed strong convergence on RPE-centered communication.

---

## Important interpretation

These analyses identify **associations**, not causal relationships.

For example, an association:

`L1FLnI-2p24.3g → Müller glia → RPE → APP`

indicates that, among wet AMD donors, higher L1FLnI-2p24.3g expression is associated with stronger inferred APP pathway signaling from Müller glia to RPE.

It does not establish that L1FLnI-2p24.3g activates APP signaling.

Likewise, significance within wet AMD does not by itself demonstrate that an association is specific to wet AMD. Formal locus-expression × disease-condition interaction analyses are required to test disease-dependent differences.

---

## Scripts

### `10_Locus_Edge_Correlations.R`

Performs donor-level locus-specific correlations between L1FLnI expression and CellChat sender-receiver edge weights.

### `11_Figure_LocusSpecific_Pearson_Networks.R`

Generates initial locus-specific CellChat network visualizations.

### `08_Locus_Edge_Pathway_Correlations.R`

Performs hierarchical pathway-level decomposition of the FDR-significant wet AMD locus-edge associations.

### `09_Wet_LocusSpecific_PathwayNetworks.R`

Generates the final wet AMD locus-specific networks annotated with the top associated CellChat signaling pathway.

---

## Figures

The `figures/` directory contains:

- Combined four-locus pathway network
- Individual locus-specific pathway networks
- PNG versions for visualization
- PDF versions for publication-quality output

---

## Tables

The `tables/` directory contains:

### `Pathways_Shown_In_Figure.csv`

Exact pathway annotation displayed for each locus-edge in the final figure.

### `Wet_Locus_Pathway_Network_DATA.csv`

Complete data used to generate the pathway-annotated networks.

### `Locus_Edge_Pathway_PEARSON_FDR005.csv`

All pathway-level Pearson associations passing FDR < 0.05.

### `TOP_Pathway_Per_LocusEdge.csv`

Top statistically supported pathway for each significant locus-edge association.

### `Pathway_Summary_By_Locus.csv`

Summary of significant signaling pathways associated with each L1FLnI locus.

---

## Statistical considerations

- Unit of inference: donor
- Wet AMD donors: n = 7
- Primary correlation: Pearson
- Sensitivity analysis: Spearman
- Multiple testing: Benjamini-Hochberg FDR
- Cell-level observations are not treated as independent biological replicates
- Pathway analysis is a hierarchical follow-up of previously identified locus-edge associations

