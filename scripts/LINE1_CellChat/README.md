# Donor-level L1FLnI–CellChat analysis in AMD

This directory contains the scripts used to investigate the relationship between
LINE-1 (L1FLnI) expression and retinal intercellular communication in
age-related macular degeneration (AMD).

The analysis integrates donor-level L1FLnI expression derived from the retinal
snRNA-seq dataset with donor-level CellChat network metrics.

## Study design

Retinal samples were classified as:

- Healthy Control: n = 6 donors
- Dry AMD: n = 3 donors
- Wet AMD: n = 7 donors
- Total: n = 16 donors

The statistical unit throughout the analysis is the **donor**, not the
individual cell.

Cells are used to estimate donor-level expression and communication
measurements. Individual cells are therefore not treated as independent
biological replicates.

---

## Main CellChat outcomes

Two donor-level CellChat metrics were selected as the primary network outcomes.

### Network Density

Network density represents the fraction of possible directed cell-cell
connections that are active:

    Network Density = Active Edges / Possible Directed Edges

where:

    Possible Directed Edges = Number of Cell Types^2

Autocrine interactions are included.

Density is used as a measure of the overall architecture/connectivity of the
retinal communication network.

### Total Communication Strength

Total Communication Strength is calculated as:

    Total Communication Strength = sum(CellChat edge weights)

It represents the aggregate predicted communication strength across the retinal
cell-cell interaction network.

---

## L1FLnI measurements

L1FLnI activity is evaluated at several complementary levels.

### 1. Retina-wide L1FLnI burden

All annotated L1FLnI loci are used to calculate a donor-level retina-wide
expression burden.

The current dataset contains approximately 13,292 L1FLnI features.

For each donor, normalized expression is averaged across L1FLnI loci and cells.

This analysis is independent of the candidate loci identified by differential
expression.

### 2. Cell-type-specific L1FLnI burden

Family-wide L1FLnI expression is calculated independently within retinal cell
types.

Examples include:

- Photoreceptor
- RGC
- Bipolar
- Amacrine
- Astrocyte
- Müller glia
- RPE

Each observation remains a donor-level measurement.

### 3. Candidate L1FLnI burden

Candidate loci identified in the differential-expression analysis are
aggregated within their corresponding retinal cell type.

The Figure 4 candidate set contains 12 unique L1FLnI loci distributed across
Astrocytes, Müller glia, and RPE.

### 4. Locus-specific analysis

Each candidate L1FLnI locus is also analyzed independently.

Two measurements are retained:

- Mean normalized expression
- Detection rate

Detection-based measurements are considered complementary/exploratory because
they may be sensitive to the number of cells available for a donor/cell type.

---

## Analysis workflow

Run the scripts in approximately the following order:

    CellChat donor objects
            |
            v
    extract_donor_network_metrics.R
            |
            v
    QC_donor_network_metrics.R
            |
            v
    LINE1_donor_burden.R
            |
            v
    LINE1_CellChat_Pearson.R
            |
            +----------------------+
            |                      |
            v                      v
    Disease-adjusted         Disease-stratified
    analysis                 correlations
            |                      |
            v                      v
    LINE1_CellChat_          LINE1_CellChat_
    Disease_Adjusted.R       ByDisease.R
                                   |
                                   v
                        Figure_LINE1_ByDisease_Paper.R

---

## Scripts

### `extract_donor_network_metrics.R`

Extracts donor-level communication metrics from the completed CellChat objects.

Major outputs include:

- Number of cells
- Number of represented cell types
- Total ligand-receptor interactions
- Active edges
- Possible edges
- Network density
- Total communication strength
- Mean edge strength
- Mean ligand-receptor interactions per edge
- Number of signaling pathways

The principal output is:

    network_metrics/donor_network_metrics.csv

---

### `QC_donor_network_metrics.R`

Performs quality-control analyses of the donor-level CellChat measurements.

It evaluates relationships between the number of cells and network metrics and
generates diagnostic plots.

In the complete 16-donor dataset, no clear evidence was observed that total
cell number systematically explained donor-level Network Density or Total
Communication Strength.

This QC does not establish absence of sampling effects; cell-type-specific
coverage should also be considered.

---

### `LINE1_donor_burden.R`

Calculates donor-level expression of the candidate L1FLnI loci.

Expression is aggregated by:

    Donor × Cell Type × Locus

Outputs include:

- Mean normalized expression
- Median expression
- Detection rate
- Summed expression
- Candidate L1FLnI burden

Candidate burden is calculated at the donor × cell-type level.

---

### `LINE1_CellChat_Pearson.R`

Performs the primary donor-level Pearson correlation analysis between L1FLnI
measurements and CellChat network properties.

Analyses include:

1. Retina-wide family-level L1FLnI
2. Cell-type-specific family-level L1FLnI
3. Candidate L1FLnI burden
4. Candidate locus-specific analyses

Primary network outcomes:

- Network Density
- Total Communication Strength

Pearson correlation is used as the primary correlation statistic.

Spearman correlation may be used as a sensitivity analysis.

Multiple testing is evaluated using the Benjamini-Hochberg false discovery
rate (FDR).

---

### `LINE1_CellChat_Disease_Adjusted.R`

Tests whether associations between L1FLnI measurements and CellChat network
metrics remain after accounting for disease state.

The general model is:

    Network Metric ~ L1FLnI Measurement + Disease State

Healthy Control is used as the reference condition.

Outputs include:

- Regression coefficient
- Standard error
- 95% confidence interval
- Standardized beta
- Disease-adjusted P value
- Model R-squared
- Adjusted R-squared
- Benjamini-Hochberg FDR

These models help distinguish an L1FLnI-network association from differences
that are explained primarily by disease group.

---

### `LINE1_CellChat_ByDisease.R`

Performs donor-level Pearson correlations independently within:

- Healthy Control
- Dry AMD
- Wet AMD

This analysis evaluates whether L1FLnI-network relationships may differ across
disease states.

Important sample-size consideration:

    Healthy Control: n = 6
    Wet AMD: n = 7
    Dry AMD: n = 3

Because Dry AMD contains only three donors, Dry AMD correlations should be
considered **descriptive/exploratory** rather than strong inferential evidence.

A significant correlation in one condition and a non-significant correlation
in another does not by itself establish that the correlations differ between
conditions. Formal comparisons require an interaction model.

---

### `Figure_LINE1_ByDisease_Paper.R`

Generates publication-oriented figures showing donor-level L1FLnI–CellChat
relationships stratified by disease state.

Current panels include examples of:

- Retina-wide L1FLnI burden
- Photoreceptor L1FLnI burden
- Bipolar L1FLnI detection breadth
- Astrocyte locus-specific L1FLnI
- RPE locus-specific L1FLnI

Healthy and Wet AMD regression lines are shown as inferential/exploratory
within-group associations.

Dry AMD is displayed separately and should be interpreted descriptively because
of its small sample size.

---

## Statistical interpretation

The following principles should be maintained when interpreting these analyses.

### Donor is the biological replicate

Correlations are performed across donors.

Individual retinal cells must not be treated as independent replicates.

### Correlation does not establish causality

An association between L1FLnI expression and CellChat communication does not
demonstrate that LINE-1 activity causes network remodeling.

### Disease stratification

Within-condition analyses have small sample sizes:

    Healthy n = 6
    Wet AMD n = 7
    Dry AMD n = 3

Large Pearson correlations can therefore be strongly affected by individual
donors.

Leave-one-donor-out sensitivity analyses are recommended for the strongest
associations.

### RPE coverage

RPE is sparsely represented in this dataset.

Some donor-level RPE measurements are derived from very few cells.

RPE locus-specific results should therefore be interpreted as exploratory and
evaluated using cell-number sensitivity analyses.

### Multiple testing

Raw P values and Benjamini-Hochberg FDR should both be reported.

Nominal P < 0.05 should not be described as FDR-significant when the
corresponding adjusted P value does not meet the FDR threshold.

---

## Reproducibility

Before running these scripts, verify:

- Seurat object path
- CellChat object path
- Candidate L1FLnI table path
- Metadata column names
- R/Seurat/CellChat versions

The analysis was developed for the SRP413248 retinal snRNA-seq dataset.

