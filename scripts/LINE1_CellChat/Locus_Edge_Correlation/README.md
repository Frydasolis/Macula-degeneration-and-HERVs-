# Locus-specific LINE-1 × CellChat edge correlation analysis

## Overview

This analysis tests whether donor-to-donor variation in the expression of individual **L1FLnI loci** is associated with variation in specific retinal cell-cell communication edges inferred with CellChat.

The analysis moves beyond condition-level comparisons of CellChat networks to determine whether individual LINE-1 loci show locus-specific relationships with retinal intercellular communication.

The biological replicate and unit of inference is the **donor**, not the individual cell.

---

## Analysis design

For each donor, candidate L1FLnI locus expression was quantified across retinal cells and integrated with the corresponding donor-level CellChat network.

For every combination of:

- L1FLnI locus
- disease condition
- sender cell type
- receiver cell type

we tested the association between donor-level locus expression and donor-level CellChat communication weight.

In other words, for each CellChat sender-to-receiver edge, we asked:

**Do donors with higher expression of a particular L1FLnI locus also show stronger or weaker communication in that edge?**

---

## Candidate loci

The analysis includes 12 unique L1FLnI loci:

- L1FLnI-11p15.5b
- L1FLnI-15q11.2n
- L1FLnI-2p24.3g
- L1FLnI-21q22.13b
- L1FLnI-3q13.13ea
- L1FLnI-3q13.13ja
- L1FLnI-3q13.13ka
- L1FLnI-3q13.13ma
- L1FLnI-4p15.2b
- L1FLnI-6q16.2d
- L1FLnI-Xq13.2ia
- L1FLnI-Xq23ab

These loci were selected from the locus-specific differential expression analysis.

---

## Disease groups

The donor-level analysis includes:

- Healthy Control: n = 6 donors
- Dry AMD: n = 3 donors
- Wet AMD: n = 7 donors

Because only three Dry AMD donors were available, correlations within Dry AMD are considered **descriptive/exploratory** and are not interpreted as robust inferential evidence.

---

## Statistical analysis

### Primary analysis: Pearson correlation

Pearson correlation was used as the primary analysis.

For each:

**locus × condition × sender × receiver**

the analysis correlates:

**donor-level L1FLnI expression**

with:

**donor-level CellChat communication weight**

The analysis reports:

- number of donors
- Pearson correlation coefficient (r)
- nominal P value
- Benjamini-Hochberg FDR

Multiple-testing correction was performed using the Benjamini-Hochberg procedure.

### Sensitivity analysis

Spearman correlation was additionally calculated as a sensitivity analysis.

---

## Interpretation of the locus-specific networks

Each network represents associations between expression of one L1FLnI locus and donor-level CellChat communication strength.

### Edge color

- **Red:** positive Pearson correlation
- **Blue:** negative Pearson correlation

A positive correlation means that donors with higher expression of the locus tend to have higher CellChat communication weight for that sender-receiver interaction.

A negative correlation means that donors with higher locus expression tend to have lower communication weight.

### Edge appearance

- **Thick / opaque:** FDR < 0.05
- **Thin / transparent:** nominal P < 0.05 but FDR >= 0.05
- **P >= 0.05:** not displayed

Edge width reflects the absolute Pearson correlation coefficient, |r|.

These associations represent covariation and should **not** be interpreted as evidence that LINE-1 expression causally changes cell-cell communication.

---

## Main results

A total of **35 locus-edge associations passed FDR < 0.05** in the primary Pearson analysis.

All 35 FDR-significant associations had positive Pearson correlation coefficients.

The significant associations were concentrated in four L1FLnI loci:

| Condition | L1FLnI locus | FDR-significant edges |
|---|---|---:|
| Healthy Control | L1FLnI-3q13.13ka | 11 |
| Wet AMD | L1FLnI-3q13.13ja | 7 |
| Wet AMD | L1FLnI-2p24.3g | 6 |
| Wet AMD | L1FLnI-3q13.13ea | 6 |
| Wet AMD | L1FLnI-3q13.13ka | 5 |

Thus, Healthy Control contained 11 FDR-significant locus-edge associations, whereas Wet AMD contained 24.

No Dry AMD association survived FDR correction. Because the Dry AMD group contained only three donors, this absence should not be interpreted as evidence that no biological association exists.

---

## Biological interpretation

Several loci showing significant communication associations were previously identified as differentially expressed in RPE, including:

- L1FLnI-3q13.13ka
- L1FLnI-3q13.13ja
- L1FLnI-3q13.13ea
- L1FLnI-2p24.3g

The results therefore indicate that donor-level variation in selected RPE-associated L1FLnI loci covaries with the strength of specific retinal intercellular communication edges.

Importantly, the locus expression used in the correlation analysis is quantified at the donor level across retinal cells. Therefore, these results do not establish that RPE-specific LINE-1 expression directly controls communication between particular cell populations.

The analysis identifies associations between locus expression and communication architecture, not causal effects.

---

## Difference from condition-level CellChat analysis

The condition-level CellChat analysis asks:

**How does retinal cell-cell communication differ between Healthy Control, Dry AMD, and Wet AMD?**

The locus-specific analysis asks a different question:

**Within a disease state, do donors with higher expression of a particular L1FLnI locus also show stronger or weaker communication in specific retinal cell-cell interactions?**

Therefore, each locus produces its own communication-association network.

---

## Scripts

### `10_Locus_Edge_Correlations.R`

Performs the donor-level locus-specific correlation analysis.

Main steps:

1. Loads the annotated retinal Seurat object.
2. Extracts the 12 candidate L1FLnI loci.
3. Calculates donor-level global locus expression and detection rate.
4. Integrates locus expression with donor-level CellChat sender-receiver weights.
5. Calculates Pearson correlations.
6. Calculates Spearman correlations as a sensitivity analysis.
7. Applies Benjamini-Hochberg FDR correction.
8. Exports complete, nominally significant, and FDR-significant results.

### `11_Figure_LocusSpecific_Pearson_Networks.R`

Generates locus-specific Pearson correlation networks.

For each of the 12 loci, a three-panel figure is generated:

**Healthy Control | Dry AMD | Wet AMD**

Each network displays sender-receiver CellChat edges whose donor-level communication weights are associated with donor-level expression of that locus.

---

## Key conclusion

Higher donor-level expression of selected L1FLnI loci was associated with greater communication strength across specific retinal cell-cell interactions.

These associations were particularly evident in Wet AMD and were concentrated among a subset of RPE-associated L1FLnI loci.

The results support a relationship between locus-specific LINE-1 activity and retinal communication architecture while remaining consistent with an associative, rather than causal, interpretation.
