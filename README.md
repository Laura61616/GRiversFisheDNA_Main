# GRiversFisheDNA_Main

**Data and codes for:**  
*Globally unified analysis of riverine fish eDNA reveals common associations of multi-faceted biodiversity with drainage characteristics*

**DOI:** to be updated

---

## Main Authors
Yan Zhang (Nanjing University, eawag, University of Zurich), Xiaowei Zhang (Yunnan University), Florian Altermatt (eawag, University of Zurich)

---

## Abstract
Freshwater biodiversity is declining at a pace that outstrips the capacity of existing monitoring approaches both in temporal and spatial dimensions, highlighting the urgent need for more scalable frameworks to track and understand biodiversity changes. Here, we present one of the first global assessments and unified analyses of riverine fish biodiversity using environmental DNA (eDNA) collected from 1818 sites across 113 river systems. We quantified species richness, functional redundancy, phylogenetic diversity, and genetic sequence diversity, and related them to drainage characteristics and human activities. Our results showed that eDNA effectively captured global patterns of biodiversity across these multiple facets. Biodiversity-area relationships were shaped by both climate and human activities. Catchments in warmer climates consistently enhanced biodiversity accumulation with area, while higher human activity intensity weakened this scaling. Species richness, functional, and genetic sequence diversity exhibited stronger negative responses to human activities in larger catchments. In contrast, phylogenetic diversity showed the strongest negative effects in smaller catchments with these impacts diminishing as catchment area increased, highlighting the facet-dependent nature of biodiversity responses to environmental gradients. Our findings demonstrate the power of eDNA-based datasets for harmonized, multi-faceted biodiversity assessments, offering a scalable approach for detecting and attributing biodiversity change and informing conservation strategies under accelerating global change.

---

## Keywords
riverine fish, eDNA metabarcoding, multi-faceted biodiversity, biodiversity-area relationship

---

## Software, Data, and Instructions for Use

**Overview**  
This repository contains all data and R codes used in the global analysis of riverine fish eDNA.  
The complete workflow, including system setup, installation, and analysis steps, is described below.  
Running the scripts in the listed order will reproduce the analyses and figures presented in the manuscript.

**System Requirements**  
- R (≥ 4.2.0)  
- Standard desktop or laptop computer  
- No GPU or non-standard hardware required

**Use Guide**  
Executing the scripts sequentially reproduces the full workflow and figures.  
All required **R packages** and all **custom functions** are contained and loaded in `00_functions.R`.

**Run the scripts in this order to complete one full analysis:**
1. `01_Modelling_alpha_HG.R` — Main alpha-diversity modeling  
2. `02_Modelling_alpha_HG_downsampling.R` — Downsampling-based modeling  
3. `02_Modelling_alpha_HG_reshuffle.R` — Reshuffling (randomization) analysis  
4. `03_Modelling_alpha_HG_predict.R` — Generates prediction outputs  
5. `Figure_*.R` — Generates main figures

**Notes**
- Basic modeling + plotting: ~1–2 hours  
- `02_Modelling_alpha_HG_downsampling.R` and `02_Modelling_alpha_HG_reshuffle.R` can be set to **1000 iterations**. This is computationally intensive and may take **days to weeks**, depending on your machine. For quick tests, lower the iteration count.
