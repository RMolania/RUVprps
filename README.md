
# RUVprps

<!-- badges: start -->
<!-- badges: end -->

**RUVprps** implements *RUV-III with pseudo-replicates of pseudo-samples (PRPS)* — a novel strategy for transcriptomics data normalization when technical replicates are unavailable or poorly designed.  

This user-friendly R package provides an end-to-end workflow for removing unwanted variation from large-scale transcriptomic datasets, whether derived from a single study or multiple studies. RUV-III effectively corrects for sources of variation such as library size, batch effects, and tumor purity. The package accommodates technical replicates, pseudo-replicates (PR), and pseudo-replicates of pseudo-samples (PRPS).  

---

## ✨ Key Features

- Comprehensive diagnostic and assessment tools to evaluate both biological and unwanted variation in RNA-seq data  
- Robust strategies to identify unknown sources of unwanted variation  
- Fast and scalable implementation of RUV-III for efficient normalization of large datasets  
- Novel unsupervised methods for identifying PRPS and negative control genes (NCGs)  
- Normalization performance summary tables, helping users select the most appropriate strategy  

---

## 📊 Overview

![RUVprps workflow](figures/RUVprps_MainFigure2.png)

---

## 📖 Citation

If you use **RUVprps**, please cite:  

Molania R, Foroutan M, Gagnon-Bartsch JA, Gandolfo LC, Jain A, Sinha A, Olshansky G, Dobrovic A, Papenfuss AT, Speed TP.  
*Removing unwanted variation from large-scale RNA sequencing data with PRPS.* **Nat Biotechnol.** 2023;41(1):82–95.  
[doi:10.1038/s41587-022-01440-w](https://doi.org/10.1038/s41587-022-01440-w)  

---

## ⚙️ Installation

After installing the required dependencies, install **RUVprps** from GitHub with:  

```r
library(devtools)
devtools::install_github(
  repo = "RMolania/RUVIIIPRPS",
  force = TRUE,
  build_vignettes = FALSE
)
