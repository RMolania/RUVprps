
# RUVprps software package

<!-- badges: start -->
<!-- badges: end -->

RUVprps is a comprehensive and user-friendly software package that streamlines the use of the RUV-III method for 
end-to-end trasncriptomics data normalization, irrespective of whether data comes from a single study or multiple studies. 
RUV-III is a powerful method that effectively removes sources of unwanted variation such as library size, batch effects,
and tumor purity, ... from complex transcriptomic datasets. RUVprps can accommodate technical replicates, pseudo-replicate (PR) and 
pseudo-replicates of pseudo-samples (PRPS). 


The key and novel contributions of RUVprps: 
\vspace{-1em} 
\begin{itemize} 
\item Comprehensive diagnostic and assessment tools for evaluating both biological and unwanted variation in RNA-seq data.
\item Robust strategies to identify sources of unwanted variation, even when such factors are unknown.
\item A fast and scalable implementation of RUV-III, enabling efficient normalization of large-scale datasets.
\item Novel unsupervised methods for identifying PRPS and negative control genes (NCGs), allowing greater flexibility and applicability to real-world datasets.
\item A normalization performance summary table, providing clear guidance to help users select the most appropriate normalization strategy.     
\end{itemize}


The RUVIII-PRPS software package itself can be cited as:

Molania R, Foroutan M, Gagnon-Bartsch JA, Gandolfo LC, Jain A, Sinha A, Olshansky G, Dobrovic A, Papenfuss AT, Speed TP. Removing unwanted variation from large-scale RNA sequencing data with PRPS. Nat Biotechnol. 2023 Jan;41(1):82-95. doi: 10.1038/s41587-022-01440-w. Epub 2022 Sep 15. PMID: 36109686; PMCID: PMC9849124.

##  RUVprps Installation

After installing the dependent libraries, RUVIII-PRPS can be installed by running the following lines

``` r
library(devtools)
devtools::install_github(
    repo = 'RMolania/RUVIIIPRPS',
    force = TRUE,
    build_vignettes = FALSE)
```

## Using RUVIII-PRPS to remove unwanted variation from large-scale RNA sequencing data

RUVIII-PRPS is a novel strategy using pseudo-replicates of pseudo-samples (PRPS) to normalize RNA-seq data in situations when technical replicate are not available or well-designed. 
We provided a vignette Introduction_to_RUVIII-PRPS.Rmd that explains step by step how to load and normalise datasets and also how to visualise the diagnostic plots before and after normalisation.


![Description of figure](figures/RUVprps_MainFigure2.png)

##  Note
The RUVIIIPRPS is under final preparation, a manuscript and comprehensive vignettes will be available soon.
