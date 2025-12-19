#' Selects genes for partial correlation and gene scoring analysis.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function selects genes from a `SummarizedExperiment` object for partial correlation and gene scoring
#' analysis, based on correlation thresholds, ranking strategies, or group-specific variables.
#'
#' @param se.obj A `SummarizedExperiment` object containing gene expression data.
#' @param assay.names Character. The name(s) of the assay(s) in the `SummarizedExperiment` object to be used.
#' @param variables Character. The variable(s) in the sample annotation to be considered for correlation analysis.
#' @param groups Character. A grouping variable in the sample annotation used to stratify the analysis.
#' @param cor.cutoff Numeric. A correlation threshold for selecting genes (e.g., 0.3).
#' @param abs.cor Logical. If `TRUE`, absolute correlation values are used when applying the cutoff. The default is `TRUE`.
#' @param top.genes Numeric. The number of top-ranked genes to select based on correlation or scoring. If `NULL`, all
#' genes above the cutoff are retained.
#' @param method Character. The method to use for correlation or scoring (e.g., `"pearson"`, `"spearman"`, `"kendall"`).
#' @param a Numeric. A tuning parameter (e.g., weighting factor) used in the gene scoring function.
#' @param rho Numeric. A correlation parameter used in scoring or thresholding (e.g., shrinkage factor).
#' @param apply.log Logical. Indicates whether to log-transform the data prior to correlation analysis. The default is `TRUE`.
#' @param pseudo.count Numeric. A pseudo-count to add before log-transformation. The default is 1.
#' @param check.se.object Logical. If `TRUE`, validates the structure of the `SummarizedExperiment` object before running
#' the analysis. The default is `TRUE`.
#' @param remove.na Character. Specifies how to handle missing values. Options are: `"genes"`, `"samples"`, `"both"`,
#' or `"none"`. The default is `"genes"`.
#' @param verbose Logical. If `TRUE`, messages describing the progress of the function will be displayed. The default is `TRUE`.
#'
#' @return A list containing the selected genes and their associated correlation or scoring statistics.
#'
#' @export

selectGenesForPPcorr <- function(
        se.obj,
        assay.names,
        variables,
        groups = NULL,
        cor.cutoff = 0.5,
        abs.cor = TRUE,
        top.genes = NULL,
        method = "spearman",
        a = 0.05,
        rho = 0,
        apply.log = TRUE,
        pseudo.count = 1,
        check.se.object = TRUE,
        remove.na = 'both',
        verbose = TRUE
        ){
    names(variables) <- assay.names
    if (!is.null(groups)){
        groups <- se.obj[[groups]]
        all.corr <- list()
        for(i in 1:length(variables)){
            corr.genes <- sapply(
                unique(groups),
                function(x){
                    gene.var.corr <- computeGenesVariableCorrelation(
                        se.obj = se.obj[ , groups == x],
                        assay.names = names(variables)[i],
                        variable = variables[i],
                        method = method,
                        a = a, rho = rho,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        check.se.obj = FALSE,
                        save.se.obj = FALSE
                        )[[names(variables)[i]]]
                    if (isTRUE(abs.cor)){
                           abs(gene.var.corr[ , 'correlation']) > cor.cutoff
                    } else gene.var.corr[ , 'correlation'] > cor.cutoff
                })
            all.corr[[variables[i]]] <- row.names(corr.genes)[rowSums(corr.genes)!= 0]
        }
    }
    return(all.corr)
}



