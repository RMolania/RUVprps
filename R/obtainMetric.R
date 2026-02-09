#' Computes ANOVA between individual gene expression and a categorical variable
#'
#' @author Ramyar Molania
#'
#' @description
#' This function calculates the ANOVA between individual gene expression of the assay(s) in a SummarizedExperiment object
#' and a categorical variable as factor.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param slot TTTT
#' @param metric.group TTTT
#' @param assay.name TTTT
#' @param metric.name TTTT
#' @param test.name TTTT
#' @param variable TTTT
#' @param show.path TTTT
#'
#' @return Either a SummarizedExperiment object containing the log2 F-statistics and p-values of ANOVA for the continuous
#' variable or a list of these results.
#'
#' @export

obtainMetric <- function(
        se.obj,
        slot = 'Metrics',
        metric.group = 'gene.level',
        assay.name,
        metric.name,
        test.name,
        variable,
        show.path = TRUE
        ){
    metric <- se.obj@'metadata'[[slot]][[assay.name]][[metric.group]][[metric.name]][[test.name]][[variable]][['fstatistics.pvalues']]
    return(metric)
}
