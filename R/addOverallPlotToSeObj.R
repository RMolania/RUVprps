#' Add overall plot to the "metadata" of SummarizedExperiment object.
#'
#' @author Ramyar Molania
#'
#' @param se.obj A SummarizedExperiment object.
#' @param slot Character A character string  indicating the name of the slots in the 'metadata' of the SummarizedExperiment object.
#' @param assessment.type Character A character string  indicating the type of assessment. Options are `gene.level` or
#' `global.level`.
#' @param assessment Character A character string  indicating the name of the metric to be checked.
#' @param method Character A character string  indicating the method used to calculate the metric.
#' @param variables Character A character string or a vector of symbols indicating the variables used to calculate the metric.
#' @param file.name Character A character string indicating the file name to which the results of the metric are assigned.
#' @param plot.data Character A character string indicating the name of the plots file for the metric.

addOverallPlotToSeObj <- function(
        se.obj,
        slot = 'Plots' ,
        assessment.type = 'gene.level',
        assessment,
        method,
        variables,
        file.name,
        plot.data
        ){
    if (!slot %in%  names(se.obj@metadata)) {
        se.obj@metadata[[slot]] <- list()
    }
    if (!assessment.type %in%  names(se.obj@metadata[[slot]]) ) {
        se.obj@metadata[[slot]][[assessment.type]] <- list()
    }
    if (!assessment %in%  names(se.obj@metadata[[slot]][[assessment.type]]) ) {
        se.obj@metadata[[slot]][[assessment.type]][[assessment]] <- list()
    }
    if (!method %in%  names(se.obj@metadata[[slot]][[assessment.type]][[assessment]] )) {
        se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]] <- list()
    }
    if (!variables %in%  names(se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]] )) {
        se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]][[variables]] <- list()
    }
    if (!file.name %in%  names(se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]][[variables]] )) {
        se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]] <- list()
    }
    se.obj@metadata[[slot]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]] <- plot.data
    return(se.obj)
}


