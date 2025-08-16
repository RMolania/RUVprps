#' Order SummarizedExperiment object based on variables
#'
#' @author Ramyar Molania
#'
#' @param se.obj A SummarizedExperiment object.
#' @param factors.to.order TTT
#' @return The function returns a log2 transformed of all specified data sets as a list object.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr arrange
#'
#' @export

orderSeObj <- function(
        se.obj,
        factors.to.order
        ){
    if (sum(factors.to.order %in% colnames(colData(se.obj))) != length(factors.to.order)){
        stop('All or some of the "factors.to.order" cannot be found in the SummarizedExperiment object.')
    }
    sample.annotation <- as.data.frame(colData(se.obj))
    sample.annotation <- dplyr::arrange(sample.annotation, !!!syms(factors.to.order))
    se.obj <- se.obj[ , row.names(sample.annotation)]
    return(se.obj)
}
