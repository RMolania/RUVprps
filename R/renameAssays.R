#' Removes datasets (assays) from a SummarizedExperiment object
#'
#' @param se.obj A SummarizedExperiment object.
#' @param new.names Character. A character to a vector of characters that specifies the name of data (assays) that
#' will be removed from the SummarizedExperiment object.
#'
#' @importFrom SummarizedExperiment assays
#'
renameAssays <- function(
        se.obj,
        new.names
){
    # Checking the input ####
    names(assays(se.obj)) <- new.names
    return(se.obj)
}
