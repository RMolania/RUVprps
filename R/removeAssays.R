#' Removes datasets (assays) from a SummarizedExperiment object
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assays.to.remove Character. A character to a vector of characters that specifies the name of data (assays) that
#' will be removed from the SummarizedExperiment object.
#'
#' @importFrom SummarizedExperiment assays
#'
#' @export

removeAssays <- function(
        se.obj,
        assays.to.remove
        ){
    # Checking the input ####
    if (sum(assays.to.remove %in% names(assays(se.obj)))!= length(assays.to.remove)){
        stop('All or some of the "assays.to.remove" cannot be found in the SummarizedExperiment object.')
    }
    for(i in assays.to.remove){
        se.obj@assays@data[[i]] <- NULL
    }
    return(se.obj)
}
