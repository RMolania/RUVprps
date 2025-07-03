#' Renames variables in SummarizedExperiment
#'
#'
#' @param se.obj description
#' @param current.names description
#' @param new.names description
#' @param keep.old.names description

renameVariables <- function(
        se.obj,
        current.names,
        new.names,
        keep.old.names = FALSE
        ){
    # Checking the input ####
    sample.annotation <- colData(se.obj)
    for(i in 1:length(current.names)){
        colnames(sample.annotation)[colnames(sample.annotation) == current.names[i] ] <- new.names[i]
    }
    colData(se.obj) <- DataFrame(sample.annotation, check.names = FALSE)
    return(se.obj)
}
