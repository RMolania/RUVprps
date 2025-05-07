#' Adds pre-selected sets of NCGs to SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function adds pre-selected sets of negative control genes (NCGs) to a SummarizedExperiment object.
#' These genes can be used for various analyses, including identifying unknown sources of variation, assessing variation,
#' performing RUV normalization, and evaluating normalization steps.

#' @details
#' A pre-selected set of negative control genes (NCGs) will be stored in the following location:
#' se.obj->metadata->NCG->pre.selected->subset.name.


#' @param se.obj A `SummarizedExperiment` object.
#' @param ncg A logical or vector of gene names of pre-selected genes as NCGs.
#' @param subset.name Character. Specifies the name of the NCG set in the metadata of the `SummarizedExperiment` object.
#' The default is 'NULL', in which case the function will generate a name  based on the number of NCGs using:
#' paste0(sum(ncg), '_genes').
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @importFrom SummarizedExperiment colData
#' @importFrom BiocSingular bsparam
#' @importFrom Matrix rowSums colSums
#' @importFrom ruv replicate.matrix
#' @export

addNCGs <- function(
        se.obj,
        ncg,
        subset.name = NULL,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The addNCGs function starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking the function inputs ####
    if (is.logical(ncg)){
        if (length(ncg) > nrow(se.obj)){
            stop('The length of the "ncg" is larger than the number of rows in the SummarizedExperiment object.')
        } else if (sum(ncg) == 0){
            stop('The "ncg" does not contain any "TRUE" value.')
        }
        printColoredMessage(
            message = paste0(
                '- ',
                sum(ncg),
                ' genes are provided as NCGs.'),
            color = 'blue',
            verbose = verbose
            )
    }
    if (is.character(ncg) | is.factor(ncg)){
        if (length(ncg) != length(unique(ncg))){
            printColoredMessage(
                message = 'The names/ids of provided gene set are not all unique.',
                color = 'yellow',
                verbose = verbose
            )
        }
        ncg <- intersect(unique(ncg), row.names(se.obj))
        if (length(ncg) == 0){
            stop('None of the genes specified in the "ncg" can be found in the SummarizedExperiment object. ')
        }
        printColoredMessage(
            message = paste0(
                '- ',
                length(ncg),
                ' genes are provided as NCGs.'),
            color = 'blue',
            verbose = verbose
        )
        ncg <- row.names(se.obj) %in% ncg
    }
    if (is.numeric(ncg)){
        if (max(ncg) > nrow(se.obj)){
            stop('- The "ncg" contains some number(s) that is larger than the number of rows in the SummarizedExperiment object.')
        }
        printColoredMessage(
            message = paste0(
                '- ',
                length(ncg),
                ' genes are provided as NCGs.'),
            color = 'blue',
            verbose = verbose
            )
        ncg.log <- rep(FALSE, nrow(se.obj))
        ncg.log[ncg] <- TRUE
        ncg <- ncg.log
    }
    # Checking the number of genes ####
    if (sum(ncg.log) < .01*nrow(se.obj)){
        printColoredMessage(
            message = '* The number of genes provided may be too few for RUV-III normalization.',
            color = 'blue',
            verbose = verbose
        )
    }
    # Adding the gene set to metadata of the SummarizedExperiment object ####
    if (is.null(subset.name)){
        subset.name <- paste0(sum(ncg), '_genes')
    }
    if (save.se.obj == TRUE){
        printColoredMessage(
            message = '-- Saving the provided NCG to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose)
        ## Check if metadata NCG already exists
        if (length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        if (!'pre.selected' %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][['pre.selected']] <- list()
        }
        if (!output.name %in% names(se.obj@metadata[['NCG']][['pre.selected']])){
            se.obj@metadata[['NCG']][['pre.selected']][[output.name]] <- list()
        }
        se.obj@metadata[['NCG']][['pre.selected']][[output.name]] <- ncg
        printColoredMessage(
            message = 'The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The addNCGs finished',
            color = 'white',
            verbose = verbose
        )
        return(se.obj)
    }
}
