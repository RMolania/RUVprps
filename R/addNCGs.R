#' Add pre-selected sets of NCGs to SummarizedExperiment object.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function adds pre-selected sets of negative control genes (NCGs) to a SummarizedExperiment object.
#'
#' @details
#' A pre-selected set of negative control genes (NCGs) will be stored in the following location:
#' se.obj->metadata->NCG->pre.selected->subset.name. These genes can be used for various analyses, including identifying
#' unknown sources of variation, assessing variation, performing RUV normalization, and evaluating normalization steps. The
#' gene set will be stored in: metadata->NCG->pre.selected->subset.name->gene.set
#'
#' @param se.obj A SummarizedExperiment object.
#' @param ncg A logical value or a vector of gene names or IDs representing pre-selected NCGs. If gene names or IDs are
#' provided, the row names of the SummarizedExperiment object must match these names or IDs.
#' @param subset.name Character. A character that specifies the name of the NCG set in the metadata of the SummarizedExperiment
#' object. The default is set to `NULL`, in which case the function will generate a name as following:
#' `paste0(sum(ncg), '_psg')`.
#' @param verbose Logical. Indicates whether to display output messages during function execution. The default is set to
#' `TRUE`.
#'
#' @return A SummarizedExperiment object with a metadata that contains the NCGs.
#'
#' @importFrom SummarizedExperiment colData
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
            stop('None of the genes specified in the "ncg" can be found in the SummarizedExperiment object.')
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
            message = '- The number of genes provided may be too few for RUV-III normalization.',
            color = 'red',
            verbose = verbose
        )
    }
    # Adding the gene set to metadata of the SummarizedExperiment object ####
    if (is.null(subset.name)){
        subset.name <- paste0(sum(ncg), '_genes')
    }
    printColoredMessage(
        message = '-- Saving the provided NCGs to the metadata of the SummarizedExperiment object.',
        color = 'magenta',
        verbose = verbose
        )
    ## Check if metadata NCG already exists
    if (length(se.obj@metadata$NCG) == 0 ) {
        se.obj@metadata[['NCG']] <- list()
    }
    if (!'pre.selected' %in% names(se.obj@metadata[['NCG']])){
        se.obj@metadata[['NCG']][['pre.selected']] <- list()
    }
    if (!subset.name %in% names(se.obj@metadata[['NCG']][['pre.selected']])){
        se.obj@metadata[['NCG']][['pre.selected']][[subset.name]] <- list()
    }
    if (!'gene.set' %in% names(se.obj@metadata[['NCG']][['pre.selected']][[subset.name]])){
        se.obj@metadata[['NCG']][['pre.selected']][[subset.name]][['gene.set']] <- list()
    }
    se.obj@metadata[['NCG']][['pre.selected']][[subset.name]][['gene.set']] <- ncg
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
