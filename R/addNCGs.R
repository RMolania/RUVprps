#' Add pre-selected set of negative control genes to SummarizedExperiment object.

#' @author Ramyar Molania

#' @description
#' This function adds a pre-selected set of negative control genes to SummarizedExperiment object.

#' @param se.obj A SummarizedExperiment object.
#' @param ncg A logical or vector of gene names of pre-selected genes as NCGs.
#' @param save.se.obj Symbol. Indicates
#' @param output.name Symbol. Indicates the name of the output in the meta data of the SummarizedExperiment object. The
#' default is 'NULL'. This means the function will create a name based the specified argument.
#' @param verbose Logical. If TRUE, displaying process messages is enabled.

#' @importFrom SummarizedExperiment colData
#' @importFrom BiocSingular bsparam
#' @importFrom Matrix rowSums colSums
#' @importFrom ruv replicate.matrix
#' @export

addNCGs <- function(
        se.obj,
        ncg,
        save.se.obj = TRUE,
        output.name = NULL,
        verbose = TRUE
        ){
    # Check inputs ####
    if (is.logical(ncg)){
        if (length(ncg) > nrow(se.obj)){
            stop('The length of the "ncg" is larger than the number of rows in the SummarizedExperiment object.')
        } else if (sum(ncg) == 0){
            stop('The "ncg" does not contain any "TRUE" value.')
        }
        printColoredMessage(
            message = paste0('- ', sum(ncg), 'NCGs genes are found in the SummarizedExperiment object.' ),
            color = 'blue',
            verbose = verbose
            )
    } else if (is.character(ncg) | is.factor(ncg)){
        ncg <- intersect(unique(ncg), row.names(se.obj))
        if(length(ncg) == 0){
            stop('None of the genes specified in the "ncg" can be found in the SummarizedExperiment object. ')
        }
        printColoredMessage(
            message = paste0(length(ncg), ' NCGs genes are found in the SummarizedExperiment object.' ),
            color = 'blue',
            verbose = verbose
            )
        ncg <- row.names(se.obj) %in% ncg
    } else if (is.numeric(ncg)){
        if(max(ncg) > nrow(se.obj)){
            stop('- The "ncg" contains some number(s) that is larger than the number of rows in the SummarizedExperiment object.')
        }
        printColoredMessage(
            message = paste0(length(ncg), ' NCGs genes are found in the SummarizedExperiment object.' ),
            color = 'blue',
            verbose = verbose)
        ncg.log <- rep(FALSE, nrow(se.obj))
        ncg.log[ncg] <- TRUE
        ncg <- ncg.log
    }

    if (isTRUE(assess.ncg)){
        if(is.null(variables.to.assess.ncg)){
            stop('The "variables.to.assess.ncg" must be provided to assess the performance of the NCG')
        }
        if(!variables.to.assess.ncg %in% colnames(colData(se.obj)) ){
            stop('Some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
        }
        if(is.null(assay.name)){
            stop('The "assay.name" must be provided to assess the performance of the NCG')
        }
    }
    # Save the results ####
    if(is.null(output.name)){
        output.name <- paste0(sum(ncg), '_genes')
    }
    if(save.se.obj == TRUE){
        printColoredMessage(
            message = '-- Saving the selected NCG to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose)
        ## Check if metadata NCG already exists
        if(length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        if(!'pre.selected' %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][['pre.selected']] <- list()
        }
        if(!output.name %in% names(se.obj@metadata[['NCG']][['pre.selected']])){
            se.obj@metadata[['NCG']][['pre.selected']][[output.name]] <- list()
        }
        if(!'gene.list' %in% names(se.obj@metadata[['NCG']][['pre.selected']][[output.name]])){
            se.obj@metadata[['NCG']][['pre.selected']][[output.name]][['gene.list']] <- list()
        }
        se.obj@metadata[['NCG']][['pre.selected']][[output.name]][['gene.list']] <- ncg

        if(isTRUE(assess.ncg)){
            if(!'assessment.plot' %in% names(se.obj@metadata[['NCG']][['pre.selected']][[output.name]])){
                se.obj@metadata[['NCG']][['pre.selected']][[output.name]][['assessment.plot']] <- list()
            }

        }
        printColoredMessage(
            message = 'The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The supervisedFindNcgTWAnova function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
}
