#' Get a metric or plot from the `metadata` of a SummarizedExperiment object.
#'
#' @author Ramyar Molania
#'
#' @param se.obj A SummarizedExperiment object.
#' @param slot Character string. Indicates the name of the slot in the `metadata` of the SummarizedExperiment object.
#' @param assay.names Character string or vector. Specifies the name(s) of the assay(s) in the SummarizedExperiment
#' object. Default is `all`.
#' @param assessment.type Character string. Specifies the type of assessment. Options are `gene.level` or `global.level`.
#' @param assessment Character string. Name of the metric to retrieve.
#' @param method Character string. Method used to calculate the metric.
#' @param variables Character string or vector. Specifies the variable(s) used to calculate the metric.
#' @param file.name Character string. Name of the file to which the metric results were saved.
#' @param sub.file.name Character string. Name of the sub-file (if applicable) to which the results were saved.
#' @param required.function Character string. Name of the function that must have been run prior to retrieving this metric.
#' @param message.to.print Character string. A message to print when retrieving the metric.
#' @param verbose Logical. If `TRUE`, prints messages describing the function's progress.
getMetricFromSeObj <- function(
        se.obj,
        slot = 'Metrics',
        assay.names,
        assessment.type = 'gene.level',
        assessment,
        method = NULL,
        variables,
        file.name,
        sub.file.name = NULL,
        required.function,
        message.to.print,
        verbose = TRUE
        ){
    if(is.null(required.function)){
        required.function <- 'required function'
    }
    # Create an empty list ####
    all.outputs <- list()

    # Check the metadata ####
    if (length(se.obj@metadata) == 0) {
        stop('The current SummarizedExperiment object does not contain "metadat".')
    }
    # Check the slot ####
    if (!slot %in% names(se.obj@metadata)) {
        stop(paste0('The "metadat" of the current SummarizedExperiment object does not contain ', slot ,' slot.'))
    }
    # Check the all ####
    for(x in assay.names){
        if (!x %in% names(se.obj@metadata[[slot]])) {
            stop(paste0('The "', slot, '" of in the "metadata" of the current SummarizedExperiment object does not contain any metrics for "', x ,'" data.'))
        }
        if (!assessment.type %in% names(se.obj@metadata[[slot]][[x]])) {
            stop(paste0('The metrics of for "', x , '" data does not contain any "', assessment.type ,'" assessments.'))
        }
        if (!assessment %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]])) {
            stop(paste0('The ', assessment.type , ' assessments of the "', x, '" data does not contain any ', assessment ,' data. ',
                        'Please run the "', required.function, '" function first.'))
        }
        if (!method %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]])) {
            stop(paste0('The ', assessment , ' metric of the "', x, '" data does not contain any data for the "', assessment ,'" method.',
                        'Please check the parameters of the "', required.function, '" function.'))
        }
        if (!variables %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]])) {
            stop(paste0('The ', assessment , ' metric of the "', x, '" data with "', method , '" method ',
                        'does not contain any data for the "', variables ,'" variable.',
                        'Please check the parameters of the "', required.function, '" function.')
                 )
        }
        if (!file.name %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]]) ) {
            stop(paste0('The ', assessment , ' metric of the "', x, '" data with "', method , '" method ', 'for the "',
                        variables ,'" variable.', 'does not contain any ', file.name, 'file.',
                        'Please check the parameters of the "', required.function, '" function.')
                 )
        }
        if(!is.null(sub.file.name)){
            if (!sub.file.name %in% names(se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]]) ) {
                stop(paste0('The ', assessment , ' metric of the "', x, '" data with "', method , '" method ', 'for the "',
                            variables ,'" variable.', ' with the ', file.name, 'file ', 'does not contain ', sub.file.name, 'file.',
                            'Please check the parameters of the "', required.function, '" function.')
                )
            }

            all.outputs[[x]] <- se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]][[sub.file.name]]
        } else {
            all.outputs[[x]] <- se.obj@metadata[[slot]][[x]][[assessment.type]][[assessment]][[method]][[variables]][[file.name]]
        }
        printColoredMessage(
            message = paste0('- Obtaining the ', message.to.print, ' the "', x , '" data.'),
            color = 'blue',
            verbose = verbose
        )
    }
    return(all.outputs)
}
