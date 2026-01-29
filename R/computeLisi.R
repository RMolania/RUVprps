#' Compute local inverse Simpson's index (LISI).
#'
#' @author Ramyar Molania
#'
#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @description
#' This functions computes local inverse Simpson's index (LISI) for any given variable using the first PCs of the specified
#' data(assay) in a SummarizedExperiment object.
#'
#' @details
#' The ARI is the corrected-for-chance version of the Rand index. The ARI measures the percentage of matches between
#' two label lists. We use the ARI to assess the performance of normalization methods in terms of sample subtype
#' separation and batch mixing. We first calculate PCs and use the first PCs to perform ARI.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names character or character vector. A character or a vector of characters of name(s) of the data(assay)
#' in the SummarizedExperiment object to use for ARI computation. The default is set to `all`, meaning all assays in the
#' object will be selected.
#' @param variable character. The name of the column containing a variable in the SummarizedExperiment object.
#' This variable can represent either a biological or an unwanted factor.
#' @param perplexity Numeric. The effective number of each cell's neighbors. The default is set to 10.
#' @param nn.eps Numeric. Error bound for nearest neighbor search with RANN:nn2(). The default is set to 0.0, implies exact
#' nearest neighbor search.
#' @param fast.pca Logical. Indicates whether to use principal components computed via fast PCA by the `comptePCA` function
#' The default is set to `TRUE`. Note that using fast PCA or standard PCA does not affect silhouette coefficient calculation.
#' See details for more information.
#' @param nb.pcs Numeric. Number of first principal components to use when calculating distances between samples. The efault
#' is to 3.
#' @param save.se.obj Logical. If `TRUE`, saves the results in the metadata of the SummarizedExperiment objet; otherwise,
#' returns the result as list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, displays messages for each step of the function.
#'
#' @return A SummarizedExperiment object or a list containing the computed ARI for specified categorical variable.
#'
#' @importFrom SummarizedExperiment assays assay
#' @importFrom lisi compute_lisi
#' @export
#'
computeLisi <- function(
        se.obj,
        assay.names = 'all',
        variable,
        perplexity = 10,
        nn.eps = 0,
        fast.pca = TRUE,
        nb.pcs = 3,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeLisi function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check the function inputs ####
    if (is.null(assay.names) | is.logical(assay.names)) {
        stop('The "assay.names" cannot be NULL or logical.')
    }
    if (is.list(assay.names)){
        stop('The "assay.names" must be a vector of the data ste names(s) or set to "all".')
    }
    if (is.null(variable)) {
        stop('The "variable" must be provided.')
    }
    if (!is.character(variable)){
        stop('The "variable" must be a character.')
    }
    if (length(variable) > 1){
        stop('The "variable" must contain only the name of a single variable in the SummarizedExperiment object.')
    }
    if (!variable %in% colnames(se.obj@colData)){
        stop('The "variable" cannot be found in the SummarizedExperiment object.')
    }
    if (length(unique(se.obj[[variable]])) == 1) {
        stop('The "variable" must have at least two levels.')
    }
    if (is.null(nb.pcs)) {
        stop('The "nb.pcs" must be specified.')
    }
    if (!is.numeric(nb.pcs) | nb.pcs < 0){
        stop('The "nb.pcs" must be a positive numerical value.')
    }
    if (sum(is.na(se.obj@colData[[variable]])) > 0){
        stop(paste0('The "', variable, '" contains NA.',
                    ' Run the checkSeObj function with "remove.na = both"',
                    ', then "computePCA"-->"computeARI".'))
    }

    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Computing the LISI on all the data set(s) ####
    printColoredMessage(
        message = paste0(
            '-- Computing local inverse Simpsons index (lisi) using the first ',
            nb.pcs,
            ' PCS for the "',
            variable,
            '" variable.') ,
        color = 'magenta',
        verbose = verbose
        )
    if (isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'

    ## Retrieving the PCA data from SummarizedExperiment ####
    all.pca.data <- getMetricFromSeObj(
        se.obj = se.obj,
        assay.names = levels(assay.names),
        slot = 'Metrics',
        assessment = 'PCA',
        assessment.type = 'global.level',
        method = method,
        variables = 'general',
        file.name = 'data',
        sub.file.name = 'svd',
        required.function = 'computePCA',
        message.to.print = 'PCs'
        )

    ## Computing LISI for all assay(s) ####
    all.lisi <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0(
                    '- Computing the LISI for the "',
                    x,
                    '" data:'),
                color = 'blue',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    '- Obtaining the first ',
                    nb.pcs,
                    ' computed PCs.'),
                color = 'blue',
                verbose = verbose
                )
            # Applying a sanity check ####
            printColoredMessage(
                message = '- Applying a sanity check.',
                color = 'blue',
                verbose = verbose
                )
            pca.data <- all.pca.data[[x]]$u
            if (ncol(pca.data) < nb.pcs){
                printColoredMessage(
                    message = paste0(
                        '- The number of PCs of the assay',
                        x,
                        'are
                        ', ncol(pca.data),
                        '.'),
                    color = 'blue',
                    verbose = verbose
                )
                stop(paste0(
                    'The number of PCs of the assay ',
                    x,
                    ' are less than',
                    nb.pcs,
                    '.',
                    'Re-run the computePCA function with nb.pcs = ',
                    nb.pcs,
                    '.'))
            }
            pca.data <- pca.data[ , seq_len(nb.pcs)]
            if (!all.equal(row.names(pca.data), colnames(se.obj))){
                stop('The column names of the SummarizedExperiment object is not the same as row names of the PCA data.')
            }
            # Clustering the data using the mclust method ####
            metadata <- data.frame(selected.var = se.obj[[variable]])
            lisi.score <- compute_lisi(
                X = pca.data,
                meta_data = metadata,
                label_colnames = "selected.var",
                perplexity = perplexity,
                nn_eps = nn.eps
                )
            lisi.score <- lisi.score[["selected.var"]]
            return(lisi.score)
        })
    names(all.lisi) <- levels(assay.names)

    # Saving the results ####
    printColoredMessage(
        message = '-- Saving all the lisi results:',
        color = 'magenta',
        verbose = verbose
        )
    ## Adding results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving the lisi of each data set(s) to the metadata in the SummarizedExperiment object:',
            color = 'blue',
            verbose = verbose
            )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'global.level',
            assessment = 'LISI',
            method = method,
            variables = variable,
            file.name = 'lisi',
            results.data = all.lisi
            )
        printColoredMessage(
            message = paste0(
                '- The lisi of induvial dataset(s) is saved to the .',
                ' ".se.obj@metadata$metric$RawCount$ARI" in the SummarizedExperiment objec.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The computeLisi function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj = se.obj)
    }
    ## Returning only the adjusted rand index results ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- Returning only the LISI results',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The computeLisi function finished.',
            color = 'white',
            verbose = verbose
            )
        return(all.lisi = all.lisi)
    }
}
