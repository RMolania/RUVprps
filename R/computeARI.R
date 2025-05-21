#' Computes Adjusted Rand Index (ARI).
#'
#' @author Ramyar Molania
#'
#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @description
#' This functions computes the adjusted rand index(ARI) for given a categorical variable using the first PCs of the specified
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
#' @param variable character. The name of the column containing the categorical variable in the SummarizedExperiment object.
#' This variable can represent either a biological or an unwanted factor.
#' @param clustering.method character. A character that indicates which clustering method to be applied on the principal
#' components to calculate the ARI. Options are `mclust` or `hclust`. The default is se to the `hclust` method.
#' @param hclust.method Character. Agglomeration method to use when `clustering.method` is set to `hclust`. Options include:
#' `ward.D`, `ward.D2`, `single`, `complete`, `average` (= UPGMA), `mcquitty` (= WPGMA), `median` (= WPGMC), or `centroid`
#' (= UPGMC). See the `hclust` function in the **stats** package for more details.
#' @param hclust.dist.measure Character. A character specifying which distance measure to be used in the `dist` function
#' when applying hierarchical clustering. Options are: `euclidean`, `maximum`, `manhattan`, `canberra`, `binary`, or
#' `minkowski`. See the `dist` function in the **stats** package for more details.
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
#' @importFrom mclust mclustBIC Mclust adjustedRandIndex
#' @importFrom SummarizedExperiment assays assay
#' @importFrom stats cutree hclust dist
#' @import ggplot2
#' @export

computeARI <- function(
        se.obj,
        assay.names = 'all',
        variable,
        clustering.method = 'hclust',
        hclust.method = 'complete',
        hclust.dist.measure = 'euclidian',
        fast.pca = TRUE,
        nb.pcs = 3,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeARI function starts:',
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
    if (class(se.obj@colData[, variable]) %in% c('numeric', 'integer')) {
        stop('The "variable" must be a categorical varible.')
    }
    if (is.null(nb.pcs)) {
        stop('The "nb.pcs" must be specified.')
    }
    if (!is.numeric(nb.pcs) | nb.pcs < 0){
        stop('The "nb.pcs" must be a positive numerical value.')
    }
    if (length(clustering.method) > 1){
        stop('The "clustering.method" must be only one of the "mclust" or "hclust".')
    }
    if (is.null(clustering.method)){
        stop('The "clustering.method" must be provided.')
    }
    if (!clustering.method %in% c('mclust', 'hclust')){
        stop('The "clustering.method" method must be one of the "mclust" or "hclust".')
    }
    if (clustering.method == 'hclust'){
        if (is.null(hclust.dist.measure)){
            stop('The "hclust.dist.measure" cannot be empty when the "clustering.method = hclust".')
        } else if (!hclust.dist.measure %in% c('euclidian',
                                               'maximum',
                                               'manhattan',
                                               'canberra',
                                               'binary',
                                               'minkowski')) {
            stop('The "hclust.dist.measure" should be one of the:"euclidean","maximum","manhattan","canberra","binary" or "minkowski".')
        }
        if (is.null(hclust.method)){
            stop('The "hclust.method" cannot be when the "clustering.method = hclust".')
        } else if (!hclust.method %in% c('complete',
                                         'ward.D',
                                         'ward.D2',
                                         'single',
                                         'average',
                                         'mcquitty',
                                         'median',
                                         'centroid')) {
            stop('The hclust.method should be one of the:"complete","ward.D","ward.D2","single","average", "mcquitty", "median" or "centroid".')
        }
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

    # Computing the ARI on all the data set(s) ####
    printColoredMessage(
        message = paste0(
            '-- Computing adjusted rand index (ARI) using the first ',
            nb.pcs,
            ' PCS for the ',
            variable,
            ' variable.') ,
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

    ## Computing ARI for all assay(s) ####
    all.ari <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0(
                    '- Computing the ARI for the "',
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
                        'The number of PCs of the assay',
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
            if (clustering.method == 'mclust'){
                printColoredMessage(
                    message = '- Clusteingr the PCs using the mclust function.',
                    color = 'blue',
                    verbose = verbose
                    )
                bic <- mclustBIC(data = pca.data)
                mod <- Mclust(
                    data = pca.data,
                    x = bic,
                    G = length(unique(se.obj@colData[[variable]]))
                    )
                printColoredMessage(
                    message = '- Calculating the adjusted rand index.',
                    color = 'blue',
                    verbose = verbose
                    )
                # Applying ARI ####
                ari <- adjustedRandIndex(mod$classification, se.obj@colData[, variable])
            }
            ### Clustering the data using the hclust method ####
            if (clustering.method == 'hclust'){
                printColoredMessage(
                    message = '- Clustering the PCs using the hclust function.',
                    color = 'blue',
                    verbose = verbose
                    )
                clusters <- cutree(
                    tree = hclust(d = dist(x = pca.data, method = hclust.dist.measure), method = hclust.method),
                    k = length(unique(se.obj@colData[, variable]))
                    )
                printColoredMessage(
                    message = '- Calculating the adjusted rand index.',
                    color = 'blue',
                    verbose = verbose
                    )
                #### Applying ARI ####
                ari <- adjustedRandIndex(clusters, se.obj@colData[[variable]])
            }
            return(ari)
        })
    names(all.ari) <- levels(assay.names)

    # Saving the results ####
    printColoredMessage(
        message = '-- Saving all the ARI results:',
        color = 'magenta',
        verbose = verbose
        )
    ## Adding results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving the ARIs of each data set(s) to the metadata in the SummarizedExperiment object:',
            color = 'blue',
            verbose = verbose
            )
        if (clustering.method == 'mclust'){
            method <- 'mclust'
        } else method <- paste0('hclust.', hclust.method, '.', hclust.dist.measure)
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'global.level',
            assessment = 'ARI',
            method = method,
            variables = variable,
            file.name = 'ari',
            results.data = all.ari
            )
        printColoredMessage(
            message = paste0(
                '- The ARI of induvial assay(s) is saved to the .',
                ' ".se.obj@metadata$metric$RawCount$ARI" in the SummarizedExperiment objec.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The computeARI function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)
    }
    ## Returning only the adjusted rand index results ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- Returning only the adjusted rand index results',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The computeARI function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.ari = all.ari)
    }
}
