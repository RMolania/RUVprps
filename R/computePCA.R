#' Performs principal component analysis (PCA)

#' @author Ramyar Molania

#' @description
#' This function uses singular value decomposition to perform principal component on the dataset(s) (assay(s) in a
#' SummarizedExperiment object. The function provides fast singular value decomposition using the BiocSingular R package.

#' @details
#' The PCs (in this context also called singular vectors) of the sample × transcript array of log counts are the linear
#' combinations of the transcript measurements having the largest, second largest, third largest, etc., variation,
#' standardized to be of unit length and orthogonal to the preceding components. Each will give a single value for
#' each sample.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or a vector of character strings for selecting the name(s) of the
#' assay(s) in the SummarizedExperiment object to compute PCA. The default is set to "all", which indicates all the
#' assays in the SummarizedExperiment object will be selected.
#' @param fast.pca Logical. Indicates whether to calculate a specific number of left singular vectors instead of the
#' full set of vectors to speed up the process. The default is set to `TRUE`.
#' @param nb.pcs Numeric. The number of first left singular vectors to be calculated for the fast PCA process.
#' The default is set to 10. If the `fast.pca = FALSE`, no need a value to for `nb.pcs`.
#' @param center Logical. Indicates whether to center the data before applying SVD. If center is `TRUE`, centering is
#' performed by subtracting the column means of the data from their corresponding columns. The default is set to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data before applying SVD. If scale is set to `TRUE`, scaling is
#' done by dividing the (centered) columns of the assays by their standard deviations,  if center is `TRUE`, and the root
#' mean square otherwise. The default is set to `FALSE`.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before computing the SVD. The
#' default is set to `TRUE`.
#' @param pseudo.count Numeric. A value used as a pseudo-count to be added to all measurements before applying the log
#' transformation. The default is set to 1.
#' @param svd.bsparam A BiocParallelParam object specifying how parallelization should be performed. The default is set
#' to `bsparam()`. See the `runSVD()` function from the BiocSingular R package for more details.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. The default is set to `TRUE`.
#' See the `checkSeObj()` function for more details.
#' @param remove.na Character. Specifies whether to remove NA or missing values from the data sets (assays). Options are
#' 'assays' and 'none'. The default is set to `assays`.
#' @param override.check Logical. When set to `TRUE`, the function checks whether PCA has already been computed for the
#' current parameters. If so, the PCA will not be recalculated. The default is set to `FALSE`.
#' @param save.se.obj Logical. Indicates whether to save the SVD results in the metadata of the SummarizedExperiment object
#' or to output the results as a list. The default is set to `TRUE`. The results can be found:
#' "se.obj->metadata->metric->AssayName->global.level->PCA$"
#' @param verbose Logical. If `TRUE`, displays messages for the different steps of the function.

#' @return A SummarizedExperiment object or a list containing the singular value decomposition results and the
#' percentage variation of each principal component (PC).

#' @importFrom SummarizedExperiment assay
#' @importFrom BiocSingular runSVD bsparam
#' @import ggplot2
#' @export

computePCA <- function(
        se.obj,
        assay.names = 'all',
        fast.pca = TRUE,
        nb.pcs = 10,
        center = TRUE,
        scale = FALSE,
        apply.log = TRUE,
        pseudo.count = 1,
        svd.bsparam = bsparam(),
        check.se.obj = TRUE,
        remove.na = 'assays',
        override.check = FALSE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computePCA function starts:',
                        color = 'white',
                        verbose = verbose)

    # Check to override or not ####
    if (isTRUE(override.check)){
        if(isTRUE(fast.pca)){
            method <- 'fast.svd'
        } else method <- 'ordinary.svd'
        override.check <- overrideCheck(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment.type = 'global.level',
            assessment = 'PCA',
            method = method,
            variable = 'general',
            file.name = 'data',
            verbose = verbose
        )
        if (is.logical(override.check)){
            compute.metric <- FALSE
        } else if (is.list(override.check)) {
            compute.metric <- TRUE
            assay.names <- override.check$selected.assays
        }
    } else if (isFALSE(override.check)) compute.metric <- TRUE

    if (isTRUE(compute.metric)){
        # Checking the function inputs ####
        if (is.null(assay.names) | is.logical(assay.names)) {
            stop('The "assay.names" cannot be empty or lgical.')
        }
        if (!is.logical(fast.pca)){
            stop('The "fast.pca" must be logical.')
        }
        if (isTRUE(fast.pca)) {
            if (is.null(nb.pcs)){
                stop('To perform fast PCA, the "nb.pcs" must be specified.')
            }
            if (nb.pcs < 0 & nb.pcs == 0) {
                stop('To perform fast PCA, the "nb.pcs" must be a postive numeric value.')
            }
        }
        if (!is.logical(center)){
            stop('The "center" must be logical.')
        }
        if (!is.logical(scale)){
            stop('The "scale" must be logical.')
        }
        if (isTRUE(scale)) {
            printColoredMessage(
                message = 'Note: highly recommend not to scale the data before computing the PCA.',
                color = 'red',
                verbose = verbose)
        }
        if (!is.logical(apply.log)){
            stop('The "apply.log" must be logical.')
        }
        if (isTRUE(apply.log)){
            if (pseudo.count < 0){
                stop('The value of "pseudo.count" cannot be negative.')
            }
            if (is.null(pseudo.count)){
                stop('A value for the "pseudo.count" must be specified.')
            }
        }
        if (!is.logical(check.se.obj)){
            stop('The "check.se.obj" must be logical.')
        }
        if (!remove.na %in% c('assays','none')){
            stop('The "remove.na" must be on of the "assays" or "none"')
        }
        if (!is.logical(override.check)){
            stop('The "override.check" must be logical.')
        }
        if (!is.logical(save.se.obj)){
            stop('The "save.se.obj" must be logical.')
        }
        if (!is.logical(verbose)){
            stop('The "verbose" must be logical.')
        }

        # Checking the assays ####
        if (length(assay.names) == 1 && assay.names == 'all') {
            assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
        } else  assay.names <- factor(x = assay.names , levels = assay.names)
        if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
        }
        # Assessing the SummarizedExperiment object ####
        if (isTRUE(check.se.obj)) {
            se.obj <- checkSeObj(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = NULL,
                remove.na = remove.na,
                verbose = verbose)
        }
        # Data transformation ####
        if (isTRUE(apply.log)){
            printColoredMessage(
                message = '-- Applying log transformation on all the specified assay(s):',
                color = 'magenta',
                verbose = verbose
            )
            all.assays <- applyLog(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                pseudo.count = pseudo.count,
                verbose = verbose
            )
        }
        if (isFALSE(apply.log)){
            printColoredMessage(
                message = '-- The specified assay(s) will be used for PCA without applying log transformation.',
                color = 'blue',
                verbose = verbose
            )
            all.assays <- lapply(
                levels(assay.names),
                function(x) assay(x = se.obj, i = x))
            names(all.assays) <- levels(assay.names)
        }

        # Compute SVD ####
        printColoredMessage(
            message = '-- Computing PCA usibg singular value decomposition (SVD):',
            color = 'magenta',
            verbose = verbose
        )
        ## compute fast SVD ####
        if (isTRUE(fast.pca)) {
            printColoredMessage(
                message = paste0(
                    '- Computing "fast" singular value decomposition (SVD) with scale = ',
                    scale,
                    ' and center = ',
                    center, '.'),
                color = 'orange',
                verbose = verbose
                )
            if (is.null(svd.bsparam))
                svd.bsparam <- bsparam()
            all.sv.decomposition <- lapply(
                levels(assay.names),
                function(x) {
                    printColoredMessage(
                        message = paste0('- Performing fast SVD on the "', x , '" data.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    sv.dec <- BiocSingular::runSVD(
                        x = t(all.assays[[x]]),
                        k = nb.pcs,
                        BSPARAM = svd.bsparam,
                        center = center,
                        scale = scale
                        )
                    rownames(sv.dec$u) <- colnames(se.obj)
                    rownames(sv.dec$v) <- row.names(se.obj)
                    percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
                    percentage <- sapply(
                        seq_along(percentage),
                        function(i) round(percentage [i], 1))
                    return(list(svd = sv.dec, percentage.variation = percentage))
                })
            printColoredMessage(
                message = paste0(
                    '- Note: in the fast svd analysis, the percentage of variation of PCs will be ',
                    'computed proportional to the highest selected number of PCs (left singular vectors), not on all the PCs.'),
                color = 'red',
                verbose = verbose)
            names(all.sv.decomposition) <- levels(assay.names)
        }
        ## compute ordinary SVD ####
        if (isFALSE(fast.pca)) {
            printColoredMessage(
                message = paste0(
                    '- Performing singular value decomposition (SVD) with scale = ',
                    scale,
                    ' and center = ',
                    center,
                    '.'),
                color = 'orange',
                verbose = verbose
            )
            all.sv.decomposition <- lapply(
                levels(assay.names),
                function(x) {
                    printColoredMessage(
                        message = paste0(
                            '- Performing SVD on the "',
                            x ,
                            '" data.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    sv.dec <- svd(scale(
                        x = t(all.assays[[x]]),
                        center = center,
                        scale = scale)
                    )
                    rownames(sv.dec$u) <- colnames(se.obj)
                    rownames(sv.dec$v) <- row.names(se.obj)
                    percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
                    percentage <- sapply(
                        seq_along(percentage),
                        function(i) round(percentage [i], 1))
                    return(list(svd = sv.dec, percentage.variation = percentage))
                })
            names(all.sv.decomposition) <- levels(assay.names)
        }

        # Save all the results ####
        printColoredMessage(
            message = '-- Saving the SVD results:',
            color = 'magenta',
            verbose = verbose
            )
        ## add the results to the SummarizedExperiment object ####
        if (isTRUE(save.se.obj)) {
            printColoredMessage(
                message = '- Saving all the SVD results to the "metadata" of the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            if(isTRUE(fast.pca)){
                method <- 'fast.svd'
            } else method <- 'ordinary.svd'
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = levels(assay.names),
                assessment.type = 'global.level',
                assessment = 'PCA',
                method = method,
                file.name = 'data',
                variables = 'general',
                results.data = all.sv.decomposition
            )
            printColoredMessage(
                message = paste0(
                    '- The SVD results of individual assay (s) are saved to the',
                    ' "se.obj@metadata$metric$AssayName$global.level$PCA$',
                    method,
                    '$data" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
                )
            printColoredMessage(message = '------------The computePCA function finished.',
                                color = 'white',
                                verbose = verbose)
            return(se.obj)
        }

        if (isFALSE(save.se.obj)) {
            ## return a list ####
            printColoredMessage(
                message = '- The SVD results of individual assays are outputed as a list.',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(message = '------------The computePCA function finished.',
                                color = 'white',
                                verbose = verbose)
            return(all.sv.decompositions = all.sv.decomposition)
        }
    } else {
        printColoredMessage(message = '------------The computePCA function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)

    }
}
