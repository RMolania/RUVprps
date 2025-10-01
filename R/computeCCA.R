#' Perform canonical correlation coordinates in RNA-seq data
#'
#' @description
#' The function performs canonical correlation analysis between each pair of groups in order to obtain a CCA dimensional
#' space. This space can be used by other functions in RUVprps to identify similar samples and create PRPS.
#'
#' @author Ramyar Molania
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character specifying a data (assay) name within the SummarizedExperiment object.
#' @param variable Character. A character specifying a variable name in the sample annotation of the SummarizedExperiment
#' object. The variable must be a categorical variable with at least two levels.
#' @param reference.group Character. The name of a group provided in the `variable` to be used as the reference group to
#' perform all possible pair-wise CCA with the reference. If it is set to `NULL`, all possible pairs of groups will be
#' selected. The default is set to `NULL`.
#' @param nb.cca Numeric. A numeric value specifying the number of CCA to be calculated. The default is set to 10.
#' @param hvg A logical vector or character vector of highly variable genes to be used to calculate CCA. The default
#' is set to `NULL`, in which case all genes will be used.
#' @param scale Logical. This specifies whether to scale the data before applying CCA. The default is set to `TRUE`.
#' @param normalization Character. A character that specifies which normalization should be performed before applying CCA.
#' The default is set to `CPM`. Refer to the `applyOtherNormalization` function for more details.
#' @param regress.out.variables Character. A character or a vector of characters that indicate the column name(s) in the sample
#' annotation in the SummarizedExperiment object. These variables will be regressed out from the data before performing CCA.
#' The default is set to `NULL`, indicating that regression will not be applied.
#' @param regress.out.rle.med Logical. This indicates whether to regress out the median of RLE data after applying either
#' normalization or regression. The default is set to `FALSE`.
#' @param cosine.norm Logical. This indicates whether to apply a cosine normalization or not. The default is set to `FALSE`.
#' This normalization will be applied after all other normalization and regression steps.
#' @param samples.to.use Logical. A logical vector indicating which samples to use for computing CCA. The default is
#' set to `all`, meaning all samples will be selected.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data. The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value as a pseudo count to be added to all measurements before log transformation.
#' The default is set to 1.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. If `TRUE`, the `checkSeObj()`
#' function will be applied inside the function. The default is set to `TRUE`.
#' @param remove.na Character. A character indicating whether to remove NA or missing values from either the 'assays',
#' the `sample.annotation`, `both`, or `none`. If `assays` is selected, the genes that contain NA or missing values will
#' be excluded. If `sample.annotation` is selected, the samples that contain NA or missing values for any `bio.variables`
#' and `uv.variables` will be excluded. The default is set to `both`.
#' @param cca.set.name Character. A character specifying the name of the output file to be saved in the metadata
#' of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0('Cca:', nb.cca, '|', 'HVG:', length(hvg), '|Scale:', scale, '|Norm:', normalization)`.
#' @param save.se.obj Logical. Indicates whether to save the results in the metadata of the SummarizedExperiment object
#' or to output the result as a list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, shows the messages of different steps of the function.
#'
#' @importFrom irlba irlba

computeCCA <- function(
        se.obj,
        assay.name,
        variable,
        reference.group = NULL,
        nb.cca = 10,
        hvg = NULL,
        scale = TRUE,
        normalization = 'CPM',
        regress.out.variables = NULL,
        regress.out.rle.med = FALSE,
        cosine.norm = FALSE,
        samples.to.use = 'all',
        apply.log = TRUE,
        pseudo.count = 1,
        check.se.obj = TRUE,
        remove.na = 'none',
        cca.set.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The computeCCA function starts',
        color = 'white',
        verbose = verbose
        )
    # Specifying the number of sample to use ####
    if (is.logical(samples.to.use)){
        if (length(samples.to.use) != ncol(se.obj)){
            stop('The "samples.to.use" must be the same length as the sample numbers in the SummarizedExperiment object.')
        }
    }
    if (is.logical(samples.to.use)){
        initial.se.obj <- se.obj
        se.obj <- se.obj[ , samples.to.use]
    }

    # Specifying highly variable genes ####
    if (is.null(hvg)){
        hvg = row.names(se.obj)
    }
    # Specifying reference data ####
    variable.names <- as.character(unique(se.obj[[variable]]))
    if (is.null(reference.group)){
        all.combinations <- combn(x = variable.names, m = 2)
    }
    if (!is.null(reference.group)){
        if (reference.group != 'random'){
            variable.names <- c(
                variable.names[variable.names %in% reference.group],
                variable.names[!variable.names %in% reference.group]
            )
            all.combinations <- combn(x = variable.names, m = 2)
            all.combinations <- all.combinations[ , 1:c(length(variable.names) -1) ]
        }
        if (reference.group == 'random'){
            reference.group <- variable.names[1]
            variable.names <- c(
                variable.names[variable.names %in% reference.group],
                variable.names[!variable.names %in% reference.group]
            )
            all.combinations <- combn(x = variable.names, m = 2)
            all.combinations <- all.combinations[ , 1:c(length(variable.names) -1) ]
        }
    }
    # Applying data normalization
    all.norm.data <- lapply(
        variable.names,
        function(x){
            index.samples <- se.obj[[variable]] == x
            norm.data <- preProcessData(
                se.obj = se.obj[ , index.samples],
                assay.name = assay.name,
                normalization = normalization,
                regress.out.variables = regress.out.variables,
                regress.out.rle.med = regress.out.rle.med,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose
                )
            norm.data
        })
    names(all.norm.data) <- variable.names

    # Applying cosine normalization ####
    if (isTRUE(cosine.norm)){
        all.norm.data <- lapply(
            unique(se.obj[[variable]]),
            function(x){
                cosineNorm(x = all.norm.data[[x]], mode = 'matrix')
            })
        names(all.norm.data) <- unique(se.obj[[variable]])
    }
    # Computing all CCA ####
    nb.cca.initial <- nb.cca
    all.cca <- lapply(
        1:ncol(all.combinations),
        function(x){
            nb.cca <- nb.cca.initial
            ## Getting the data ####
            data.1 <- all.norm.data[[all.combinations[1 , x]]][hvg , ]
            data.2 <- all.norm.data[[all.combinations[2 , x]]][hvg , ]

            ## Checking gene level expression ###
            remove.genes <- c(
                unname((which(rowSums(data.1) == 0))) ,
                unname((which(rowSums(data.2) == 0)))
                )
            if (length(remove.genes) > 0){
                data.1 <- data.1[-remove.genes , ]
                data.2 <- data.2[-remove.genes , ]
            }

            ## Scaling the data ####
            if (isTRUE(scale)){
                data.1 <- t(scale(t(data.1)))
                data.2 <-  t(scale(t(data.2)))
            }
            ## Scaling the data ####
            data.1 <- scale(data.1)
            data.2 <- scale(data.2)

            ## Computing CAA
            cross.pro.data <- crossprod(x = data.1, y = data.2)
            if (nb.cca >= min(nrow(cross.pro.data), ncol(cross.pro.data))){
                printColoredMessage(
                    message = '- Note, the "nb.cca" is capped',
                    color = 'red',
                    verbose = verbose
                    )
                nb.cca <- min(nrow(cross.pro.data), ncol(cross.pro.data)) - 1
                if (nb.cca < 1){
                    stop('TTTT')
                }
            }
            cca.svd <- irlba(A = cross.pro.data, nv = nb.cca)
            cca.data <- rbind(cca.svd$u, cca.svd$v)
            colnames(x = cca.data) <- paste0("CC", 1:nb.cca)
            rownames(cca.data) <- c(colnames(data.1), colnames(data.2))
            cca.data <- apply(
                X = cca.data,
                MARGIN = 2,
                FUN = function(x) {
                    if (sign(x[1]) == -1) {
                        x <- x * -1
                    }
                    return(x)
                })
            return(cca.data)
        })
    all.combinations <- as.matrix(all.combinations)
    names(all.cca) <- sapply(
        1:ncol(all.combinations),
        function(x) paste0(all.combinations[, x] , collapse = '_'))

    if (is.null(cca.set.name)){
        cca.set.name <- paste0(
            'Cca:',
            nb.cca,
            '|',
            'HVG:',
            length(hvg),
            '|Scale:',
            scale,
            '|Norm:',
            normalization
            )
    }

    # Specifying the number of sample to use ####
    if (is.logical(samples.to.use)){
        se.obj <- initial.se.obj
    }
    #### Saving all the results ####
    if (isTRUE(save.se.obj)){
        if (length(se.obj@metadata$CCA) == 0 ) {
            se.obj@metadata[['CCA']] <- list()
        }
        if (!assay.name %in% names(se.obj@metadata[['CCA']])){
            se.obj@metadata[['CCA']][[assay.name]] <- list()
        }
        if (!variable %in% names(se.obj@metadata[['CCA']][[assay.name]] )){
            se.obj@metadata[['CCA']][[assay.name]][[variable]] <- list()
        }
        if (!cca.set.name %in% names(se.obj@metadata[['CCA']][[assay.name]][[variable]])){
            se.obj@metadata[['CCA']][[assay.name]][[variable]][[cca.set.name]] <- list()
        }
        se.obj@metadata[['CCA']][[assay.name]][[variable]][[cca.set.name]] <- all.cca

        printColoredMessage(
            message = '- The hvg are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = '------------The findNcgByTwoWayAnova function finished.',
            color = 'white',
            verbose = verbose
        )
        return(se.obj)
    }
    ### Export results as logical vector ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = '-- All the CCA results are outputed as a list.',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The computeCCA function finished.',
            color = 'white',
            verbose = verbose
            )
        return(all.cca)
    }
}



