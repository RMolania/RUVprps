#' Performs several normalization methods for RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function provides several normalization methods for RNA-seq data, including CPM, TMM,  upper-quartile, median,
#' full quantile (nonlinear), and VST. It utilizes the **EDASeq**, **DESeq2**, and **edgeR** R packages for normalization.
#' Note that the default parameters of these underlying functions are used.


#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character that specifies the name of the data (assay) in the SummarizedExperiment object.
#' The selected data should contain raw count data.
#' @param method Character. A character that indicates the normalization method to apply to the assay specified by `assay.name`.
#' Options include:
#' - `'CPM'`: Counts Per Million, from the **edgeR** package.
#' - `'TMM'`: Trimmed Mean of M-values, from the **edgeR** package.
#' - `'upper'`: A scaling normalization that adjusts the upper quartile of each lane, from the **EDASeq** package.
#' - `'median'`: A scaling normalization that adjusts the median of each lane, from the **EDASeq** package.
#' - `'full'`: Full quantile normalization (nonlinear), from the **EDASeq** package.
#' - `'VST'`: Variance Stabilizing Transformation, from the **DESeq2** package.
#' The default is `'CPM'`.
#' @param apply.log Logical. Indicates whether to apply a log transformation to the data after normalization.
#' The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value as pseudo count to be added to all measurements after normalization and
#' before log transformation to avoid `-Inf` for zero values. The default is set to `1`.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not.
#' The default is set to `TRUE`. We refer to the `checkSeObj()` function for more details.
#' @param remove.na Character. A character that indicates whether to remove NA or missing values from the data sets or
#' not. The options are `'assays'` or `'none'`. The default is set to `'assays'`. Refer to the `checkSeObj()` function
#' for more details.
#' @param new.assay.name Character. Specifies the name of the newly normalized assay. If set to `NULL` (default),
#' the function will automatically generate a name using the format: paste0(method, ' on ', assay.name).
#' @param save.se.obj Logical. Indicates whether to save the normalized data as new data (assay) within the
#' `SummarizedExperiment` object. If set to `FALSE`, the normalized data will be returned as a matrix.
#' The default is set to `TRUE`. Then, each normalized data will be added a new data (assay) to `SummarizedExperiment`
#' object.
#' @param verbose Logical. Indicates whether to display messages and output during function execution.
#' The default is to `TRUE`.

#' @return A SummarizedExperiment object containing the newly normalized data, or a normalized expression matrix
#' if `save.se.obj` is set to `FALSE`.

#' @importFrom EDASeq betweenLaneNormalization
#' @importFrom SummarizedExperiment assay
#' @importFrom edgeR cpm normLibSizes
#' @importFrom DESeq2 vst
#' @export

applyOtherNormalizations <- function(
        se.obj,
        assay.name,
        method = 'CPM',
        apply.log = TRUE,
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'assays',
        new.assay.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The applyOtherNormalizations function starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking the function inputs ####
    if (is.logical(assay.name) | is.null(assay.name)){
        stop ('The "assay.name" cannot be NULL or logical.')
    }
    if (length(assay.name) > 1 | assay.name == 'all') {
        stop('The "assay.name" must be a name of a data (assay) in the SummarizedExperiment object.')
    }
    if (is.logical(method) | is.null(method) | is.numeric(method) ){
        stop ('The "method" cannot be NULL, logical and numeric.')
    }
    if (!method %in% c('CPM', 'TMM', 'upper', 'median', 'full', 'VST')){
        stop ('The method must be one of the "CPM", "TMM", "upper", "median", "full", or "VST".')
    }
    if (is.logical(apply.log)){
        stop('The "apply.log" must be logical.')
    }
    if (isTRUE(apply.log)){
        if (!is.numeric(pseudo.count)){
            stop('The "pseudo.count" must be a postive numeric value.')
        }
        if (pseudo.count < 0){
            stop('The value for "pseudo.count" should be postive.')
        }
    }
    if (!is.logical(assess.se.obj)){
        stop('The "assess.se.obj" must be logical.')
    }
    if (is.logical(remove.na) | is.null(remove.na)){
        stop('The "remove.na" cannot be NULL or logical.')
    }
    if (!remove.na %in% c('assays', 'none')){
        stop('The "remove.na" must be one of the "assays" or "none".')
    }
    if (is.logical(new.assay.name) | is.numeric(new.assay.name)){
        stop('The "new.assay.name" cannot be logical or numeric.')
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical.')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical.')
    }
    # Checking SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose
            )
    }
    # Applying normalization ####
    printColoredMessage(
        message = '-- Normalizing the data for mainly library size:',
        color = 'magenta',
        verbose = verbose
        )
    if (isTRUE(apply.log)){
        if (is.null(pseudo.count)) pseudo.count <- 0
    }
    ## CPM normalization ####
    if (method == 'CPM' & isTRUE(apply.log)) {
        printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method , and then performing log2 transformation.'),
                color = 'blue',
                verbose = verbose
                )
        norm.data <- edgeR::cpm(y = assay(se.obj, i = assay.name))
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method == "CPM" & isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'Applying the ',
                method,
                ' method.'),
            color = 'blue',
            verbose = verbose
            )
        norm.data <- cpm(y = assay(se.obj, i = assay.name))
        norm.data
    }
    ## TMM method ####
    if (method == "TMM" & isTRUE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'Applying the ',
                method,' method , and then performing log2 transformation.'),
            color = 'blue',
            verbose = verbose
            )
        norm.data <- edgeR::normLibSizes(object = assay(x = se.obj, i = assay.name))
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method == "TMM" & isFALSE(apply.log)) {
        printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method.'),
                color = 'blue',
                verbose = verbose
                )
        norm.data <- normLibSizes(object = assay(x = se.obj, i = assay.name))
        norm.data

    }
    ## Median, upper or full quartile Methods ####
    if (method %in% c("median", "upper", "full") && isTRUE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'Applying the ',
                method,
                ' method and then performing log2 transformation.'),
            color = 'blue',
            verbose = verbose
            )
        norm.data <- EDASeq::betweenLaneNormalization(
            x = assay(x = se.obj, i = assay.name),
            which = method
            )
        norm.data <- log2(norm.data + pseudo.count)
        norm.data
    } else if (method %in% c("median", "upper", "full") && isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'Applying the ',
                method,
                ' method.'),
            color = 'blue',
            verbose = verbose
            )
        norm.data <- EDASeq::betweenLaneNormalization(
            x = assay(x = se.obj, i = assay.name),
            which = method
            )
        norm.data
    }
    ## VST method ####
    if (method == 'VST') {
        printColoredMessage(
                message = paste0(
                    'Applying the ',
                    method,
                    ' method.'),
                color = 'blue',
                verbose = verbose
                )
        norm.data <- DESeq2::vst(
            object = assay(x = se.obj, i = assay.name),
            blind = TRUE,
            nsub = 1000,
            fitType = "parametric"
            )
        norm.data
    }
    # Save the data ####
    printColoredMessage(
        message = '-- Saving the normalized data:',
        color = 'magenta',
        verbose = verbose
        )
    if (is.null(new.assay.name)){
        new.assay.name <- paste0(method, ' on ', assay.name)
    }

    ## save the data into  SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        if (!new.assay.name %in% names(assays(se.obj))) {
            se.obj@assays@data[[new.assay.name]] <- norm.data
        }
        printColoredMessage(
            message = paste0(
                'The normalized data ',
                new.assay.name,
                ' is saved to SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The applyOtherNormalizations function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ## out the data ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = paste0(
                'The normalized data ',
                new.assay.name,
                ' is outputed as data marix.'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = '------------The applyOtherNormalizations function finished.',
            color = 'white',
            verbose = verbose
        )
        return(norm.data)
    }
}
