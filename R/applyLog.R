#' Applies log2 with a pseudo count on data sets in SummarizedExperiment object.

#' @author Ramyar Molania

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character or character vectos specifying the name(s) of the data (assays) in the
#' `SummarizedExperiment` object to be selected. These assays will be log2-transformed with a pseudo count. The default
#' is `all`, which indicates that all assays in the `SummarizedExperiment` object will be selected.
#' @param pseudo.count Numeric. A pseudo count value to be added to all measurements in the selected assay(s) before
#' applying the log2 transformation, to avoid `-Inf` values for zero measurements. The default is 1.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default it is
#'  set to `TRUE`.
#' @param remove.na Character. A Character. that indicates whether to remove NA or missing values from the data sets or
#' not. The options are `assays` or `none`. The default is set to `assays`.  Refer to the `checkSeObj()` function for more
#' details.
#' @param verbose Logical. If it is set to `TRUE`, displays messages describing the steps of the function.

#' @return The function returns a log2 transformed of all specified data sets as a list object.

applyLog <- function(
        se.obj,
        assay.names = 'all',
        pseudo.count = 1,
        assess.se.obj = TRUE,
        remove.na = 'assays',
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The applyLog starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking the function inputs ####
    if (is.logical(assay.names) | is.null(assay.names)){
        stop('The "assay.names" cannot be NULL or logical.')
    }
    if (!is.numeric(pseudo.count)){
        stop('The "pseudo.count" must be a postive numeric value.')
    }
    if (pseudo.count < 0){
        stop('The "pseudo.count" must be a postive numeric value.')
    }
    if (!is.logical(assess.se.obj)){
        stop('The "assess.se.obj" must be logical.')
    }
    if (!remove.na %in% c('assays', 'none')){
        stop('The "remove.na" must be on of the "assays" or "none".')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical.')
    }
    # Checking SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose
            )
    }
    # Checking the data sets names ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Applying log transformation ####
    all.assays.loged <- lapply(
        levels(assay.names),
        function(x){
            if (!is.null(pseudo.count)) {
                printColoredMessage(
                    message = paste0(
                        '- Applying log2 on the "',
                        x,
                        '" + ',
                        pseudo.count,
                        ' (pseudo.count) data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                expr <- log2(assay(x = se.obj, i = x) + pseudo.count)
            } else if (is.null(pseudo.count)){
                printColoredMessage(
                    message = paste0(
                        '- Applying log2 on the "',
                        x,
                        '" data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                expr <- log2(assay(x = se.obj, i = x))
            }
            return(expr)
        })
    names(all.assays.loged) <- levels(assay.names)
    printColoredMessage(
        message = '------------The applyLog finished',
        color = 'white',
        verbose = verbose
        )
    return(all.assays.loged)
}


