#' Regress out variables from RNA-seq data
#'
#' @author Ramyar Molania
#'
#' @description
#' This functions uses a linear model to regresses out any specified variables from data in SummarizedExperiment object.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character specifying the name of the data (assay) in the SummarizedExperiment object
#' to be selected.
#' @param variables Character. A character string or vector specifying column names in the sample annotation of the
#' `SummarizedExperiment` object. These columns can be categorical, continuous, or a combination of both. These variable
#' will be regressed out from the data.
#' @param apply.log Logical. If `TRUE`, a log2 transformation + `pseudo.count` will be applied to the data before applying
#' regression analysis.
#' @param pseudo.count Numeric. A numerical value to be addedd to all measurements before applying a log transformation.
#' The default is set to 1
#' @param check.se.obj Logical. Indicates whether to validate the SummarizedExperiment object before analysis. The default
#' is set to `TRUE`.
#' @param remove.na Character. Indicates whether to remove NA or missing values from the data sets. Options are: `assays`
#' or `none`. Default is `assays`. Refer to the `checkSeObj()` function for more details.
#' @param verbose Logical. Indicates whether to display messages during function execution. The default is set `TRUE`.
#'
#' @importFrom stats model.matrix residuals
#' @importFrom limma lmFit
#'
#' @export

regressOutVariables <- function(
        se.obj,
        assay.name,
        variables,
        apply.log = TRUE,
        pseudo.count = 1,
        check.se.obj = TRUE,
        remove.na = 'none',
        verbose = TRUE
        ){
    if (isTRUE(apply.log)){
        expr.data <- applyLog(
            se.obj = se.obj,
            assay.names = assay.name,
            pseudo.count = pseudo.count,
            check.se.obj = check.se.obj,
            remove.na = remove.na,
            verbose = verbose
            )[[assay.name]]
    } else expr.data <- assay(x = se.obj, i = assay.name)
    sample.info <- as.data.frame(colData(se.obj))
    adjustment.form <- as.formula(paste0('~', paste0(variables, collapse = '+')))
    lm.fit <- limma::lmFit(
        object = expr.data,
        design =  model.matrix(adjustment.form, sample.info)
    )
    res.data <- residuals(object = lm.fit, expr.data)
    return(res.data)
}
