#' Regress out variables from RNA-seq data
#'
#' @author Ramyar Molania
#'
#' @description
#' This functions uses a linear model to regresses out any specified variables from data in SummarizedExperiment object.
#'
#' @param se.obj description
#' @param assay.name description
#' @param variables description
#' @param apply.log description
#' @param pseudo.count TTTT
#' @param check.se.obj description
#' @param remove.na description
#' @param verbose description
#'
#' @importFrom stats model.matrix residuals
#' @importFrom limma lmFit
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
