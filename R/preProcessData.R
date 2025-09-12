#' Pre-process RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function applies ....
#'
#' @param se.obj description
#' @param assay.name description
#' @param normalization description
#' @param regress.out.variables description
#' @param regress.out.rle.med description
#' @param apply.log description
#' @param pseudo.count description
#' @param check.se.obj description
#' @param remove.na description
#' @param verbose description
#' @importFrom SummarizedExperiment assays colData
#' @importFrom CMScaller CMScaller
#' @importFrom parallel mclapply
#' @export

preProcessData <- function(
        se.obj ,
        assay.name,
        normalization = 'CPM',
        regress.out.variables = NULL,
        regress.out.rle.med = FALSE,
        apply.log = TRUE,
        pseudo.count = 1,
        check.se.obj = TRUE,
        remove.na = 'assay',
        verbose = TRUE
        ){
    printColoredMessage(
        message =  '------------The preProcessData function starts:',
        color = 'white',
        verbose = verbose
        )

    # checking the SummarizedExperiment object ####
    if (isTRUE(check.se.obj)){
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = NULL,
            remove.na = remove.na,
            verbose = verbose
            )
    }
    # Applying library size normalization ####
    if (!is.null(normalization)){
        expr.data <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            pseudo.count = pseudo.count,
            apply.log = apply.log,
            check.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
        )
    }
    # Regressing out unwanted variables ####
    if (!is.null(regress.out.variables) & !is.null(normalization)){
        printColoredMessage(
            message = paste0(
                '- The ',
                paste0(regress.out.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- Note, we do not recommend regressing out the ',
                paste0(regress.out.variables, collapse = ' & '),
                ' if they are largely associated with the biological variables.'),
            color = 'red',
            verbose = verbose
            )
        expr.data <- t(expr.data)
        uv.variables.all <- paste('se.obj', regress.out.variables, sep = '$')
        expr.data <- lm(as.formula(paste(
            'expr.data',
            paste0(uv.variables.all, collapse = '+') ,
            sep = '~')))
        expr.data <- t(expr.data$residuals)
        colnames(expr.data) <- colnames(se.obj)
        row.names(expr.data) <- row.names(se.obj)
    }
    if (!is.null(regress.out.variables) & is.null(normalization)){
        if (isTRUE(apply.log)){
            expr.data <- applyLog(
                se.obj = se.obj,
                assay.names = assay.name,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose
                )[[assay.name]]
        }
        if (isFALSE(apply.log)){
            expr.data <- assay(x = se.obj, i = assay.name)
        }
        printColoredMessage(
            message = paste0(
                '- The',
                paste0(regress.out.variables, collapse = ' & '),
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- Note, we do not recommend regressing out the ',
                paste0(regress.out.variables, collapse = ' & '),
                ' if they are largely associated with the biological variables.'),
            color = 'red',
            verbose = verbose
            )
        expr.data <- t(expr.data)
        uv.variables.all <- paste('se.obj', regress.out.variables, sep = '$')
        expr.data <- lm(as.formula(paste(
            'expr.data',
            paste0(uv.variables.all, collapse = '+') ,
            sep = '~'
            )))
        expr.data <- t(expr.data$residuals)
        colnames(expr.data) <- colnames(se.obj)
        row.names(expr.data) <- row.names(se.obj)
    }
    if (is.null(regress.out.variables) & is.null(normalization)){
        if (isTRUE(apply.log)){
            expr.data <- applyLog(
                se.obj = se.obj,
                assay.names = assay.name,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose
                )[[assay.name]]
        }
        if (isFALSE(apply.log)){
            expr.data <- assay(x = se.obj, i = assay.name)
        }
    }
    if (isTRUE(regress.out.rle.med)){
        rle.med <- matrixStats::colMedians(expr.data - matrixStats::rowMedians(expr.data))
        expr.data <- lm(t(expr.data) ~ rle.med)
        expr.data <- t(expr.data$residuals)
        colnames(expr.data) <- colnames(se.obj)
        row.names(expr.data) <- row.names(se.obj)
    }
    printColoredMessage(
        message =  '------------The preProcessData function finished.',
        color = 'white',
        verbose = verbose
    )
    return(expr.data)
}










