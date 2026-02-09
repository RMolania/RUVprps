#' Compute gene-level variance composition across variables in RNA-seq.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function computes the proportion of variance in gene expression explained by different biological
#' and technical variables in RNA-seq data, using variance partitioning methods. It can adjust for
#' unwanted variation, apply linear mixed models, and optionally log-transform data.
#'
#' @param se.obj A `SummarizedExperiment` object containing gene expression data and sample annotations.
#' @param assay.names Character. A character string or vector of character strings specifying the name(s) of the assay(s)
#' in the `SummarizedExperiment` object to compute gene-level variance composition. By default, set to `"all"`, which means all
#' assays in the object will be selected.
#' @param form Formula. A formula specifying the model used for variance partitioning (e.g., `~ Batch + Condition + Age`).
#' @param adjust.data Logical. If `TRUE`, adjusts the data before variance partitioning using the model defined in
#' `adjustment.form`. The default is `FALSE`.
#' @param adjustment.form Formula. A formula specifying the model to adjust the data (e.g., unwanted variation factors).
#' Used only if `adjust.data = TRUE`.
#' @param adjustment.method Character. The method used to adjust for unwanted variation (e.g., `"lm"`, `"combat"`).
#' @param samples.to.use Character or NULL. A vector of sample identifiers to include in the analysis. If `NULL`,
#' all samples are used.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing variance
#' partitioning. The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value added as a pseudo-count to avoid log-transformation of zeros.
#' The default is 1.
#' @param nb.cores Numeric. Number of CPU cores to use for parallel computation. If `NULL`, maximum available cores - 1
#' will be used.
#' @param plot.output Logical. If `TRUE`, generates plots summarizing variance composition across variables. The default
#' is `TRUE`.
#' @param output.name Character. A label to assign to the results when stored in the metadata of the `SummarizedExperiment`
#' object.
#' @param check.se.obj Logical. Indicates whether to validate the `SummarizedExperiment` object before running the
#' analysis. The default is `TRUE`.
#' @param remove.na Character. Specifies how to handle missing values. Options are `"assays"`, `"sample.annotation"`,
#' `"both"`, or `"none"`. The default is `"both"`.
#' @param save.se.obj Logical. If `TRUE`, saves the results (variance components, F-statistics, p-values) in the metadata
#' of the `SummarizedExperiment` object. If `FALSE`, returns them as a list. The default is `TRUE`.
#' @param verbose Logical. If `TRUE`, displays progress messages during execution. The default is `TRUE`.
#'
#' @return Either a `SummarizedExperiment` object containing variance partitioning results in the metadata,
#' or a list of variance components, depending on `save.se.obj`.
#'
#' @importFrom variancePartition fitExtractVarPartModel
#' @importFrom variancePartition fitVarPartModel
#' @importFrom BiocParallel MulticoreParam
#' @importFrom limma lmFit
#'
#' @export

computeGenesVarianceComposition <- function(
        se.obj,
        assay.names,
        form,
        adjust.data = FALSE,
        adjustment.form = NULL,
        adjustment.method = 'lm',
        samples.to.use = 'all',
        apply.log = TRUE,
        pseudo.count = 1,
        nb.cores = NULL,
        plot.output = TRUE,
        output.name = NULL,
        check.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The computeGeneVarianceComposition function starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking the function inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be NULL.')
    }
    if (is.list(assay.names) | is.logical(assay.names) | is.numeric(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s) or "assay.names = all".')
    }
    if (!is.logical(adjust.data)){
        stop('The "adjust.data" must be logical.')
    }
    if (isTRUE(adjust.data)){
        if (is.null(adjustment.form)){
            stop('The "adjustment.form" cannot be NULL when the "adjust.data" is TRUE')
        }
        if (!is.character(adjustment.method)){
            stop('The "adjustment.method" must be a character string.')
        }
        if (!adjustment.method %in% c('lm', 'lmm')){
            stop('The "adjustment.method" must be one of the "lm" or "lmm".')
        }
    }
    if (!is.logical(samples.to.use)){
        if (!is.character(samples.to.use) | length(samples.to.use) > 1){
            stop('The "samples.to.use" must be logical or set to "all".')
        }
        if (samples.to.use != 'all'){
            stop('The "samples.to.use" must be logical or set to "all".')
        }
    }
    if (is.logical(samples.to.use)){
        if (length(samples.to.use) != ncol(se.obj)){
            stop('The length of the "samples.to.use" must be the same as the number of samples in the data.')
        }
        if (sum(samples.to.use) == 0){
            stop('All elements of the "samples.to.use" are FALSE.')
        }
    }
    if (!is.logical(apply.log)){
        stop('The "apply.log" must be logical.')
    }
    if (isTRUE(apply.log)){
        if (!is.numeric(pseudo.count)){
            stop('The "pseudo.count" must be a numeric value.')
        }
        if (pseudo.count < 0){
            stop('The "pseudo.count" must be a postive numeric value.')
        }
    }
    if (!is.logical(plot.output)){
        stop('The "plot.output" must be logical.')
    }
    if (!is.logical(check.se.obj)){
        stop('The "check.se.obj" must be logical.')
    }
    if (!is.logical(save.se.obj)){
        stop('The "check.se.obj" must be logical.')
    }
    if (!is.logical(verbose)){
        stop('The "check.se.obj" must be logical.')
    }

    # Selecting samples to use ####
    if (is.logical(samples.to.use)){
        se.obj.initial <- se.obj
        se.obj <- se.obj[ , samples.to.use]
    }
    # Specifying cores
    if (is.null(nb.cores)){
        if (.Platform$OS.type == "windows") {
            nb.cores <- as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS", unset = 1))
        } else {
            # macOS or Unix
            nb.cores <- as.numeric(system("sysctl -n hw.ncpu", intern = TRUE)) - 1
            if (is.na(nb.cores) || length(nb.cores) == 0) {
                nb.cores <- 1
            }
        }
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
        variables <- c(
            changeLmmFormula(form),
            changeLmmFormula(adjustment.form)
            )
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            variables = variables[variables!= 'NULL'],
            remove.na = remove.na,
            verbose = verbose
            )
    }
    sample.info <- droplevels(as.data.frame(colData(se.obj)))

    # Applying data log transformation ####
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
            check.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
    }
    if (isFALSE(apply.log)){
        printColoredMessage(
            message = '- The specified assay(s) will be used for LMM, without applying log transformation.',
            color = 'blue',
            verbose = verbose
            )
        all.assays <- lapply(
            levels(assay.names),
            function(x) assay(x = se.obj, i = x)
            )
        names(all.assays) <- levels(assay.names)
    }
    # Applying linear mixed model ####
    printColoredMessage(
        message = '-- Applying linear mixed model on the specified dataset(s):',
        color = 'magenta',
        verbose = verbose
        )
    all.gene.var.decop <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0(
                    '- Performing linear mixed model on the "',
                    x,
                    '" data.'),
                color = 'blue',
                verbose = verbose
                )
            if (isTRUE(adjust.data)){
                if (adjustment.method == 'lm'){
                    printColoredMessage(
                        message = paste0(
                            '- Adjusting the data for ',
                            adjustment.form,
                            ' using linear model and then use the residuals for gene level variance decomposition analysis.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    adjustment.form <- as.formula(
                        paste0(
                            gsub("\\(1 \\| ([^)]+)\\)", "\\1", adjustment.form),
                            collapse = '')
                        )
                    lm.fit <- limma::lmFit(
                        object = all.assays[[x]],
                        design =  model.matrix(adjustment.form, sample.info)
                        )
                    res.data <- residuals(object = lm.fit, all.assays[[x]])
                    printColoredMessage(
                        message = paste0(
                            '- Fitting linear mixed model for ',
                            paste(deparse(form), collapse = ""),
                            ' to obtain gene level variance composition.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    gene.var.decomposition <- fitExtractVarPartModel(
                        exprObj = res.data,
                        formula = form,
                        data = sample.info,
                        BPPARAM = MulticoreParam(workers = nb.cores)
                        )
                }
                if (adjustment.method == 'lmm'){
                    printColoredMessage(
                        message = paste0(
                            '- Adjusting the data for ',
                            deparse(adjustment.form),
                            ' using linear mixed model and then use the residuals for gene level variance decomposition.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    lmm.fit <- fitVarPartModel(
                        exprObj = all.assays[[x]],
                        formula = adjustment.form,
                        data = sample.info,
                        BPPARAM = MulticoreParam(workers = nb.cores)
                        )
                    res.data <- residuals(lmm.fit)
                    printColoredMessage(
                        message = paste0(
                            '- Fitting linear mixed model for ',
                            paste(deparse(form), collapse = ""),
                            ' to obtain gene level variance decomposition.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    gene.var.decomposition <- fitExtractVarPartModel(
                        exprObj = res.data,
                        formula = form,
                        data = sample.info,
                        BPPARAM = MulticoreParam(workers = nb.cores)
                        )
                    row.names(gene.var.decomposition) <- row.names(se.obj)
                }
            }
            if (isFALSE(adjust.data)){
                printColoredMessage(
                    message = paste0(
                        '- Fitting a linear mixed model for ',
                        paste(deparse(form), collapse = ""),
                        ' to obtain gene level variance composition.'),
                    color = 'blue',
                    verbose = verbose
                    )
                gene.var.decomposition <- fitExtractVarPartModel(
                    exprObj = all.assays[[x]],
                    formula = form,
                    data = sample.info,
                    BPPARAM = MulticoreParam(workers = nb.cores)
                    )
                gene.var.decomposition
            }
            return(gene.var.decomposition)
        })
    names(all.gene.var.decop) <- levels(assay.names)
    # Selecting samples to use ####
    if (is.logical(samples.to.use)) se.obj <- se.obj.initial

    # Saving the results
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- The gene variance decomposition results for the indiviaul assay(s) are saved to the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        if (isFALSE(adjust.data)){
            form <- paste0(
                    gsub("\\(1 \\| ([^)]+)\\)", "\\1", form),
                    collapse = ''
                    )
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = levels(assay.names),
                assessment.type = 'gene.level',
                assessment = 'GeneVarDecompostion',
                method = 'LMM',
                variables = form,
                file.name = 'gene.variance',
                results.data = all.gene.var.decop
                )
        }
        if (isTRUE(adjust.data)){
            form <- paste0(
                    gsub("\\(1 \\| ([^)]+)\\)", "\\1", form),
                    collapse = ''
                    )
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = levels(assay.names),
                assessment.type = 'gene.level',
                assessment = 'GeneVarDecompostion',
                method = paste0('Adjusted.by.', adjustment.method),
                variables = form,
                file.name = 'gene.variance',
                results.data = all.gene.var.decop
            )
        }
        printColoredMessage(
            message = '------------The computeGeneVarianceComposition function finished.',
            color = 'white',
            verbose = verbose
        )
        return(se.obj = se.obj)
    }
    if (isFALSE(save.se.obj)){
        return(all.gene.var.decop)
    }
}
