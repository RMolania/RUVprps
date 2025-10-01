#' Estimates the CMS subtypes of colorectal cancer RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function employs the `CMScaller` function from the `CMScaller` R package to estimate the CMS subtypes of
#' colorectal cancer RNA-seq data. The function allows flexibility in choosing assays, applying transformations,
#' handling missing values, and parallelizing computations.
#'
#' @param se.obj A `SummarizedExperiment` object containing the RNA-seq data.
#' @param assay.names Character or character vector. Specifies the name(s) of the assay(s) in the `SummarizedExperiment`
#' object to be used for CMS subtype estimation.
#' @param raw.count.assay.name Character. The name of the assay containing raw counts, used if normalization or log
#' transformation is applied.
#' @param subtype Character. The subtype classification method to use. Options are those supported by `CMScaller`.
#' @param tissue.type Character. A character string specifying the tissue type (e.g., `colon`, `rectum`) for subtype assignment.
#' @param groups Character. A vector of group labels to be used for stratification or annotation in the subtype analysis.
#' @param templates Character. A template set for CMS classification, as required by `CMScaller`.
#' @param row.names Character or logical. Indicates whether to use row names as gene identifiers.
#' @param count.data Logical. If `TRUE`, the input assay data are treated as raw counts. The default is `FALSE`.
#' @param nb.perm Integer. Number of permutations for the subtype classification. The default is as in `CMScaller`.
#' @param fdr Numeric. False discovery rate threshold to call CMS subtypes.
#' @param generate.plot Logical. If `TRUE`, diagnostic plots of subtype classification will be generated. The default is `TRUE`.
#' @param nb.cores Integer. Number of CPU cores to use for parallel processing. The default is 1.
#' @param apply.log Logical. Indicates whether to apply log2 transformation to the input assay data. The default is `TRUE`.
#' @param pseudo.count Numeric. A pseudo count value to be added prior to log2 transformation to avoid `-Inf` values.
#' The default is set to 1.
#' @param out.put.name Character. The name prefix for saving the output results or plots.
#' @param check.se.obj Logical. If `TRUE`, the input `SummarizedExperiment` object will be checked for validity before
#' processing. The default is `TRUE`.
#' @param remove.na Character. Indicates how missing values should be handled. Options are `assays` or `none`. The
#' default is `assays`.
#' @param seed Integer. A random seed to ensure reproducibility of the results.
#' @param verbose Logical. If `TRUE`, messages and progress updates will be displayed during function execution. The
#' default is `TRUE`.
#'
#' @return A `SummarizedExperiment` object with CMS subtype classification results added to the metadata, along with
#' optional plots if `generate.plot = TRUE`.
#'
#' @importFrom SummarizedExperiment assays colData colData<-
#' @importFrom CMScaller CMScaller
#' @importFrom parallel mclapply
#' @export

estimateCMS <- function(
        se.obj,
        assay.names = 'all',
        raw.count.assay.name = NULL,
        subtype = 'CMS',
        tissue.type = NULL,
        groups = NULL,
        templates = CMScaller::templates.CMS,
        row.names = 'symbol',
        count.data = FALSE,
        nb.perm  = 1000,
        fdr = 0.05,
        generate.plot = FALSE,
        nb.cores = 1,
        apply.log = TRUE,
        pseudo.count = 1,
        out.put.name = NULL,
        check.se.obj = TRUE,
        remove.na = 'none',
        seed = NULL,
        verbose = TRUE
        ){
    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
        } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
        }

    # Assessing the SummarizedExperiment object ####
    if (isTRUE(check.se.obj)) {
        if (is.null(tissue.type) & is.null(groups)){
            variables = NULL
        }
        if (is.null(tissue.type) & !is.null(groups)){
            variables = groups
        }
        if (!is.null(tissue.type) & is.null(groups)){
            variables = tissue.type
        }
        if (!is.null(tissue.type) & !is.null(groups)){
            variables = c(tissue.type, groups)
        }
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.names,
            variables = variables,
            remove.na = remove.na,
            verbose = verbose
        )
    }
    # Checking the tissue.type variable ####
    if (!is.null(tissue.type)){
        tissue.type.name <- tissue.type
        if (!is.character(tissue.type)){
            stop('The "groups" must be a single name of a column in the SummarizedExperiment object.')
        }
        if (length(tissue.type) > 1){
            stop('The "tissue.type" must be a single name of a column in the SummarizedExperiment object.')
        }
        if (!tissue.type %in% colnames(SummarizedExperiment::colData(se.obj))){
            stop('The "tissue.type" cannot be found in the SummarizedExperiment object.')
        }
        tissue.type <- factor(x = se.obj[[tissue.type]], levels = unique(se.obj[[tissue.type]]))
        if (sum(tissue.type %in% c('Tumor', 'Normal')) != length(tissue.type)){
            stop('The "tissue.type" must contain "Tumor" and "Normal" factors.')
        }
        if (!'Tumor' %in% tissue.type){
            stop('The "tissue.type" has no any "Tumor" factor.')
        }
    }
    if (!is.null(groups)){
        groups.name <- groups
        groups <- factor(x = se.obj[[groups]], levels = unique(se.obj[[groups]]))
    }

    if (is.null(raw.count.assay.name)){
        names(assay.names) <- assay.names
    }
    if (!is.null(raw.count.assay.name)){
        assay.names.raw <- assay.names[assay.names %in% raw.count.assay.name]
        assay.names.norm <- assay.names[!assay.names %in% raw.count.assay.name]
        assay.names <- c(assay.names.raw, assay.names.norm)
        names(assay.names) <- c('Quantile', levels(assay.names)[-1])
    }
    if (is.null(tissue.type) & is.null(groups)){
        cms.calls <- mclapply(
            1:length(assay.names),
            function(x){
                expr.matrix <- as.matrix(assay(x = se.obj, i = assay.names[x]))
                if (names(assay.names[x]) == 'Quantile'){
                    rna.seq = TRUE
                } else {
                    rna.seq = FALSE
                    if (isTRUE(apply.log)){
                        expr.matrix <- log2(expr.matrix + pseudo.count)
                        }
                    }
                cms.cluster <- CMScaller(
                    emat = expr.matrix,
                    templates = templates,
                    rowNames = row.names,
                    RNAseq = rna.seq,
                    nPerm = nb.perm,
                    seed = seed,
                    FDR = fdr,
                    doPlot = generate.plot,
                    verbose = verbose
                    )
                return(cms.cluster)
            }, mc.cores = nb.cores)
        names(cms.calls) <- names(assay.names)
        for(i in names(cms.calls)){
            if (!is.null(out.put.name)){
                col.name <- paste0(i, '.CMS.', out.put.name)
            } else col.name <- paste0(i, '.CMS')
            col.name <- paste0(i, '.CMS')
            SummarizedExperiment::colData(se.obj)[ , col.name] <- 'Normal'
            index <- match(row.names(cms.calls[[col.name]]), colnames(se.obj))
            SummarizedExperiment::colData(se.obj)[ , col.name][index] <- as.character(cms.calls[[col.name]]$prediction)
            SummarizedExperiment::colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not.classified'
            SummarizedExperiment::colData(se.obj)[ , col.name] <- factor(
                x = SummarizedExperiment::colData(se.obj)[ , col.name],
                levels = c('CMS1','CMS2','CMS3','CMS4','Not.classified'))
        }
    }
    if (!is.null(tissue.type) & is.null(groups)){
        cms.calls <- mclapply(
            1:length(assay.names),
            function(x){
                cancer.samples <- se.obj[[tissue.type.name]] == 'Tumor'
                if (names(assay.names[x]) == 'Quantile'){
                    rna.seq = TRUE
                    expr.matrix <- as.matrix(assay(x = se.obj[ , cancer.samples], i = assay.names[x]))
                }
                if (names(assay.names[x]) != 'Quantile'){
                    expr.matrix <- as.matrix(assay(x = se.obj[ , cancer.samples], i = assay.names[x]))
                    rna.seq = FALSE
                    if (isTRUE(apply.log)){
                        expr.matrix <- log2(expr.matrix + pseudo.count)
                    }
                }
                cms.cluster <- CMScaller(
                    emat = expr.matrix,
                    templates = templates,
                    rowNames = row.names,
                    RNAseq = rna.seq,
                    nPerm = nb.perm,
                    seed = seed,
                    FDR = fdr,
                    doPlot = generate.plot,
                    verbose = verbose
                )
                return(cms.cluster)
            }, mc.cores = nb.cores
            )
        names(cms.calls) <- names(assay.names)
        if (subtype == 'CMS'){
            for(i in names(cms.calls)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.CMS.', out.put.name)
                } else col.name <- paste0(i, '.CMS')
                SummarizedExperiment::colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls[[i]]), colnames(se.obj))
                SummarizedExperiment::colData(se.obj)[ , col.name][index] <- as.character(cms.calls[[i]]$prediction)
                SummarizedExperiment::colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not.classified'
                SummarizedExperiment::colData(se.obj)[ , col.name] <- factor(
                    x = SummarizedExperiment::colData(se.obj)[ , col.name],
                    levels = c('CMS1','CMS2','CMS3','CMS4', 'Normal', 'Not.classified'))
            }
        }
        if (subtype == 'MSI'){
            for(i in names(cms.calls)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.MSI.', out.put.name)
                } else col.name <- paste0(i, '.MSI')
                SummarizedExperiment::colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls[[i]]), colnames(se.obj))
                SummarizedExperiment::colData(se.obj)[ , col.name][index] <- as.character(cms.calls[[i]]$prediction)
                SummarizedExperiment::colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not.classified'
                SummarizedExperiment::colData(se.obj)[ , col.name] <- factor(
                    x = SummarizedExperiment::colData(se.obj)[ , col.name],
                    levels = c('MSS','MSI', 'Normal', 'Not.classified'))
            }
        }
        if (subtype == 'CRIS'){
            for(i in names(cms.calls)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.CRIS.', out.put.name)
                } else col.name <- paste0(i, '.CRIS')
                colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls[[i]]), colnames(se.obj))
                SummarizedExperiment::colData(se.obj)[ , col.name][index] <- as.character(cms.calls[[i]]$prediction)
                SummarizedExperiment::colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not.classified'
                SummarizedExperiment::colData(se.obj)[ , col.name] <- factor(
                    x = SummarizedExperiment::colData(se.obj)[ , col.name],
                    levels = c('CRISA','CRISB', 'CRISC', 'CRISD', 'CRISE', 'Normal', 'Not.classified'))
            }
        }
    }
    if (!is.null(tissue.type) & !is.null(groups)){
        cms.calls.per.group <- mclapply(
            assay.names,
            function(x){
                cms.calls <- lapply(
                    levels(groups),
                    function(y){
                        cancer.time.samples <- se.obj[[groups.name]] == y & se.obj[[tissue.type.name]] == 'Tumor'
                        if(names(assay.names[x]) == 'Quantile'){
                            expr.matrix = as.matrix(assay(x = se.obj[ , cancer.time.samples], i = x))
                            rna.seq = TRUE
                        }
                        if (names(assay.names[x]) != 'Quantile'){
                            expr.matrix = as.matrix(assay(x = se.obj[ , cancer.time.samples], i = x))
                            rna.seq = FALSE
                            if (isTRUE(apply.log)){
                                expr.matrix <- log2(expr.matrix + pseudo.count)
                            }
                        }
                        cms.cluster <- CMScaller(
                            emat = expr.matrix,
                            templates = templates,
                            rowNames = row.names,
                            RNAseq = rna.seq,
                            nPerm = nb.perm,
                            seed = seed,
                            FDR = fdr,
                            doPlot = generate.plot,
                            verbose = verbose
                        )
                        return(cms.cluster)
                    })
                cms.calls <- do.call(rbind, cms.calls)
                return(cms.calls)
            }, mc.cores = nb.cores)
        names(cms.calls.per.group) <- names(assay.names)

        if (subtype == 'CMS'){
            for(i in names(cms.calls.per.group)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.CMS.', out.put.name)
                } else col.name <- paste0(i, '.CMS')
                SummarizedExperiment::colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls.per.group[[i]]), colnames(se.obj))
                SummarizedExperiment::colData(se.obj)[ , col.name][index] <- as.character(cms.calls.per.group[[i]]$prediction)
                SummarizedExperiment::colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not.classified'
                SummarizedExperiment::colData(se.obj)[ , col.name] <- factor(
                    x = SummarizedExperiment::colData(se.obj)[ , col.name],
                    levels = c('CMS1','CMS2','CMS3','CMS4','Normal','Not.classified'))
            }
        }
        if (subtype == 'MSI'){
            for(i in names(cms.calls.per.group)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.MSI.', out.put.name)
                } else col.name <- paste0(i, '.MSI')
                SummarizedExperiment::colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls.per.group[[i]]), colnames(se.obj))
                SummarizedExperiment::colData(se.obj)[ , col.name][index] <- as.character(cms.calls.per.group[[i]]$prediction)
                SummarizedExperiment::colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not.classified'
                SummarizedExperiment::colData(se.obj)[ , col.name] <- factor(
                    x = SummarizedExperiment::colData(se.obj)[ , col.name],
                    levels = c('MSS','MSI','Normal','Not.classified'))
            }
        }
        if (subtype == 'CRIS'){
            for(i in names(cms.calls.per.group)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.CRIS.', out.put.name)
                } else col.name <- paste0(i, '.CRIS')
                SummarizedExperiment::colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls.per.group[[i]]), colnames(se.obj))
                SummarizedExperiment::colData(se.obj)[ , col.name][index] <- as.character(cms.calls.per.group[[i]]$prediction)
                SummarizedExperiment::colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not.classified'
                SummarizedExperiment::colData(se.obj)[ , col.name] <- factor(
                    x = SummarizedExperiment::colData(se.obj)[ , col.name],
                    levels = c('CRISA','CRISB', 'CRISC', 'CRISD', 'CRISE', 'Normal', 'Not.classified'))
            }
        }
    }
    if (is.null(tissue.type) & !is.null(groups)){
        cms.calls.per.group <- lapply(
            assay.names,
            function(x){
                cms.calls <- lapply(
                    levels(groups),
                    function(y){
                        cancer.time.samples <- se.obj[[groups.name]] == y & se.obj[[tissue.type.name]] == 'Tumor'
                        expr.matrix <- as.matrix(assay(x = se.obj[ , cancer.time.samples], i = x))
                        if (isTRUE(apply.log)){
                            expr.matrix <- log2(expr.matrix + pseudo.count)
                        }
                        if(x == 'Quantile'){
                            rna.seq = TRUE
                        } else rna.seq = FALSE
                        cms.cluster <- CMScaller(
                            emat = expr.matrix,
                            templates = templates,
                            rowNames = row.names,
                            RNAseq = rna.seq,
                            nPerm = nb.perm,
                            seed = seed,
                            FDR = fdr,
                            doPlot = generate.plot,
                            verbose = verbose
                        )
                        return(cms.cluster)
                    })
                cms.calls <- do.call(rbind, cms.calls)
                return(cms.calls)
            })
        names(cms.calls.per.group) <- names(assay.names)
        for(i in names(cms.calls.per.group)){
            if (!is.null(out.put.name)){
                col.name <- paste0(i, '.CMS.', out.put.name)
            } else col.name <- paste0(i, '.CMS')
            SummarizedExperiment::colData(se.obj)[ , col.name] <- 'Normal'
            index <- match(row.names(cms.calls.per.group[[i]]), colnames(se.obj))
            SummarizedExperiment::colData(se.obj)[ , col.name][index] <- as.character(cms.calls.per.group[[i]]$prediction)
            SummarizedExperiment::colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not.classified'
            SummarizedExperiment::colData(se.obj)[ , col.name] <- factor(
                x = SummarizedExperiment::colData(se.obj)[ , col.name],
                levels = c('CMS1','CMS2','CMS3','CMS4','Normal','Not.classified'))
        }
    }
    return(se.obj)
}
