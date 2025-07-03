#' Estimates the CMS subtypes of colorectal cancer RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function employees the `CMScaller` function from the `CMScaller` R package to estimates the CMS subtypes of
#' colorectal cancer RNA-seq data
#'
#' @param se.obj description
#' @param assay.names description
#' @param raw.count.assay.name description
#' @param subtype description
#' @param tissue.type description
#' @param templates description
#' @param row.names description
#' @param count.data description
#' @param nb.perm description
#' @param fdr description
#' @param generate.plot description
#' @param nb.cores description
#' @param apply.log description
#' @param pseudo.count description
#' @param out.put.name description
#' @param check.se.obj description
#' @param remove.na description
#' @param seed description
#' @param verbose description
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom CMScaller CMScaller
#' @importFrom parallel mclapply
#' @export

computeCMS <- function(
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
        if (!tissue.type %in% colnames(colData(se.obj))){
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
                print(x)
                expr.matrix <- as.matrix(assay(x = se.obj, i = assay.names[x]))
                if (isTRUE(apply.log)){
                    expr.matrix <- log2(expr.matrix + pseudo.count)
                }
                if (names(assay.names[x]) == 'Quantile'){
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
            }, mc.cores = nb.cores)
        names(cms.calls) <- names(assay.names)
        for(i in names(cms.calls)){
            if (!is.null(out.put.name)){
                col.name <- paste0(i, '.CMS.', out.put.name)
            } else col.name <- paste0(i, '.CMS')
            col.name <- paste0(i, '.CMS')
            colData(se.obj)[ , col.name] <- 'Normal'
            index <- match(row.names(cms.calls[[col.name]]), colnames(se.obj))
            colData(se.obj)[ , col.name][index] <- as.character(cms.calls[[col.name]]$prediction)
            colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not classified'
            colData(se.obj)[ , col.name] <- factor(
                x = colData(se.obj)[ , col.name],
                levels = c('CMS1','CMS2','CMS3','CMS4','Not classified'))
        }
    }
    if (!is.null(tissue.type) & is.null(groups)){
        cms.calls <- mclapply(
            1:length(assay.names),
            function(x){
                cancer.samples <- se.obj[[tissue.type.name]] == 'Tumor'
                expr.matrix <- as.matrix(assay(x = se.obj[ , cancer.samples], i = assay.names[x]))
                if (isTRUE(apply.log)){
                    expr.matrix <- log2(expr.matrix + pseudo.count)
                }
                if (names(assay.names[x]) == 'Quantile'){
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
            }, mc.cores = nb.cores)
        names(cms.calls) <- names(assay.names)
        if (subtype == 'CMS'){
            for(i in names(cms.calls)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.CMS.', out.put.name)
                } else col.name <- paste0(i, '.CMS')
                colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls[[i]]), colnames(se.obj))
                colData(se.obj)[ , col.name][index] <- as.character(cms.calls[[i]]$prediction)
                colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not classified'
                colData(se.obj)[ , col.name] <- factor(
                    x = colData(se.obj)[ , col.name],
                    levels = c('CMS1','CMS2','CMS3','CMS4', 'Normal', 'Not classified'))
            }
        }
        if (subtype == 'MSI'){
            for(i in names(cms.calls)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.MSI.', out.put.name)
                } else col.name <- paste0(i, '.MSI')
                colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls[[i]]), colnames(se.obj))
                colData(se.obj)[ , col.name][index] <- as.character(cms.calls[[i]]$prediction)
                colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not classified'
                colData(se.obj)[ , col.name] <- factor(
                    x = colData(se.obj)[ , col.name],
                    levels = c('MSS','MSI', 'Normal', 'Not classified'))
            }
        }
        if (subtype == 'CRIS'){
            for(i in names(cms.calls)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.CRIS.', out.put.name)
                } else col.name <- paste0(i, '.CRIS')
                colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls[[i]]), colnames(se.obj))
                colData(se.obj)[ , col.name][index] <- as.character(cms.calls[[i]]$prediction)
                colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not classified'
                colData(se.obj)[ , col.name] <- factor(
                    x = colData(se.obj)[ , col.name],
                    levels = c('CRISA','CRISB', 'CRISC', 'CRISD', 'CRISE', 'Normal', 'Not classified'))
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
                        expr.matrix <- as.matrix(assay(x = se.obj[ , cancer.time.samples], i = x))
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
            }, mc.cores = nb.cores)
        names(cms.calls.per.group) <- names(assay.names)

        if (subtype == 'CMS'){
            for(i in names(cms.calls.per.group)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.CMS.', out.put.name)
                } else col.name <- paste0(i, '.CMS')
                colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls.per.group[[i]]), colnames(se.obj))
                colData(se.obj)[ , col.name][index] <- as.character(cms.calls.per.group[[i]]$prediction)
                colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not classified'
                colData(se.obj)[ , col.name] <- factor(
                    x = colData(se.obj)[ , col.name],
                    levels = c('CMS1','CMS2','CMS3','CMS4','Normal','Not classified'))
            }
        }
        if (subtype == 'MSI'){
            for(i in names(cms.calls.per.group)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.MSI.', out.put.name)
                } else col.name <- paste0(i, '.MSI')
                colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls.per.group[[i]]), colnames(se.obj))
                colData(se.obj)[ , col.name][index] <- as.character(cms.calls.per.group[[i]]$prediction)
                colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not classified'
                colData(se.obj)[ , col.name] <- factor(
                    x = colData(se.obj)[ , col.name],
                    levels = c('MSS','MSI','Normal','Not classified'))
            }
        }
        if (subtype == 'CRIS'){
            for(i in names(cms.calls.per.group)){
                if (!is.null(out.put.name)){
                    col.name <- paste0(i, '.CRIS.', out.put.name)
                } else col.name <- paste0(i, '.CRIS')
                colData(se.obj)[ , col.name] <- 'Normal'
                index <- match(row.names(cms.calls.per.group[[i]]), colnames(se.obj))
                colData(se.obj)[ , col.name][index] <- as.character(cms.calls.per.group[[i]]$prediction)
                colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not classified'
                colData(se.obj)[ , col.name] <- factor(
                    x = colData(se.obj)[ , col.name],
                    levels = c('CRISA','CRISB', 'CRISC', 'CRISD', 'CRISE', 'Normal', 'Not classified'))
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
            colData(se.obj)[ , col.name] <- 'Normal'
            index <- match(row.names(cms.calls.per.group[[i]]), colnames(se.obj))
            colData(se.obj)[ , col.name][index] <- as.character(cms.calls.per.group[[i]]$prediction)
            colData(se.obj)[ , col.name][is.na(colData(se.obj)[ , col.name])] <- 'Not classified'
            colData(se.obj)[ , col.name] <- factor(
                x = colData(se.obj)[ , col.name],
                levels = c('CMS1','CMS2','CMS3','CMS4','Normal','Not classified'))
        }
    }
    return(se.obj)
}
