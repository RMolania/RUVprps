#' Find canonical correlation coordinates in RNA-seq data

#' @author Ramyar Molania

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name TTTTT
#' @param variable TTTTT
#' @param reference.group TTTTT
#' @param nb.cca TTTTT
#' @param hvg TTTTT
#' @param scale TTTTT
#' @param normalization TTTTT
#' @param regress.out.variables TTTTT
#' @param regress.out.rle.med TTTTT
#' @param cosine.norm TTTTT
#' @param samples.to.use TTTTT
#' @param apply.log TTTTT
#' @param pseudo.count TTTTT
#' @param check.se.obj TTTTT
#' @param remove.na TTTTT
#' @param cca.set.name TTTTT
#' @param save.se.obj TTTTT
#' @param verbose TTTTT
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
            message = '-- The set of HVG is outpputed as a logical vector.',
            color = 'magenta',
            verbose = verbose
        )
        printColoredMessage(
            message = '------------The findHVG function finished.',
            color = 'white',
            verbose = verbose
        )
        return(all.cca)
    }

}

# pairs(all.cca$Batch1_Batch2[,1:4], col = factor(read.se.obj$Time.interval))
# pairs(all.cca$Batch1_Batch2[,1:4], col = factor(read.se.obj$Tissues))
# pairs(all.cca$Batch1_Batch2[,1:4], col = factor(read.se.obj$Quantile.CMS.TI))
#
# pairs(cca.data[,1:4], col = factor(read.se.obj$Time.interval))
# pairs(cca.data[,1:4], col = factor(read.se.obj$Tissues))
# pairs(cca.data[,1:4], col = factor(read.se.obj$Quantile.CMS.TI))
#
#
#
# data.0 <- assay(read.se.obj, 'RawCount')
# data.1 <- data.0[ , read.se.obj$Time.interval == 'Batch1']
# data.2 <- data.0[ , read.se.obj$Time.interval == 'Batch2']
#
#
#
# data.0 <- cpm(assay(read.se.obj, 'RawCount'), log = TRUE)
# data.a <- data.0[ , read.se.obj$Time.interval == 'Batch1']
# data.b <- data.0[ , read.se.obj$Time.interval == 'Batch2']
#
# data.0 <- assay(read.se.obj, 'RawCount')
# data.a <- cpm(data.0[ , read.se.obj$Time.interval == 'Batch1'], log = T)
# data.b <- cpm(data.0[ , read.se.obj$Time.interval == 'Batch2'], log = T)
#
#
# library(Seurat)
# s1 <- CreateSeuratObject(counts = data.1[ , ])
# s2 <- CreateSeuratObject(counts = data.2[, ])
# s1 <- NormalizeData(s1)
# s2 <- NormalizeData(s2)
#
#
# s1 <- SetAssayData(s1, layer = "data", new.data = as.matrix(data.a))
# s2 <- SetAssayData(s2, layer = "data", new.data = as.matrix(data.b))
#
#
# s1 <- FindVariableFeatures(s1)
# s2 <- FindVariableFeatures(s2)
# s1 <- ScaleData(object = s1 , features = row.names(data.1))
# s2 <- ScaleData(object = s2 , features = row.names(data.2))
# pp <- intersect(VariableFeatures(s1), VariableFeatures(s2))
# cca_out <- RunCCA(s1, s2, features = pp, num.cc = 10)
# cca_embeddings <- Embeddings(cca_out, reduction = "cca")
# pairs(cca_embeddings[,1:4], col = factor(read.se.obj$Time.interval))
# pairs(cca_embeddings[,1:4], col = factor(read.se.obj$Tissues))
# pairs(cca_embeddings[,1:4], col = factor(read.se.obj$Quantile.CMS.TI))
#
# plot(cca_embeddings[,1], cca.data[,1])
# plot(cca_embeddings[,2], cca.data[,2])
# plot(cca_embeddings[,3], cca.data[,3])
#
# pairs(cca_embeddings[, 1:3], col = factor(brca.se.obj$paper_BRCA_Subtype_PAM50))
# pairs(cca_embeddings[, 1:3], col = factor(all.sample$studies))
# cca.mnn <- findMutualNN(data1 = cca_embeddings[colnames(data.a) , ], data2 = cca_embeddings[colnames(data.b),  ], k1  = 3)
# cca.mnn <- data.frame(left = cca.mnn$first, right = cca.mnn$second + 221)
# cca.mnn$bio1 <- all.sample$pam50[cca.mnn$left]
# cca.mnn$bio2 <- all.sample$pam50[cca.mnn$right]
# table(apply(cca.mnn[ , 3:4], 1, function(x) length(unique(x))))
#
# pp <- as.data.frame(cca.mnn@metadata$merge.info$pairs@listData)
# pp$bio1 <- all.sample$pam50[pp$left]
# pp$bio2 <- all.sample$pam50[pp$right]
# table(apply(pp[ , 3:4], 1, function(x) length(unique(x))))
# 368/108

# data.1 <- log2(assay(read.se.obj[ , read.se.obj$Years == '2013'], 'RawCount') + 1)
# data.2 <- log2(assay(read.se.obj[ , read.se.obj$Years == '2014'], 'RawCount') + 1)
#
# cross.pro.data <- crossprod(x = data.1, y = data.2)
# cca.svd <- irlba(A = cross.pro.data, nv = nb.cca)

