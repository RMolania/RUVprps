#' Creates PRPS sets using mutual nearest neighbors in RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function uses mutual nearest neighbors approach to create PRPS in the RNA-seq data. This function can be used in
#' situation where the biological variation are entirely unknown.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character indicating the name of the data (assay) in the SummarizedExperiment object.
#' This data will be used to create PRPS data for RUV-III normalization. This data must be the one that will be
#' used for the RUV-III normalization.
#' @param main.uv.variable Character. Indicates the name of a column in the sample annotation of the SummarizedExperiment
#' object. The `uv.variable` can be either categorical or continuous. If `uv.variable` is a continuous variable, this will
#' be divided into `nb.clusters` groups using the `clustering.method`.
#' @param clustering.method Character. A character indicating the choice of clustering method for grouping the
#' `uv.variable` if a continuous variable is provided. Options include `kmeans`, `cut`, and `quantile`. The default is set
#' to `kmeans`.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the `uv.variable` is a
#' continuous variable. The default is set to 3.
#' @param other.uv.variables Character. A character or character vector representing the name(s) of the columns of
#' unwanted variable(s) within the sample annotation (colData) of the SummarizedExperiment object. These can be categorical,
#' continuous, or a combination. These variables will be considered when generating PRPS sets for the `main.uv.variable`
#' to help avoid potential contamination. The default is set to `NULL`
#' @param other.uv.clustering.method Character. A character indicating which clustering method should be used to
#' group each continuous unwanted variable, if specified in `other.uv.variables`. Options include `kmeans`, `cut`,
#' and `quantile`. The default is set to `kmeans`. See createHomogeneousUVGroups() for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous
#' unwanted variable specified in the `other.uv.variables`. The default is set to 3.
#' @param min.sample.for.ps Numeric. Minimum number of samples required for pseudo-replicate creation. The default is set
#' to 3.
#' @param select.extreme.groups Logical. Indicates whether to select only the extreme groups e.g., highest and lowest
#' clusters, when the `uv.variable` is continuous. Default is set to `TRUE`. This will increase the variation between
#' PR sets in order to better capture the unwanted variation.
#' @param filter.prps.sets Logical. If `TRUE`, the number of PRPS sets across each pair of batches will be filtered if
#' they are higher than the `max.prps.sets` value. The default is set to `TRUE`. A high number of PRPS sets will increase
#' the computational time for the RUV-III normalization.
#' @param max.prps.sets Numeric. A numeric value specifying the maximum number for PRPS sets across each pair of batches.
#' The default is set to 10.
#' @param min.batches.to.cover Numeric. Minimum number of batches that must be covered by PRPS set. The default is set to
#' `all`, indicating all possible batch must have enough samples to create PRPS, otherwise the function gives error.
#' @param check.prps.connectedness Logical. Indicates whether to assess the `connectedness` between the PRPS sets across
#' all batches. Default is set to `TRUE`, indicating if there is not connections between all PRPS sets across all batches,
#' the function will give error. We refer to the checkPRPSconnectedness() function for more details.
#' @param nb.mnn Numeric. A numeric value specifying the maximum number of mutual nearest neighbors to compute. The
#' default is set to 1.
#' @param hvg Vector. A logical vector or a vector of the names (feature ids) of the highly variable genes. These genes
#' will be used to prepare the input data for knn and mnn analysis. The default is set to `NULL`, this means all genes
#' will be used.
#' @param normalization Character. A character that indicates which normalization method should be applied on the
#' data before finding the knn. Options are: `CPM`, `TMM`, `upper`, `median`, `full`, and `VST`. The default is set to
#' `cpm`.
#' If set to `NULL`, no normalization will be applied. See the applyOtherNormalizations() function for more details.
#' @param apply.cosine.norm Logical. Indicates whether cosine normalization should be applied before finding MNN. Default
#' is set to `TRUE`.
#' @param regress.out.variables Character. A character or a vector of character that indicate the column name(s) in the
#' sample annotation in the SummarizedExperiment object. These variables will be regressed out from the data before
#' finding MNN. The default is set to `NULL`, indicating that regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not for down-stream analysis.
#' The default is set to `TRUE`.
#' @param pseudo.count Numeric. A positive numeric value as a pseudo count to be added to all measurements of the specified
#' assay(data) before applying log transformation to avoid -Inf for measurements that are equal to 0. The default is set
#' to 1.
#' @param mnn.bpparam Character. A BiocParallelParam object specifying how palatalization should be performed to find MNN.
#' The default is set to SerialParam(). We refer to the `findMutualNN()` function from the **BiocNeighbors** R package.
#' @param mnn.nbparam Character. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' The default is KmknnParam(). We refer to the `findMutualNN()` function from the `BiocNeighbors` R package.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default is set
#' to `TRUE`. See the checkSeObj() function for more details.
#' @param remove.na Character. To remove NA or missing values from the assay (data) or not. The options are `assays` and
#' `none`. The default is set to `assays`, so all the NA or missing values from the assay(s) will be removed before computing
#' performing any down-stream analysis. See the `checkSeObj()` function for more details.
#' @param plot.output Logical. If `TRUE`, the function plots the distribution of MNN across the batches and PRPS sets
#' across the `main.uv.variable`.
#' @param prps.group.name Character. A character specifying the name of the prps.group.name to which the current KNN belong.
#' If set to `NULL`, the function will automatically assign a name using  `main.uv.variable`.
#' @param prps.sets.name Character. A character specifying the name of the output file to be saved in the metadata
#' of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0(uv.variable, '|', assay.name)`.
#' @param save.se.obj Logical. Indicates whether to save the KNN results in the metadata of the SummarizedExperiment object
#' or to output the result as a list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, shows the messages of different steps of the function.
#'
#' @return The SummarizedExperiment object containing PRPS data and plot results in the metadata, or a list of
#' these results.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocNeighbors findMutualNN
#' @importFrom batchelor cosineNorm
#' @importFrom RANN nn2
#' @export

# qs::qsave(read.se.obj, 'tes.read.se.obj.qs')
# read.se.obj <- qs::qread('tes.read.se.obj.qs')
# optimized_function <- function(A, B) {
#     # Convert B to a sparse matrix for efficiency if it is large
#     B <- as(B, "CsparseMatrix")
#     # Compute intermediate matrices
#     BtB <- t(B) %*% B
#     BtB_inv <- solve(BtB)
#     BtA <- t(B) %*% A
#     # Perform the matrix operations
#     result <- A - B %*% BtB_inv %*% BtA
#     return(result)
# }
#
# read.se.obj <- findHVG(
#     se.obj = read.se.obj,
#     assay.name = 'RawCount',
#     approach = 'lmm',
#     group.name = 'Time.interval',
#     uv.variables = c('Time.interval', 'Library.size'),
#     form = ~ (1 | Time.interval) + Library.size,
#     nb.hvg = .1,
#     nb.cores = 14,
#     hvg.selection = 'intersect'
#     )
#
#
# se.obj = read.se.obj
# assay.name = 'RawCount'
# main.uv.variable = 'Library.size'
# other.uv.variables = 'Tumour.purity'
# nb.cca = 5
# nb.pca = 10
# min.sample.for.ps = 3
# filter.prps.sets = TRUE
# select.extreme.groups = FALSE
# max.prps.sets = 3
# min.batches.to.cover = 'all'
# reference.group = NULL
# hvg = read.se.obj@metadata$HVG$hvg_lmm$`1479|Bio:Bio.pc1&Bio.pc2&Bio.pc3|Uv:Time.interval&Library.size|RawCount`$hvg.set
# ncg = NULL
# scale = TRUE
# regress.out.rle.med = FALSE
# sample.to.use = 'all'
# check.prps.connectedness = TRUE
# apply.ruviii.norm = TRUE
# use.ruviii.for.mnn = TRUE
# k = 3
# return.ruviii.w = TRUE
# nb.mnn = 3
# min.ps = 4
# min.nb.for.mnn =  1
# similarity.approach = 'euclidean'
# clustering.method = 'kmeans'
# nb.clusters = 3
# other.uv.clustering.method = 'kmeans'
# nb.other.uv.clusters = 2
# cover.all.batches = FALSE
# nb.batches.to.cover = 2
# normalization = 'CPM'
# cosine.norm = FALSE
# regress.out.variables = NULL
# apply.log = TRUE
# apply.log.for.prps = TRUE
# pseudo.count = 1
# use.ruviii = 'similarity'
# assess.variables.association = TRUE
# mnn.bpparam = SerialParam()
# mnn.nbparam = KmknnParam()
# check.se.obj = TRUE
# remove.na = 'both'
# plot.output = TRUE
# prps.group.name = NULL
# prps.sets.name = NULL
# cca.set.name = NULL
# save.se.obj = TRUE
# verbose = TRUE
# coordinates.to.use = 'both'


createPrPsUnSupervisedByCcaPCA <- function(
        se.obj,
        assay.name,
        main.uv.variable,
        other.uv.variables = NULL,
        nb.cca = 5,
        nb.pca = 10,
        min.sample.for.ps = 3,
        filter.prps.sets = TRUE,
        select.extreme.groups = FALSE,
        max.prps.sets = 3,
        min.batches.to.cover = 'all',
        reference.group = NULL,
        hvg = NULL,
        ncg = NULL,
        scale = TRUE,
        regress.out.rle.med = FALSE,
        sample.to.use = 'all',
        check.prps.connectedness = TRUE,
        apply.ruviii.norm = TRUE,
        use.ruviii.for.mnn = TRUE,
        k = 5,
        coordinates.to.use = 'both',
        return.ruviii.w = TRUE,
        nb.mnn = 3,
        min.ps = 10,
        min.nb.for.mnn =  1,
        similarity.approach = 'euclidean',
        clustering.method = 'kmeans',
        nb.clusters = 3,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        cover.all.batches = FALSE,
        nb.batches.to.cover = 2,
        normalization = 'CPM',
        cosine.norm = FALSE,
        regress.out.variables = NULL,
        apply.log = TRUE,
        apply.log.for.prps = TRUE,
        pseudo.count = 1,
        assess.variables.association = TRUE,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        check.se.obj = TRUE,
        remove.na = 'both',
        plot.output = TRUE,
        prps.group.name = NULL,
        prps.sets.name = NULL,
        cca.set.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    # Assessing and grouping the main unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the main unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    ## Considering only main variable ####
    if (is.null(other.uv.variables)){
        all.varibles <- prepareVariableForPrPs(
            se.obj = se.obj,
            main.variable = main.uv.variable,
            other.variables = other.uv.variables,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            other.uv.clustering.method = other.uv.clustering.method,
            nb.other.uv.clusters = nb.other.uv.clusters,
            min.sample.for.ps = min.sample.for.ps,
            select.extreme.groups = select.extreme.groups,
            cover.all.batches = cover.all.batches,
            nb.batches.to.cover = nb.batches.to.cover,
            assess.variables.association = assess.variables.association,
            plot.output = plot.output,
            nb.mnn = nb.mnn,
            verbose = verbose)
    }
    ## Considering other variables ####
    if (!is.null(other.uv.variables)){
        all.varibles <- prepareVariableForPrPs(
            se.obj = se.obj,
            main.variable = main.uv.variable,
            other.variables = other.uv.variables,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            other.uv.clustering.method = other.uv.clustering.method,
            nb.other.uv.clusters = nb.other.uv.clusters,
            min.sample.for.ps = min.sample.for.ps,
            select.extreme.groups = select.extreme.groups,
            cover.all.batches = cover.all.batches,
            nb.batches.to.cover = nb.batches.to.cover,
            assess.variables.association = assess.variables.association,
            plot.output = plot.output,
            nb.mnn = nb.mnn,
            verbose = verbose
            )
        colnames(all.varibles)[3] <- 'main.uv.variable'
    }
    ## Applying a sanity check ####
    if (isFALSE(all.equal(colnames(se.obj), row.names(all.varibles)))){
        stop('There are issues with the order of samples.')
        }
    if (sum(all.varibles$selected == 'TRUE') == 0){
        stop('There are issues with assessing and grouping the main unwanted variable.')
        }
    ## Updating the Summarized Experiment object  ####
    se.obj.initial <- se.obj
    se.obj.initial[[main.uv.variable]] <- all.varibles$main.uv.variable
    selected.samples <- all.varibles$selected == 'TRUE'
    sub.all.varibles <- droplevels(all.varibles[selected.samples , ])
    other.uv.variables <- unique(sub.all.varibles$other.variables)

    all.prps.sets <- lapply(
        other.uv.variables,
        function(d){
            # Performing CCA between batches ####
            printColoredMessage(
                message = '-- Computing CCA:',
                color = 'magenta',
                verbose = verbose
                )
            printColoredMessage(
                message = '- Computing CCA between all specified pairs of batches',
                color = 'orange',
                verbose = verbose
                )
            sample.index <- all.varibles[all.varibles$other.variables == d , ]
            se.obj <- se.obj.initial[ , sample.index$sampl.ids]
            all.cca <- computeCCA(
                se.obj = se.obj,
                assay.name = assay.name,
                variable = main.uv.variable,
                reference.group = reference.group,
                nb.cca = nb.cca,
                hvg = hvg,
                scale = scale,
                normalization = normalization,
                regress.out.variables = regress.out.variables,
                regress.out.rle.med = regress.out.rle.med,
                cosine.norm = cosine.norm,
                sample.to.use = sample.to.use,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                cca.set.name = cca.set.name,
                save.se.obj = FALSE,
                verbose = verbose
                )

            ## Retrieving initial variable names ####
            all.cca.data <- lapply(
                names(all.cca),
                function(name) {
                    matches <- unique(se.obj[[main.uv.variable]])[sapply(
                        unique(se.obj[[main.uv.variable]]),
                        function(pattern) grepl(pattern, name))]
                    matches
                })
            ## Applying a sanity check ####
            printColoredMessage(
                message = '- Applying a sanity check on the calculated CCA.',
                color = 'orange',
                verbose = verbose
            )
            for(i in 1:length(all.cca.data)){
                sample.ids <- c(
                    colnames(se.obj)[se.obj[[main.uv.variable]] == all.cca.data[[i]][1]],
                    colnames(se.obj)[se.obj[[main.uv.variable]] == all.cca.data[[i]][2]]
                )
                if (isFALSE(all.equal(row.names(all.cca[[i]]), sample.ids))){
                    stop('The order of samples in the calculated CCA have issues.')
                }
                rm(sample.ids)
            }

            # Creating PRPS sets ####
            printColoredMessage(
                message = '- Creating all PRPS sets between all pairs of batches:',
                color = 'magenta',
                verbose = verbose
            )
            all.sample.annot <- as.data.frame(colData(se.obj))
            if (is.null(hvg)){
                hvg <- rep(TRUE, nrow(se.obj))
            }
            min.nb.for.mnn.initial <- min.nb.for.mnn
            all.prps.data <- lapply(
                1:length(all.cca.data),
                function(x){
                    printColoredMessage(
                        message = paste0(
                            '- Creating PRPS sets between the ',
                            all.cca.data[[x]][1],
                            ' and ',
                            all.cca.data[[x]][2],
                            ' groups:'),
                        color = 'orange',
                        verbose = verbose
                    )
                    ### Sub-setting the CCA data for each Batch ####
                    temp.cca <- all.cca[[x]]
                    sample.annot.a <- all.sample.annot[all.sample.annot[[main.uv.variable]] == all.cca.data[[x]][1], ]
                    sample.annot.b <- all.sample.annot[all.sample.annot[[main.uv.variable]] == all.cca.data[[x]][2], ]
                    all.samples <- c(row.names(sample.annot.a), row.names(sample.annot.b))
                    if (coordinates.to.use == 'cca'){
                        ## Finding MNN using the CCA coordinates ####
                        printColoredMessage(
                            message = '- Finding MNN using the CCA coordinates.',
                            color = 'orange',
                            verbose = verbose
                        )
                        if (!identical(row.names(temp.cca), all.samples)){
                            stop('The order of samples in the calculated CCA have issues.')
                        }
                        ### Applying MNN ####
                        cca.mnn <- findMutualNN(
                            data1 = temp.cca[row.names(sample.annot.a) , ],
                            data2 = temp.cca[row.names(sample.annot.b) , ],
                            k1 = min.nb.for.mnn
                        )
                        ### Creating MNN - CCA data and adding samples ids ####
                        all.pp <- data.frame(
                            left.index = cca.mnn$first,
                            right.index = cca.mnn$second,
                            left.samples = row.names(sample.annot.a)[cca.mnn$first],
                            right.samples = row.names(sample.annot.b)[cca.mnn$second]
                        )
                    }
                    if (coordinates.to.use == 'pca'){
                        ## Finding MNN using the PCA coordinates ####
                        printColoredMessage(
                            message = '- Finding MNN using the PCA coordinates.',
                            color = 'orange',
                            verbose = verbose
                        )
                        ### Applying normalization, regression an log on the data ####
                        temp.data <- preProcessData(
                            se.obj = se.obj[ , all.samples],
                            assay.name = assay.name,
                            normalization = normalization,
                            regress.out.variables = regress.out.variables,
                            apply.log = apply.log,
                            pseudo.count = pseudo.count,
                            check.se.obj = FALSE,
                            remove.na = 'none',
                            verbose = verbose
                        )
                        if (!identical(colnames(temp.data), all.samples)){
                            stop('The order of samples in the calculated CCA have issues.')
                        }
                        if (isTRUE(cosine.norm)){
                            temp.data <- cosineNorm(x = temp.data, mode = "matrix")
                        }
                        ### Performing PCA on the data ####
                        pca.data <- irlba::prcomp_irlba(
                            x = t(temp.data[hvg , ]),
                            center = TRUE,
                            scale = FALSE,
                            n = nb.pca
                        )
                        ### Applying MNN ####
                        row.names(pca.data$x) <- colnames(temp.data)
                        pca.mnn <- findMutualNN(
                            data1 = pca.data$x[row.names(sample.annot.a) , ],
                            data2 = pca.data$x[row.names(sample.annot.b) , ],
                            k1  = min.nb.for.mnn
                        )
                        ### Creating MNN - PCA data and adding samples ids ####
                        pca.mnn <- data.frame(
                            left.index = pca.mnn$first,
                            right.index = pca.mnn$second,
                            left.samples = row.names(sample.annot.a)[pca.mnn$first],
                            right.samples = row.names(sample.annot.b)[pca.mnn$second]
                        )
                    }
                    if (coordinates.to.use == 'both'){
                        printColoredMessage(
                            message = '- Finding MNN using both CCA and PCA coordinates.',
                            color = 'orange',
                            verbose = verbose
                        )
                        all.pp <- matrix(1)
                        while(nrow(all.pp) < min.ps){
                            ### Applying MNN ####
                            cca.mnn <- findMutualNN(
                                data1 = temp.cca[row.names(sample.annot.a) , ],
                                data2 = temp.cca[row.names(sample.annot.b) , ],
                                k1 = min.nb.for.mnn
                            )
                            ### Creating MNN - CCA data and adding samples ids ####
                            cca.mnn <- data.frame(
                                left.index = cca.mnn$first,
                                right.index = cca.mnn$second,
                                left.samples = row.names(sample.annot.a)[cca.mnn$first],
                                right.samples = row.names(sample.annot.b)[cca.mnn$second]
                            )
                            ## Finding MNN using the PCA coordinates ####
                            printColoredMessage(
                                message = '- Finding MNN using the PCA coordinates.',
                                color = 'orange',
                                verbose = verbose
                            )
                            ### Applying normalization, regression an log on the data ####
                            temp.data <- preProcessData(
                                se.obj = se.obj[ , all.samples],
                                assay.name = assay.name,
                                normalization = normalization,
                                regress.out.variables = regress.out.variables,
                                apply.log = apply.log,
                                pseudo.count = pseudo.count,
                                check.se.obj = FALSE,
                                remove.na = 'none',
                                verbose = verbose
                            )
                            if (!identical(colnames(temp.data), all.samples)){
                                stop('The order of samples in the calculated CCA have issues.')
                            }
                            if (isTRUE(cosine.norm)){
                                temp.data <- cosineNorm(x = temp.data, mode = "matrix")
                            }
                            ### Performing PCA on the data ####
                            pca.data <- irlba::prcomp_irlba(
                                x = t(temp.data[hvg , ]),
                                center = TRUE,
                                scale = FALSE,
                                n = nb.pca
                            )
                            ### Applying MNN ####
                            row.names(pca.data$x) <- colnames(temp.data)
                            pca.mnn <- findMutualNN(
                                data1 = pca.data$x[row.names(sample.annot.a) , ],
                                data2 = pca.data$x[row.names(sample.annot.b) , ],
                                k1  = min.nb.for.mnn
                            )
                            ### Creating MNN - PCA data and adding samples ids ####
                            pca.mnn <- data.frame(
                                left.index = pca.mnn$first,
                                right.index = pca.mnn$second,
                                left.samples = row.names(sample.annot.a)[pca.mnn$first],
                                right.samples = row.names(sample.annot.b)[pca.mnn$second]
                            )
                            ## Performing RUV-III normalization ####
                            ### Finding common PS across batches ####
                            all.pp <- c(
                                paste(pca.mnn$left.index, pca.mnn$right.index, sep = '_'),
                                paste(cca.mnn$left.index, cca.mnn$right.index, sep = '_')
                            )
                            all.pp <- names(which(table(all.pp) > 1))
                            all.pp <- do.call(rbind, strsplit(all.pp, "_"))
                            all.pp <- data.frame(
                                left.index = as.numeric(all.pp[,1]),
                                right.index = as.numeric(all.pp[,2])
                            )
                            all.pp$left.samples <- row.names(sample.annot.a)[all.pp$left.index]
                            all.pp$right.samples <- row.names(sample.annot.b)[all.pp$right.index]
                            min.nb.for.mnn <- min.nb.for.mnn + 1
                        }
                    }
                    ### Applying RUV-III
                    if (isTRUE(apply.ruviii.norm)){
                        ### Creating M matrix ####
                        tr <- all.samples
                        for(j in 1:nrow(all.pp)){
                            tr[tr == all.pp$left.samples[j]] <- tr[tr == all.pp$right.samples[j]]
                        }
                        if (sum(table(tr)> 1) == 0){
                            stop('There are some issues with the applying RUV-III normalization.')
                        }
                        if (is.null(ncg)){
                            ncg <- row.names(se.obj) %in% row.names(se.obj)
                        }
                        m.matrix <- ruv::replicate.matrix(a = tr)
                        temp.data <- applyLog(
                            se.obj = se.obj[ , all.samples],
                            assay.names = assay.name,
                            pseudo.count = pseudo.count,
                            check.se.obj = FALSE,
                            remove.na = 'none',
                            verbose = verbose
                        )[[assay.name]]
                        Y <- t(temp.data)
                        Y.stand <- scale(
                            x = Y,
                            center = TRUE,
                            scale = FALSE
                        )
                        Y0 <- optimized_function(Y, m.matrix)
                        left.sing.value <- BiocSingular::runSVD(
                            x = Y0,
                            k = k,
                            BSPARAM = bsparam(),
                            center = FALSE,
                            scale = FALSE
                        )$u
                        alpha <- t(left.sing.value[, 1:k, drop = FALSE]) %*% Y
                        ac <- alpha[, ncg, drop = FALSE]
                        W <- Y.stand[, ncg] %*% t(ac) %*% solve(ac %*% t(ac))
                        ruv.adj <- Y - W %*% alpha
                    }
                    if (isTRUE(use.ruviii.for.mnn)){
                        min.nb.for.mnn <- min.nb.for.mnn.initial
                        all.pp.new <- matrix(1)
                        while(nrow(all.pp.new) < min.ps){
                            ruv.mnn <- findMutualNN(
                                data1 = ruv.adj[row.names(sample.annot.a) , ],
                                data2 = ruv.adj[row.names(sample.annot.b) , ],
                                k1  = min.nb.for.mnn
                            )
                            ### Creating MNN - CCA data and adding samples ids ####
                            ruv.mnn <- data.frame(
                                left.index = ruv.mnn$first,
                                right.index = ruv.mnn$second,
                                left.samples = row.names(sample.annot.a)[ruv.mnn$first],
                                right.samples = row.names(sample.annot.b)[ruv.mnn$second]
                            )
                            all.pp.new <- c(
                                paste(ruv.mnn$left.index, ruv.mnn$right.index, sep = '_'),
                                paste(all.pp$left.index, all.pp$right.index, sep = '_')
                            )
                            all.pp.new <- names(which(table(all.pp.new) > 1))
                            all.pp.new <- do.call(rbind, strsplit(all.pp.new, "_"))
                            all.pp.new <- data.frame(
                                left.index = as.numeric(all.pp.new[,1]),
                                right.index = as.numeric(all.pp.new[,2])
                            )
                            all.pp.new$left.samples <- row.names(sample.annot.a)[all.pp.new$left.index]
                            all.pp.new$right.samples <- row.names(sample.annot.b)[all.pp.new$right.index]
                            min.nb.for.mnn <- min.nb.for.mnn + 1
                        }
                        all.pp <- all.pp.new
                    }
                    if (use.ruviii == 'similarity'){
                        similarity.data <- t(ruv.adj)
                    } else {
                        similarity.data <- preProcessData(
                            se.obj = se.obj[, all.samples],
                            assay.name = assay.name,
                            normalization = normalization,
                            regress.out.variables = regress.out.variables,
                            apply.log = apply.log,
                            pseudo.count = pseudo.count,
                            check.se.obj = FALSE,
                            remove.na = 'none',
                            verbose = verbose)
                    }

                    ## Finding the most similar samples ####
                    ### Applying KNN ####
                    #### Applying KNN for first batch ####
                    min.sample.for.ps.initial <- min.sample.for.ps
                    if (similarity.approach == 'euclidean'){
                        min.sample.for.ps <- min.sample.for.ps - 1
                        distance.a <- BiocNeighbors::findKNN(
                            X = t(similarity.data[hvg, row.names(sample.annot.a) , drop = FALSE]),
                            k = min.sample.for.ps
                        )
                        distance.a.index <- as.data.frame(distance.a$index)
                        colnames(distance.a.index) <- paste0('V', 2:(min.sample.for.ps + 1))
                        distance.a.index$V1 <- as.numeric(row.names(distance.a.index))
                        distance.a.index <- distance.a.index[order(colnames(distance.a.index))]
                        for(i in 1:ncol(distance.a.index)){
                            col.name <- paste0('Sample', i)
                            distance.a.index[col.name] <- row.names(sample.annot.a)[distance.a.index[ , i]]
                        }
                        #### Applying KNN for second batch ####
                        distance.b <- BiocNeighbors::findKNN(
                            X = t(similarity.data[hvg, row.names(sample.annot.b) , drop = FALSE]),
                            k = min.sample.for.ps
                        )
                        distance.b.index <- as.data.frame(distance.b$index)
                        colnames(distance.b.index) <- paste0('V', 2:(min.sample.for.ps + 1))
                        distance.b.index$V1 <- as.numeric(row.names(distance.b.index))
                        distance.b.index <- distance.b.index[order(colnames(distance.b.index))]
                        for(i in 1:ncol(distance.b.index)){
                            col.name <- paste0('Sample', i)
                            distance.b.index[col.name] <- row.names(sample.annot.b)[distance.b.index[ , i]]
                        }
                    }
                    ### Applying correlation analysis ####
                    if (similarity.approach == 'correlation'){
                        #### Applying correlation for first batch ####
                        min.sample.for.ps <- min.sample.for.ps - 1
                        cor.matrix.a <- cor(
                            x = similarity.data[hvg, row.names(sample.annot.a), drop = FALSE ],
                            method = "pearson"
                        )
                        cor.matrix.a <- round(x = cor.matrix.a, digits = 3)
                        distance.a.index <- lapply(
                            1:ncol(cor.matrix.a),
                            function(i) {
                                sample_name <- colnames(cor.matrix.a)[i]
                                cor.values <- cor.matrix.a[ , i]
                                cor.values <- cor.values[-i]
                                sorted <- sort(cor.values, decreasing = TRUE)
                                top.samples <- head(sorted, min.sample.for.ps)
                                temp.df <- data.frame(
                                    Sample1 = sample_name,
                                    TopMatch = names(top.samples),
                                    Correlation = mean(as.numeric(top.samples))
                                )
                                temp.df$TopMatchID <- c(2:c(min.sample.for.ps + 1))
                                temp.df <- temp.df %>%
                                    tidyr::pivot_wider(
                                        names_from = TopMatchID,
                                        values_from = TopMatch,
                                        names_prefix = "Sample") %>%
                                    data.frame()
                            })
                        distance.a.index <- do.call(rbind, distance.a.index)
                        distance.a.index <- distance.a.index[ , c(1, c(3:c(min.sample.for.ps + 2)), 2)]
                        if (isFALSE(all.equal(distance.a.index$Sample1, row.names(sample.annot.a)))){
                            stop('There are something wrong with the ...')
                        }
                        #### Applying correlation for first batch ####
                        cor.matrix.b <- cor(
                            x = similarity.data[hvg, row.names(sample.annot.b), drop = FALSE ],
                            method = "pearson"
                        )
                        cor.matrix.b <- round(cor.matrix.b, 3)
                        distance.b.index <- lapply(
                            1:ncol(cor.matrix.b),
                            function(i) {
                                sample_name <- colnames(cor.matrix.b)[i]
                                cor.values <- cor.matrix.b[ , i]
                                cor.values <- cor.values[-i]
                                sorted <- sort(cor.values, decreasing = TRUE)
                                top.samples <- head(sorted, min.sample.for.ps)
                                temp.df <- data.frame(
                                    Sample1 = sample_name,
                                    TopMatch = names(top.samples),
                                    Correlation = mean(as.numeric(top.samples))
                                )
                                temp.df$TopMatchID <-  c(2:c(min.sample.for.ps + 1))
                                temp.df <- temp.df %>%
                                    tidyr::pivot_wider(
                                        names_from = TopMatchID,
                                        values_from = TopMatch,
                                        names_prefix = "Sample") %>%
                                    data.frame()
                            })
                        distance.b.index <- do.call(rbind, distance.b.index)
                        distance.b.index <- distance.b.index[ , c(1, c(3:c(min.sample.for.ps + 2)), 2)]
                        if (isFALSE(all.equal(distance.b.index$Sample1, row.names(sample.annot.b)))){
                            stop('There are something wrong with the ...')
                        }
                    }
                    ### Applying cosine analysis ####
                    if (similarity.approach == 'cosine'){
                        #### Applying correlation for first batch ####
                        min.sample.for.ps <- min.sample.for.ps - 1
                        cor.matrix.a <- lsa::cosine(
                            x = similarity.data[hvg, row.names(sample.annot.a), drop = FALSE]
                        )
                        cor.matrix.a <- round(x = cor.matrix.a, digits = 3)
                        distance.a.index <- lapply(
                            1:ncol(cor.matrix.a),
                            function(i) {
                                sample_name <- colnames(cor.matrix.a)[i]
                                cor.values <- cor.matrix.a[ , i]
                                cor.values <- cor.values[-i]
                                sorted <- sort(cor.values, decreasing = TRUE)
                                top.samples <- head(sorted, min.sample.for.ps)
                                temp.df <- data.frame(
                                    Sample1 = sample_name,
                                    TopMatch = names(top.samples),
                                    Correlation = mean(as.numeric(top.samples))
                                )
                                temp.df$TopMatchID <- c(2:c(min.sample.for.ps + 1))
                                temp.df <- temp.df %>%
                                    tidyr::pivot_wider(
                                        names_from = TopMatchID,
                                        values_from = TopMatch,
                                        names_prefix = "Sample") %>%
                                    data.frame()
                            })
                        distance.a.index <- do.call(rbind, distance.a.index)
                        distance.a.index <- distance.a.index[ , c(1, c(3:c(min.sample.for.ps + 2)), 2)]
                        if (isFALSE(all.equal(distance.a.index$Sample1, row.names(sample.annot.a)))){
                            stop('There are something wrong with the ...')
                        }
                        #### Applying correlation for first batch ####
                        cor.matrix.b <- lsa::cosine(
                            x = similarity.data[hvg, row.names(sample.annot.b), drop = FALSE]
                        )
                        cor.matrix.b <- round(cor.matrix.b, 3)
                        distance.b.index <- lapply(
                            1:ncol(cor.matrix.b),
                            function(i) {
                                sample_name <- colnames(cor.matrix.b)[i]
                                cor.values <- cor.matrix.b[ , i]
                                cor.values <- cor.values[-i]
                                sorted <- sort(cor.values, decreasing = TRUE)
                                top.samples <- head(sorted, min.sample.for.ps)
                                temp.df <- data.frame(
                                    Sample1 = sample_name,
                                    TopMatch = names(top.samples),
                                    Correlation = mean(as.numeric(top.samples))
                                )
                                temp.df$TopMatchID <-  c(2:c(min.sample.for.ps + 1))
                                temp.df <- temp.df %>%
                                    tidyr::pivot_wider(
                                        names_from = TopMatchID,
                                        values_from = TopMatch,
                                        names_prefix = "Sample") %>%
                                    data.frame()
                            })
                        distance.b.index <- do.call(rbind, distance.b.index)
                        distance.b.index <- distance.b.index[ , c(1, c(3:c(min.sample.for.ps + 2)), 2)]
                        if (isFALSE(all.equal(distance.b.index$Sample1, row.names(sample.annot.b)))){
                            stop('There are something wrong with the ...')
                        }
                    }
                    # Creating PRPS data ####
                    if (isTRUE(apply.log.for.prps)){
                        prps.data <- applyLog(
                            se.obj = se.obj,
                            assay.names = assay.name,
                            pseudo.count = pseudo.count,
                            check.se.obj = check.se.obj,
                            remove.na = 'none',
                            verbose = verbose
                        )[[assay.name]]
                    } else {
                        prps.data <- assay(x = se.obj, i = assay.name)
                    }
                    prps.data <- lapply(
                        1:nrow(all.pp),
                        function(ps){
                            selected.cols <- paste0('Sample', c(1:min.sample.for.ps.initial))
                            ps1 <- unlist(unname(distance.a.index[distance.a.index$Sample1 == all.pp$left.samples[ps] , selected.cols ]))
                            ps1 <- rowMeans(prps.data[ , ps1])
                            ps2 <- unlist(unname(distance.b.index[distance.b.index$Sample1 == all.pp$right.samples[ps] , selected.cols ]))
                            ps2 <- rowMeans(prps.data[ , ps2])
                            prps.set <- data.frame(ps1 = ps1 , ps2 = ps2)
                            colnames(prps.set) <- rep(x = paste0(d, '.sample', ps, x), each = 2)
                            prps.set
                        })
                    prps.data <- do.call(cbind , prps.data)
                    prps.data
                })
            all.prps.data <- do.call(cbind, all.prps.data)
        })
    all.prps.sets <- do.call(cbind, all.prps.sets)
    return(all.prps.sets)
}



# ruv.data <- cbind(log2(assay(read.se.obj, 'RawCount') + 1), all.prps.sets)
# k.vals = 2
# M <- ruv::replicate.matrix(a = colnames(ruv.data))
# dim(M)
# Y <- t(log2(ruv.data + 1))
# # Y <- t(ruv.data)
# Y.stand <- scale(Y, center = TRUE, scale = FALSE)
# Y0 <- optimized_function(Y, M)
# left.sing.value <- BiocSingular::runSVD(
#     x = Y0,
#     k = max(k.vals),
#     BSPARAM = bsparam(),
#     center = FALSE,
#     scale = FALSE
# )$u
# hvg <- singscore::getStableGenes(n_stable = 2000)
# hvg <- row.names(se.obj) %in% hvg
# alpha <- t(left.sing.value[, 1:k.vals, drop = FALSE]) %*% Y
# ac <- alpha[, hvg, drop = FALSE]
# W <- Y.stand[, hvg] %*% t(ac) %*% solve(ac %*% t(ac))
# newY <- Y - W %*% alpha
# newY <- t(newY[1:ncol(se.obj) ,])
#
# p <- prcomp(t(newY[ , ]))
# pairs(p$x[,1:3], col = factor(read.se.obj$Quantile.CMS.TI))
# pairs(p$x[,1:3], col = factor(read.se.obj$Time.interval))
# plot(W[ 1:ncol(se.obj) , 1])
# plot(W[ 1:ncol(se.obj) , 2])
# plot(W[ 1:ncol(se.obj) , 3])
# plot(W[ 1:ncol(se.obj) , 4])
# plot(p$x[,1], read.se.obj$Tumour.purity)
#
# newY <- p$x[,1:40]
# cca.mnn <- findMutualNN(
#     data1 = newY[read.se.obj$Time.interval == 'Batch1' , ],
#     data2 = newY[read.se.obj$Time.interval == 'Batch2' , ],
#     k1  = 1
# )
# cca.mnn <- data.frame(left = cca.mnn$first, right = cca.mnn$second + 72)
# cca.mnn$bio1 <- read.se.obj$Quantile.CMS.TI[cca.mnn$left]
# cca.mnn$bio2 <- read.se.obj$Quantile.CMS.TI[cca.mnn$right]
#
# cca.mnn <- cca.mnn[cca.mnn$bio1!= 'Not.classified', ]
# cca.mnn <- cca.mnn[cca.mnn$bio2!= 'Not.classified', ]
# cca.mnn <- cca.mnn[cca.mnn$bio1!= 'Normal', ]
# cca.mnn <- cca.mnn[cca.mnn$bio2!= 'Normal', ]
#
# table(apply(cca.mnn[ , 3:4], 1, function(x) length(unique(x))))
#
#
# knn.data <- BiocNeighbors::findKNN(X = t(newY[hvg , ]), k = 3)
# knn.data <- as.data.frame(knn.data$index)
# knn.data$bio1 <- read.se.obj$Quantile.CMS.TI
# knn.data$bio2 <- read.se.obj$Quantile.CMS.TI[knn.data$V1]
# knn.data$bio3 <- read.se.obj$Quantile.CMS.TI[knn.data$V2]
# knn.data$bio4 <- read.se.obj$Quantile.CMS.TI[knn.data$V3]
# table(apply(knn.data[ , 4:6], 1, function(x) length(unique(x))))
#
#
# ll <- mclust::Mclust(data = t(newY[ , ]))
# pairs(p$x[,1:3], col = factor(ll$classification))
# table(ll$classification, read.se.obj$Quantile.CMS.TI)
#
#
#
# kk <- kmeans(x = all.cca, centers = 6)
# table(kk$cluster, read.se.obj$Quantile.CMS.TI)
#
# cca.data.1 = all.cca[read.se.obj$Time.interval == 'Batch1' , ]
# cca.data.2 = all.cca[read.se.obj$Time.interval == 'Batch2' , ]
# # Normalize each row vector to unit length
# normalize <- function(x) x / sqrt(sum(x^2))
#
# cca.data.1.norm <- t(apply(cca.data.1, 1, normalize))
# cca.data.2.norm <- t(apply(cca.data.2, 1, normalize))
#
# # Compute cosine similarity matrix
# similarity.matrix <- cca.data.1.norm %*% t(cca.data.2.norm)  # Rows = data.1 samples, Cols = data.2 samples
# # For each sample in data.1, find the most similar in data.2
# most.similar.pairs <- apply(similarity.matrix, 1, function(x) {
#     row.names(cca.data.2)[which.max(x)]
# })
#
# # Optionally, get similarity scores
# similarity.scores <- apply(similarity.matrix, 1, max)
#
# # Combine into a data.frame
# similar.sample.pairs <- data.frame(
#     sample.1 = rownames(cca.data.1),
#     sample.2 = most.similar.pairs,
#     similarity = similarity.scores,
#     stringsAsFactors = FALSE
# )
# similar.sample.pairs <- similar.sample.pairs[order(-similar.sample.pairs$similarity) , ]
# for(i in 1:nrow(similar.sample.pairs)){
#     a <- read.se.obj$Quantile.CMS.TI[colnames(read.se.obj) == similar.sample.pairs$sample.1[i]]
#     b <- read.se.obj$Quantile.CMS.TI[colnames(read.se.obj) == similar.sample.pairs$sample.2[i]]
#     print(c(a, b))
# }
#
#
#
#
#
#
#
#
#
#
# head(all.cca@metadata$CCA$RawCount$Time.interval$`Cca:10|HVG:14799|Scale:TRUE|Norm:CPM`$Batch1_Batch2)
# mm <- computeGenesVariableAnova(se.obj = se.obj, assay.names = 'FPKM.UQ', variable = 'Quantile.CMS.TI')
# mm <- mm@metadata$Metrics$FPKM.UQ$gene.level$ANOVA$aov$Quantile.CMS.TI$fstatistics.pvalues
# mm <- mm[order(-mm$statistic) , ]
# all.paired.batches <- combn(x = unique(se.obj[[main.uv.variable]]), m = 2)
# all.mnn <- lapply(
#     1:ncol(all.paired.batches),
#     function(x){
#         index.a <- se.obj[[main.uv.variable]] == all.paired.batches[1, x]
#         data.a.a <- log2(assay(x = se.obj[ , index.a], i = assay.name) + 1)
#         data.a <- CreateSeuratObject(
#             counts = data.a.a
#         )
#         data.a <- SetAssayData(
#             object = data.a,
#             layer = "data",
#             new.data = as.matrix(data.a.a)
#         )
#         data.a <- ScaleData(object = data.a)
#         data.a <- FindVariableFeatures(object = data.a, nfeatures = 1000)
#
#         index.b <- se.obj[[main.uv.variable]] == all.paired.batches[2, x]
#         data.b.b <- log2(assay(x = se.obj[ , index.b], i = assay.name) + 1)
#         data.b <- CreateSeuratObject(
#             counts = data.b.b
#         )
#         data.b <- SetAssayData(
#             object = data.b,
#             layer = "data",
#             new.data = as.matrix(data.b.b)
#         )
#         data.b <- ScaleData(data.b)
#         data.b <- FindVariableFeatures(object = data.b, nfeatures = 1000)
#
#         cca.out <- RunCCA(
#             object1 = data.a,
#             object2 = data.b,
#             features = NULL,
#             num.cc = 71
#         )
#         cca.embeddings <- Embeddings(object = cca.out, reduction = "cca")
#         pairs(cca.embeddings[, 1:3], col = factor(read.se.obj$Time.interval))
#         pairs(cca.embeddings[, 1:3], col = factor(read.se.obj$Quantile.CMS.TI))
#
#         unsuper.cc <- mclust::Mclust(data = cca.embeddings[, 1:2], G = 5:15)
#         pairs(cca.embeddings[, 1:3], col = factor(unsuper.cc$classification))
#
#         pairs(cca.embeddings[, 1:3], col = factor(ff))
#         table(ff, read.se.obj$Quantile.CMS.TI)
#         pairs(cca.embeddings[, 1:3], col = factor(rowSums(mm)))
#
#
#
#         all.anchors <- Seurat::FindIntegrationAnchors(
#             object.list = list(data.a = data.b, data.b = data.b) , reduction = 'cca',dims = 1:10
#         )
#         anchor <- all.anchors@anchors
#         anchor.mnn <- data.frame(left = anchor$cell1, right = anchor$cell2 + ncol(data.a))
#         anchor.mnn$bio1 <- read.se.obj$Quantile.CMS.TI[anchor.mnn$left]
#         anchor.mnn$bio2 <- read.se.obj$Quantile.CMS.TI[anchor.mnn$right]
#         table(apply(anchor.mnn[ , 3:4], 1, function(x) length(unique(x))))
#
#
#         library(cluster)
#         data <- as.danchorsdata <- as.data.frame(cca.embeddings[, 1:2])
#         fanny.result <- fanny(data, k = 6)
#         mm <- fanny.result$membership > .3
#         rowSums(mm)
#
#         ff <- max.col(fanny.result$membership)
#
#         library(ggplot2)
#
#         data$Cluster1 <- fanny.result$membership[,1]
#         data$Cluster2 <- fanny.result$membership[,2]
#         data$Cluster3 <- fanny.result$membership[,3]
#
#         # Plot membership to Cluster 1 vs. Cluster 2
#         ggplot(data, aes(x = Cluster1, y = Cluster2)) +
#             geom_point() +
#             labs(title = "Soft Cluster Memberships")
#
#
#         cca.mnn <- findMutualNN(
#             data1 = cca.embeddings[colnames(data.a) , ],
#             data2 = cca.embeddings[colnames(data.b),  ],
#             k1  = 3
#         )
#         cca.mnn <- data.frame(left = cca.mnn$first, right = cca.mnn$second + ncol(data.a))
#         cca.mnn$bio1 <- read.se.obj$Quantile.CMS.TI[cca.mnn$left]
#         cca.mnn$bio2 <- read.se.obj$Quantile.CMS.TI[cca.mnn$right]
#         table(apply(cca.mnn[ , 3:4], 1, function(x) length(unique(x))))
#
#
#         kk <- kmeans(x = cca.embeddings[, 1:3], centers = 5)
#         table(read.se.obj$Quantile.CMS.TI[kk$cluster == 1])
#         table(read.se.obj$Quantile.CMS.TI[kk$cluster == 2])
#         table(read.se.obj$Quantile.CMS.TI[kk$cluster == 3])
#         table(read.se.obj$Quantile.CMS.TI[kk$cluster == 4])
#
#
#         cca.mnn$vv <- vv
#         cca.mnn <- cca.mnn[order(cca.mnn$vv) , ]
#
#         full.data <- log2(assay(x = se.obj, i = assay.name) + 1)
#         full.sample <- as.data.frame(colData(se.obj))
#         all <- lapply(
#             1:nrow(cca.mnn),
#             function(x){
#                 temp.data <- data.frame(
#                     s1 = full.data[ , cca.mnn$left[x]],
#                     s2 =  full.data[ , cca.mnn$right[x]]
#                 )
#                 no <- paste0(cca.mnn$left[x], '_', cca.mnn$right[x])
#                 colnames(temp.data) <- rep(paste0('sample', no), 2)
#                 return(temp.data)
#             })
#         sample.infor <- unlist(lapply(
#             1:nrow(cca.mnn),
#             function(x){
#                 c(cca.mnn$left[x], cca.mnn$right[x])
#             }))
#         sample.infor <- full.sample[sample.infor , ]
#         dim(sample.infor)
#
#         all <- do.call(cbind, all)
#         M <- replicate.matrix(colnames(all))
#         rr <- ruv::residop(A = t(all), B = M)
#         vv <- rowVars(rr)
#         vv <- vv[seq(1, nrow(rr), 2)]
#         p <- prcomp(x = rr)
#         pairs(p$x[,1:5], col = factor(sample.infor$Time.interval))
#         pairs(p$x[,1:5], col = factor(sample.infor$Quantile.CMS.TI))
#         boxplot(p$x[,7]~ sample.infor$Quantile.CMS.TI)
#
#         install.packages('factoextra')
#         rn <- row.names(rr)
#         row.names(rr) <- paste0('sample', 1:nrow(rr))
#         p <- prcomp(x = rr)
#         ccp.ca <- fviz_pca_ind(p, contrib = "yes")
#         fviz_pca_ind(p, axes = c(1, 3), repel = TRUE)
#         cca.mnn
#         library(factoextra)
#
#         # Extract contributions of individuals (samples)
#         sample_contrib <- get_pca_ind(p)$contrib  # rows = samples, columns = PCs
#
#         # Order by contribution to PC1
#         ordered_samples <- rownames(sample_contrib)[order(sample_contrib[, 1], decreasing = TRUE)]
#         ordered_samples[1:10]
#         rn[163]
#
#
#     })
#
#
#
#
#
#
# #
# # ### CCA analysis
# # library(Seurat)
# # data.a.cos <- cosineNorm(x = data.a)
# # data.b.cos <- cosineNorm(x = data.b)
# #
# # s1 <- CreateSeuratObject(counts = data.a[ , ])
# # s2 <- CreateSeuratObject(counts = data.b[, ])
# # s1 <- SetAssayData(s1, layer = "data", new.data = as.matrix(data.a))
# # s2 <- SetAssayData(s2, layer = "data", new.data = as.matrix(data.b))
# # s1 <- FindVariableFeatures(s1)
# # s2 <- FindVariableFeatures(s2)
# # s1 <- ScaleData(s1)
# # s2 <- ScaleData(s2)
# # pp <- intersect(VariableFeatures(s1), VariableFeatures(s2))
# # cca_out <- RunCCA(s1, s2, features = row.names(data.a), num.cc = 100)
# # cca_embeddings <- Embeddings(cca_out, reduction = "cca")
# # pairs(cca_embeddings[, 1:3], col = factor(all.sample$pam50))
# # pairs(cca_embeddings[, 1:3], col = factor(all.sample$studies))
# # cca.mnn <- findMutualNN(data1 = cca_embeddings[colnames(data.a) , ], data2 = cca_embeddings[colnames(data.b),  ], k1  = 3)
# # cca.mnn <- data.frame(left = cca.mnn$first, right = cca.mnn$second + 221)
# # cca.mnn$bio1 <- all.sample$pam50[cca.mnn$left]
# # cca.mnn$bio2 <- all.sample$pam50[cca.mnn$right]
# # table(apply(cca.mnn[ , 3:4], 1, function(x) length(unique(x))))
# #
# # pp <- as.data.frame(cca.mnn@metadata$merge.info$pairs@listData)
# # pp$bio1 <- all.sample$pam50[pp$left]
# # pp$bio2 <- all.sample$pam50[pp$right]
# # table(apply(pp[ , 3:4], 1, function(x) length(unique(x))))
#
# library(limma)
# count.data <- assay(x = read.se.obj, i = 'RawCount')
# logCPM <- cpm(count.data, log = TRUE, prior.count = 1)
# logCPM_corrected <- removeBatchEffect(logCPM, batch = read.se.obj$Time.interval)
# p <- prcomp(t(logCPM_corrected))
# pairs(p$x[,1:3], col = factor(read.se.obj$Quantile.CMS.TI))
#
#
# cor_matrix <- cor(logCPM_corrected[ 1:20, ], method = "pearson")
# dist_matrix <- as.dist(1 - cor_matrix)
# samples_batch1 <- colnames(logCPM_corrected)[read.se.obj$Time.interval == "Batch1"]
# samples_batch2 <- colnames(logCPM_corrected)[read.se.obj$Time.interval == "Batch2"]
#
# # For each sample in batch1, find the most similar sample in batch2
# similar_pairs <- lapply(samples_batch1, function(s1) {
#     sims <- cor_matrix[s1, samples_batch2]
#     best_match <- names(sort(sims, decreasing = TRUE))[1]  # most similar
#     return(data.frame(sample1 = s1, sample2 = best_match, correlation = sims[best_match]))
# })
# similar_pairs_df <- do.call(rbind, similar_pairs)
# similar_pairs_df$bio1 <- read.se.obj$Quantile.CMS.TI[match( similar_pairs_df$sample1, colnames(read.se.obj))]
# similar_pairs_df$bio2 <- read.se.obj$Quantile.CMS.TI[match( similar_pairs_df$sample2, colnames(read.se.obj))]
# similar_pairs_df <- similar_pairs_df[order(-similar_pairs_df$correlation) , ]
#
# table(apply(similar_pairs_df[ , 4:5], 1, function(x) length(unique(x))))
#
# eread.se.obj$Quantile.CMS.TI[colnames(read.se.obj) == 'TCGA-AG-3592-01A-02R-1736-07']
#
# # For batch1
# top3_batch1 <- lapply(samples_batch1, function(s) {
#     sims <- cor_matrix[s, samples_batch1]
#     sims <- sims[names(sims) != s]  # exclude self
#     top3 <- sort(sims, decreasing = TRUE)[1:3]
#     data.frame(reference = s, top_similar = names(top3), correlation = top3)
# })
#
# # For batch2
# top3_batch2 <- lapply(samples_batch2, function(s) {
#     sims <- cor_matrix[s, samples_batch2]
#     sims <- sims[names(sims) != s]
#     top3 <- sort(sims, decreasing = TRUE)[1:3]
#     data.frame(reference = s, top_similar = names(top3), correlation = top3)
# })
#
# top3_batch1_df <- do.call(rbind, top3_batch1)
# top3_batch2_df <- do.call(rbind, top3_batch2)
#
#
# library(ggplot2)
# library(uwot)
#
# pca <- prcomp(t(logCPM_corrected), scale. = TRUE)
# df_pca <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], batch = read.se.obj$Time.interval)
#
# ggplot(df_pca, aes(PC1, PC2, color = batch)) +
#     geom_point(size = 3) +
#     theme_minimal()
#
#
#
# library(batchelor)
#
# # Assume: logCPM matrix or log-normalized matrix per batch
# batch1 <- logCPM[hvg, read.se.obj$Time.interval == "Batch1"]
# batch2 <- logCPM[hvg, read.se.obj$Time.interval == "Batch2"]
#
# # Apply MNN correction
# mnn_out <- mnnCorrect(batch1, batch2, k = 2, cos.norm.in = T, cos.norm.out = T)
#
# m <- as.data.frame(mnn_out@metadata$merge.info@listData$pairs@listData)
# m$bio1 <- read.se.obj$Quantile.CMS.TI[m$left]
# m$bio2 <- read.se.obj$Quantile.CMS.TI[m$right]
# table(apply(m[ , 3:4], 1, function(x) length(unique(x))))
#
# # Combined corrected matrix
# logCPM_mnn <- assay(mnn_out)
#
# p <- prcomp(t(logCPM_mnn))
# pairs(p$x[,1:3], col = factor(read.se.obj$Quantile.CMS.TI))
#
# read.se.obj$Quantile.CMS.TI[158]
#
#
# all.data <- cbind(batch1, batch2)
# all.data.cos <- cosineNorm(x = all.data, mode = "matrix")
# pca.data <- irlba::prcomp_irlba(x = t(all.data.cos), n = 177, center = T, scale. = F)
# row.names(pca.data$x) <- colnames(all.data)
# mnn.mnn <- findMutualNN(data1 = pca.data$x[colnames(data.a) , ] , pca.data$x[colnames(data.b) , ], k1 = 3)
# mnn.mnn <- data.frame(left = mnn.mnn$first, right = mnn.mnn$second + 72)
# mnn.mnn$bio1 <- all.sample$pam50[mnn.mnn$left]
# mnn.mnn$bio2 <- all.sample$pam50[mnn.mnn$right]
# table(apply(mnn.mnn[ , 3:4], 1, function(x) length(unique(x))))
#
#
#
#
#
#
#
#
# if (is.numeric(se.obj[[main.uv.variable]])){
#     ## Keeping original values of the main unwanted variable ####
#     initial.variable <- se.obj[[main.uv.variable]]
#     initial.variable2 <- se.obj[[main.uv.variable]]
#     initial.variable3 <- se.obj[[main.uv.variable]]
#
#     ## Grouping main the main unwanted variable ####
#     se.obj[[main.uv.variable]] <- groupContinuousVariable(
#         se.obj = se.obj,
#         variable = main.uv.variable,
#         nb.clusters = nb.clusters,
#         clustering.method = clustering.method,
#         perfix = '.',
#         verbose = verbose
#     )
#     ## Selecting subgroups of the main unwanted variable with highest and lowest values ####
#     if (isTRUE(select.extreme.groups)){
#         printColoredMessage(
#             message = paste0(
#                 '- Selecting the two subgroups of the ',
#                 main.uv.variable,
#                 ' variable with highest and lowest values.'),
#             color = 'blue',
#             verbose = verbose
#         )
#         initial.se.obj <- se.obj
#         max.group <- se.obj[[main.uv.variable]][initial.variable == max(initial.variable)]
#         min.group <- se.obj[[main.uv.variable]][initial.variable == min(initial.variable)]
#         selected.samples <- se.obj[[main.uv.variable]] %in% c(max.group, min.group)
#         se.obj <- se.obj[ , selected.samples]
#         se.obj[[main.uv.variable]] <- droplevels(se.obj[[main.uv.variable]])
#         initial.variable <- initial.variable[selected.samples]
#     }
#     if (isFALSE(select.extreme.groups)) initial.se.obj <- se.obj
# }
# if (!is.numeric(se.obj[[main.uv.variable]])){
#     if (length(unique(initial.variable)) == 1){
#         stop('To create MNN, the "main.uv.variable" must have at least two groups/levels.')
#     }
#     if (length(unique(initial.variable)) > 1){
#         printColoredMessage(
#             message = paste0(
#                 '- The "',
#                 main.uv.variable,
#                 '" is a categorical variable with ',
#                 length(unique(se.obj[[main.uv.variable]])),
#                 ' levels.'),
#             color = 'blue',
#             verbose = verbose
#         )
#         se.obj[[main.uv.variable]] <- factor(x = se.obj[[main.uv.variable]])
#     }
# }
#
# ## Checking the sample size of each group in the variable ####
# subgroups.size <- findRepeatingPatterns(
#     vec = se.obj[[main.uv.variable]],
#     n.repeat = max(min.sample.for.ps, nb.mnn)
# )
# if (min.batches.to.cover == 'all') {
#     if (length(subgroups.size) != length(unique(se.obj[[main.uv.variable]])) ){
#         stop(paste0(
#             'Some sub-groups of the variable "',
#             main.uv.variable,
#             '" have less than ',
#             max(min.sample.for.ps, nb.mnn),
#             ' samples. Then, MNN cannot be created across all batches.')
#         )
#     }
# }
# if (is.numeric(min.batches.to.cover)) {
#     if (length(subgroups.size) >= min.batches.to.cover){
#         printColoredMessage(
#             message = paste0(
#                 '- At least ',
#                 min.batches.to.cover,
#                 ' sub-groups of the variable ',
#                 main.uv.variable,
#                 ' have at least ',
#                 max(min.sample.for.ps, nb.mnn),
#                 ' samples.'),
#             color = 'blue',
#             verbose = verbose
#         )
#     } else {
#         stop(paste0(
#             'Some sub-groups of the variable "',
#             main.uv.variable,
#             '" have less than ',
#             max(min.sample.for.ps, nb.mnn),
#             ' samples. Then, MNN cannot be created across all batches.')
#         )
#     }
# }
