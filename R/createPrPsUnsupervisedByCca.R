#' Create PRPS sets using CCA and MNN.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function uses canonical correlation analysis (CCA) and mutual nearest neighbors (MNN) approach to create PRPS data
#' in the RNA-seq data. This function can be used in situation where the biological variation are entirely unknown.
#'
#' @param se.obj A `SummarizedExperiment` object.
#' @param assay.name Character. A character indicating the name of the data (assay) in the SummarizedExperiment object.
#' This data will be used to create PRPS data for RUV-III normalization. This data must be the one that will be
#' used for the RUV-III normalization.
#' @param main.uv.variable Character. Indicates the name of a column in the sample annotation of the SummarizedExperiment
#' object. The `uv.variable` can be either categorical or continuous. If `uv.variable` is a continuous variable, this will
#' be divided into `nb.clusters` groups using the `clustering.method`.
#' @param reference.group Character. A character string specifying which group in the `main.uv.variable` should be used
#' as the reference group when constructing PRPS. Default is set to `NULL`, meaning no explicit reference group is used.
#' @param other.uv.variables Character. A character or character vector representing the name(s) of the columns of
#' unwanted variable(s) within the sample annotation (colData) of the SummarizedExperiment object. These can be categorical,
#' continuous, or a combination. These variables will be considered when generating PRPS sets for the `main.uv.variable`
#' to help avoid potential contamination. The default is set to `NULL`
#' @param coordinates.to.use Character. Indicates which coordinates (e.g., PCA or CCA embedding) should be used to
#' identify mutual nearest neighbors. Default is set to `NULL`, meaning the raw assay data is used.
#' @param nb.cca Numeric. Number of canonical correlation components to use for constructing PRPS. Default is set to 20.
#' @param nb.pcs Numeric. Number of principal components to use for constructing PRPS. Default is set to 20.
#' @param samples.to.use Character. A vector of sample names to restrict the analysis to. Default is `NULL`, meaning
#' all samples in the SummarizedExperiment object are used.
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
#' @param cover.all.batches Logical. Indicates whether PRPS must be created such that all batches are covered. If `TRUE`,
#' the function ensures coverage across all batches; otherwise, partial coverage is allowed. Default is set to `TRUE`.
#' @param check.prps.connectedness Logical. Indicates whether to assess the `connectedness` between the PRPS sets across
#' all batches. Default is set to `TRUE`, indicating if there is not connections between all PRPS sets across all batches,
#' the function will give error. We refer to the checkPRPSconnectedness() function for more details.
#' @param hvg Vector. A logical vector or a vector of the names (feature ids) of the highly variable genes. These genes
#' will be used to prepare the input data for knn and mnn analysis. The default is set to `NULL`, this means all genes
#' will be used.
#' @param scale.cca Logical. Indicates whether to scale the CCA components before constructing PRPS. Default is set to `TRUE`.
#' @param apply.ruviii.norm Logical. Indicates whether to apply RUV-III normalization on the constructed PRPS before
#' downstream analysis. Default is set to `FALSE`.
#' @param use.ruviii.norm.for.mnn Logical. Indicates whether RUV-III normalized data should be used to identify mutual
#' nearest neighbors. Default is set to `FALSE`.
#' @param ncg Character. A vector of negative control genes to be used in the RUV-III normalization when
#' `apply.ruviii.norm = TRUE`. Default is set to `NULL`.
#' @param k Numeric. The number of nearest neighbors to use when finding mutual nearest neighbors. Default is set to 20.
#' @param nb.mnn Numeric. A numeric value specifying the maximum number of mutual nearest neighbors to compute. The
#' default is set to 1.
#' @param min.ps Numeric. Minimum number of pseudo-samples required for each PRPS set to be considered valid. Default is set to 2.
#' @param min.nb.for.mnn Numeric. Minimum number of neighbors required to define a mutual nearest neighbor pair. Default is set to 5.
#' @param similarity.approach Character. A character string specifying how similarity between samples is measured
#' (e.g., `cosine`, `correlation`, or `euclidean`). Default is set to `cosine`.
#' @param data.for.similarity Character. Specifies which type of data should be used for computing similarity
#' (e.g., `raw`, `normalized`, or `reduced`). Default is set to `normalized`.
#' @param clustering.method Character. A character indicating the choice of clustering method for grouping the
#' `uv.variable` if a continuous variable is provided. Options include `kmeans`, `cut`, and `quantile`. The default is set
#' to `kmeans`.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the `uv.variable` is a
#' continuous variable. The default is set to 3.
#' @param other.uv.clustering.method Character. A character indicating which clustering method should be used to
#' group each continuous unwanted variable, if specified in `other.uv.variables`. Options include `kmeans`, `cut`,
#' and `quantile`. The default is set to `kmeans`. See createHomogeneousUVGroups() for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous
#' unwanted variable specified in the `other.uv.variables`. The default is set to 3.
#' @param nb.batches.to.cover Numeric. A numeric value specifying the minimum number of batches that each PRPS set must
#' span to be considered valid. Default is set to 2.
#' @param normalization Character. A character that indicates which normalization method should be applied on the
#' data before finding the knn. Options are: `CPM`, `TMM`, `upper`, `median`, `full`, and `VST`. The default is set to
#' `cpm`.
#' If set to `NULL`, no normalization will be applied. See the applyOtherNormalizations() function for more details.
#' @param cosine.norm Logical. Indicates whether cosine normalization should be applied before finding MNN. Default
#' is set to `TRUE`.
#' @param regress.out.variables Character. A character or a vector of character that indicate the column name(s) in the
#' sample annotation in the SummarizedExperiment object. These variables will be regressed out from the data before
#' finding MNN. The default is set to `NULL`, indicating that regression will not be applied.
#' @param regress.out.rle.med Logical. Indicates whether to regress out the relative log expression (RLE) median
#' from the data before finding MNN. Default is set to `FALSE`.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not for down-stream analysis.
#' The default is set to `TRUE`.
#' @param apply.log.for.prps Logical. Indicates whether log-transformation should be applied only for the construction of
#' PRPS sets, regardless of global log-transformation. Default is set to `FALSE`.
#' @param pseudo.count Numeric. A positive numeric value as a pseudo count to be added to all measurements of the specified
#' assay(data) before applying log transformation to avoid -Inf for measurements that are equal to 0. The default is set
#' to 1.
#' @param assess.variables.association Logical. Indicates whether to assess the association between constructed PRPS and
#' biological/unwanted variables. Default is set to `TRUE`.
#' @param create.prps.map Logical. Indicates whether to generate a visualization (map) of the constructed PRPS sets across
#' batches and unwanted variables. Default is set to `FALSE`.
#' @param plot.output Logical. If `TRUE`, the function plots the distribution of MNN across the batches and PRPS sets
#' across the `main.uv.variable`.
#' @param mnn.bpparam Character. A BiocParallelParam object specifying how palatalization should be performed to find MNN.
#' The default is set to SerialParam(). We refer to the `findMutualNN()` function from the **BiocNeighbors** R package.
#' @param mnn.nbparam Character. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' The default is KmknnParam(). We refer to the `findMutualNN()` function from the `BiocNeighbors` R package.
#' @param residop.fun Character. A character indicating which function to use to calculate residuals. The options are
#' `c1`, `c2`, `lqr`, `r1` and `r2`. The default is set to `r2`.
#' @param nb.cores Numeric. The number of CPU cores to use for computation. The default is `NULL`, which automatically
#' uses the maximum number of available cores minus one.
#' @param use.annoy Logical. If `TRUE`, the Annoy algorithm is used for approximate nearest-neighbor (kNN) search. The
#' default is `FALSE`.
#' @param annot.nb.trees Numeric. The number of trees used in the Annoy index, corresponding to \code{AnnoyParam(ntrees = annot.nb.trees)}.
#' The default is set to 50.
#' @param max.iter Numeric. The maximum number of iterations allowed for the optimization or iterative procedure.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default is set
#' to `TRUE`. See the checkSeObj() function for more details.
#' @param remove.na Character. To remove NA or missing values from the assay (data) or not. The options are `assays` and
#' `none`. The default is set to `assays`, so all the NA or missing values from the assay(s) will be removed before computing
#' performing any down-stream analysis. See the `checkSeObj()` function for more details.
#' @param cca.set.name Character. A character string specifying the name of the CCA set used for PRPS construction.
#' Default is set to `NULL`, meaning a name will be automatically generated.
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
#' @importFrom BiocNeighbors findMutualNN AnnoyParam
#' @importFrom SummarizedExperiment assay colData
#' @importFrom batchelor cosineNorm
#' @importFrom Matrix solve
#' @importFrom methods as
#' @importFrom utils head
#' @importFrom RANN nn2
#' @export

createPrPsUnSupervisedByCca <- function(
        se.obj,
        assay.name,
        main.uv.variable,
        reference.group = NULL,
        other.uv.variables = NULL,
        coordinates.to.use = 'both',
        nb.cca = 2,
        nb.pcs = 5,
        samples.to.use = 'all',
        min.sample.for.ps = 3,
        select.extreme.groups = FALSE,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        min.batches.to.cover = 'all',
        cover.all.batches = FALSE,
        check.prps.connectedness = FALSE,
        hvg = NULL,
        scale.cca = TRUE,
        apply.ruviii.norm = TRUE,
        use.ruviii.norm.for.mnn = TRUE,
        ncg = NULL,
        k = 2,
        nb.mnn = 3,
        min.ps = 10,
        min.nb.for.mnn =  1,
        similarity.approach = 'euclidean',
        data.for.similarity = 'ruv',
        clustering.method = 'cut',
        nb.clusters = 3,
        other.uv.clustering.method = 'cut',
        nb.other.uv.clusters = 2,
        nb.batches.to.cover = 2,
        normalization = 'CPM',
        cosine.norm = FALSE,
        regress.out.variables = NULL,
        regress.out.rle.med = FALSE,
        apply.log = TRUE,
        apply.log.for.prps = TRUE,
        pseudo.count = 1,
        assess.variables.association = TRUE,
        create.prps.map = FALSE,
        plot.output = TRUE,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        residop.fun = 'r2',
        nb.cores = NULL,
        use.annoy = FALSE,
        annot.nb.trees = 50,
        max.iter = 200,
        check.se.obj = TRUE,
        remove.na = 'both',
        cca.set.name = NULL,
        prps.group.name = NULL,
        prps.sets.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The createPrPsUnSupervisedByCca function starts:',
        color = 'white',
        verbose = verbose
        )
    # Assessing and grouping the main unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the main unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    # Assessing which samples to use ####
    if (is.logical(samples.to.use)){
        initial.se.obj <- se.obj
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
    # Finding the class of the variable ####
    main.uv.variable.class <- class(se.obj[[main.uv.variable]])

    # Considering only main variable ####
    if (is.null(other.uv.variables)){
        ## Pre-processing variables #####
        all.variables <- prepareVariableForPrPs(
            se.obj = se.obj,
            main.variable = main.uv.variable,
            other.variables = other.uv.variables,
            nb.mnn = nb.mnn,
            min.sample.for.ps = min.sample.for.ps,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            other.uv.clustering.method = other.uv.clustering.method,
            nb.other.uv.clusters = nb.other.uv.clusters,
            select.extreme.groups = select.extreme.groups,
            cover.all.batches = cover.all.batches,
            nb.batches.to.cover = nb.batches.to.cover,
            assess.variables.association = assess.variables.association,
            plot.output = plot.output,
            verbose = verbose
            )
        # Applying a sanity check ####
        if (isFALSE(identical(colnames(se.obj), all.variables$sampl.ids ))){
            stop('There are issues with the order of samples.')
        }
        if (sum(all.variables$selected == 'TRUE') == 0){
            stop('There are issues with assessing and grouping the main unwanted variable.')
        }
        # Updating the SE object  ####
        if (samples.to.use == 'all'){
            initial.se.obj <- se.obj
        }
        se.obj[[main.uv.variable]] <- as.character(all.variables$groups)
        selected.samples <- all.variables$selected == 'TRUE'
        se.obj <- se.obj[ , selected.samples]

        # Finding CCA between all possible pairs of batches #####
        all.cca <- computeCCA(
            se.obj = se.obj,
            assay.name = assay.name,
            variable = main.uv.variable,
            reference.group = reference.group,
            nb.cca = nb.cca,
            hvg = hvg,
            scale = scale.cca,
            normalization = normalization,
            regress.out.variables = regress.out.variables,
            regress.out.rle.med = regress.out.rle.med,
            cosine.norm = cosine.norm,
            samples.to.use = 'all',
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            check.se.obj = FALSE,
            remove.na = 'none',
            cca.set.name = cca.set.name,
            save.se.obj = FALSE,
            verbose = verbose
            )
        # Retrieving initial variable names ####
        all.cca.data <- lapply(
            names(all.cca),
            function(name) {
                matches <- unique(se.obj[[main.uv.variable]])[sapply(
                    unique(se.obj[[main.uv.variable]]),
                    function(pattern) grepl(pattern, name))]
                as.character(matches)
            })
        # Applying a sanity check ####
        printColoredMessage(
            message = '- Applying a sanity check on the calculated CCA.',
            color = 'orange',
            verbose = verbose
            )
        for(i in 1:length(all.cca.data)){
            sample.ids <- c(
                colnames(se.obj)[se.obj[[main.uv.variable]] == all.cca.data[[i]][1] ],
                colnames(se.obj)[se.obj[[main.uv.variable]] == all.cca.data[[i]][2] ]
                )
            if (isFALSE(identical(row.names(all.cca[[i]]), sample.ids))){
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
        if (is.null(hvg))  hvg <- rep(TRUE, nrow(se.obj))
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
                ## Sub-setting the CCA data for each pairs of batches ####
                temp.cca <- all.cca[[x]]
                sample.annot.a <- all.sample.annot[all.sample.annot[[main.uv.variable]] == all.cca.data[[x]][1], , drop = FALSE ]
                sample.annot.b <- all.sample.annot[all.sample.annot[[main.uv.variable]] == all.cca.data[[x]][2], , drop = FALSE]
                all.samples <- c(row.names(sample.annot.a), row.names(sample.annot.b))
                ## Applying a sanity check ####
                if (isFALSE(identical(row.names(temp.cca), all.samples))){
                    stop('The order of samples in the calculated CCA have issues.')
                }
                ## Finding MNN using the CCA coordinates ####
                if (coordinates.to.use == 'cca'){
                    printColoredMessage(
                        message = '- Finding MNN using the CCA coordinates.',
                        color = 'orange',
                        verbose = verbose
                        )
                    ### Applying MNN ####
                    cca.mnn <- findMutualNN(
                        data1 = temp.cca[row.names(sample.annot.a) , ],
                        data2 = temp.cca[row.names(sample.annot.b) , ],
                        k1 = min.nb.for.mnn,
                        BPPARAM = mnn.bpparam,
                        nbparam = mnn.nbparam
                        )
                    ## Creating MNN - CCA data and adding samples ids ####
                    all.pp <- data.frame(
                        left.index = cca.mnn$first,
                        right.index = cca.mnn$second,
                        left.samples = row.names(sample.annot.a)[cca.mnn$first],
                        right.samples = row.names(sample.annot.b)[cca.mnn$second]
                        )
                }
                ## Finding MNN using the PCA coordinates ####
                if (coordinates.to.use == 'pca'){
                    printColoredMessage(
                        message = '- Finding MNN using the PCA coordinates.',
                        color = 'orange',
                        verbose = verbose
                        )
                    ### Applying normalization, regression an  the data ####
                    temp.data <- preProcessData(
                        se.obj = se.obj[ , all.samples],
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
                    ## Applying a sanity check ####
                    if (!identical(colnames(temp.data), all.samples)){
                        stop('The order of samples in the calculated CCA have issues.')
                    }
                    if (isTRUE(cosine.norm)){
                        temp.data <- cosineNorm(x = temp.data, mode = "matrix")
                    }
                    ## Performing PCA on the data ####
                    pca.data <- irlba::prcomp_irlba(
                        x = t(temp.data[hvg , ]),
                        center = TRUE,
                        scale = FALSE,
                        n = nb.pcs
                        )
                    ## Applying MNN ####
                    row.names(pca.data$x) <- colnames(temp.data)
                    pca.mnn <- findMutualNN(
                        data1 = pca.data$x[row.names(sample.annot.a) , ],
                        data2 = pca.data$x[row.names(sample.annot.b) , ],
                        k1  = min.nb.for.mnn,
                        BPPARAM = mnn.bpparam,
                        nbparam = mnn.nbparam
                        )
                    ## Creating MNN - PCA data and adding samples ids ####
                    pca.mnn <- data.frame(
                        left.index = pca.mnn$first,
                        right.index = pca.mnn$second,
                        left.samples = row.names(sample.annot.a)[pca.mnn$first],
                        right.samples = row.names(sample.annot.b)[pca.mnn$second]
                        )
                }
                ## Finding MNN using both CCA and PCA coordinates ####
                if (coordinates.to.use == 'both'){
                    printColoredMessage(
                        message = '- Finding MNN using both CCA and PCA coordinates.',
                        color = 'orange',
                        verbose = verbose
                        )
                    all.pp <- matrix(1)
                    iter <- 0
                    while(nrow(all.pp) < min.ps){
                        if (iter > max.iter) {
                            message("Stopped because maximum iterations reached.")
                            break
                        }
                        ### Finding MNN using the CCA coordinates ####
                        cca.mnn <- findMutualNN(
                            data1 = temp.cca[row.names(sample.annot.a) , ],
                            data2 = temp.cca[row.names(sample.annot.b) , ],
                            k1 = min.nb.for.mnn,
                            BPPARAM = mnn.bpparam,
                            nbparam = mnn.nbparam
                            )
                        ## Creating MNN - CCA data and adding samples ids ####
                        cca.mnn <- data.frame(
                            left.index = cca.mnn$first,
                            right.index = cca.mnn$second,
                            left.samples = row.names(sample.annot.a)[cca.mnn$first],
                            right.samples = row.names(sample.annot.b)[cca.mnn$second]
                            )
                        if (sum(is.na(cca.mnn)) > 0 | nrow(cca.mnn) == 0){
                            stop('There are some issues with finding MNN using CAA.')
                        }

                        ## Finding MNN using the PCA coordinates ####
                        # printColoredMessage(
                        #     message = '- Finding MNN using the PCA coordinates.',
                        #     color = 'orange',
                        #     verbose = verbose
                        #     )
                        ### Applying normalization, regression an log on the data ####
                        temp.data <- preProcessData(
                            se.obj = se.obj[ , all.samples],
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
                        if (!identical(colnames(temp.data), all.samples)){
                            stop('The order of samples in the calculated CCA have issues.')
                        }
                        if (isTRUE(cosine.norm)){
                            temp.data <- cosineNorm(x = temp.data, mode = "matrix")
                        }
                        ### Performing PCA on the data ####
                        pca.data <- irlba::prcomp_irlba(
                            x = t(temp.data[hvg , ]),
                            n = nb.pcs,
                            center = TRUE,
                            scale = FALSE
                            )
                        ### Applying MNN ####
                        row.names(pca.data$x) <- colnames(temp.data)
                        pca.mnn <- findMutualNN(
                            data1 = pca.data$x[row.names(sample.annot.a) , ],
                            data2 = pca.data$x[row.names(sample.annot.b) , ],
                            k1  = min.nb.for.mnn,
                            BPPARAM = mnn.bpparam,
                            nbparam = mnn.nbparam
                            )
                        ### Creating MNN - PCA data and adding samples ids ####
                        pca.mnn <- data.frame(
                            left.index = pca.mnn$first,
                            right.index = pca.mnn$second,
                            left.samples = row.names(sample.annot.a)[pca.mnn$first],
                            right.samples = row.names(sample.annot.b)[pca.mnn$second]
                            )
                        if (sum(is.na(pca.mnn)) > 0 | nrow(pca.mnn) == 0){
                            stop('There are some issues with finding MNN using CAA.')
                        }
                        ### Finding common PS across batches ####
                        all.pp <- c(
                            paste(cca.mnn$left.index, cca.mnn$right.index, sep = '_'),
                            paste(pca.mnn$left.index, pca.mnn$right.index, sep = '_')
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

                # Applying RUV-III normalization ####
                if (isTRUE(apply.ruviii.norm)){
                    ps.tr <- all.samples
                    for(j in 1:nrow(all.pp)){
                        ps.tr[ps.tr == all.pp$left.samples[j]] <- paste0('sample', j)
                        ps.tr[ps.tr == all.pp$right.samples[j]] <- paste0('sample', j)
                    }
                    if (sum(table(ps.tr)> 1) == 0){
                        stop('There are some issues with the applying RUV-III normalization.')
                    }
                    if (is.null(ncg)){
                        ncg <- row.names(se.obj) %in% row.names(se.obj)
                    }
                    m.matrix <- ruv::replicate.matrix(a = ps.tr)
                    if (isTRUE(apply.log)){
                        temp.data <- applyLog(
                            se.obj = se.obj[ , all.samples],
                            assay.names = assay.name,
                            pseudo.count = pseudo.count,
                            check.se.obj = FALSE,
                            remove.na = 'none',
                            verbose = verbose
                        )[[assay.name]]
                    } else temp.data <- assay(se.obj[, all.samples], i = assay.name)
                    Y <- t(temp.data)
                    Y.stand <- scale(
                        x = Y,
                        center = TRUE,
                        scale = FALSE
                        )
                    if (residop.fun == 'c1'){
                        Y0 <- fastResidopC1(Y, m.matrix)
                    } else if (residop.fun == 'c2'){
                        Y0 <- fastResidopC2(Y, m.matrix)
                    } else if (residop.fun == 'lqr'){
                        Y0 <- fastResidoplQR(Y, m.matrix)
                    } else if (residop.fun == 'r1'){
                        Y0 <- ruv::residop(Y, m.matrix)
                    } else if (residop.fun == 'r2'){
                        fastResidopR <- function(A, B) {
                            B <- Matrix::Matrix(B, sparse = TRUE)
                            BtB <- Matrix::t(B) %*% B
                            BtB.inv <- Matrix::solve(BtB)
                            BtA <- Matrix::t(B) %*% A
                            result <- A - B %*% BtB.inv %*% BtA
                            return(result)
                        }
                        Y0 <- fastResidopR(Y, m.matrix)
                    }
                    left.sing.value <- BiocSingular::runSVD(
                        x = Y0,
                        k = k,
                        BSPARAM = bsparam(),
                        center = FALSE,
                        scale = FALSE)$u
                    alpha <- t(left.sing.value[, 1:k, drop = FALSE]) %*% Y
                    ac <- alpha[, ncg, drop = FALSE]
                    W <- Y.stand[, ncg] %*% t(ac) %*% solve(ac %*% t(ac))
                    ruv.adj.data <- Y - W %*% alpha
                    pca.ruv <- irlba::prcomp_irlba(
                        x = ruv.adj.data,
                        n = nb.pcs,
                        center = TRUE,
                        scale = FALSE
                    )
                    ### Applying MNN ####
                    row.names(pca.ruv$x) <- colnames(temp.data)
                    ruv.adj <- pca.ruv$x
                }
                # Finding MNN using RUV-III normalized data ####
                if (isTRUE(use.ruviii.norm.for.mnn)){
                    printColoredMessage(
                        message = '- Applying RUV-III norm222.',
                        color = 'orange',
                        verbose = verbose
                    )
                    min.nb.for.mnn <- min.nb.for.mnn.initial
                    all.pp.new <- matrix(1)
                    iter <- 0
                    while (nrow(all.pp.new) < min.ps){
                        iter <- iter + 1
                        if (iter > max.iter) {
                            message("Stopped because maximum iterations reached.")
                            break
                        }
                        ### Applying MNN on RUV-III  ####
                        ruv.mnn <- findMutualNN(
                            data1 = ruv.adj[row.names(sample.annot.a) , ],
                            data2 = ruv.adj[row.names(sample.annot.b) , ],
                            k1  = min.nb.for.mnn,
                            BPPARAM = mnn.bpparam,
                            nbparam = mnn.nbparam
                            )
                        ### Creating MNN - RUV-III data and adding samples ids ####
                        ruv.mnn <- data.frame(
                            left.index = ruv.mnn$first,
                            right.index = ruv.mnn$second,
                            left.samples = row.names(sample.annot.a)[ruv.mnn$first],
                            right.samples = row.names(sample.annot.b)[ruv.mnn$second]
                            )
                        ### Applying a sanity check ####
                        if (sum(is.na(ruv.mnn)) > 0 | nrow(ruv.mnn) == 0){
                            stop('There are some issues with finding MNN using CAA.')
                        }
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
                        print(dim(all.pp.new))
                        min.nb.for.mnn <- min.nb.for.mnn + 1

                    }
                    all.pp <- all.pp.new
                }

                ## Finding the most similar samples ####
                if (data.for.similarity == 'ruv'){
                    similarity.data <- t(ruv.adj.data)
                }
                if (data.for.similarity == 'other'){
                    similarity.data <- preProcessData(
                        se.obj = se.obj[, all.samples],
                        assay.name = assay.name,
                        normalization = normalization,
                        regress.out.variables = regress.out.variables,
                        regress.out.rle.med = regress.out.rle.med,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        check.se.obj = FALSE,
                        remove.na = 'none',
                        verbose = verbose)
                }
                ### Using KNN approach  ####
                min.sample.for.ps.initial <- min.sample.for.ps
                if (similarity.approach == 'euclidean'){
                    #### Applying KNN for first batch ####
                    min.sample.for.ps <- min.sample.for.ps - 1
                    if (isTRUE(use.annoy)){
                        distance.a <- BiocNeighbors::findKNN(
                            X = t(similarity.data[hvg, row.names(sample.annot.a) , drop = FALSE]),
                            k = min.sample.for.ps,
                            num.threads = nb.cores,
                            BNPARAM = AnnoyParam(ntrees = annot.nb.trees)
                        )
                    } else {
                        distance.a <- BiocNeighbors::findKNN(
                            X = t(similarity.data[hvg, row.names(sample.annot.a) , drop = FALSE]),
                            k = min.sample.for.ps,
                            num.threads = nb.cores
                        )
                    }
                    distance.a.index <- as.data.frame(distance.a$index)
                    colnames(distance.a.index) <- paste0('V', 2:(min.sample.for.ps + 1))
                    distance.a.index$V1 <- as.numeric(row.names(distance.a.index))
                    distance.a.index <- distance.a.index[order(colnames(distance.a.index))]

                    for(i in 1:ncol(distance.a.index)){
                        col.name <- paste0('Sample', i)
                        distance.a.index[col.name] <- row.names(sample.annot.a)[distance.a.index[ , i]]
                    }
                    ### Applying KNN for second batch ####
                    print(dim(t(similarity.data[hvg, row.names(sample.annot.b) , drop = FALSE])))
                    if (isTRUE(use.annoy)){
                        distance.b <- BiocNeighbors::findKNN(
                            X = t(similarity.data[hvg, row.names(sample.annot.b) , drop = FALSE]),
                            k = min.sample.for.ps,
                            num.threads = nb.cores,
                            BNPARAM = AnnoyParam(ntrees = annot.nb.trees)
                        )
                    } else{
                        distance.b <- BiocNeighbors::findKNN(
                            X = t(similarity.data[hvg, row.names(sample.annot.b) , drop = FALSE]),
                            k = min.sample.for.ps,
                            num.threads = nb.cores
                        )
                    }
                    distance.b.index <- as.data.frame(distance.b$index)
                    colnames(distance.b.index) <- paste0('V', 2:(min.sample.for.ps + 1))
                    distance.b.index$V1 <- as.numeric(row.names(distance.b.index))
                    distance.b.index <- distance.b.index[order(colnames(distance.b.index))]
                    for(i in 1:ncol(distance.b.index)){
                        col.name <- paste0('Sample', i)
                        distance.b.index[col.name] <- row.names(sample.annot.b)[distance.b.index[ , i]]
                    }
                }
                ### Using correlation approach ####
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
                            top.samples <- utils::head(x = sorted, n = min.sample.for.ps)
                            temp.df <- data.frame(
                                Sample1 = sample_name,
                                Topmatch = names(top.samples),
                                Correlation = mean(as.numeric(top.samples))
                            )
                            temp.df$Topmatchid <- c(2:c(min.sample.for.ps + 1))
                            temp.df <- temp.df %>%
                                tidyr::pivot_wider(
                                    names_from = Topmatchid,
                                    values_from = Topmatch,
                                    names_prefix = "Sample") %>%
                                data.frame()
                        })
                    distance.a.index <- do.call(rbind, distance.a.index)
                    distance.a.index <- distance.a.index[ , c(1, c(3:c(min.sample.for.ps + 2)), 2)]
                    if (isFALSE(all.equal(distance.a.index$Sample1, row.names(sample.annot.a)))){
                        stop('There are something wrong with the ...')
                    }
                    ### Applying correlation for first batch ####
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
                            top.samples <- utils::head(sorted, min.sample.for.ps)
                            temp.df <- data.frame(
                                Sample1 = sample_name,
                                Topmatch = names(top.samples),
                                Correlation = mean(as.numeric(top.samples))
                            )
                            temp.df$Topmatchid <-  c(2:c(min.sample.for.ps + 1))
                            temp.df <- temp.df %>%
                                tidyr::pivot_wider(
                                    names_from = Topmatchid,
                                    values_from = Topmatch,
                                    names_prefix = "Sample") %>%
                                data.frame()
                        })
                    distance.b.index <- do.call(rbind, distance.b.index)
                    distance.b.index <- distance.b.index[ , c(1, c(3:c(min.sample.for.ps + 2)), 2)]
                    if (isFALSE(all.equal(distance.b.index$Sample1, row.names(sample.annot.b)))){
                        stop('There are something wrong with the ...')
                    }
                }
                ### Using cosine approach ####
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
                            top.samples <- utils::head(sorted, min.sample.for.ps)
                            temp.df <- data.frame(
                                Sample1 = sample_name,
                                Topmatch = names(top.samples),
                                Correlation = mean(as.numeric(top.samples))
                            )
                            temp.df$Topmatchid <- c(2:c(min.sample.for.ps + 1))
                            temp.df <- temp.df %>%
                                tidyr::pivot_wider(
                                    names_from = Topmatchid,
                                    values_from = Topmatch,
                                    names_prefix = "Sample") %>%
                                data.frame()
                        })
                    distance.a.index <- do.call(rbind, distance.a.index)
                    distance.a.index <- distance.a.index[ , c(1, c(3:c(min.sample.for.ps + 2)), 2)]
                    if (isFALSE(all.equal(distance.a.index$Sample1, row.names(sample.annot.a)))){
                        stop('There are something wrong with the ...')
                    }
                    ### Applying correlation for first batch ####
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
                            top.samples <- utils::head(sorted, min.sample.for.ps)
                            temp.df <- data.frame(
                                Sample1 = sample_name,
                                Topmatch = names(top.samples),
                                Correlation = mean(as.numeric(top.samples))
                            )
                            temp.df$Topmatchid <-  c(2:c(min.sample.for.ps + 1))
                            temp.df <- temp.df %>%
                                tidyr::pivot_wider(
                                    names_from = Topmatchid,
                                    values_from = Topmatch,
                                    names_prefix = "Sample") %>%
                                data.frame()
                        })
                    distance.b.index <- do.call(rbind, distance.b.index)
                    distance.b.index <- distance.b.index[ , c(1, c(3:c(min.sample.for.ps + 2)), 2)]
                    if (isFALSE(all.equal(distance.b.index$Sample1, row.names(sample.annot.b)))){
                        stop('There are something wrong with the ...')
                    }
                }
                ## Creating PRPS data ####
                ### Pre-processing the data ####
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
                ### Filtering the PRPS sets ####
                if (isTRUE(filter.prps.sets)){
                    if (nrow(all.pp) > max.prps.sets){
                        all.pp <- all.pp[1:max.prps.sets, ]
                    }
                }
                ### Generating the PRPS matrix ####
                prps.data <- lapply(
                    1:nrow(all.pp),
                    function(ps){
                        selected.cols <- paste0('Sample', c(1:min.sample.for.ps.initial))
                        ps1 <- unlist(unname(distance.a.index[distance.a.index$Sample1 == all.pp$left.samples[ps] , selected.cols ]))
                        ps1 <- rowMeans(prps.data[ , ps1])
                        ps2 <- unlist(unname(distance.b.index[distance.b.index$Sample1 == all.pp$right.samples[ps] , selected.cols ]))
                        ps2 <- rowMeans(prps.data[ , ps2])
                        prps.set <- data.frame(
                            ps1 = ps1 ,
                            ps2 = ps2
                            )
                        colnames(prps.set) <- rep(
                            x = paste0('sample', ps, x, main.uv.variable),
                            each = 2
                            )
                        prps.set
                    })
                prps.data <- do.call(cbind , prps.data)
                ### Generating the PRPS map ####
                if (isTRUE(create.prps.map)){
                    prps.map <- lapply(
                        1:nrow(all.pp),
                        function(ps){
                            selected.cols <- paste0('Sample', c(1:min.sample.for.ps.initial))
                            ps1 <- unlist(unname(distance.a.index[distance.a.index$Sample1 == all.pp$left.samples[ps] , selected.cols ]))
                            ps2 <- unlist(unname(distance.b.index[distance.b.index$Sample1 == all.pp$right.samples[ps] , selected.cols ]))
                            prps.set <- data.frame(
                                samples.1 = ps1,
                                samples.2 = ps2,
                                group = paste0('set', ps, x)
                                )
                            prps.set
                    })
                    prps.map
                }
                if (isTRUE(create.prps.map)){
                    return(list(prps.data = prps.data, prps.map = prps.map, all.mnn = all.pp))
                } else return(list(prps.data = prps.data, all.mnn = all.pp))
            })

        if (isTRUE(create.prps.map)){
            all.prps.sets <- lapply(seq_len(length(all.prps.data)), function(x) all.prps.data[[x]]$prps.data)
            all.prps.sets <- do.call(cbind, all.prps.sets)
            all.mnn <- lapply(seq_len(length(all.prps.data)), function(x) all.prps.data[[x]]$all.mnn)
            all.mnn <- do.call(rbind, all.mnn)
        } else {
            all.prps.sets <- do.call(cbind, all.prps.data)
            all.mnn <- lapply(seq_len(length(all.cca.data)), function(x) all.prps.data[[x]]$all.mnn)
            all.mnn <- do.call(rbind, all.mnn)
        }

        # Plotting the PRPS map ####
        if (isTRUE(create.prps.map) & main.uv.variable.class %in% c('numeric', 'integer')){
            prps.map <- lapply(
                seq_len(length(all.cca.data)),
                function(x) {
                    do.call(rbind, all.prps.data[[x]]$prps.map)
                })
            prps.map <- do.call(rbind, prps.map)
            df <- tibble(
                col = all.variables$sampl.ids,
                variable = all.variables$variable,
                groups = all.variables$groups
                )
            prps.map <- prps.map %>%
                left_join(df, by = c("samples.1" = "col")) %>%
                dplyr::rename(group.set.1 = variable, group.1 = groups) %>%
                left_join(df, by = c("samples.2" = "col")) %>%
                dplyr::rename(group.set.2 = variable, group.2 = groups)
            prps.map <- prps.map %>%
                pivot_longer(-c('group', 'samples.1', 'samples.2', 'group.1', 'group.2'), values_to = 'val', names_to = 'set') %>%
                group_by(group) %>%
                mutate(
                    mean.group.set1 = mean(val[set == "group.set.1"], na.rm = TRUE),
                    mean.group.set2 = mean(val[set == "group.set.2"], na.rm = TRUE),
                    dominant.set = dplyr::if_else(mean.group.set1 >= mean.group.set2, "group.set.1", "group.set.2"),
                    set = dplyr::if_else(set == dominant.set, "top", "bottom")
                ) %>%
                dplyr::select(-mean.group.set1, -mean.group.set2, -dominant.set) %>%
                dplyr::ungroup() %>%
                dplyr::select(group, val, set) %>%
                mutate(category = 'PRPS sets')
            all.sample.info <- data.frame(
                group = 'uv',
                val = all.variables$variable,
                set = main.uv.variable,
                category = 'UV'
                )
            all.info <- rbind(all.sample.info, prps.map)
            all.info$set <- factor(all.info$set, levels = c(main.uv.variable, 'top', 'bottom'))
            all.info$category <- factor(all.info$category, levels = c(main.uv.variable, 'UV', 'PRPS sets'))
            prps.map <- ggplot(all.info, aes(x = val, y = group, color = set)) +
                geom_boxplot() +
                geom_point() +
                facet_grid(category~., scales = 'free', space = 'free') +
                scale_color_manual(values = c('darkgreen', 'orange', 'navy'), name = 'Groups') +
                xlab(main.uv.variable) +
                ylab('Homogeneous groups') +
                theme_bw() +
                theme(
                    legend.key = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text.y = element_text(size = 12),
                    axis.text.x = element_text(size = 12, angle = 35, hjust = 1, vjust = 1),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    strip.text.y = element_text(size = 15)
                )
            if (isTRUE(plot.output)) print(prps.map)
        }
        if (isTRUE(create.prps.map) & !main.uv.variable.class %in% c('numeric', 'integer')){
            prps.map <- lapply(
                seq_len(length(all.cca.data)),
                function(x) {
                    do.call(rbind, all.prps.data[[x]]$prps.map)
                })
            prps.map <- do.call(rbind, prps.map)
            df <- tibble(
                col = all.variables$sampl.ids,
                variable = all.variables$variable,
                groups = all.variables$groups
            )
            prps.map <- prps.map %>%
                dplyr::left_join(df, by = c("samples.1" = "col")) %>%
                dplyr::rename(group.set.1 = variable, group.1 = groups) %>%
                dplyr::left_join(df, by = c("samples.2" = "col")) %>%
                dplyr::rename(group.set.2 = variable, group.2 = groups)
            all.info <- prps.map[ , c('group', 'group.1', 'group.2')] %>%
                pivot_longer(-group, values_to = 'batch', names_to = 'cate')
            prps.map <- ggplot(all.info, aes(x = batch, y = group, group = group)) +
                geom_line() +
                geom_point() +
                xlab(main.uv.variable) +
                ylab('Homogeneous groups') +
                theme_bw() +
                theme(
                    legend.key = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text.y = element_text(size = 12),
                    axis.text.x = element_text(size = 12, angle = 35, hjust = 1, vjust = 1),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    strip.text.y = element_text(size = 15)
                )
            if (isTRUE(plot.output)) print(prps.map)
        }
    }
    # Considering other uv variables ####
    if (!is.null(other.uv.variables)){
        ## Pre-processing variables #####
        all.variables <- prepareVariableForPrPs(
            se.obj = se.obj,
            main.variable = main.uv.variable,
            other.variables = other.uv.variables,
            nb.mnn = nb.mnn,
            min.sample.for.ps = min.sample.for.ps,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            other.uv.clustering.method = other.uv.clustering.method,
            nb.other.uv.clusters = nb.other.uv.clusters,
            select.extreme.groups = select.extreme.groups,
            cover.all.batches = cover.all.batches,
            nb.batches.to.cover = nb.batches.to.cover,
            assess.variables.association = assess.variables.association,
            plot.output = plot.output,
            verbose = verbose
            )
        colnames(all.variables)[3] <- 'main.uv.variable'
        # Applying a sanity check ####
        if (isFALSE(identical(colnames(se.obj), all.variables$sampl.ids))){
            stop('There are issues with the order of samples.')
        }
        if (sum(all.variables$selected == 'TRUE') == 0){
            stop('There are issues with assessing and grouping the main unwanted variable.')
        }

        # Updating the SummarizedExperiment object  ####
        initial.se.obj <- se.obj
        initial.se.obj[[main.uv.variable]] <- all.variables$main.uv.variable
        selected.samples <- all.variables$selected == 'TRUE'
        sub.all.variables <- droplevels(all.variables[selected.samples , ])
        other.uv.variables <- unique(sub.all.variables$other.variables)

        # Creating PRPS per each group of other unwanted variables ####
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
                sample.index <- all.variables[all.variables$other.variables == d , ]
                se.obj <- initial.se.obj[ , sample.index$sampl.ids]
                all.cca <- computeCCA(
                    se.obj = se.obj,
                    assay.name = assay.name,
                    variable = main.uv.variable,
                    reference.group = reference.group,
                    nb.cca = nb.cca,
                    hvg = hvg,
                    scale = scale.cca,
                    normalization = normalization,
                    regress.out.variables = regress.out.variables,
                    regress.out.rle.med = regress.out.rle.med,
                    cosine.norm = cosine.norm,
                    samples.to.use = 'all',
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
                    if (isFALSE(identical(row.names(all.cca[[i]]), sample.ids))){
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
                if (is.null(hvg)) hvg <- rep(TRUE, nrow(se.obj))
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
                        if (!identical(row.names(temp.cca), all.samples)){
                            stop('The order of samples in the calculated CCA have issues.')
                        }
                        if (coordinates.to.use == 'cca'){
                            ## Finding MNN using the CCA coordinates ####
                            printColoredMessage(
                                message = '- Finding MNN using the CCA coordinates.',
                                color = 'orange',
                                verbose = verbose
                                )
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
                                n = nb.pcs,
                                center = TRUE,
                                scale = FALSE
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
                            iter <- 0
                            while(nrow(all.pp) < min.ps){
                                ### Applying MNN ####
                                if (iter > max.iter) {
                                    message("Stopped because maximum iterations reached.")
                                    break
                                }
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
                                    n = nb.pcs,
                                    center = TRUE,
                                    scale = FALSE
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
                            if (residop.fun == 'c1'){
                                Y0 <- fastResidopC1(Y, m.matrix)
                            } else if (residop.fun == 'c2'){
                                Y0 <- fastResidopC2(Y, m.matrix)
                            } else if (residop.fun == 'lqr'){
                                Y0 <- fastResidoplQR(Y, m.matrix)
                            } else if (residop.fun == 'r1'){
                                Y0 <- ruv::residop(Y, m.matrix)
                            } else if (residop.fun == 'r2'){
                                fastResidopR <- function(A, B) {
                                    B <- Matrix::Matrix(B, sparse = TRUE)
                                    BtB <- Matrix::t(B) %*% B
                                    BtB.inv <- Matrix::solve(BtB)
                                    BtA <- Matrix::t(B) %*% A
                                    result <- A - B %*% BtB.inv %*% BtA
                                    return(result)
                                }
                                Y0 <- fastResidopR(Y, m.matrix)
                            }
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
                        if (isTRUE(use.ruviii.norm.for.mnn)){
                            min.nb.for.mnn <- min.nb.for.mnn.initial
                            all.pp.new <- matrix(1)
                            iter <- 0
                            while(nrow(all.pp.new) < min.ps){
                                if (iter > max.iter) {
                                    message("Stopped because maximum iterations reached.")
                                    break
                                }
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
                        if (data.for.similarity == 'ruv'){
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
                            if (isTRUE(use.annoy)){
                                distance.a <- BiocNeighbors::findKNN(
                                    X = t(similarity.data[hvg, row.names(sample.annot.a) , drop = FALSE]),
                                    k = min.sample.for.ps,
                                    num.threads = nb.cores,
                                    BNPARAM = AnnoyParam(ntrees = annot.nb.trees)
                                )
                            } else{
                                distance.a <- BiocNeighbors::findKNN(
                                    X = t(similarity.data[hvg, row.names(sample.annot.a) , drop = FALSE]),
                                    k = min.sample.for.ps,
                                    num.threads = nb.cores
                                )
                            }

                            distance.a.index <- as.data.frame(distance.a$index)
                            colnames(distance.a.index) <- paste0('V', 2:(min.sample.for.ps + 1))
                            distance.a.index$V1 <- as.numeric(row.names(distance.a.index))
                            distance.a.index <- distance.a.index[order(colnames(distance.a.index))]
                            for(i in 1:ncol(distance.a.index)){
                                col.name <- paste0('Sample', i)
                                distance.a.index[col.name] <- row.names(sample.annot.a)[distance.a.index[ , i]]
                            }
                            #### Applying KNN for second batch ####
                            if (isTRUE(use.annoy)){
                                distance.b <- BiocNeighbors::findKNN(
                                    X = t(similarity.data[hvg, row.names(sample.annot.b) , drop = FALSE]),
                                    k = min.sample.for.ps,
                                    num.threads = nb.cores,
                                    BNPARAM = AnnoyParam(ntrees = annot.nb.trees)
                                )
                            } else{
                                distance.b <- BiocNeighbors::findKNN(
                                    X = t(similarity.data[hvg, row.names(sample.annot.b) , drop = FALSE]),
                                    k = min.sample.for.ps,
                                    num.threads = nb.cores
                                )
                            }
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
                                    top.samples <- utils::head(x = sorted, n = min.sample.for.ps)
                                    temp.df <- data.frame(
                                        Sample1 = sample_name,
                                        Topmatch = names(top.samples),
                                        Correlation = mean(as.numeric(top.samples))
                                    )
                                    temp.df$Topmatchid <- c(2:c(min.sample.for.ps + 1))
                                    temp.df <- temp.df %>%
                                        tidyr::pivot_wider(
                                            names_from = Topmatchid,
                                            values_from = Topmatch,
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
                                    top.samples <- utils::head(sorted, min.sample.for.ps)
                                    temp.df <- data.frame(
                                        Sample1 = sample_name,
                                        Topmatch = names(top.samples),
                                        Correlation = mean(as.numeric(top.samples))
                                    )
                                    temp.df$Topmatchid <-  c(2:c(min.sample.for.ps + 1))
                                    temp.df <- temp.df %>%
                                        tidyr::pivot_wider(
                                            names_from = Topmatchid,
                                            values_from = Topmatch,
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
                                    top.samples <- utils::head(sorted, min.sample.for.ps)
                                    temp.df <- data.frame(
                                        Sample1 = sample_name,
                                        Topmatch = names(top.samples),
                                        Correlation = mean(as.numeric(top.samples))
                                    )
                                    temp.df$Topmatchid <- c(2:c(min.sample.for.ps + 1))
                                    temp.df <- temp.df %>%
                                        tidyr::pivot_wider(
                                            names_from = Topmatchid,
                                            values_from = Topmatch,
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
                                    top.samples <- utils::head(sorted, min.sample.for.ps)
                                    temp.df <- data.frame(
                                        Sample1 = sample_name,
                                        Topmatch = names(top.samples),
                                        Correlation = mean(as.numeric(top.samples))
                                    )
                                    temp.df$Topmatchid <-  c(2:c(min.sample.for.ps + 1))
                                    temp.df <- temp.df %>%
                                        tidyr::pivot_wider(
                                            names_from = Topmatchid,
                                            values_from = Topmatch,
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

    }
    # Saving the results ####
    ## Selecting prps.sets.name ####
    prps.sets.name <- paste0(main.uv.variable, '|', 'CcaPca', '|', assay.name)
    if (is.null(prps.group.name)) {
        prps.group.name <- main.uv.variable
    }
    printColoredMessage(
        message = '-- Saving the PRPS data',
        color = 'magenta',
        verbose = verbose
        )
    se.obj <- initial.se.obj
    # se.obj <- initial.se.obj.a
    ## Saving the PRPS data in the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = 'Save all the PRPS data into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        if (!'PRPS' %in% names(se.obj@metadata)) {
            se.obj@metadata[['PRPS']] <- list()
        }
        if (!'un.supervised' %in% names(se.obj@metadata[['PRPS']])) {
            se.obj@metadata[['PRPS']][['un.supervised']] <- list()
        }
        if (!prps.group.name %in% names(se.obj@metadata[['PRPS']][['un.supervised']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]] <- list()
        }
        if (!prps.sets.name %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]] <- list()
        }
        if (!'prps.data' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.data']] <- list()
        }
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.data']] <- all.prps.sets
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['mnn.samples']] <- all.mnn

        if (isTRUE(create.prps.map)){
            if (!'prps.map.plot' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]])) {
                se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.map.plot']] <- list()
            }
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.map.plot']]  <- prps.map
        }

        printColoredMessage(
            message = '------------The createPrPsUnSupervisedByCca function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ## Outputting the PRPS data as matrix ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '------------The createPrPsUnSupervisedByCca function finished.',
            color = 'white',
            verbose = verbose
            )
        return(list(prps.data = all.prps.sets))
    }
}
