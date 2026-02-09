#' Creates PRPS sets using k and mutual nearest neighbors in RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function uses the k and mutual nearest neighbors approaches to create PRPS data for RUV-III normalization in
#' RNA-seq data. This function can be used in situation where biological variation are entirely unknown.
#'
#' @details
#' The function applies the `findKnn` function to find similar samples per batch and then average them to create
#' pseudo-samples. Then, function uses the `findMnn` to match up pseudo samples across batches to create pseudo-replicates.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character indicating the name of the data(assay) in the SummarizedExperiment
#' object. This data will be used to create PRPS data for RUV-III normalization. This data must be the one that will be
#' used for the RUV-III normalization.
#' @param main.uv.variable Character. Indicates the name of a column in the sample annotation of the SummarizedExperiment
#' object. The `main.uv.variable` can be either categorical or continuous. If `main.uv.variable` is a continuous variable,
#' this will be divided into `nb.clusters` groups using the `clustering.method`.
#' @param clustering.method Character. A character indicating the choice of clustering method for grouping the
#' `main.uv.variable` if a continuous variable is provided. Options include `kmeans`, `cut`, and `quantile`. The default
#' is set to `kmeans`.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the `main.uv.variable` is
#' a continuous variable. The default is set to 3.
#' @param other.uv.variables Character. A character or character vector representing the name(s) of the columns of
#' unwanted variable(s) within the sample annotation (colData) of the SummarizedExperiment object. These can be categorical,
#' continuous, or a combination in PRPS data . These variables will be considered when generating PRPS sets for the
#' `main.uv.variable` to help avoid potential contamination. The default is set to `NULL`
#' @param other.uv.clustering.method Character. A character indicating which clustering method should be used to
#' group each continuous unwanted variable, if specified in `other.uv.variables`. Options include `kmeans`, `cut`,
#' and `quantile`. The default is `kmeans`. See `createHomogeneousUVGroups()` for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous
#' unwanted variable specified in the `other.uv.variables`. The default is set to 3.
#' @param min.sample.for.ps Numeric. Minimum number of samples required for pseudo-replicate creation. The default is set
#' to 3.
#' @param select.extreme.groups Logical. Indicates whether to select only the extreme groups e.g., highest and lowest
#' clusters, when the `main.uv.variable` is a continuous variable. Default is set to `TRUE`. This will increase the
#' variation between PR sets in order to better capture the unwanted variation.
#' @param filter.prps.sets Logical. If `TRUE`, the number of PRPS sets across each pair of batches will be filtered if
#' they are higher than the `max.prps.sets` value. A high number of PRPS sets will increase the computational time for
#' the RUV-III normalization. The default is set to `TRUE`.
#' @param max.prps.sets Numeric. A numeric value specifying the maximum number for PRPS sets across each pair of batches(
#' subgrouos of the `main.uv.variable`. The default is set to 10.
#' @param min.batches.to.cover Numeric. Minimum number of batches that must be covered by PRPS set. The default is set to
#' `all`, indicating all possible batch must have enough samples to create PRPS, otherwise the function gives error.
#' @param check.prps.connectedness Logical. Indicates whether to assess the `connectedness` between the PRPS sets across
#' all batches. Default is set to `TRUE`, indicating if there is not connections between all PRPS sets across all batches,
#' the function will give error.We refer to the checkPRPSconnectedness() function for more details.
#' @param data.input Character. A character that indicates which data should be used as input for finding the k
#' and mutual nearest neighbors. Options: `expr` and `pcs`. If `pcs` is selected, the first `nb.pcs` principal components
#' will be used. Default is set to `expr`.
#' @param nb.pcs Numeric. Number of principal components to be calculated and used when data.input = `pcs`. Default is
#' set to 2.
#' @param center Logical. Indicates whether to center the data or not before calculating PCs. If center is `TRUE`, then
#' centering is done by subtracting the column means of the assay from their corresponding columns. The default is set
#' to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data or not before calculating PCs. If scale is set to `TRUE`, then
#' scaling is done by dividing the (centered) columns of the assays by their standard deviations if center is `TRUE`, and
#' the root mean square otherwise. The default is set to `FALSE`.
#' @param svd.bsparam Character. A BiocParallelParam object specifying how palatalization should be performed. The default
#' is set to bsparam(). We refer to the `runSVD` function from the BiocSingular R package for further details.
#' @param nb.knn Numeric. A numeric number that indicates the maximum number of k nearest neighbors to compute for each
#' sample. The default is set to 3.
#' @param nb.mnn Numeric. A numeric value specifying the maximum number of mutual nearest neighbors to compute. The
#' default is set to 1.
#' @param hvg Vector. A logical vector or a vector of the names (feature ids) of the highly variable genes. These genes
#' will be used to prepare the input data for knn and mnn analysis. The default is set to `NULL`, this means all genes
#' will be used.
#' @param normalization Character. A character that indicates which normalization method should be applied on the
#' data before finding the knn. Options are: `CPM`, `TMM`, `upper`, `median`, `full` and `VST`. The default is set to `cpm`.
#' If set to `NULL`, no normalization will be applied. See the `applyOtherNormalizations()` function for more details.
#' @param apply.cosine.norm Logical. Indicates whether cosine normalization should be applied before finding MNN. Default
#' is set to `TRUE`.
#' @param regress.out.variables Character. A character or strings that indicate the column name(s) in the sample
#' annotation in the SummarizedExperiment object. These variables will be regressed out from the data before
#' finding KNN and MNN. The default is set to `NULL`, indicating that regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not for down-stream analysis.
#' The default is set to `TRUE`.
#' @param apply.log.for.prps Logical. Indicates whether to apply a log-transformation to the data before averaging samples
#' to create PS or not The default is set to `TRUE`.
#' @param pseudo.count Numeric. A positive numeric value as a pseudo count to be added to all measurements of the specified
#' assay(data) before applying log transformation to avoid -Inf for measurements that are equal to 0. The default is set
#' to 1.
#' @param mnn.bpparam Character. A BiocParallelParam object specifying how palatalization should be performed to find MNN.
#' The default is SerialParam(). We refer to the **`findMutualNN()`** function from the **BiocNeighbors** R package.
#' @param mnn.nbparam Character. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' The default is KmknnParam(). We refer to the `findMutualNN()` function from the **BiocNeighbors** R package.
#' @param samples.to.use Logical. A logifcal vector specifiyung wich samples should be used for the analysis. The default
#' is set to `all`, then all samples will be used.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default is set
#' to TRUE. See the checkSeObj() function for more details.
#' @param remove.na Character. To remove NA or missing values from the assay (data) or not. The options are `assays` and
#' `none`. The default is set to `assays`, so all the NA or missing values from the assay(s) will be removed before
#' computing performing any down-stream analysis. See the checkSeObj() function for more details.
#' @param plot.output Logical. If `TRUE`, the function plots the distribution of MNN across the batches and PRPS sets
#' across the `main.uv.variable`.
#' @param knn.group.name Character. A character specifying the name of the knn  to which the current KNN belong.
#' If set to `NULL`, the function will automatically assign a name using  `main.uv.variable`.
#' @param knn.sets.name  Character. A character specifying the name of the knn set names to be saved in the metadata
#' of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0(uv.variable, '|', assay.name)`.
#' @param mnn.group.name Character. A character specifying the name of the mnn to which the current MNN belong.
#' If set to `NULL`, the function will automatically assign a name using  `main.uv.variable`.
#' @param mnn.sets.name A character specifying the name of the mnn set names to be saved in the metadata of the
#' SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0(uv.variable, '|', assay.name)`.
#' @param prps.group.name Character. A character specifying the name of the prps.group.name to which the current PRPS belong.
#' If set to `NULL`, the function will automatically assign a name using  `main.uv.variable`.
#' @param prps.sets.name Character. A character specifying the name of the output file to be saved in the metadata
#' of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0(uv.variable, '|', assay.name)`.
#' @param save.se.obj Logical. Indicates whether to save the KNN results in the metadata of the SummarizedExperiment object
#' or to output the result as a list. By default, it is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, shows the messages of different steps of the function.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom BiocNeighbors findMutualNN KmknnParam
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocParallel SerialParam
#' @importFrom stats dist
#' @importFrom RANN nn2
#'
#' @export

createPrPsUnsupervisedByKnnMnn <- function(
        se.obj,
        assay.name,
        main.uv.variable,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        other.uv.variables = NULL,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        min.sample.for.ps = 3,
        select.extreme.groups = FALSE,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        min.batches.to.cover = 'all',
        check.prps.connectedness = TRUE,
        data.input = 'expr',
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        nb.knn = 2,
        nb.mnn = 1,
        hvg = NULL,
        normalization = 'CPM',
        apply.cosine.norm = FALSE,
        regress.out.variables = NULL,
        apply.log = TRUE,
        apply.log.for.prps = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        samples.to.use = 'all',
        check.se.obj = TRUE,
        remove.na = 'both',
        plot.output = TRUE,
        knn.group.name = NULL,
        knn.sets.name = NULL,
        mnn.group.name = NULL,
        mnn.sets.name = NULL,
        prps.group.name = NULL,
        prps.sets.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The createPrPsByKnnMnn function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the function inputs ####
    if (is.null(assay.name) | is.logical(assay.name)) {
        stop('The "assay.name" cannot be empty or logical.')
    }
    if (length(assay.name) > 1 | assay.name == 'all') {
        stop('The "assay.name" must be an assay name in the SummarizedExperiment object.')
    }
    if (isFALSE(check.se.obj)){
        if (!assay.name %in% names(assays(se.obj))){
            stop('The "assay.name" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (length(main.uv.variable) > 1) {
        stop('The "main.uv.variable" must a categorical or continuous variable in the SummarizedExperiment object.')
    }
    if (is.null(main.uv.variable) | is.logical(main.uv.variable)) {
        stop('The "main.uv.variable" cannot be empty or logical(TRUE or FALSE).')
    }
    if (isFALSE(check.se.obj)){
        if (!main.uv.variable %in% colnames(colData(se.obj))){
            stop('The "main.uv.variable" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (is.numeric(colData(se.obj)[[main.uv.variable]])){
        if (var(colData(se.obj)[[main.uv.variable]]) == 0){
            stop('The variance of the "main.uv.variable" is 0. No need to create PRPS for this variable.')
        }
    }
    if (!is.null(other.uv.variables)){
        if (main.uv.variable %in% other.uv.variables){
            stop('The "main.uv.variable" must not be in the "other.uv.variables".')
        }
        if (isFALSE(check.se.obj)){
            if (sum(other.uv.variables %in% colnames(colData(se.obj))) != length(other.uv.variables)){
                stop('All or some of the "other.uv.variables" cannot be found in the SummarizedExperiment object.')
            }
        }
    }
    if (!is.logical(filter.prps.sets)){
        stop('The "filter.prps.sets" must be logical (TRUE or FALSE)')
    }
    if (isTRUE(filter.prps.sets)){
        if (!is.numeric(max.prps.sets) | max.prps.sets < 0){
            stop('The "max.prps.sets" must be postive numeric value.')
        }
    }
    if (min.sample.for.ps <= 1) {
        stop('The minimum value for the "min.sample.for.ps" is 2.')
    }
    if (!is.logical(apply.log)){
        stop('The "apply.log" must be logical (TRUE or FALSE).')
    }
    if (!is.logical(check.prps.connectedness)){
        stop('The "check.prps.connectedness" must be logical.')
    }
    if (!is.logical(check.se.obj)){
        stop('The "check.se.obj" must be logical (TRUE or FALSE).')
    }
    if (isTRUE(apply.log)){
        if (pseudo.count < 0){
            stop('The value for "pseudo.count" can not be negative.')
        }
    }
    if (!is.null(regress.out.variables)){
        if (isFALSE(check.se.obj)){
            if (sum(regress.out.variables %in% colnames(colData(se.obj))) != length(regress.out.variables)){
                stop('All or some of the "regress.out.variables" cannot be found in the SummarizedExperiment object.')
            }
        }
        if (main.uv.variable %in% regress.out.variables){
            stop('The "main.uv.variable" can not be in the "regress.out.variables" variables.')
        }
    }
    if (!is.logical(plot.output)){
        stop('The "plot.output" must be logical (TRUE or FALSE).')
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical (TRUE or FALSE).')
    }
    if (is.logical(knn.group.name)){
        stop('The "knn.group.name" must be a character or NULL.')
    }
    if (is.logical(knn.sets.name)){
        stop('The "knn.sets.name" must be a character or NULL.')
    }
    if (is.logical(mnn.group.name)){
        stop('The "mnn.group.name" must be a character or NULL.')
    }
    if (is.logical(mnn.sets.name)){
        stop('The "mnn.sets.name" must be a character or NULL.')
    }
    if (is.logical(prps.group.name)){
        stop('The "prps.group.name" must be a character or NULL.')
    }
    if (is.logical(prps.sets.name)){
        stop('The "prps.sets.name" must be a character or NULL.')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical (TRUE or FALSE).')
    }

    # Assessing the SummarizedExperiment object ####
    if (isTRUE(check.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(main.uv.variable, other.uv.variables, regress.out.variables),
            remove.na = remove.na,
            verbose = verbose
        )
    }
    #
    if (is.logical(samples.to.use)){
        initial.se.obj.a <- se.obj
        se.obj <- se.obj[ , samples.to.use]
    }

    # Assessing and grouping the main unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    initial.variable <- se.obj[[main.uv.variable]]
    initial.variable2 <- se.obj[[main.uv.variable]]
    if (is.numeric(initial.variable)){
        se.obj[[main.uv.variable]] <- groupContinuousVariable(
            se.obj = se.obj,
            variable = main.uv.variable,
            nb.clusters = nb.clusters,
            clustering.method = clustering.method,
            perfix = '_',
            verbose = verbose
        )
        if (isTRUE(select.extreme.groups)){
            printColoredMessage(
                message = paste0(
                    '- Selecting the two subgroups of the ',
                    main.uv.variable,
                    ' variable with highest and lowest values.'),
                color = 'blue',
                verbose = verbose
            )
            initial.se.obj <- se.obj
            max.group <- se.obj[[main.uv.variable]][initial.variable == max(initial.variable)]
            min.group <- se.obj[[main.uv.variable]][initial.variable == min(initial.variable)]
            selected.samples <- se.obj[[main.uv.variable]] %in% c(max.group, min.group)
            se.obj <- se.obj[ , selected.samples]
            se.obj[[main.uv.variable]] <- droplevels(se.obj[[main.uv.variable]])
            initial.variable <- initial.variable[selected.samples]
        }
        if (isFALSE(select.extreme.groups)){
            initial.se.obj <- se.obj
        }
    }
    if (!is.numeric(initial.variable)){
        initial.se.obj <- se.obj
        length.variable <- length(unique(initial.variable))
        if (length.variable == 1){
            stop('To create PRPS, the "main.uv.variable" must have at least two groups/levels.')
        } else if (length.variable > 1){
            printColoredMessage(
                message = paste0(
                    '- The "',
                    main.uv.variable,
                    '" is a categorical variable with ',
                    length(unique(se.obj[[main.uv.variable]])),
                    ' levels.'),
                color = 'blue',
                verbose = verbose
            )
            se.obj[[main.uv.variable]] <- factor(x = se.obj[[main.uv.variable]])
        }
    }

    # Creating PRPS data with KNN and MNN ####
    if (!is.null(other.uv.variables)){
        ## Considering other unwanted variables ####
        printColoredMessage(
            message = '- Creating PRPS data by considering other unwanted variables.',
            color = 'magenta',
            verbose = verbose
            )
        ## Grouping the other unwanted variables ####
        printColoredMessage(
            message = '- Assessing and grouping the other unwanted variable(s):',
            color = 'blue',
            verbose = verbose
            )
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = other.uv.variables,
            nb.clusters = nb.other.uv.clusters,
            clustering.method = other.uv.clustering.method,
            check.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
        all.uv.groups <- data.frame(
            main.uv = se.obj[[main.uv.variable]],
            other.uv = homo.uv.groups
            )
        min.sample.size <- max(
            min.sample.for.ps,
            nb.mnn,
            nb.knn
            )
        ## Finding and plotting covered batches ####
        covered.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                subgroups.size <- findRepeatingPatterns(
                    vec = all.uv.groups$main.uv[all.uv.groups$other.uv == x],
                    n.repeat = min.sample.size
                )
            })
        names(covered.batches) <- unique(all.uv.groups$other.uv)
        covered.batches.plot <- table(all.uv.groups$main.uv, all.uv.groups$other.uv) %>%
            data.frame(.) %>%
            dplyr::mutate(selected = Freq >= min.sample.size) %>%
            ggplot(data = ., aes(x = Var2, y = Var1, color = selected)) +
            geom_point(size = 4) +
            geom_text(aes(label = Freq , hjust = 0.5, vjust = 0.5), color = 'black', size = 5) +
            xlab('Homogeneous groups (other unwanted variables)') +
            ylab('Main unwanted variable') +
            theme_bw() +
            theme(
                legend.key = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.y = element_text(size = 14),
                axis.text.x = element_text(size = 14, angle = 90, vjust = 1, hjust = 1),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 18),
                strip.text.y = element_text(size = 0)
            )
        if (isTRUE(plot.output)) print(covered.batches.plot)

        ## Checking covered batches ####
        printColoredMessage(
            message = '-- Checking the distribution of the "main.uv.variable" across the "other.uv.variables":',
            color = 'blue',
            verbose = verbose
            )
        selected.covered.batches <- lapply(
            1:length(covered.batches),
            function(x) length(covered.batches[[x]])
            )
        if (sum(selected.covered.batches == 1) == length(selected.covered.batches)){
            stop(paste0(
                ' Non of the sample groups with respect to the other unwanted variables that have at least ',
                max(min.sample.for.ps, nb.mnn),
                ' samples across at least two sub-groups of the ',
                main.uv.variable,
                ' variable.'))
        }
        if (sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) == 0){
            printColoredMessage(
                message = paste0(
                    '- Non of the sample groups with respect to the other unwanted variables have at least ',
                    max(min.sample.for.ps, nb.mnn),
                    ' samples across all the sub-groups of the "',
                    main.uv.variable,
                    '" variable.'),
                color = 'blue',
                verbose = verbose
            )
            if (isFALSE(check.prps.connectedness)){
                printColoredMessage(
                    message = '-- We recommend applying the "check.prps.connectedness"',
                    color = 'red',
                    verbose = verbose
                )
            } else if (isTRUE(check.prps.connectedness)){
                checkPRPSconnectedness(
                    data.input = table(all.uv.groups$main.uv, all.uv.groups$other.uv),
                    min.samples = c(nb.mnn, min.sample.for.ps),
                    batch.name = main.uv.variable,
                    verbose = verbose
                )
            }
        }
        if (sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) > 0  ){
            printColoredMessage(
                message = paste0(
                    '-- There are ',
                    sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) ,
                    ' groups with respect to the other unwanted variables that have at least ',
                    min.sample.size,
                    ' samples across all sub-groups of the main unwanted variable.'),
                color = 'blue',
                verbose = verbose
            )
        }
        ## Finding K nearest neighbors ####
        printColoredMessage(
            message = '-- Finding k nearest neighbor by applying the findKnn function:',
            color = 'magenta',
            verbose = verbose
            )
        all.possible.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                possible.batch <- findRepeatingPatterns(
                    vec = all.uv.groups[all.uv.groups$other.uv == x, ]$main.uv,
                    n.repeat = max(min.sample.for.ps, nb.mnn)
                )
                if (length(possible.batch) > 1){
                    combn(x = possible.batch , m = 2)
                } else NA

            })
        names(all.possible.batches) <- unique(all.uv.groups$other.uv)
        all.possible.batches <- all.possible.batches[!is.na(all.possible.batches)]
        all.knn <- lapply(
            1:length(all.possible.batches),
            function(x){
                all.knn.samples <- findKnn(
                    se.obj = se.obj[ , homo.uv.groups == names(all.possible.batches)[x] ],
                    assay.name = assay.name,
                    uv.variable = main.uv.variable,
                    data.input = data.input,
                    nb.pcs = nb.pcs,
                    center = center,
                    scale = scale,
                    svd.bsparam = svd.bsparam,
                    clustering.method = clustering.method,
                    nb.clusters = nb.clusters,
                    nb.knn = nb.knn,
                    hvg = hvg,
                    normalization = normalization,
                    regress.out.variables = regress.out.variables,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    check.se.obj = FALSE,
                    remove.na = remove.na,
                    knn.group.name = knn.group.name,
                    knn.sets.name = knn.sets.name,
                    save.se.obj = FALSE,
                    verbose = verbose
                )
                all.knn.samples$other.group <- names(all.possible.batches[x])
                all.knn.samples
            })
        names(all.knn) <- names(all.possible.batches)

        ## Finding mutual nearest neighbors ####
        printColoredMessage(
            message = '-- Finding mutual nearest neighbors by applying the findMnn function:',
            color = 'magenta',
            verbose = verbose
            )
        all.mnn <- lapply(
            1:length(all.possible.batches),
            function(x){
                all.mnn.samples <- findMnn(
                    se.obj = se.obj[ , homo.uv.groups == names(all.possible.batches)[x] ],
                    assay.name = assay.name,
                    uv.variable = main.uv.variable,
                    nb.mnn = nb.mnn,
                    clustering.method = clustering.method,
                    nb.clusters = nb.clusters,
                    data.input = data.input,
                    nb.pcs = nb.pcs,
                    center = center,
                    scale = scale,
                    svd.bsparam = svd.bsparam,
                    normalization = normalization,
                    apply.cosine.norm = apply.cosine.norm,
                    regress.out.variables = regress.out.variables,
                    hvg = hvg,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    mnn.bpparam = mnn.bpparam,
                    mnn.nbparam = mnn.nbparam,
                    check.se.obj = check.se.obj,
                    remove.na = remove.na,
                    plot.output = FALSE,
                    mnn.group.name = mnn.group.name,
                    mnn.sets.name = mnn.sets.name,
                    save.se.obj = FALSE,
                    verbose = verbose
                    )
                all.mnn.samples$mnn$other.group <- names(all.possible.batches[x])
                all.mnn.samples
            })
        ## Plotting all mnn ####
        all.mnn.plots <- lapply(
            1:length(all.possible.batches),
            function(x){
                all.mnn[[x]]$mnn.plot + ggtitle(names(all.possible.batches[x]))
            })
        all.mnn.plots <- ggarrange(plotlist = all.mnn.plots)
        if (isTRUE(plot.output)) print(all.mnn.plots)
        ## Obtaining mnn data ####
        all.mnn <- lapply(
            1:length(all.possible.batches),
            function(x){
                all.mnn[[x]]$mnn
            })
        names(all.mnn) <- names(all.possible.batches)

        ## Matching mnn and knn data ####
        printColoredMessage(
            message = '- Matching the MNN sets with the corresponding KNN sets:',
            color = 'magenta',
            verbose = verbose
            )
        all.prps.sets <- lapply(
            names(all.possible.batches),
            function(i){
                sub.all.knn <- all.knn[[i]]
                sub.all.mnn <- all.mnn[[i]]
                # Applying a sanity check ####
                if (length(unique(sub.all.knn$other.group)) > 1){
                    stop('There something wrong with knn and mnn.')
                }
                if (length(unique(sub.all.mnn$other.group)) > 1){
                    stop('There something wrong with knn and mnn.')
                }
                if (!unique(sub.all.knn$other.group) == unique(sub.all.mnn$other.group)){
                    stop('There something wrong with knn and mnn.')
                }
                sub.prps.sets <- lapply(
                    1:nrow(sub.all.mnn),
                    function(x) {
                        # ps set 1
                        ps.set.1 <- sub.all.knn[ ,   grep('sample.ids', colnames(sub.all.knn))] == sub.all.mnn$sample.ids.1[x]
                        ps.set.1 <- sub.all.knn[rowSums(ps.set.1) > 0 , ]
                        ps.set.1$mnn.sets <- paste0(paste0(
                            sort(c(sub.all.mnn[x , 3], sub.all.mnn[x , 4])),
                            collapse = '_'),
                            '_',
                            x,
                            '_',
                            i
                            )
                        ps.set.1$mnn.sets.data <- paste0(sort(c(sub.all.mnn[x , 1], sub.all.mnn[x , 2])), collapse = '_')
                        if (nrow(ps.set.1) > 1) {
                            ps.set.1 <- ps.set.1[ps.set.1$rank.aver.dist == min(ps.set.1$rank.aver.dist) , ]
                        }
                        # ps set 2
                        ps.set.2 <- sub.all.knn[ ,   grep('sample.ids', colnames(sub.all.knn))] == sub.all.mnn$sample.ids.2[x]
                        ps.set.2 <- sub.all.knn[rowSums(ps.set.2) > 0 , ]
                        ps.set.2$mnn.sets <- paste0(paste0(
                            sort(c(sub.all.mnn[x , 3], sub.all.mnn[x , 4])),
                            collapse = '_'),
                            '_',
                            x,
                            '_',
                            i
                            )
                        ps.set.2$mnn.sets.data <- paste0(sort(c(sub.all.mnn[x , 1], sub.all.mnn[x , 2])), collapse = '_')
                        if (nrow(ps.set.2) > 1) {
                            ps.set.2 <- ps.set.2[ps.set.2$rank.aver.dist == min(ps.set.2$rank.aver.dist) , ]
                        }
                        prps.set <- rbind(ps.set.1, ps.set.2)
                        # Applying a sanity check ####
                        if ( sum(prps.set[1 , grep('sample.ids', colnames(sub.all.knn))] %in% prps.set[2 , grep('sample.ids', colnames(sub.all.knn))]) > 1 ){
                            stop('There something wrong with knn and mnn.')
                        }
                        prps.set
                    })
                sub.prps.sets <- do.call(rbind, sub.prps.sets)
            })
        names(all.prps.sets) <- names(all.possible.batches)

        if (is.null(all.prps.sets)) {
            stop('PRPS cannot be created. You may want to increase the value of the mnn.')
        }
        ## Applying a sanity check ####
        printColoredMessage(
            message =  '- Applying a sanity check on the mnn and knn.',
            color = 'blue',
            verbose = verbose
        )
        sanity.check <- unlist(lapply(
            names(all.prps.sets),
            function(x){
                nrow(all.prps.sets[[x]]) != 2*nrow(all.mnn[[x]])
            }))
        if (sum(sanity.check) == 0){
            printColoredMessage(
                message =  '- The rows of the of the matched MNN and KNN is correct.',
                color = 'blue',
                verbose = verbose
            )
        } else {
            stop('For individual MNN set, the corresponding KNN sets cannot be found. Check the the input.')
        }

        ## Adding the average distances of each knn sets for each PRPS set and then rank them ####
        printColoredMessage(
            message = '- Averaging the distances of each knn sets for each MNN set and then rank them:',
            color = 'blue',
            verbose = verbose
            )
        all.prps.sets <- do.call(rbind, all.prps.sets)
        all.prps.sets$aver.mnn.sets <- unlist(lapply(
            seq(1, nrow(all.prps.sets), 2),
            function(x)
                rep(mean(all.prps.sets$aver.dist[x:(x + 1)]), 2))
            )
        set.seed(2233)
        all.prps.sets$rank.aver.mnn.sets <- rank(
            x = all.prps.sets$aver.mnn.sets,
            ties.method = 'random'
            )
        ## Filtering PRPS sets ####
        if (isTRUE(filter.prps.sets)) {
            printColoredMessage(
                message = '- Filtering the PRPS sets across each pair of batches:',
                color = 'orange',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    '- The maximum number of PRPS sets across each pairs of batches is ',
                    max.prps.sets,
                    '.'),
                color = 'blue',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    '- The PRPS sets will be filtered based on the distances between each knn sets.'),
                color = 'blue',
                verbose = verbose
                )
            all.prps.sets <- lapply(
                unique(all.prps.sets$mnn.sets.data),
                function(x) {
                    temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                    if (length(unique(temp.prps.set$mnn.sets)) >= max.prps.sets) {
                        printColoredMessage(
                            message = paste0(
                                '* The number of PRPS sets across the batches "',
                                x,
                                '" is ',
                                length(unique(temp.prps.set$mnn.sets)), '.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- arrange(temp.prps.set, aver.mnn.sets, mnn.sets)
                        printColoredMessage(
                            message = paste0(
                                '* ',
                                length(unique(temp.prps.set$mnn.sets)) - max.prps.sets,
                                ' PRPS sets are removed.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- temp.prps.set[1:c(2 * max.prps.sets) , ]
                    } else if (length(unique(temp.prps.set$mnn.sets))  < max.prps.sets) {
                        printColoredMessage(
                            message = paste0(
                                '* The number of PRPS sets across the batches "',
                                x,
                                '" is ',
                                length(unique(temp.prps.set$mnn.sets)) , '.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                    }
                    return(temp.prps.set)
                })
            all.prps.sets <- do.call(rbind, all.prps.sets)
        }
        printColoredMessage(
            message = paste0(
                '- ',
                length(unique(all.prps.sets$mnn.sets)),
                ' PRPS stes are found in total.'),
            color = 'blue',
            verbose = verbose
            )
        ## Creating PRPS data ####
        printColoredMessage(
            message = '-- Creating PRPS data matrix:',
            color = 'magenta',
            verbose = verbose
            )
        ## Applying log ####
        printColoredMessage(
            message = '- Applying data log transformation before creating the PRPS expression data:',
            color = 'blue',
            verbose = verbose
            )
        if (isTRUE(apply.log.for.prps) & !is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0(
                    '- Applying log2 on the "',
                    assay.name,
                    '" + ',
                    pseudo.count,
                    ' (pseudo.count)  data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        } else if (isTRUE(apply.log.for.prps) & is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0(
                    'Applying log2 on the "',
                    assay.name,
                    '" data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log.for.prps)) {
            printColoredMessage(
                message = paste0(
                    'The "',
                    assay.name,
                    '" data will be used without any log transformation.' ),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- assay(x = se.obj, i = assay.name)
        }
        printColoredMessage(
            message = '- Aeveraging samples to create pseudo samples:',
            color = 'blue',
            verbose = verbose
        )
        prps.data <- lapply(
            unique(all.prps.sets$mnn.sets),
            function(i){
                tep.mnn.sets <- all.prps.sets[all.prps.sets$mnn.sets == i , ]
                set.a <- tep.mnn.sets[1 , grep('sample.ids', colnames(all.prps.sets))]
                set.a <- rowMeans(expr.data[ , unlist(as.vector(set.a[1 , ])) ])
                set.b <- tep.mnn.sets[2 , grep('sample.ids', colnames(all.prps.sets))]
                set.b <- rowMeans(expr.data[ , unlist(as.vector(set.b[1 , ])) ])
                prps <- cbind(set.a, set.b)
                colnames(prps) <- rep(paste(
                    main.uv.variable,
                    i,
                    sep = '_'),  2)
                prps
            })
        prps.data <- do.call(cbind, prps.data)

        ## Applying a sanity check ####
        if (!sum(table(colnames(prps.data)) == 2) == ncol(prps.data) / 2) {
            stop('There someting wrong with PRPS sets.')
        }
        se.obj[[main.uv.variable]] <- initial.variable

        ## Plotting the PRPS map ####
        prps.map <- lapply(
            unique(all.prps.sets$mnn.sets),
            function(i){
                tep.mnn.sets <- all.prps.sets[all.prps.sets$mnn.sets == i , ]
                set.a <- tep.mnn.sets[1 , grep('sample.ids', colnames(all.prps.sets))]
                set.a <- unlist(as.vector(set.a[1 , ]))
                set.b <- tep.mnn.sets[2 , grep('sample.ids', colnames(all.prps.sets))]
                set.b <- unlist(as.vector(set.b[1 , ]))
                initial.variable.set.a <- initial.variable[colnames(se.obj) %in% set.a]
                initial.variable.set.b <- initial.variable[colnames(se.obj) %in% set.b]
                data.frame(
                    group1 = initial.variable.set.a,
                    group2 = initial.variable.set.b,
                    set = rep(i, min.sample.for.ps)
                )
            })
        prps.map <- do.call(rbind, prps.map)
        prps.map <- pivot_longer(prps.map, -set, names_to = 'group', values_to = 'var')
        prps.map$group2 <- 'PRPS sets'
        all.uv.variable <- data.frame(
            set = main.uv.variable,
            group = initial.se.obj[[main.uv.variable]],
            var = initial.variable2,
            group2 = 'UV'
        )
        prps.map <- rbind(prps.map, all.uv.variable)
        prps.map$group2 <- factor(x = prps.map$group2 , levels = c('UV', 'PRPS sets'))
        if (is.numeric(initial.variable)){
            prps.map.plot <- ggplot(prps.map, aes(x = set, y = var, color = group)) +
                geom_boxplot() +
                geom_point(size = 2) +
                xlab('') +
                ylab(main.uv.variable) +
                scale_color_manual(
                    values = c('darkgreen', 'tomato', c(viridis(nb.clusters)))
                ) +
                facet_grid(.~group2, scales = 'free', space = 'free') +
                scale_x_discrete(expand = c(0, 0.5)) +
                theme_bw() +
                theme(
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    axis.line = element_line(colour = 'black', linewidth = .85),
                    axis.title.x = element_text(size = 16),
                    strip.text.x.top = element_text(size = 20),
                    axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5 , vjust = 0.5),
                    axis.text.y = element_text(size = 12),
                    legend.position = 'right') +
                guides(color = guide_legend(title = "Groups"))
            if (isTRUE(verbose)) print(prps.map.plot)
            se.obj[[main.uv.variable]] <- initial.variable
        }
        if ((!is.numeric(initial.variable))){
            prps.map$group <- factor(
                x = prps.map$group ,
                levels = c('group1', 'group2', names(table(initial.variable)))
            )
            prps.map.plot <- ggplot(prps.map, aes(x = set, y = var, color = group)) +
                geom_point(size = 2) +
                xlab('') +
                ylab(main.uv.variable) +
                scale_color_manual(
                    values = c('darkgreen', 'tomato', c(viridis(n = length(unique(se.obj[[main.uv.variable]])))))
                ) +
                facet_grid(.~group2, scales = 'free', space = 'free') +
                scale_x_discrete(expand = c(0, 0.5)) +
                theme_bw() +
                theme(
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    axis.line = element_line(colour = 'black', linewidth = .85),
                    axis.title.x = element_text(size = 16),
                    strip.text.x.top = element_text(size = 20),
                    axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5 , vjust = 0.5),
                    axis.text.y = element_text(size = 12),
                    legend.position = 'right') +
                guides(color = guide_legend(title = "Groups"))
            if (isTRUE(verbose)) print(prps.map.plot)
            se.obj[[main.uv.variable]] <- initial.variable
        }
    }

    ## Considering only main unwanted variable ####
    if (is.null(other.uv.variables)){
        ### Checking sample sizes of each sub group ####
        ### KNN ####
        sub.group.sample.size.knn <- findRepeatingPatterns(
            vec = se.obj[[main.uv.variable]],
            n.repeat = nb.knn + 1
            )
        if (length(sub.group.sample.size.knn) == 0){
            stop(paste0(
                'All subgroups of the unwanted variable have less than ',
                nb.knn + 1,
                ' (nb.knn + 1) samples. KNN cannot be found.')
                )
        } else if (length(sub.group.sample.size.knn) != length(unique(se.obj[[main.uv.variable]])) ){
            printColoredMessage(
                message = paste0(
                    'All or some subgroups of the unwanted variable have less than ',
                    nb.knn + 1,
                    ' (nb.knn + 1) samples. Then KNN for those sub-groups cannot be created.'),
                color = 'red',
                verbose = verbose
                )
        } else {
            printColoredMessage(
                message = paste0(
                    '- All the sub-groups of the unwanted variable have at least ',
                    nb.knn + 1,
                    ' (nb.knn + 1) samples.'),
                color = 'blue',
                verbose = verbose
                )
        }
        ### MNN ####
        sub.group.sample.size.mnn <- findRepeatingPatterns(
            vec = se.obj[[main.uv.variable]],
            n.repeat = nb.mnn + 1
            )
        if (length(sub.group.sample.size.mnn) == 0){
            stop(paste0(
                'All subgroups of the unwanted variable have less than ',
                nb.mnn + 1,
                ' (nb.mnn + 1) samples. MNN cannot be found.')
                )
        } else if (length(sub.group.sample.size.mnn) != length(unique(se.obj[[main.uv.variable]])) ){
            printColoredMessage(
                message = paste0(
                    'All or some subgroups of the unwanted variable have less than ',
                    nb.mnn + 1,
                    ' (nb.mnn + 1) samples.'),
                color = 'red',
                verbose = verbose
            )
        } else {
            printColoredMessage(
                message = paste0(
                    '- All the subgroups of the unwanted variable have at least, ',
                    nb.mnn + 1,
                    ' (nb.mnn + 1) samples.'),
                color = 'blue',
                verbose = verbose
            )
        }
        ### Finding k nearest neighbor ####
        printColoredMessage(
            message = '-- Finding k nearest neighbor by applying the findKnn function:',
            color = 'magenta',
            verbose = verbose
        )
        all.knn <- findKnn(
            se.obj = se.obj,
            assay.name = assay.name,
            uv.variable = main.uv.variable,
            data.input = data.input,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            svd.bsparam = svd.bsparam,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            nb.knn = nb.knn,
            hvg = hvg,
            normalization = normalization,
            regress.out.variables = regress.out.variables,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            check.se.obj = check.se.obj,
            remove.na = remove.na,
            knn.group.name = knn.group.name,
            knn.sets.name = knn.sets.name,
            save.se.obj = FALSE,
            verbose = verbose
            )
        if (isTRUE(save.se.obj)){
            if(is.null(knn.group.name)){
                knn.group.name <- main.uv.variable
            }
            if (is.null(knn.sets.name)){
                if (is.numeric(se.obj[[main.uv.variable]])){
                    knn.sets.name <- paste0(assay.name, '|', nb.clusters ,'groups|', nb.knn, 'knn')
                } else knn.sets.name <- paste0(assay.name, '|', length(unique(se.obj[[main.uv.variable]])) ,'groups|', nb.knn, 'knn')
            }
            if (!'KnnMnn' %in% names(se.obj@metadata)) {
                se.obj@metadata[['KnnMnn']] <- list()
            }
            if (!'Knn' %in% names(se.obj@metadata[['KnnMnn']])) {
                se.obj@metadata[['KnnMnn']][['knn']] <- list()
            }
            if (!knn.group.name %in% names(se.obj@metadata[['KnnMnn']][['knn']])) {
                se.obj@metadata[['KnnMnn']][['knn']][[knn.group.name]] <- list()
            }
            if (!knn.sets.name %in% names(se.obj@metadata[['KnnMnn']][['knn']][[knn.group.name]])) {
                se.obj@metadata[['KnnMnn']][['knn']][[knn.group.name]][[knn.sets.name]] <- list()
            }
            se.obj@metadata[['KnnMnn']][['knn']][[knn.group.name]][[knn.sets.name]] <- all.knn
        }
        ### Finding mutual nearest neighbor ####
        printColoredMessage(
            message = '-- Finding mutual nearest neighbors by applying the findMnn function:',
            color = 'magenta',
            verbose = verbose
            )
        all.mnn.res <- findMnn(
            se.obj = se.obj,
            assay.name = assay.name,
            uv.variable = main.uv.variable,
            nb.mnn = nb.mnn,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            data.input = data.input,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            svd.bsparam = svd.bsparam,
            normalization = normalization,
            apply.cosine.norm = apply.cosine.norm,
            regress.out.variables = regress.out.variables,
            hvg = hvg,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            mnn.bpparam = mnn.bpparam,
            mnn.nbparam = mnn.nbparam,
            check.se.obj = check.se.obj,
            remove.na = remove.na,
            plot.output = plot.output,
            mnn.group.name = mnn.group.name,
            mnn.sets.name = mnn.sets.name,
            save.se.obj = FALSE,
            verbose = verbose
            )
        all.mnn <- all.mnn.res$mnn
        if (isTRUE(save.se.obj)) {
            if(is.null(mnn.group.name)){
                mnn.group.name <- main.uv.variable
            }
            if (is.null(mnn.sets.name)){
                if (is.numeric(se.obj[[main.uv.variable]])){
                    mnn.sets.name <- paste0(assay.name, '|', nb.clusters ,'groups|', nb.mnn, 'mnn')
                } else mnn.sets.name <- paste0(assay.name, '|', length(unique(se.obj[[main.uv.variable]])) ,'groups|', nb.mnn, 'mnn')
            }
            if (!'KnnMnn' %in% names(se.obj@metadata)) {
                se.obj@metadata[['KnnMnn']] <- list()
            }
            if (!'Mnn' %in% names(se.obj@metadata[['KnnMnn']])) {
                se.obj@metadata[['KnnMnn']][['Mnn']] <- list()
            }
            if (!mnn.group.name %in% names(se.obj@metadata[['KnnMnn']][['Mnn']])) {
                se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]] <- list()
            }
            if (!mnn.sets.name %in% names(se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.sets.name]])) {
                se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]][[mnn.sets.name]] <- list()
            }
            if (!'data' %in% names(se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.sets.name]])) {
                se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]][[mnn.sets.name]][['data']] <- list()
            }
            se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]][[mnn.sets.name]][['data']]  <- all.mnn
            if (!'plot' %in% names(se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.sets.name]])) {
                se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]][[mnn.sets.name]][['plot']] <- list()
            }
            se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]][[mnn.sets.name]][['plot']] <- all.mnn.res$mnn.plot
        }

        ###  Matching KNN and MNN  ####
        printColoredMessage(
            message = '-- Finding all possible similar samples across batches:',
            color = 'magenta',
            verbose = verbose
            )
        ### Finding the knn for each mnn set ####
        printColoredMessage(
            message = '- Matching the results of the KNN and MNN data:',
            color = 'orange',
            verbose = verbose
            )
        printColoredMessage(
            message = '* Matching the MNN sets with the corresponding KNN sets:',
            color = 'blue',
            verbose = verbose
            )
        ### Finding all the prps sets ####
        all.prps.sets <- lapply(
            1:nrow(all.mnn),
            function(x) {
                # ps set 1
                ps.set.1 <- all.knn[ ,   grep('sample.ids', colnames(all.knn))] == all.mnn$sample.ids.1[x]
                ps.set.1 <- all.knn[rowSums(ps.set.1) > 0 , ]
                ps.set.1$mnn.sets <- paste0(paste0(
                    sort(c(all.mnn[x , 3], all.mnn[x , 4])),
                    collapse = '_'),
                    '||',
                    x
                    )
                ps.set.1$mnn.sets.data <- paste0(sort(c(all.mnn[x , 1], all.mnn[x , 2])), collapse = '||')
                if (nrow(ps.set.1) > 1) {
                    ps.set.1 <- ps.set.1[ps.set.1$rank.aver.dist == min(ps.set.1$rank.aver.dist) , ]
                }
                # ps set 2
                ps.set.2 <- all.knn[ ,   grep('sample.ids', colnames(all.knn))] == all.mnn$sample.ids.2[x]
                ps.set.2 <- all.knn[rowSums(ps.set.2) > 0 , ]
                ps.set.2$mnn.sets <- paste0(paste0(
                    sort(c(all.mnn[x , 3], all.mnn[x , 4])),
                    collapse = '_'),
                    '||',
                    x
                )
                ps.set.2$mnn.sets.data <- paste0(sort(c(all.mnn[x , 1], all.mnn[x , 2])), collapse = '||')
                if (nrow(ps.set.2) > 1) {
                    ps.set.2 <- ps.set.2[ps.set.2$rank.aver.dist == min(ps.set.2$rank.aver.dist) , ]
                }
                prps.set <- rbind(ps.set.1, ps.set.2)
                # Applying a sanity check ####
                if ( sum(prps.set[1 , grep('sample.ids', colnames(all.knn))] %in% prps.set[2 , grep('sample.ids', colnames(all.knn))]) > 1 ){
                    stop('There something wrong with knn and mnn findings, please check.')
                }
                prps.set
            })
        all.prps.sets <- do.call(rbind, all.prps.sets)

        if (is.null(all.prps.sets)) {
            stop('PRPS cannot be created. You may want to increase the value of the mnn.')
        }
        ### Applying a sanity check ####
        if (nrow(all.prps.sets) == 2*nrow(all.mnn)){
            printColoredMessage(
                message = paste0(
                    '* The nrow of the matched MNN and KNN is ',
                    nrow(all.prps.sets),
                    '.'),
                color = 'blue',
                verbose = verbose
            )
        } else {
            stop('For individual MNN set, the corresponding KNN sets cannot be found. Check the the input.')
        }

        ### Adding the average of the knn sets for each PRPS set and then rank them ####
        printColoredMessage(
            message = '- Averaging the average distances of each knn for individual MNN set and then rank them:',
            color = 'blue',
            verbose = verbose
        )
        all.prps.sets$aver.mnn.sets <- unlist(lapply(
            seq(1, nrow(all.prps.sets), 2),
            function(x)
                rep(mean(all.prps.sets$aver.dist[x:(x + 1)]), 2))
        )
        set.seed(2233)
        all.prps.sets$rank.aver.mnn.sets <- rank(
            x = all.prps.sets$aver.mnn.sets,
            ties.method = 'random'
        )

        ### Filtering PRPS sets ####
        if (isTRUE(filter.prps.sets)) {
            printColoredMessage(
                message = '- Filtering the PRPS sets across each pair of batches:',
                color = 'orange',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    '- The maximum number of PRPS sets across each pairs of batches is ',
                    max.prps.sets,
                    '.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    '- The PRPS sets will be filtered based on the average distances between each knn sets.'),
                color = 'blue',
                verbose = verbose
            )
            all.prps.sets <- lapply(
                unique(all.prps.sets$mnn.sets.data),
                function(x) {
                    temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                    if (length(unique(temp.prps.set$mnn.sets)) >= max.prps.sets) {
                        printColoredMessage(
                            message = paste0(
                                '* The number of PRPS sets across the batches "',
                                x,
                                '" is ',
                                length(unique(temp.prps.set$mnn.sets)), '.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- arrange(temp.prps.set, aver.mnn.sets, mnn.sets)
                        printColoredMessage(
                            message = paste0(
                                '* ',
                                length(unique(temp.prps.set$mnn.sets)) - max.prps.sets,
                                ' PRPS sets are removed.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- temp.prps.set[1:c(2 * max.prps.sets) , ]
                    } else if (length(unique(temp.prps.set$mnn.sets))  < max.prps.sets) {
                        printColoredMessage(
                            message = paste0(
                                '* The number of PRPS sets across the batches "',
                                x,
                                '" is ',
                                length(unique(temp.prps.set$mnn.sets)) , '.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                    }
                    return(temp.prps.set)
                })
            all.prps.sets <- do.call(rbind, all.prps.sets)
        }
        printColoredMessage(
            message = paste0(
                '- ',
                length(unique(all.prps.sets$mnn.sets)),
                ' PRPS stes are kept in total.'),
            color = 'blue',
            verbose = verbose
        )

        ### Ccreating PRPS data ####
        printColoredMessage(
            message = '-- Creating PRPS data:',
            color = 'magenta',
            verbose = verbose
        )
        ### Applying log ####
        printColoredMessage(
            message = '- Applying data log transformation before creating the PRPS expression data:',
            color = 'blue',
            verbose = verbose
        )
        if (isTRUE(apply.log.for.prps) & !is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0(
                    '- Applying log2 on the "',
                    assay.name,
                    '" + ',
                    pseudo.count,
                    ' (pseudo.count)  data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        } else if (isTRUE(apply.log.for.prps) & is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0(
                    'Applying log2 on the "',
                    assay.name,
                    '" data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log.for.prps)) {
            printColoredMessage(
                message = paste0(
                    'The "',
                    assay.name,
                    '" data will be used without any log transformation.' ),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- assay(x = se.obj, i = assay.name)
        }
        printColoredMessage(
            message = '- Aeveraging samples to create pseudo samples:',
            color = 'blue',
            verbose = verbose
        )
        prps.data <- lapply(
            unique(all.prps.sets$mnn.sets),
            function(x) {
                temp.prps <- all.prps.sets[all.prps.sets$mnn.sets == x, ]
                index.a <- unlist(unname(temp.prps[1, grep('sample.ids', colnames(temp.prps))]))
                index.b <- unlist(unname(temp.prps[2, grep('sample.ids', colnames(temp.prps))]))
                prps.a <- rowMeans(expr.data[, index.a])
                prps.b <- rowMeans(expr.data[, index.b])
                prps <- cbind(prps.a, prps.b)
                colnames(prps) <- paste(main.uv.variable, temp.prps$mnn.sets, sep = '_')
                return(prps)
            })
        prps.data <- do.call(cbind, prps.data)
        ### Applying a sanity check ####
        if (!sum(table(colnames(prps.data)) == 2) == ncol(prps.data) / 2) {
            stop('There someting wrong with PRPS sets.')
        }
        ## Plotting the PRPS map ####
        prps.map <- lapply(
            unique(all.prps.sets$mnn.sets),
            function(i){
                tep.mnn.sets <- all.prps.sets[all.prps.sets$mnn.sets == i , ]
                set.a <- tep.mnn.sets[1 , grep('sample.ids', colnames(all.prps.sets))]
                set.a <- unlist(as.vector(set.a[1 , ]))
                set.b <- tep.mnn.sets[2 , grep('sample.ids', colnames(all.prps.sets))]
                set.b <- unlist(as.vector(set.b[1 , ]))
                initial.variable.set.a <- initial.variable[colnames(se.obj) %in% set.a]
                initial.variable.set.b <- initial.variable[colnames(se.obj) %in% set.b]
                data.frame(
                    group1 = initial.variable.set.a,
                    group2 = initial.variable.set.b,
                    set = rep(i, min.sample.for.ps)
                )
            })
        prps.map <- do.call(rbind, prps.map)
        prps.map <- pivot_longer(prps.map, -set, names_to = 'group', values_to = 'var')
        prps.map$group2 <- 'PRPS sets'
        all.uv.variable <- data.frame(
            set = main.uv.variable,
            group = initial.se.obj[[main.uv.variable]],
            var = initial.variable2,
            group2 = 'UV'
            )
        prps.map <- rbind(prps.map, all.uv.variable)
        prps.map$group2 <- factor(x = prps.map$group2 , levels = c('UV', 'PRPS sets'))
        if (is.numeric(initial.variable)){
            prps.map.plot <- ggplot(prps.map, aes(x = set, y = var, color = group)) +
                geom_boxplot() +
                geom_point(size = 2) +
                xlab('') +
                ylab(main.uv.variable) +
                scale_color_manual(
                    values = c('darkgreen', 'tomato', c(viridis(nb.clusters)))
                ) +
                facet_grid(.~group2, scales = 'free', space = 'free') +
                scale_x_discrete(expand = c(0, 0.5)) +
                theme_bw() +
                theme(
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    axis.line = element_line(colour = 'black', linewidth = .85),
                    axis.title.x = element_text(size = 16),
                    strip.text.x.top = element_text(size = 20),
                    axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5 , vjust = 0.5),
                    axis.text.y = element_text(size = 12),
                    legend.position = 'right') +
                guides(color = guide_legend(title = "Groups"))
            if (isTRUE(verbose)) print(prps.map.plot)
            se.obj[[main.uv.variable]] <- initial.variable
        }
        if (!is.numeric(initial.variable) ){
            prps.map$group <- factor(
                x = prps.map$group ,
                levels = c('group1', 'group2', names(table(initial.variable)))
                )
            prps.map.plot <- ggplot(prps.map, aes(x = set, y = var, color = group)) +
                geom_point(size = 2) +
                xlab('') +
                ylab(main.uv.variable) +
                scale_color_manual(
                    values = c('darkgreen', 'tomato', c(viridis(n = length(unique(se.obj[[main.uv.variable]])))))
                ) +
                facet_grid(.~group2, scales = 'free', space = 'free') +
                scale_x_discrete(expand = c(0, 0.5)) +
                theme_bw() +
                theme(
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    axis.line = element_line(colour = 'black', linewidth = .85),
                    axis.title.x = element_text(size = 16),
                    strip.text.x.top = element_text(size = 20),
                    axis.title.y = element_text(size = 16),
                    axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5 , vjust = 0.5),
                    axis.text.y = element_text(size = 12),
                    legend.position = 'right') +
                guides(color = guide_legend(title = "Groups"))
            if (isTRUE(verbose)) print(prps.map.plot)
            se.obj[[main.uv.variable]] <- initial.variable
        }
    }
    # Saving the results ####
    ## Selecting prps.sets.name ####
    prps.sets.name <- paste0(main.uv.variable, '|', 'KnnMnn', '|', assay.name)
    if (is.null(prps.group.name)) {
        prps.group.name <- main.uv.variable
    }
    printColoredMessage(message = '-- Saving the PRPS data',
                        color = 'magenta',
                        verbose = verbose)
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
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.data']] <- prps.data

        # plot
        if (!'prps.map.plot' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.map.plot']] <- list()
        }
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.map.plot']]  <- prps.map.plot

        printColoredMessage(message = '------------The createPrPsByKnnMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## Outputting the PRPS data as matrix ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The createPrPsByKnnMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(list(prps.data = prps.data, mnn.data = all.mnn, knn.data = all.knn, prps.map.plot = prps.map.plot))
    }
}



