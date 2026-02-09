#' Create PRPS sets using mutual nearest neighbors in RNA-seq data.
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
#' @param other.uv.variables Character. A character or character vector representing the name(s) of the columns of
#' unwanted variable(s) within the sample annotation (colData) of the SummarizedExperiment object. These can be categorical,
#' continuous, or a combination. These variables will be considered when generating PRPS sets for the `main.uv.variable`
#' to help avoid potential contamination. The default is set to `NULL`
#' @param data.input Character. Specifies the data type used to identify nearest neighbours. The options are `pcs` and
#' `expr`. `pcs` refers to principal components, while `expr` refers to gene expression data. The default is set to `pcs`.
#' @param bio.dims Numeric. A numeric value specifying the number of leading principal components to use for identifying
#' nearest neighbours when `data.input` is set to `pcs`. The default is set to 5.
#' @param clustering.method Character. A character string indicating the choice of clustering method for grouping the
#' `uv.variable` if a continuous variable is provided. Options include `kmeans`, `cut`, and `quantile`. The default is set
#' to `kmeans`.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the `uv.variable` is a
#' continuous variable. The default is set to 3.
#' @param other.uv.clustering.method Character. A character indicating which clustering method should be used to
#' group each continuous unwanted variable, if specified in `other.uv.variables`. Options include `kmeans`, `cut`,
#' and `quantile`. The default is set to `kmeans`. See createHomogeneousUVGroups() for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous
#' unwanted variable specified in the `other.uv.variables`. The default is set to 3.
#' @param select.extreme.groups Logical. Indicates whether to select only the extreme groups e.g., highest and lowest
#' clusters, when the `uv.variable` is continuous. Default is set to `TRUE`. This will increase the variation between
#' PR sets in order to better capture the unwanted variation.
#' @param min.sample.for.ps Numeric. Minimum number of samples required for pseudo-replicate creation. The default is set
#' to 3.
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
#' @param samples.to.use Logical, Character, or Integer. Specifies which samples to use for the analysis. If `all`, all
#' samples are used.
#' @param normalization Character. A character that indicates which normalization method should be applied on the
#' data before finding the knn. Options are: `CPM`, `TMM`, `upper`, `median`, `full`, and `VST`. The default is set to
#' `cpm`. If set to `NULL`, no normalization will be applied. See the applyOtherNormalizations() function for more details.
#' @param apply.cosine.norm Logical. Indicates whether cosine normalization should be applied before finding MNN. Default
#' is set to `TRUE`.
#' @param regress.out.variables Character. A character or a vector of character that indicate the column name(s) in the
#' sample annotation in the SummarizedExperiment object. These variables will be regressed out from the data before
#' finding MNN. The default is set to `NULL`, indicating that regression will not be applied.
#' @param regress.out.rle.med Logical. Indicates whether to regress out the RLE medians of the expression data during data
#' pre-processing. The default is set to `FALSE`.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not for down-stream analysis.
#' The default is set to `TRUE`.
#' @param apply.log.for.prps Logical. Indicates whether to apply log transformation to gene expression data before averaging
#' to create pseudo-samples (PS). The default is set to `FALSE`. Note that if the selected data (assay) specified in
#' `assay.name` is already log-transformed, this argument must be set to `FALSE`.
#' @param pseudo.count Numeric. A positive numeric value as a pseudo count to be added to all measurements of the specified
#' assay(data) before applying log transformation to avoid -Inf for measurements that are equal to 0. The default is set
#' to 1.
#' @param mnn.bpparam Character. A BiocParallelParam object specifying how palatalization should be performed to find MNN.
#' The default is set to `SerialParam()`. We refer to the `findMutualNN()` function from the **BiocNeighbors** R package.
#' @param mnn.nbparam Character. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' The default is set to `KmknnParam()`. We refer to the `findMutualNN()` function from the `BiocNeighbors` R package.
#' @param create.prps.map Logical. Whether to create a PRPS map. The default is set to `TRUE`.
#' @param plot.output Logical. If `TRUE`, the function plots the distribution of MNN across the batches and PRPS sets
#' across the `main.uv.variable`.
#' @param prps.group.name Character. A character specifying the name of the prps.group.name to which the current KNN belong.
#' If set to `NULL`, the function will automatically assign a name using  `main.uv.variable`.
#' @param prps.sets.name Character. A character specifying the name of the output file to be saved in the metadata
#' of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0(uv.variable, '|', assay.name)`.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default is set
#' to `TRUE`. See the `checkSeObj()` function for more details.
#' @param remove.na Character. To remove NA or missing values from the assay (data) or not. The options are `assays` and
#' `none`. The default is set to `assays`, so all the NA or missing values from the assay(s) will be removed before computing
#' performing any down-stream analysis. See the `checkSeObj()` function for more details.
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
#' @importFrom irlba prcomp_irlba
#' @importFrom RANN nn2
#' @export

createPrPsUnsupervisedByMnn <- function(
        se.obj,
        assay.name,
        main.uv.variable,
        other.uv.variables = NULL,
        data.input = 'pcs',
        bio.dims = 5,
        clustering.method = 'cut',
        nb.clusters = 3,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        select.extreme.groups = FALSE,
        min.sample.for.ps = 3,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        min.batches.to.cover = 'all',
        check.prps.connectedness = TRUE,
        nb.mnn = 3,
        hvg = NULL,
        samples.to.use = 'all',
        normalization = 'CPM',
        apply.cosine.norm = FALSE,
        regress.out.variables = NULL,
        regress.out.rle.med = FALSE,
        apply.log = TRUE,
        apply.log.for.prps = FALSE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        create.prps.map = TRUE,
        plot.output = TRUE,
        prps.group.name = NULL,
        prps.sets.name = NULL,
        check.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The createPrPsByMnn function starts:',
        color = 'white',
        verbose = verbose
        )
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
    if (is.numeric(colData(se.obj))){
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
    if (is.logical(prps.sets.name)){
        stop('The "prps.sets.name" must be a character or NULL.')
    }
    if (is.logical(prps.group.name)){
        stop('The "prps.group.name" must be a character or NULL.')
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
            verbose = verbose)
    }

    # Selecting a subset of specified samples ####
    if (is.logical(samples.to.use)){
        initial.se.obj.a <- se.obj
        se.obj <- se.obj[ , samples.to.use]
    }

    # Applying data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Applying data normalization, transformation and regression:',
        color = 'magenta',
        verbose = verbose
        )
    norm.data <- preProcessData(
        se.obj = se.obj,
        assay.name = 'RawCount',
        normalization = normalization,
        regress.out.variables = regress.out.variables,
        regress.out.rle.med = regress.out.rle.med,
        apply.log = apply.log,
        pseudo.count = pseudo.count,
        check.se.obj = FALSE,
        remove.na = 'none',
        verbose = verbose
        )
    sample.names <- colnames(norm.data)

    # Selecting the input data ####
    printColoredMessage(
        message = '-- Selecting the data input for MNN identification:',
        color = 'magenta',
        verbose = verbose
        )
    ## Selecting PCS ####
    if (data.input == 'pcs'){
        printColoredMessage(
            message = paste0(
                '- Performing PCA and using',
                ' PCs as data input for MNN identification.'),
            color = 'blue',
            verbose = verbose
            )
        ### Using highly variable genes to perform PCS ####
        if (!is.null(hvg)){
            printColoredMessage(
                message = '- Performing PCA using the provided highly variable genes.',
                color = 'blue',
                verbose = verbose
                )
            norm.data <- irlba::prcomp_irlba(
                x = t(norm.data[hvg , ]),
                n = bio.dims,
                scale = FALSE,
                center = TRUE
                )
            ### Using provided number of PCs ####
            if (!is.null(bio.dims)){
                printColoredMessage(
                    message = paste0(
                        '- Selection the first ',
                        bio.dims,
                        ' PCs as data input for MNN identification.'),
                    color = 'blue',
                    verbose = verbose
                    )
                norm.data <- as.data.frame(t(norm.data$x[ , 1:bio.dims]))
                colnames(norm.data) <- sample.names
            }
            ### Using all PCs ####
            if (is.null(bio.dims)){
                printColoredMessage(
                    message = paste0(
                        '- Using all PCs as data input for MNN identification:'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- as.data.frame(t(norm.data$x))
                colnames(norm.data) <- sample.names
            }
        }
        ## Using all genes to perform PCS ####
        if (is.null(hvg)){
            printColoredMessage(
                message = '- Performing PCA using all genes.',
                color = 'blue',
                verbose = verbose
                )
            norm.data <- irlba::prcomp_irlba(
                x = t(norm.data),
                n = bio.dims,
                scale. = FALSE,
                center = TRUE
                )
            ##### Using provided number of PCs ####
            if (!is.null(bio.dims)){
                printColoredMessage(
                    message = paste0(
                        '- Selection the first ',
                        bio.dims,
                        ' PCs as data input for MNN identification.'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- as.data.frame(t(norm.data$x[ , 1:bio.dims]))
                colnames(norm.data) <- sample.names
            }
            ### Using all PCs ####
            if (is.null(bio.dims)){
                printColoredMessage(
                    message = paste0(
                        '- Using all PCs as data input for MNN identification:'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- as.data.frame(t(norm.data$x))
                colnames(norm.data) <- sample.names
            }
        }
    }

    ## Selecting expression data ####
    if (data.input == 'expr'){
        printColoredMessage(
            message = '- Using gene expression matrix as data input fot MNN correction.',
            color = 'blue',
            verbose = verbose
            )
        ### Using all genes ####
        if (is.null(hvg)){
            printColoredMessage(
                message = '- Using all genes in gene expression matrix.',
                color = 'blue',
                verbose = verbose
            )
        }
        ### Using highly variable genes ####
        if (!is.null(hvg)){
            printColoredMessage(
                message = '- Using the provided highly variable genes in gene expression matrix.',
                color = 'blue',
                verbose = verbose
            )
        }
    }

    # Assessing and grouping the main unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the main unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    initial.variable <- se.obj[[main.uv.variable]]
    initial.variable2 <- se.obj[[main.uv.variable]]
    ## Continuous variable ####
    if (is.numeric(initial.variable)){
        se.obj[[main.uv.variable]] <- groupContinuousVariable(
            se.obj = se.obj,
            variable = main.uv.variable,
            nb.clusters = nb.clusters,
            clustering.method = clustering.method,
            perfix = '.',
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
        if (isFALSE(select.extreme.groups)) initial.se.obj <- se.obj
    }
    ## Categorical variable ####
    if (!is.numeric(initial.variable)){
        initial.se.obj <- se.obj
        if (length(unique(initial.variable)) == 1){
            stop('To create MNN, the "main.uv.variable" must have at least two groups/levels.')
        }
        if (length(unique(initial.variable)) > 1){
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
    initial.variable3 <- se.obj[[main.uv.variable]]

    # Creating PRPS data ####
    ## Considering other unwanted variables ####
    if (!is.null(other.uv.variables)){
        ## Assessing and grouping the other unwanted variable ####
        printColoredMessage(
            message = '- Assessing and grouping the other specified unwanted variable(s):',
            color = 'magenta',
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
        ## Finding subgroups that has enough samples for generating PS ####
        printColoredMessage(
            message = '- Finding subgroups with respect to other unwanted variable(s) that have enoug samples for PRPS:',
            color = 'magenta',
            verbose = verbose
            )
        covered.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                subgroups.size <- findRepeatingPatterns(
                    vec = all.uv.groups$main.uv[all.uv.groups$other.uv == x],
                    n.repeat = max(min.sample.for.ps, nb.mnn)
                )
            })
        names(covered.batches) <- unique(all.uv.groups$other.uv)
        printColoredMessage(
            message = paste0(
                '- There are ',
                length(covered.batches),
                ' subgroups with respect to other unwanted variables that have enough samples for generating PS.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '- Checking the distribution on main unwanted variable across subgroups with respect to other unwanted variables.',
            color = 'blue',
            verbose = verbose
        )
        covered.batches.table <- as.data.frame(
            table(all.uv.groups$main.uv, all.uv.groups$other.uv)
            )
        covered.batches.table$groups <- covered.batches.table$Freq >= max(min.sample.for.ps, nb.mnn)

        ## Plotting the distribution of th main unwanted variable across the subgroups of other unwanted variable ####
        if (isTRUE(plot.output)){
            printColoredMessage(
                message = '- Plotting distribution of samples across subgroups with respect to other unwanted variable(s):',
                color = 'magenta',
                verbose = verbose
            )
        }
        covered.batches.table$groups <- ifelse(covered.batches.table$groups == 'TRUE', 'selected', 'unselected')
        covered.batches.plot <- ggplot(covered.batches.table, aes(x = Var1, y = Var2, color = groups)) +
            geom_point(size = 6) +
            geom_text(aes(label = Freq , hjust = 0.5, vjust = 0.5), color = 'black') +
            xlab(main.uv.variable) +
            ylab('Homogeneous groups\n(other unwanted variables)') +
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

        ## Finding subgroups of other unwanted variable that have enough samples for PRPS ####
        printColoredMessage(
            message = '- Finding subgroups of other unwanted variable that have enough samples for generating PRPS:',
            color = 'magenta',
            verbose = verbose
            )
        selected.covered.batches <- lapply( covered.batches, length)
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
                    ' samples across all the sub-groups of the main unwanted variable: "',
                    main.uv.variable,
                    '" variable.'),
                color = 'blue',
                verbose = verbose
            )
            if (isFALSE(check.prps.connectedness)){
                printColoredMessage(
                    message = '-- We recommend specifiying the "check.prps.connectedness"',
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
                    max(min.sample.for.ps, nb.mnn),
                    ' samples across all sub-groups of the main unwanted variable: ',
                    main.uv.variable,
                    ' variable.'),
                color = 'blue',
                verbose = verbose
            )
        }
        ## Finding mutual nearest neighbor ####
        printColoredMessage(
            message = '-- Creating PRPS data:',
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
        prps.data <- assay(x = se.obj, i = assay.name)
        sample.annotation <- as.data.frame(x = colData(x = se.obj))
        colnames(sample.annotation) <- colnames(colData(x = se.obj))
        all.prps.data <- lapply(
            1:length(all.possible.batches),
            function(x){
                sub.sample.annotation <- sample.annotation[all.uv.groups$other.uv == names(all.possible.batches)[x] , ]
                pairs.batch <- all.possible.batches[[x]]
                sub.prps.data <- lapply(
                    1:ncol(pairs.batch),
                    function(y){
                        printColoredMessage(
                            message = paste0(
                                '- Creating PRPS data between the "',
                                pairs.batch[1 , y],
                                '" and "' ,
                                pairs.batch[2 , y],
                                '" subgroups within the "',
                                names(all.possible.batches)[x],
                                '" subgroup.'),
                            color = 'orange',
                            verbose = verbose
                        )
                        ## Obtaning sample annotation ####
                        sample.annot.a <- sub.sample.annotation[sub.sample.annotation[[main.uv.variable]] == pairs.batch[1 , y] , , drop = FALSE]
                        sample.annot.b <- sub.sample.annotation[sub.sample.annotation[[main.uv.variable]] == pairs.batch[2 , y] , , drop = FALSE]

                        ## Using highly variable genes ####
                        if (is.null(hvg)){
                            printColoredMessage(
                                message = '-- Highly variable are not specified, then using all genes.',
                                color = 'blue',
                                verbose = verbose
                                )
                            data.a <- norm.data[ , row.names(sample.annot.a)]
                            data.b <- norm.data[ , row.names(sample.annot.b)]
                        }
                        if (!is.null(hvg)){
                            printColoredMessage(
                                message = '-- Using the specified highly variable genes.',
                                color = 'blue',
                                verbose = verbose
                                )
                            data.a <- norm.data[hvg , row.names(sample.annot.a)]
                            data.b <- norm.data[hvg , row.names(sample.annot.b)]
                        }

                        ## Applying cosine normalization ####
                        if (isTRUE(apply.cosine.norm)){
                            printColoredMessage(
                                message = '- Applying cosine normalization on the data:',
                                color = 'blue',
                                verbose = verbose
                                )
                            data.a <- cosineNorm(x = data.a, mode = 'matrix')
                            data.b <- cosineNorm(x = data.b, mode = 'matrix')
                        }

                        ## Finding MNN for data b ####
                        printColoredMessage(
                            message = paste0(
                                '-- Finding ',
                                min.sample.for.ps,
                                ' nearest neighbours for each sample of the "',
                                pairs.batch[2 , y],
                                '" using the "' ,
                                pairs.batch[1 , y],
                                '" subgroup:'),
                            color = 'blue',
                            verbose = verbose
                            )
                        knn.data.a.b <- RANN::nn2(
                            data = t(data.a),
                            query = t(data.b),
                            k = min.sample.for.ps
                            )
                        ## Obtaining knn index ####
                        printColoredMessage(
                            message = '-- Obtaining the knn indexs.',
                            color = 'blue',
                            verbose = verbose
                            )
                        knn.data.a.b.index <- knn.data.a.b$nn.idx
                        colnames(knn.data.a.b.index) <- paste0('knn', seq(min.sample.for.ps))
                        row.names(knn.data.a.b.index) <- c(1:ncol(data.b))

                        ## Obtaining the distance ####
                        printColoredMessage(
                            message = '-- Obtaining the distances.',
                            color = 'blue',
                            verbose = verbose
                            )
                        knn.data.a.b.distance <- as.data.frame(knn.data.a.b$nn.dists)
                        colnames(knn.data.a.b.distance) <- paste0('knn', seq(min.sample.for.ps))
                        row.names(knn.data.a.b.distance) <- c(1:ncol(data.b))
                        knn.data.a.b.distance$aver.dist <- rowMeans(knn.data.a.b.distance)

                        ## Finding KNN for data a ####
                        printColoredMessage(
                            message = paste0(
                                '-- Finding ',
                                min.sample.for.ps,
                                ' nearest neighbours for each sample of the "',
                                pairs.batch[1 , y],
                                '" using "' ,
                                pairs.batch[2 , y],
                                '" subgroups:'),
                            color = 'blue',
                            verbose = verbose
                            )
                        knn.data.b.a <- RANN::nn2(
                            data = t(data.b),
                            query = t(data.a),
                            k = min.sample.for.ps
                            )
                        printColoredMessage(
                            message = '-- Obtaining the knn indexs.',
                            color = 'blue',
                            verbose = verbose
                            )
                        knn.data.b.a.index <- knn.data.b.a$nn.idx
                        colnames(knn.data.b.a.index) <- paste0('knn', seq(min.sample.for.ps))
                        row.names(knn.data.b.a.index) <- c(1:ncol(data.a))
                        ### distance
                        printColoredMessage(
                            message = '-- Obtaining the distances.',
                            color = 'blue',
                            verbose = verbose
                            )
                        knn.data.b.a.distance <- as.data.frame(knn.data.b.a$nn.dists)
                        colnames(knn.data.b.a.distance) <- paste0('knn', seq(min.sample.for.ps))
                        row.names(knn.data.b.a.distance) <- c(1:ncol(data.a))
                        knn.data.b.a.distance$aver.dist <- rowMeans(knn.data.b.a.distance)

                        ## Finding final MNN ####
                        printColoredMessage(
                            message = paste0(
                                '- Finding MNN between the "',
                                pairs.batch[1 , y],
                                '" using "' ,
                                pairs.batch[2 , y],
                                '" subgroups:'),
                            color = 'blue',
                            verbose = verbose
                            )
                        all.mnn <- BiocNeighbors::findMutualNN(
                            data1 = t(data.a),
                            data2 = t(data.b),
                            k1 = nb.mnn,
                            BNPARAM = mnn.nbparam,
                            BPPARAM = mnn.bpparam
                            )
                        printColoredMessage(
                            message = paste0(
                                '-- ',
                                length(all.mnn$first),
                                ' MNN are found.'),
                            color = 'blue',
                            verbose = verbose
                            )
                        ## Creating prps index ####
                        printColoredMessage(
                            message = paste0('- Creating PRPS indexs and score.'),
                            color = 'blue',
                            verbose = verbose
                            )
                        all.prps.index <- lapply(
                            1:length(all.mnn$first),
                            function(z){
                                sample.to.ave.b <- knn.data.b.a.index[all.mnn$first[z] , ]
                                sample.to.ave.a <- knn.data.a.b.index[all.mnn$second[z] , ]
                                prps.index <- data.frame(
                                    set.a = sample.to.ave.a,
                                    set.b = sample.to.ave.b
                                    )
                                samples.ids <- data.frame(
                                    set.a = row.names(sample.annot.a)[sample.to.ave.a],
                                    set.b = row.names(sample.annot.b)[sample.to.ave.b]
                                    )
                                avre.dist.a <- knn.data.b.a.distance[all.mnn$first[z] , 'aver.dist' ]
                                avre.dist.b <- knn.data.a.b.distance[all.mnn$second[z] , 'aver.dist' ]
                                list(
                                    prps.index = prps.index,
                                    aver.dist = c(avre.dist.a + avre.dist.b)/2,
                                    samples.ids = samples.ids
                                    )
                            })
                        names(all.prps.index) <- paste0('prps.set', 1:length(all.mnn$first))

                        ## Filtering the number of PRPS sets ####
                        if (isTRUE(filter.prps.sets)){
                            printColoredMessage(
                                message = '- Filtering the number of PRPS sets.',
                                color = 'blue',
                                verbose = verbose
                            )
                            if (length(all.prps.index) >=  max.prps.sets){
                                printColoredMessage(
                                    message = '- The number of the PRPS set is larger than the specified "max.prps.sets", filtering PRPS sets.',
                                    color = 'blue',
                                    verbose = verbose
                                    )
                                aver.dists <- sapply(
                                    1:length(all.prps.index),
                                    function(p) all.prps.index[[p]]$aver.dist
                                    )
                                names(aver.dists) <- 1:length(all.prps.index)
                                aver.dists <- aver.dists[order(aver.dists, decreasing = FALSE)]
                                select.prps <- names(aver.dists)[1:max.prps.sets]
                                all.prps.index <- all.prps.index[as.numeric(select.prps)]
                                printColoredMessage(
                                    message = paste0(
                                        '-- ',
                                        length(all.prps.index),
                                        ' PRPS set are kept.'),
                                    color = 'blue',
                                    verbose = verbose
                                    )
                            }
                        }
                        ## prps data
                        # Finding PRPS sets ####
                        printColoredMessage(
                            message = '- Creating the PRPS data:',
                            color = 'blue',
                            verbose = verbose
                        )
                        data.a <- prps.data[ , colnames(data.a)]
                        data.b <- prps.data[ , colnames(data.b)]
                        if (isTRUE(apply.log)){
                            printColoredMessage(
                                message = '- Applying log transformation on the data before creating PRPS. ',
                                color = 'blue',
                                verbose = verbose
                                )
                            data.a <- log2(data.a + pseudo.count)
                            data.b <- log2(data.b + pseudo.count)
                        } else {
                            printColoredMessage(
                                message = '- The data will be used without any transformation. ',
                                color = 'blue',
                                verbose = verbose
                                )
                        }
                        all.prps.data <- lapply(
                            1:length(all.prps.index),
                            function(a){
                                prps.set <- all.prps.index[[a]]
                                ps.a <- rowMeans(data.a[ , prps.set$prps.index$set.a, drop = FALSE])
                                ps.b <- rowMeans(data.b[ , prps.set$prps.index$set.b, drop = FALSE])
                                prps <- data.frame(ps.a, ps.b)
                                prps
                            })
                        all.prps.data <- do.call(cbind, all.prps.data)
                        colnames(all.prps.data) <- paste0(
                            pairs.batch[1 , y],
                            '||',
                            pairs.batch[2 , y],
                            '||',
                            rep(1:c(ncol(all.prps.data)/2), each = 2),
                            '||',
                            names(all.possible.batches[x])
                            )
                        all.prps.sample.annot <- lapply(
                            1:length(all.prps.index),
                            function(f){
                                prps.set <- all.prps.index[[f]]$samples.ids
                            })
                        names(all.prps.sample.annot) <- rep(
                            paste0(pairs.batch[1 , y], '_' , pairs.batch[2 , y]),
                            length(all.prps.sample.annot)
                        )
                        list(
                            all.prps.data = all.prps.data,
                            all.prps.sample.annot = all.prps.sample.annot
                        )
                    })
                sub.prps.epxr.data <- lapply(
                    1:length(sub.prps.data),
                    function(t) sub.prps.data[[t]]$all.prps.data
                )
                sub.prps.epxr.data <- do.call(
                    cbind,
                    sub.prps.epxr.data
                )
                sub.prps.info.data <- lapply(
                    1:length(sub.prps.data),
                    function(t) sub.prps.data[[t]]$all.prps.sample.annot
                )
                return(list(
                    sub.prps.epxr.data = sub.prps.epxr.data,
                    sub.prps.info.data = sub.prps.info.data,
                    groups = names(all.possible.batches[x]))
                )
            })
        ### Putting all the PRPS data together ####
        printColoredMessage(
            message = '-- Putting all the PRPS data togather: ',
            color = 'blue',
            verbose = verbose
            )
        all.prps.expr.data <- lapply(
            1:length(all.prps.data),
            function(x){
                all.prps.data[[x]]$sub.prps.epxr.data
            })
        all.prps.expr.data <- do.call(cbind, all.prps.expr.data)
        printColoredMessage(
            message = paste0(
                '-- The totall number of ',
                ncol(all.prps.expr.data),
                ' pseudo-samples are created.'),
            color = 'blue',
            verbose = verbose
            )
        if (sum(is.na(all.prps.expr.data)) > 0){
            stop('The PRPS data contains NA, this is not supported.')
        }
        se.obj[[main.uv.variable]] <- initial.variable
        sample.annotation <- as.data.frame(colData(x = se.obj))
        colnames(sample.annotation) <- colnames(colData(x = se.obj))

        ### Plotting PRPS map ####
        printColoredMessage(
            message = '-- Plotting all the PRPS data togather: ',
            color = 'blue',
            verbose = verbose
            )
        all.prps.plot <- lapply(
            1:length(all.prps.data),
            function(x){
                prps.index <- all.prps.data[[x]]$sub.prps.info.data
                expr.data <- lapply(
                    1:length(prps.index),
                    function(y){
                        m <- lapply(
                            1:length(prps.index[[y]]),
                            function(z){
                                data.frame(
                                    group1 = sample.annotation[prps.index[[y]][[z]]$set.a , ][[ main.uv.variable]],
                                    group2 =  sample.annotation[prps.index[[y]][[z]]$set.b , ][[ main.uv.variable]],
                                    set.name = rep(
                                        unique(names(prps.index[[y]])),
                                        length(prps.index[[y]][[z]]$set.b)
                                    ),
                                    group.name = rep(
                                        all.prps.data[[x]]$groups,
                                        length(prps.index[[y]][[z]]$set.b)
                                    )
                                )
                            })
                        m <- do.call(rbind, m )
                    })
                expr.data <- do.call(rbind, expr.data)
                expr.data$prps.set <- rep(1:c(nrow(expr.data)/3), each = min.sample.for.ps)
                expr.data
            })
        all.prps.plot <- do.call(rbind, all.prps.plot) %>%
            pivot_longer(-c(set.name, group.name, prps.set), values_to = 'var', names_to = 'rep')
        all.prps.plot$new.g <- paste(
            all.prps.plot$group.name,
            all.prps.plot$set.name,
            all.prps.plot$prps.set,
            sep = '||'
            )
        all.prps.plot$new.group <- 'PRPS sets'
        all.uv <- data.frame(
            set.name = rep(main.uv.variable, ncol(se.obj)),
            group.name = rep(main.uv.variable, ncol(se.obj)),
            prps.set = c(1:ncol(se.obj)),
            rep = initial.variable3,
            var = se.obj[[main.uv.variable]],
            new.g = rep('UV', ncol(se.obj)),
            new.group = rep('UV', ncol(se.obj))
            )
        new <- rbind(all.uv , all.prps.plot)
        new$new.group <- factor(x = new$new.group, levels = c('UV', 'PRPS sets'))
        if (is.numeric(se.obj[[main.uv.variable]])){
            prps.map.plot <- ggplot(data = new, aes(x = var , y = new.g , color = rep)) +
                geom_boxplot() +
                geom_point() +
                scale_color_manual(
                    values = c('darkgreen', 'tomato', c(viridis(nb.clusters))),
                    name = 'Groups') +
                facet_grid(new.group~., scales = 'free', space = 'free') +
                scale_x_discrete(expand = c(0, 0.5)) +
                xlab(main.uv.variable) +
                ylab('Homogeneous groups') +
                xlim(c(
                    min(se.obj[[main.uv.variable]]),
                    max(se.obj[[main.uv.variable]])
                )) +
                geom_vline(xintercept = c(
                    min(se.obj[[main.uv.variable]]),
                    max(se.obj[[main.uv.variable]])), color = 'gray70') +
                theme_bw() +
                theme(
                    legend.key = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text.y = element_text(size = 14),
                    axis.text.x = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    strip.text.y = element_text(size = 15)
                    )
        }
        if (!is.numeric(se.obj[[main.uv.variable]])){
            prps.map.plot <- ggplot(data = all.prps.plot, aes(x = var , y = new.g  , color = rep)) +
                geom_point(size = 3) +
                scale_color_manual(values = c('darkgreen', 'tomato', 'navy'), name = 'Groups') +
                # facet_grid(new.group~., scales = 'free', space = 'free') +
                scale_x_discrete(expand = c(0, 0.5)) +
                xlab(main.uv.variable) +
                ylab('Homogeneous groups') +
                theme_bw() +
                theme(
                    legend.key = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    axis.text.y = element_text(size = 14),
                    axis.text.x = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 18),
                    strip.text.y = element_text(size = 15)
                )
        }
        if (isTRUE(plot.output)) print(prps.map.plot)

        ### Creating sample annotation of PRPS data ####
        group <- 1
        all.prps.sample.annot <- lapply(
            1:length(all.prps.data),
            function(x){
                sample.annot <- all.prps.data[[x]]$sub.prps.info.data
                sample.annot <- sample.annot[[1]]
                sample.annot <- do.call(rbind, sample.annot)
                sample.annot$group <- all.prps.data[[x]]$groups
                sample.annot
            })
        all.prps.sample.annot <- do.call(rbind, all.prps.sample.annot)
    }

    # Considering only main unwanted variation ####
    if (is.null(other.uv.variables)){
        ## Checking the sample size of each group in the variable ####
        subgroups.size <- findRepeatingPatterns(
            vec = se.obj[[main.uv.variable]],
            n.repeat = max(min.sample.for.ps, nb.mnn)
            )
        if (min.batches.to.cover == 'all') {
            if (length(subgroups.size) != length(unique(se.obj[[main.uv.variable]])) ){
                stop(paste0(
                    'Some sub-groups of the variable "',
                    main.uv.variable,
                    '" have less than ',
                    max(min.sample.for.ps, nb.mnn),
                    ' samples. Then, MNN cannot be created across all batches.')
                    )
            }
        }
        if (is.numeric(min.batches.to.cover)) {
            if (length(subgroups.size) >= min.batches.to.cover){
                printColoredMessage(
                    message = paste0(
                        '- At least ',
                        min.batches.to.cover,
                        ' sub-groups of the variable ',
                        main.uv.variable,
                        ' have at least ',
                        max(min.sample.for.ps, nb.mnn),
                        ' samples.'),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                stop(paste0(
                    'Some sub-groups of the variable "',
                    main.uv.variable,
                    '" have less than ',
                    max(min.sample.for.ps, nb.mnn),
                    ' samples. Then, MNN cannot be created across all batches.')
                )
            }
        }
        ### Finding mutual nearest neighbor ####
        printColoredMessage(
            message = '-- Creating PRPS data across all pairs of the subgroups:',
            color = 'magenta',
            verbose = verbose
            )
        pairs.batch <- combn(x = subgroups.size, m = 2)
        prps.data <- assay(x = se.obj, i = assay.name)
        sample.annotation <- as.data.frame(colData(x = se.obj))
        colnames(sample.annotation) <- colnames(colData(x = se.obj))
        all.prps.data <- lapply(
            1:ncol(pairs.batch),
            function(y){
                printColoredMessage(
                    message = paste0(
                        '- Creating PRPS data between the "',
                        pairs.batch[1 , y],
                        '" and "' ,
                        pairs.batch[2 , y],
                        '" subgroups."'),
                    color = 'blue',
                    verbose = verbose
                    )
                ## Obtaining sample annotation ####
                sample.annot.a <- sample.annotation[sample.annotation[[main.uv.variable]] == pairs.batch[1 , y] , , drop = FALSE]
                sample.annot.b <- sample.annotation[sample.annotation[[main.uv.variable]] == pairs.batch[2 , y] , , drop = FALSE]

                ## Using highly variable genes ####
                if (data.input == 'pcs'){
                    data.a <- norm.data[ , row.names(sample.annot.a), drop = FALSE]
                    data.b <- norm.data[ , row.names(sample.annot.b), drop = FALSE]
                }
                if (!is.null(hvg) & data.input == 'expr'){
                    data.a <- norm.data[hvg , row.names(sample.annot.a), drop = FALSE]
                    data.b <- norm.data[hvg , row.names(sample.annot.b), drop = FALSE]
                }
                if (is.null(hvg) & data.input == 'expr'){
                    data.a <- norm.data[ , row.names(sample.annot.a), drop = FALSE]
                    data.b <- norm.data[ , row.names(sample.annot.b), drop = FALSE]
                }

                ## Applying cosine normalization ####
                if (isTRUE(apply.cosine.norm)){
                    printColoredMessage(
                        message = '- Applying cosine normalization on the data.',
                        color = 'blue',
                        verbose = verbose
                        )
                    data.a <- cosineNorm(x = data.a, mode = 'matrix')
                    data.b <- cosineNorm(x = data.b, mode = 'matrix')
                }

                ## Finding MNN for data b ####
                printColoredMessage(
                    message = paste0(
                        '- Finding ',
                        min.sample.for.ps,
                        ' nearest neighbours for each sample of the "',
                        pairs.batch[2 , y],
                        '" using the "' ,
                        pairs.batch[1 , y],
                        '" subgroup:'),
                    color = 'blue',
                    verbose = verbose
                    )
                knn.data.a.b <- RANN::nn2(
                    data = t(data.a),
                    query = t(data.b),
                    k = min.sample.for.ps
                    )
                ## Obtaining knn index ####
                printColoredMessage(
                    message = ' - Obtaining the knn indexs.',
                    color = 'blue',
                    verbose = verbose
                    )
                knn.data.a.b.index <- knn.data.a.b$nn.idx
                colnames(knn.data.a.b.index) <- paste0('knn', seq(min.sample.for.ps))
                row.names(knn.data.a.b.index) <- c(1:ncol(data.b))

                ## Obtaining the distance ####
                printColoredMessage(
                    message = '- Obtaining the distances.',
                    color = 'blue',
                    verbose = verbose
                    )
                knn.data.a.b.distance <- as.data.frame(knn.data.a.b$nn.dists)
                colnames(knn.data.a.b.distance) <- paste0('knn', seq(min.sample.for.ps))
                row.names(knn.data.a.b.distance) <- c(1:ncol(data.b))
                knn.data.a.b.distance$aver.dist <- rowMeans(knn.data.a.b.distance)

                ## Finding knn for data a ####
                printColoredMessage(
                    message = paste0(
                        '- Finding ',
                        min.sample.for.ps,
                        ' nearest neighbours for each sample of the "',
                        pairs.batch[1 , y],
                        '" using "' ,
                        pairs.batch[2 , y],
                        '" subgroups:'),
                    color = 'blue',
                    verbose = verbose
                    )
                knn.data.b.a <- RANN::nn2(
                    data = t(data.b),
                    query = t(data.a),
                    k = min.sample.for.ps
                    )
                printColoredMessage(
                    message = '- Obtaining the knn indexs.',
                    color = 'blue',
                    verbose = verbose
                    )
                knn.data.b.a.index <- knn.data.b.a$nn.idx
                colnames(knn.data.b.a.index) <- paste0('knn', seq(min.sample.for.ps))
                row.names(knn.data.b.a.index) <- c(1:ncol(data.a))
                ### distance
                printColoredMessage(
                    message = '-- Obtaining the distances.',
                    color = 'blue',
                    verbose = verbose
                    )
                knn.data.b.a.distance <- as.data.frame(knn.data.b.a$nn.dists)
                colnames(knn.data.b.a.distance) <- paste0('knn', seq(min.sample.for.ps))
                row.names(knn.data.b.a.distance) <- c(1:ncol(data.a))
                knn.data.b.a.distance$aver.dist <- rowMeans(knn.data.b.a.distance)

                ## Finding final MNN ####
                printColoredMessage(
                    message = paste0(
                        '- Finding MNN between the "',
                        pairs.batch[1 , y],
                        '" using "' ,
                        pairs.batch[2 , y],
                        '" subgroups:'),
                    color = 'blue',
                    verbose = verbose
                    )
                all.mnn <- BiocNeighbors::findMutualNN(
                    data1 = t(data.a),
                    data2 = t(data.b),
                    k1 = nb.mnn,
                    BNPARAM = mnn.nbparam,
                    BPPARAM = mnn.bpparam
                    )
                printColoredMessage(
                    message = paste0(
                        '- Totall number of ',
                        length(all.mnn$first),
                        ' MNN are found.'),
                    color = 'blue',
                    verbose = verbose
                    )
                ## Creating prps index ####
                printColoredMessage(
                    message = paste0('- creating PRPS indexs and score.'),
                    color = 'blue',
                    verbose = verbose
                    )
                all.prps.index <- lapply(
                    1:length(all.mnn$first),
                    function(z){
                        sample.to.ave.b <- knn.data.b.a.index[all.mnn$first[z] , ]
                        sample.to.ave.a <- knn.data.a.b.index[all.mnn$second[z] , ]
                        prps.index <- data.frame(
                            set.a = sample.to.ave.a,
                            set.b = sample.to.ave.b
                            )
                        samples.ids <- data.frame(
                            set.a = row.names(sample.annot.a)[unname(sample.to.ave.a)],
                            set.b = row.names(sample.annot.b)[unname(sample.to.ave.b)]
                            )
                        avre.dist.a <- knn.data.b.a.distance[all.mnn$first[z] , 'aver.dist' ]
                        avre.dist.b <- knn.data.a.b.distance[all.mnn$second[z] , 'aver.dist' ]
                        list(
                            prps.index = prps.index,
                            aver.dist = c(avre.dist.a + avre.dist.b)/2,
                            samples.ids = samples.ids
                            )
                    })
                names(all.prps.index) <- paste0('prps.set', 1:length(all.mnn$first))
                keep.prps.unique <- sapply(
                    1:length(all.prps.index),
                    function(x){
                        target.set.a <- c(
                            all.prps.index[[x]]$samples.ids$set.a,
                            all.prps.index[[x]]$samples.ids$set.b
                            )
                        m <- sapply(
                            1:length(all.prps.index),
                            function(y){
                                target.set.b <- c(
                                    all.prps.index[[y]]$samples.ids$set.a,
                                    all.prps.index[[y]]$samples.ids$set.b
                                )
                                length(intersect(x = target.set.b, target.set.a))
                            })
                        if (sum(m == 6)>1){
                            rep.set <- which(m == 6)[!which(m == 6) %in% x]
                            dist.diff <- sapply(
                                c(x, rep.set),
                                function(s){
                                    all.prps.index[[s]]$aver.dist
                                })
                            c(x, rep.set)[which(dist.diff == min(dist.diff))]
                        } else x
                    })
                keep.prps.unique <- unique(keep.prps.unique)
                all.prps.index <- lapply(
                    keep.prps.unique,
                    function(x){
                        all.prps.index[[x]]
                    })
                names(all.prps.index) <- paste0('prps.set', keep.prps.unique)

                ## Filtering the number of PRPS sets ####
                if (isTRUE(filter.prps.sets)){
                    printColoredMessage(
                        message = '- Filtering the number of PRPS sets.',
                        color = 'blue',
                        verbose = verbose
                        )
                    if (length(all.prps.index) >=  max.prps.sets){
                        printColoredMessage(
                            message = '- The number of the PRPS set is larger than the specified "max.prps.sets", filtering PRPS sets.',
                            color = 'blue',
                            verbose = verbose
                            )
                        aver.dists <- sapply(
                            1:length(all.prps.index),
                            function(p) all.prps.index[[p]]$aver.dist
                            )
                        names(aver.dists) <- 1:length(all.prps.index)
                        aver.dists <- aver.dists[order(aver.dists, decreasing = FALSE)]
                        select.prps <- names(aver.dists)[1:max.prps.sets]
                        all.prps.index <- all.prps.index[as.numeric(select.prps)]
                        printColoredMessage(
                            message = paste0(
                                '- ',
                                length(all.prps.index),
                                ' PRPS set are kept.'),
                            color = 'blue',
                            verbose = verbose
                        )
                    }
                }
                ## prps data
                # Finding PRPS sets ####
                printColoredMessage(
                    message = '- Creating the PRPS data:',
                    color = 'blue',
                    verbose = verbose
                    )
                data.a <- prps.data[ , colnames(data.a)]
                data.b <- prps.data[ , colnames(data.b)]
                if (isTRUE(apply.log.for.prps)){
                    printColoredMessage(
                        message = paste0(
                            '- Applying log transformation on the "',
                            assay.name,
                            '" data before averaging samples to create PS.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    data.a <- log2(data.a + pseudo.count)
                    data.b <- log2(data.b + pseudo.count)
                } else {
                    printColoredMessage(
                        message = paste0(
                        '- The "',
                        assay.name,
                        '" data will be used without any transformation before averaging samples to create PS.'),
                        color = 'blue',
                        verbose = verbose
                    )
                }
                all.prps.data <- lapply(
                    1:length(all.prps.index),
                    function(a){
                        prps.set <- all.prps.index[[a]]
                        ps.a <- rowMeans(data.a[ , prps.set$prps.index$set.a, drop = FALSE])
                        ps.b <- rowMeans(data.b[ , prps.set$prps.index$set.b, drop = FALSE])
                        prps <- data.frame(ps.a, ps.b)
                        prps
                    })
                all.prps.data <- do.call(cbind, all.prps.data)
                colnames(all.prps.data) <- paste0(
                    pairs.batch[1 , y],
                    '||',
                    pairs.batch[2 , y],
                    '||',
                    rep(1:c(ncol(all.prps.data)/2), each = 2)
                    )
                all.prps.sample.annot <- lapply(
                    1:length(all.prps.index),
                    function(f){
                        prps.set <- all.prps.index[[f]]$samples.ids
                    })
                names(all.prps.sample.annot) <- rep(
                    paste0(pairs.batch[1 , y], '_' , pairs.batch[2 , y]),
                    length(all.prps.sample.annot)
                    )
                list(
                    all.prps.data = all.prps.data,
                    all.prps.sample.annot = all.prps.sample.annot
                    )
            })
        ### Putting all the PRPS data together ####
        printColoredMessage(
            message = '- Putting all the PRPS data togather: ',
            color = 'blue',
            verbose = verbose
            )
        all.prps.expr.data <- lapply(
            1:length(all.prps.data),
            function(x){
                all.prps.data[[x]]$all.prps.data
            })
        all.prps.expr.data <- do.call(cbind, all.prps.expr.data)
        se.obj[[main.uv.variable]] <- initial.variable
        sample.annotation <- as.data.frame(colData(x = se.obj))
        colnames(sample.annotation) <- colnames(colData(x = se.obj))

        ### Plotting PRPS map ####
        if (isTRUE(create.prps.map)){
            printColoredMessage(
                message = '- Plotting the PRPS map:',
                color = 'magenta',
                verbose = verbose
            )
            all.prps.sample.annot <- lapply(
                1:length(all.prps.data),
                function(x){
                    sample.annot <- all.prps.data[[x]]$all.prps.sample.annot
                    possible.max.prps.sets <- length(sample.annot)
                    sample.annot <- do.call(rbind, sample.annot)
                    sample.annot$group <- paste0(
                        unique(names(all.prps.data[[x]]$all.prps.sample.annot)),
                        rep(paste0('-PRPS', 1:possible.max.prps.sets), each = min.sample.for.ps)
                    )
                    sample.annot <- sample.annot %>%
                        pivot_longer(-c(group), values_to = 'var', names_to = 'rep')
                    sample.annot$var.obj <- unlist(lapply(
                        sample.annot$var,
                        function(g){
                            se.obj[[main.uv.variable]][colnames(se.obj) == g]
                        }))
                    sample.annot
                })
            all.prps.sample.annot <- do.call(rbind, all.prps.sample.annot)
            all.prps.sample.annot$rep <- ifelse(
                all.prps.sample.annot$rep == 'set.a', 'Set1', 'Set2'
            )
            all.prps.sample.annot$new.group <- 'PRPS sets'
            all.uv <- data.frame(
                group = rep(main.uv.variable, ncol(initial.se.obj)),
                rep = initial.se.obj[[main.uv.variable]],
                var = colnames(initial.se.obj),
                var.obj = initial.variable2,
                new.group = rep('UV', ncol(initial.se.obj))
            )
            all.uv <- rbind(all.uv , all.prps.sample.annot)
            all.uv$new.group <- factor(
                x = all.uv$new.group,
                levels = c('UV', 'PRPS sets')
            )
            all.uv$rep <- factor(
                x = all.uv$rep,
                levels = c('Set1', 'Set2', unique(as.character(initial.se.obj[[main.uv.variable]])))
            )
            if (is.numeric(initial.variable)){
                prps.sets <- all.uv$group[which(!all.uv$group %in% main.uv.variable)]
                for(i in prps.sets){
                    temp.data <- droplevels(all.uv[all.uv$group == i , , drop = FALSE])
                    g1 <- mean(temp.data$var.obj[temp.data$rep == 'Set1'])
                    g2 <- mean(temp.data$var.obj[temp.data$rep == 'Set2'])
                    if (g1 > g2 ){
                        all.uv$rep[all.uv$group == i] <- ifelse (temp.data$rep == 'Set1', 'Set2', 'Set1')
                    }
                }
            }
            if (is.numeric(se.obj[[main.uv.variable]])){
                prps.map.plot <- ggplot(data = all.uv, aes(x = var.obj , y = group , color = rep)) +
                    geom_boxplot() +
                    geom_point() +
                    scale_color_manual(
                        values = c('orange', 'navy', viridis(n = nb.clusters)),
                        name = 'Groups') +
                    facet_grid(new.group~., scales = 'free', space = 'free') +
                    #scale_x_discrete(expand = c(0, 0.5)) +
                    xlab(main.uv.variable) +
                    ylab('') +
                    xlim(c(
                        min(se.obj[[main.uv.variable]]),
                        max(se.obj[[main.uv.variable]])
                    )) +
                    geom_vline(xintercept = c(
                        min(se.obj[[main.uv.variable]]),
                        max(se.obj[[main.uv.variable]])), color = 'gray70') +
                    theme_bw() +
                    theme(
                        legend.text = element_text(size = 14),
                        legend.title = element_text(size = 18),
                        axis.line = element_line(colour = 'black', linewidth = .85),
                        axis.title.x = element_text(size = 16),
                        strip.text = element_text(size = 20),
                        axis.title.y = element_text(size = 16),
                        axis.text.x = element_text(size = 12, angle = 35, hjust = 1, vjust = 1),
                        axis.text.y = element_text(size = 12),
                        legend.position = 'right') +
                    guides(color = guide_legend(title = "Groups"))
            }
            if (!is.numeric(se.obj[[main.uv.variable]])){
                prps.map.plot <- ggplot(data = all.prps.sample.annot, aes(x = var.obj , y = group , color = rep)) +
                    geom_point(size = 3) +
                    scale_color_manual(
                        values = c('darkgreen', 'tomato', viridis(n = length(unique(initial.se.obj[[main.uv.variable]])))),
                        name = 'Groups') +
                    facet_grid(new.group~., scales = 'free', space = 'free') +
                    scale_x_discrete(expand = c(0, 0.5)) +
                    xlab(main.uv.variable) +
                    ylab('') +
                    theme_bw() +
                    theme(
                        legend.key = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 16),
                        axis.title.y = element_text(size = 16),
                        axis.text.y = element_text(size = 14),
                        axis.text.x = element_text(size = 14),
                        legend.text = element_text(size = 14),
                        legend.title = element_text(size = 18),
                        strip.text.y = element_text(size = 15)
                    )
            }
            if (isTRUE(plot.output)) print(prps.map.plot)
        }
        ### Obtaining sample ids of RPPS data ####
        all.prps.sample.annot <- lapply(
            1:length(all.prps.data),
            function(x){
                sample.annot <- all.prps.data[[x]]$all.prps.sample.annot
                sample.annot <- do.call(rbind, sample.annot)
                sample.annot$group <- all.prps.data[[x]]$groups
                sample.annot
            })
        all.prps.sample.annot <- do.call(rbind, all.prps.sample.annot)
    }

    # Saving the results ####
    ## Selecting names for the PRPS data ####
    if (is.null(prps.sets.name))
        prps.sets.name <- paste0(
            main.uv.variable,
            '|',
            'mnn',
            '|',
            assay.name
            )
    if (is.null(prps.group.name))
        prps.group.name <- main.uv.variable
    printColoredMessage(
        message = '-- Saving the PRPS data',
        color = 'magenta',
        verbose = verbose
        )
    se.obj <- initial.se.obj
    # se.obj <- initial.se.obj.a
    se.obj[[main.uv.variable]] <- initial.variable2
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
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.data']] <- all.prps.expr.data

        # plot
        if (!'prps.map.plot' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.map.plot']] <- list()
        }
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group.name]][[prps.sets.name]][['prps.map.plot']] <- prps.map.plot
        printColoredMessage(message = '------------The createPrPsByMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## output the PRPS data as matrix ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The createPrPsByMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(list(prps.data = all.prps.expr.data, prps.map.plot = prps.map.plot))
    }
}
