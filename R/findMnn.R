#' Find mutual nearest neighbors in RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function finds mutual nearest neighbors (MNN) between all pairs of specified batches in RNA-seq data. The mutual
#' nearest neighbors will be used to find and create pseudo samples (Ps) and eventually pseudo-replicates (PS )for RUV-III
#' with PRPS normalization. This function is used in the `createPRPSByKnnMnn()` function.
#'
#' @param se.obj A summarized experiment object.
#' @param assay.name Character. A character that indicates the name of the data (assay) in the SummarizedExperiment
#' object.
#' @param uv.variable Character. A character indicating the name of a column in the sample annotation in the
#' SummarizedExperiment object. The `uv.variable` can be either a categorical or continuous source of unwanted variation.
#' If `uv.variable` is a continuous variable, this will be divided into `nb.clusters` groups using the `clustering.method`
#' method.
#' @param nb.mnn Numeric. A numeric value specifying the maximum number of mutual nearest neighbors to compute across
#' batches or subgroups. The default is set to 1.
#' @param clustering.method Character. A character that indicates the choice of clustering method for grouping the
#' `uv.variable`, if a continuous variable is provided. Options include `kmeans`, `cut`, and `quantile`. The default is
#' set to `kmeans`.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the `uv.variable` is a
#' continuous variable. The default is set to 3.
#' @param data.input Character. A character that indicates which data should be used as input for finding the mutual
#' nearest neighbors. Options include: `expr` and `pcs`. If `pcs` is selected, the first `nb.pcs` of PCs of the data will
#' be computed and used as input. If `expr` is selected, the expression data will be used as input. The default is set
#' to `expr`.
#' @param nb.pcs Numeric. A numeric value that indicates the number of PCs to be calculated and then used as data input
#' for finding the mutual nearest neighbors. The default is set to 2. The `nb.pcs` must be set when `data.input = pcs`.
#' @param center Logical. Indicates whether to scale the data or not before calculating PCs. If center is set to `TRUE`,
#' then centering is done by subtracting the column means of the assay from their corresponding columns. The default is
#' se to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data or not before calculating PCs. If scale is set to `TRUE`,
#' then scaling is done by dividing the (centered) columns of the assays by their standard deviations if center is `TRUE`,
#' and the root mean square otherwise. The default is set to `FALSE`.
#' @param svd.bsparam Character. A BiocParallelParam object specifying how palatalization should be performed to compute
#' PCs. The default is set to bsparam(). We refer to the `runSVD()` function from the **BiocSingular** R package for
#' further details.
#' @param normalization Character. A character that indicates which normalization methods should be applied before finding
#' the MNN. The options are: `CPM`, `TMM`, `VST`, `full`, `upper` and `median`. The default is set to `CPM`. Refer to the
#' `applyOtherNormalization()` for more details.
#' @param apply.cosine.norm Logical. Indicates whether to apply cosine normalization on the data before finding MNN. The
#'  default is set to `TRUE`.
#' @param regress.out.variables Character. A character or a vector of characters indicating the column name(s) in sample
#' annotation in the SummarizedExperiment object. These variables will be regressed out from the data before finding MNN.
#' The default is set to `NULL`, indicating the regression will not be applied.
#' @param hvg Vector. A logical vector or a vector of the names (feature ids) of the highly variable genes. These genes
#' will be used to prepare the input data, either `pcs` or `expre`, for MNN analysis. The default is set to `NULL`, this
#' means all genes will be used.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is set to
#' `TRUE`.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid `-Inf` for measurements that are equal to 0. The default is set to 1.
#' @param mnn.bpparam Character. A BiocParallelParam object specifying how palatalization should be performed to find MNN.
#' The default is set to `SerialParam()`. We refer to the `findMutualNN()` function from the `BiocNeighbors` R package.
#' @param mnn.nbparam Character. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' The default is set to `KmknnParam()`. We refer to the 'findMutualNN()' function from the 'BiocNeighbors' R package.
#' @param check.se.obj Logical. Indicates whether to check the SummarizedExperiment object or not. The default is set
#' to `TRUE`. See the `checkSeObj()` function for more details.
#' @param remove.na Character. Specifies whether to remove NA or missing values from the assays (data) or not. The options
#' are `assays`, `sample.annotation`, `both` and `none`. The default is set to `both`, so all the NA or missing values
#' from the provided data and variables will be removed before any analysis See the `checkSeObj()` function for more details.
#' @param plot.output Logical. If `TRUE`, the function plots the distribution of MNN across the batches or subgroups. The
#' default is set to `TRUE`.
#' @param mnn.group.name Character. A character specifying the name of the MNN group to which the current MNN belong.
#' If set to `NULL`, the function will automatically assign a name using the `uv.variable`.
#' @param mnn.sets.name Character. A character specifying the name of the output file (MNN data) to be saved in the
#' metadata of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on:
#' `paste0(assay.name, '|', length(unique(se.obj[[uv.variable]])) ,'groups|', nb.mnn, 'mnn')`.
#' @param save.se.obj Logical. Indicates whether to save the MNN results in the metadata of the SummarizedExperiment object
#' or to output the result as a list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, shows the messages of different steps of the function.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocNeighbors findMutualNN KmknnParam
#' @importFrom BiocParallel SerialParam
#' @importFrom utils setTxtProgressBar
#'
#' @export

findMnn <- function(
        se.obj,
        assay.name,
        uv.variable,
        nb.mnn = 1,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        data.input = 'pcs',
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        normalization = 'CPM',
        apply.cosine.norm = FALSE,
        regress.out.variables = NULL,
        hvg = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        check.se.obj = TRUE,
        remove.na = 'both',
        plot.output = TRUE,
        mnn.group.name = NULL,
        mnn.sets.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findMnn function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the function inputs #####
    if (is.list(assay.name)) {
        stop('The "assay.name" cannot be a list.')
    }
    if (length(assay.name) > 1) {
        stop('The "assay.name" must be the name of signle assay in the SummarizedExperiment object.')
    }
    if (is.null(uv.variable)) {
        stop('The "uv.variable" variable cannot be empty.')
    }
    if (length(uv.variable) > 1) {
        stop('The "uv.variable" must contain the name of signle variable in the SummarizedExperiment object.')
    }
    if (!uv.variable %in% colnames(colData(se.obj))) {
        stop('The "uv.variable" variable cannot be found in the SummarizedExperiment object.')
    }
    if (isFALSE(check.se.obj)){
        if (!uv.variable %in% colnames(colData(se.obj))) {
            stop('The "uv.variable" cannot be found in the SummarizedExperiment object.')
        }
        if (sum(regress.out.variables %in% colnames(colData(se.obj)))!=length(regress.out.variables) ) {
            stop('All or some of the "regress.out.variables" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (!data.input %in% c('expr', 'pcs')) {
        stop('The "data.input" must be one of the "expr" or "pcs".')
    }
    if (data.input == 'pcs' & is.null(nb.pcs)) {
        stop('The valuse of "nb.pcs" must be sepcified when the data.input = pcs.')
    }
    if (!is.null(nb.pcs)){
        if (!is.numeric(nb.pcs)){
            stop('The "nb.pcs" must be a postive numeric value.')
        }
        if (nb.pcs <= 0 & nb.pcs != as.integer(nb.pcs)) {
            stop('The "nb.pcs" must be a postive whole numeric value')
        }
    }
    if (!is.numeric(nb.mnn)){
        stop('The "nb.mnn" must be a postive numeric value.')
    }
    if (nb.mnn <= 0 & nb.mnn != as.integer(nb.mnn)) {
        stop('The "nb.mnn" must be a postive whole numeric value')
    }
    if (!is.null(hvg)) {
        if(is.logical(hvg)){
            if(sum(hvg) == 0){
                stop('The "hvg" does not contain any genes.')
            } else if (sum(hvg) <= 3){
                stop('The number of "hvg" must be at least 3')
            } else if (length(hvg) > nrow(se.obj)){
                stop('The length of the "hvg" cannot be larger than the row numbers of the SammarizedExperiment object.')
            }
        }
        if (is.character(hvg)){
            if (sum(hvg %in% row.names(se.obj)) != length(hvg)){
                stop('All the hvg genes are not found in the SummarizedExperiment object.')
            } else if (length(hvg) > nrow(se.obj)){
                stop('The length of the "hvg" cannot be larger than the row numbers of the SammarizedExperiment object.')
            } else if (length(hvg) <= 3){
                stop('The number of "hvg" must be at least 3')
            }
            hvg <- row.names(se.obj) %in% hvg
        }
    } else if (is.null(hvg)){
        hvg <- rep(TRUE, nrow(se.obj))
    }
    if(!is.logical(apply.log)){
        stop('The "apply.log" must be logical.')
    }
    if (!is.null(pseudo.count)){
        if (!is.numeric(pseudo.count)){
            stop('The "pseudo.count" must be a numeric value.')
        }
        if (pseudo.count < 0){
            stop('The "pseudo.count" must be a postive numeric value.')
        }
    }
    if (!is.logical(check.se.obj)){
        stop('The "apply.log" must be logical.')
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical.')
    }
    if (!is.null(mnn.sets.name)){
        if (!is.character(mnn.sets.name)){
            stop('The "mnn.sets.name" must be a character.')
        }
    }
    if (!is.null(mnn.group.name)){
        if (!is.character(mnn.group.name)){
            stop('The "mnn.group.name" must be logical.')
        }
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical.')
    }

    # Assessing the SummarizedExperiment object #####
    if (isTRUE(check.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(uv.variable, regress.out.variables),
            remove.na = remove.na,
            verbose = verbose
        )
    }

    # Keeping the original sample orders and the unwanted variable ####
    all.samples.index <- c(1:ncol(se.obj))
    initial.variable <- se.obj[[uv.variable]]
    initial.sample.names <- colnames(se.obj)
    colnames(se.obj) <- paste0('sample_', seq(ncol(se.obj)))

    # Assessing and grouping the unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    if (is.numeric(initial.variable)){
        se.obj[[uv.variable]] <- groupContinuousVariable(
            se.obj = se.obj,
            variable = uv.variable,
            nb.clusters = nb.clusters,
            clustering.method = clustering.method,
            perfix = '_group',
            verbose = verbose
        )
    }
    if (!is.numeric(initial.variable)){
        length.variable <- length(unique(initial.variable))
        if (length.variable == 1){
            stop('To create MNN, the "uv.variable" must have at least two groups/levels.')
        } else if (length.variable > 1){
            printColoredMessage(
                message = paste0(
                    '- The "',
                    uv.variable,
                    '" is a categorical variable with ',
                    length(unique(se.obj[[uv.variable]])),
                    ' levels.'),
                color = 'blue',
                verbose = verbose
            )
            se.obj[[uv.variable]] <- factor(x = se.obj[[uv.variable]])
        }
    }

    # Checking sample sizes of each sub group ####
    printColoredMessage(
        message = '-- Checking the sample size of each sub-group of the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    sub.group.sample.size <- findRepeatingPatterns(
        vec = se.obj[[uv.variable]],
        n.repeat = nb.mnn + 1
        )
    if (length(sub.group.sample.size) == 0){
        stop(paste0(
            'All subgroups of the unwanted variable have less than ',
            nb.mnn + 1,
            ' (nb.mnn + 1) samples. MNN cannot be found.')
            )
    } else if (length(sub.group.sample.size) != length(unique(se.obj[[uv.variable]])) ){
        printColoredMessage(
            message = paste0(
                'All or some subgroups of the unwanted variable have less than ',
                nb.mnn + 1,
                ' (nb.mnn + 1) samples. '),
            color = 'red',
            verbose = verbose
        )
    } else {
        printColoredMessage(
            message = paste0(
                '- All the subgroups of the unwanted variable have at least, ',
                nb.mnn + 1,
                ' (nb.mnn + 1) samples'),
            color = 'blue',
            verbose = verbose
        )
    }

    # Data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Applying data normalization, transformation and regression:',
        color = 'magenta',
        verbose = verbose
    )
    pairs.batch <- combn(
        x = sub.group.sample.size,
        m = 2
        )
    pb <- utils::txtProgressBar(
        min = 0,
        max = ncol(pairs.batch),
        style = 3
        )
    all.mnn.sets <- lapply(
        1:ncol(pairs.batch),
        function(x){
            index.samples <- se.obj[[uv.variable]] %in% pairs.batch[ , x]
            sub.se.obj <- se.obj[, index.samples]
            sub.se.obj[[uv.variable]] <- droplevels(sub.se.obj[[uv.variable]])
            norm.data <- preProcessData(
                se.obj = sub.se.obj,
                assay.name = assay.name,
                normalization = normalization,
                regress.out.variables = regress.out.variables,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = remove.na,
                verbose = verbose
                )
            # Selecting input data for knn analysis ####
            printColoredMessage(
                message = '-- Selecting the data input for knn analysis',
                color = 'magenta',
                verbose = verbose
                )
            if (data.input == 'expr' & !is.null(hvg)) {
                printColoredMessage(
                    message = '- Selecting the gene expression matrix with the highly variable genes as the data input.',
                    color = 'blue',
                    verbose = verbose
                    )
                norm.data <- t(norm.data[hvg,])
                }
            ## data input: expression matrix with all genes #####
            if (data.input == 'expr' & is.null(hvg)) {
                printColoredMessage(
                    message = paste0(
                        '- Selecting the gene expression matrix with all genes as the data input for',
                        'the sub-group: ',
                        x,
                        '.'),
                    color = 'blue',
                    verbose = verbose
                    )
                norm.data <- t(norm.data)
                }
            ## data input: PCA with hvg #####
            if (data.input == 'pcs' & !is.null(hvg)) {
                printColoredMessage(
                    message = paste0(
                        '- Performing PCA on the gene expression matrix using highly',
                        ' variable genes and using PCs as the data input.'),
                    color = 'blue',
                    verbose = verbose
                    )
                pca.data <- irlba::prcomp_irlba(
                    x = t(norm.data[hvg,]),
                    center = center,
                    scale. = scale
                    )
                cols.names <- colnames(norm.data)
                norm.data <- pca.data$x
                row.names(norm.data) <- cols.names
                }
            ## data input: PCA all genes #####
            if (data.input == 'pcs' & is.null(hvg)) {
                printColoredMessage(
                    message = '- Performing PCA on the gene expression matrix and select PCs as the data input.',
                    color = 'blue',
                    verbose = verbose
                    )
                pca.data <- irlba::prcomp_irlba(
                    x = t(norm.data),
                    center = center,
                    scale. = scale
                    )
                cols.names <- colnames(norm.data)
                norm.data <- pca.data$x
                row.names(norm.data) <- cols.names
            }
            # Finding MNN between batches ####
            printColoredMessage(
                message = '-- Finding MNN across all possible pairs of sub-groups of the unwanted variable:',
                color = 'magenta',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    '- All MNN between all ',
                    ncol(pairs.batch),
                    ' pairs of the sub-groups of the "',
                    uv.variable,
                    '" variable will be identified:'),
                color = 'orange',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    '* Finding MNN between the "',
                    pairs.batch[1, x],
                    '" and the "',
                    pairs.batch[2, x],
                    '" sub-groups.'),
                color = 'blue',
                verbose = verbose
                )
            # Applying cosine normalization ####
            if (isTRUE(apply.cosine.norm)){
                printColoredMessage(
                    message = '-- Applying cosine normalization of the data.',
                    color = 'blue',
                    verbose = verbose
                    )
                norm.data = cosineNorm(
                    x = norm.data,
                    mode = 'matrix'
                    )
            }
            sample.index.1 <- sub.se.obj[[uv.variable]] == pairs.batch[1 , x]
            data1 <- norm.data[sample.index.1 , ]
            sample.index.2 <- sub.se.obj[[uv.variable]] == pairs.batch[2 , x]
            data2 <- norm.data[sample.index.2 , ]
            # Applying mnn function ####
            mnn.samples <- BiocNeighbors::findMutualNN(
                data1 = data1,
                data2 = data2,
                k1 = nb.mnn,
                k2 = nb.mnn,
                BPPARAM = mnn.bpparam,
                nbparam = mnn.nbparam
                )
            # Checking the results of the mnn function ####
            if (is.null(mnn.samples)){
                stop('- MNN cannot be found.')
            }
            # Putting all the results of the mnn together ####
            mnn.data <- data.frame(
                group.a = rep(pairs.batch[1, x], length(mnn.samples$first)),
                group.b = rep(pairs.batch[2, x], length(mnn.samples$second)),
                overal.index.1 = all.samples.index[se.obj[[uv.variable]] == pairs.batch[1, x]][mnn.samples$first],
                overal.index.2 = all.samples.index[se.obj[[uv.variable]] == pairs.batch[2, x]][mnn.samples$second],
                overal.index.1.1 = as.numeric(gsub('sample_', '', row.names(data1)[mnn.samples$first])),
                overal.index.2.1 = as.numeric(gsub('sample_', '', row.names(data2)[mnn.samples$second]))
                )
            # Applying a sanity check ####
            if (!all.equal(mnn.data$overal.index.1, mnn.data$overal.index.1.1)){
                stop('There something wrong with the MNN.')
            }
            if (!all.equal(mnn.data$overal.index.2, mnn.data$overal.index.2.1)){
                stop('There something wrong with the MNN.')
            }
            printColoredMessage(
                message = paste0(
                    '- ',
                    nrow(mnn.data),
                    ' MNN are found.'),
                color = 'blue',
                verbose = verbose
            )
            setTxtProgressBar(pb, x)
            return(mnn.data)
        })

    all.mnn.sets <- do.call(rbind, all.mnn.sets)
    se.obj[[uv.variable]] <- initial.variable
    colnames(se.obj) <- initial.sample.names
    all.mnn.sets$sample.ids.1 <- colnames(se.obj)[all.mnn.sets$overal.index.1]
    all.mnn.sets$sample.ids.2 <- colnames(se.obj)[all.mnn.sets$overal.index.2]
    if (is.null(all.mnn.sets)){
        stop('MNN are not found across any batches. You may want to increase the valeu of the mnn.')
    }

    # Plotting the distribution of MNN  ####
    # First, calculate counts manually
    counts <- all.mnn.sets %>%
        count(group.a, group.b)
    p.mnn <- ggplot(all.mnn.sets, aes(x = group.a, y = group.b)) +
        geom_count() +
        geom_text(data = counts, aes(x = group.a, y = group.b, label = n), vjust = -0.5) +
        ggtitle('Distribution of MNN across batches') +
        xlab('') +
        ylab('') +
        scale_size_continuous(
            breaks = function(x) pretty(x[x == floor(x)], n = 3),
            name = "Count") +
        theme_bw() +
        theme(
            axis.line = element_line(colour = 'black', linewidth = 1),
            plot.title = element_text(size = 12),
            axis.text.x = element_text(
                size = 12,
                angle = 35,
                vjust = 1,
                hjust = 1),
            axis.text.y = element_text(
                size = 12,
                angle = 35,
                vjust = 1,
                hjust = 1)
        )
    if (isTRUE(plot.output)) print(p.mnn)
    message(' ')
    printColoredMessage(
        message = paste0(
            '- in total ' ,
            nrow(all.mnn.sets),
            ' MNN are found.'),
        color = 'blue',
        verbose = verbose
        )
    if (isTRUE(sum(is.na(all.mnn.sets))))
        stop('There are NA in the MNN. This is not supported.')

    # Saving the results ####
    ## Selecting prps name ####
    if(is.null(mnn.group.name)){
        mnn.group.name <- uv.variable
    }
    ## Selecting output name ####
    ## Selecting mnn sets name ####
    if (is.null(mnn.sets.name)){
        if (is.numeric(se.obj[[uv.variable]])){
            mnn.sets.name <- paste0(assay.name, '|', nb.clusters ,'groups|', nb.mnn, 'mnn')
        } else mnn.sets.name <- paste0(assay.name, '|', length(unique(se.obj[[uv.variable]])) ,'groups|', nb.mnn, 'mnn')
    }

    ## Saving the results in the SummarizedExperiment object ####
    message(' ')
    printColoredMessage(
        message = '-- Saving the results.',
        color = 'magenta',
        verbose = verbose
        )
    if (isTRUE(save.se.obj)) {
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
        se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]][[mnn.sets.name]][['data']]  <- all.mnn.sets
        if (!'plot' %in% names(se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.sets.name]])) {
            se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]][[mnn.sets.name]][['plot']] <- list()
        }
        se.obj@metadata[['KnnMnn']][['Mnn']][[mnn.group.name]][[mnn.sets.name]][['plot']] <- p.mnn
        printColoredMessage(
            message = '- All the mnn results are saved in the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The findMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## output the results as matrix  ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = '- All the mnn results are outputed as matrix.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The findMnn function finished.',
            color = 'white',
            verbose = verbose
            )
        return(list(
            mnn = all.mnn.sets,
            mnn.plot = p.mnn)
            )
    }
}
