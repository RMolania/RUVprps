#' Finds k-nearest neighbors in RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function finds k nearest neighbors for individual samples within each groups of unwanted variation in RNA-seq
#' data. Further, the distance between individual neighbors are also calculated. The k nearest neighbors will be used to
#' create pseudo-samples within individual batches.
#'
#' @param se.obj A summarized experiment object.
#' @param assay.name Character. A character representing the name of the data (assay) in the SummarizedExperiment object
#'  to be used to find k nearest neighbors.
#' @param uv.variable Character. A character that indicates the name of the column in the sample annotation of the
#' SummarizedExperiment object. The `uv.variable` can be either categorical or continuous. If `uv.variable` is a continuous
#' variable, this will be divided into `nb.clusters` groups using the `clustering.method`.
#' @param nb.knn Numeric. A numeric number that indicates the maximum number of nearest neighbors to compute for each
#' sample. The default is set to 3.
#' @param clustering.method Character. A character indicating the choice of clustering method for grouping the
#' `uv.variable` if a continuous variable is provided. Options include `kmeans`, `cut`, and `quantile`. The default is
#' set to `kmeans`.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the `uv.variable` is a
#' continuous variable. The default is set to 3.
#' @param data.input Character. A character that indicates which data type should be used as input for finding the k nearest
#' neighbors. Options include: `expr` and `pcs`. If `pcs` is selected, the first `nb.pcs` of PCs of the data will be used
#' as input. If `expr` is selected, the expression data will be used as input. The default is set to `expr`.
#' @param nb.pcs Numeric. A numeric value that indicates the number of PCs to be calculated and then used as data input
#' for finding the k nearest neighbors. The default is set to 2. The `nb.pcs` must be set when `data.input = pcs`.
#' @param center Logical. Indicates whether to center the data or not before calculating PCs. If center is `TRUE`, then
#' centering is done by subtracting the column means of the assay from their corresponding columns. The default is set
#' to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data or not before calculating PCs. If scale is set to `TRUE`, then
#' scaling is done by dividing the (centered) columns of the assays by their standard deviations if center is `TRUE`, and
#' the root mean square otherwise. The default is set to `FALSE`.
#' @param svd.bsparam Character. A BiocParallelParam object specifying how palatalization should be performed for performing
#' PCA using SVD. The default is set to `bsparam()`. We refer to the `runSVD()` function from the **BiocSingular** R package
#' for further details.
#' @param hvg Vector. A logical vector or a vector of the names (feature ids) of the highly variable genes. These genes
#' will be used to prepare the input data for knn analysis. The default is set to `NULL`, this means all genes will be used.
#' @param normalization Character. A character that indicates which normalization method should be applied on the
#' data before finding the knn. Options are: `CPM`,`TMM`, `upper`, `median`, `full` and `VST`. The default is set to `CPM`.
#' If set to `NULL`, no normalization will be applied. See the `applyOtherNormalizations()` function for more details.
#' @param regress.out.variables Character. A character or a vector characters that indicate the column name(s) in the sample
#' annotation in the SummarizedExperiment object. These variables will be regressed out from the data before finding KNN.
#' The default is set to `NULL`, indicating that regression will not be applied.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not for down-stream analysis.
#' The default is set to `TRUE`.
#' @param pseudo.count Numeric. A positive numeric value as a pseudo count to be added to all measurements of the specified
#' assay (data) before applying log transformation to avoid `-Inf`for measurements that are equal to 0. The default is set
#' to 1.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default is set
#' to `TRUE`. See the `checkSeObj()` function for more details.
#' @param remove.na Character. To remove NA or missing values from the assay (data) or not. The options are `assays` and
#''none'. The default is set to `assays`, so all the NA or missing values from the data(assay) will be removed before
#' computing performing any down-stream analysis. See the `checkSeObj()` function for more details.
#' @param knn.group.name Character. A character specifying the name of the PRPS group to which the current KNN belong.
#' If it is set to `NULL`, the function will automatically assign a name using `uv.variable`.
#' @param knn.sets.name Character. A character specifying the name of the output file to be saved in the metadata
#' of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0(assay.name, '|', length(unique(se.obj[[uv.variable]])) ,'groups|', nb.knn, 'lnn')`.
#' @param save.se.obj Logical. Indicates whether to save the KNN results in the metadata of the SummarizedExperiment object
#' or to output the result as a list. The defaults is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, shows the messages of different steps of the function.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom utils txtProgressBar
#' @importFrom stats dist
#' @importFrom RANN nn2
#' @export

findKnn <- function(
        se.obj,
        assay.name,
        uv.variable,
        nb.knn = 3,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        data.input = 'expr',
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        hvg = NULL,
        normalization = 'CPM',
        regress.out.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        check.se.obj = TRUE,
        remove.na = 'both',
        knn.group.name = NULL,
        knn.sets.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findKnn function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the inputs #####
    if (is.list(assay.name)) {
        stop('The "assay.name" cannot be a list.')
    }
    if (length(assay.name) > 1) {
        stop('The "assay.name" must be a name of signle assay(data) in the SummarizedExperiment object.')
    }
    if (is.null(uv.variable)) {
        stop('The "uv.variable" variable cannot be empty or NULL.')
    }
    if (length(uv.variable) > 1) {
        stop('The "uv.variable" must contain a name of signle variable in the SummarizedExperiment object.')
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
    if (!is.numeric(nb.knn)){
        stop('The "nb.knn" must be a postive numeric value.')
    }
    if (nb.knn <= 0 & nb.knn != as.integer(nb.knn)) {
        stop('The "nb.knn" must be a postive whole numeric value')
    }

    if (!is.null(hvg)){
        if (is.logical(hvg)){
            if (length(hvg) != nrow(se.obj)){
                stop('The length of the "hvg" logical vector should the same as the number fo rows in the SummarizedExperiment object.')
            }
            if (sum(hvg) == 0){
                stop('The "hvg" does not contain any selected genes. All are FALSE')
            }
        }
        if (is.character(hvg)){
            if (sum(hvg %in% row.names(se.obj)) != length(hvg))
                stop('All or some of the "hvg" are not found in the SummarizedExperiment object.')
        }
        if (is.numeric(hvg)){
            stop('The "hvg" must be either logical vector or feature ids.')
        }
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
    if (!is.null(knn.sets.name)){
        if (!is.character(knn.sets.name)){
            stop('The "knn.sets.name" must be a character.')
        }
    }
    if (!is.null(knn.group.name)){
        if (!is.character(knn.group.name)){
            stop('The "knn.group.name" must be logical.')
        }
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical.')
    }
    # Assessing the SummarizedExperiment object ####
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
    colnames.seobj <- colnames(se.obj)
    colnames(se.obj) <- paste0(
        colnames.seobj,
        '_',
        all.samples.index
        )

    # Assessing and grouping the unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    initial.variable <- se.obj[[uv.variable]]
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

    # Checking sample sizes of each sub group ####
    printColoredMessage(
        message = '-- Checking the sample size of each subgroup of the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    sub.group.sample.size <- findRepeatingPatterns(
        vec = se.obj[[uv.variable]],
        n.repeat =  nb.knn + 1
        )
    if (length(sub.group.sample.size) == 0){
        stop(paste0(
            'All subgroups of the specified unwanted variable have less than ',
             nb.knn + 1,
            ' samples. KNN cannot be found.')
            )
    } else if (length(sub.group.sample.size) != length(unique(se.obj[[uv.variable]])) ){
        printColoredMessage(
            message = paste0(
                'All or some subgroups of the specified unwanted variable have less than ',
                 nb.knn + 1,
                ' (nb.knn + 1) samples. Then KNN for those sub-groups cannot be created.'),
            color = 'red',
            verbose = verbose
        )
    } else {
        printColoredMessage(
            message = paste0(
                '- All the sub-groups of the specified unwanted variable have at least ',
                 nb.knn + 1,
                ' (nb.knn + 1) samples.'),
            color = 'blue',
            verbose = verbose
            )
    }

    # Applying data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Applying data normalization, transformation and regression:',
        color = 'magenta',
        verbose = verbose
        )
    all.norm.data <- lapply(
        sub.group.sample.size,
        function(x) {
            selected.samples <- colData(se.obj)[[uv.variable]] == x
            ## Applying library size normalization ####
            if (!is.null(normalization) & is.null(regress.out.variables)) {
                printColoredMessage(
                    message = paste0(
                        '- Applying the ',
                        normalization,
                        ' on the samples from the "',
                        x,
                        '" sub-group.'),
                    color = 'blue',
                    verbose = verbose
                )
                norm.data <- applyOtherNormalizations(
                    se.obj = se.obj[, selected.samples],
                    assay.name = assay.name,
                    method = normalization,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    check.se.obj = FALSE,
                    save.se.obj = FALSE,
                    remove.na = 'none',
                    verbose = FALSE
                )
            }
            ## Applying library size normalization and regressing out variables ####
            if (!is.null(normalization) & !is.null(regress.out.variables)) {
                printColoredMessage(
                    message = paste0(
                        '- Applying the ',
                        normalization,
                        ' on the samples from "',
                        x,
                        '" group and then regressing out the ',
                        paste0(regress.out.variables, collapse = '&'),
                        ' variable(s) from the data.'),
                    color = 'blue',
                    verbose = verbose
                )
                ### normalization ####
                norm.data <- applyOtherNormalizations(
                    se.obj = se.obj[, selected.samples],
                    assay.name = assay.name,
                    method = normalization,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    check.se.obj = FALSE,
                    save.se.obj = FALSE,
                    remove.na = 'none',
                    verbose = FALSE
                )
                ## regression ####
                sample.info <- as.data.frame(colData(se.obj[, selected.samples]))
                norm.data <- t(norm.data)
                lm.formua <- paste('sample.info', regress.out.variables, sep = '$')
                norm.data <- lm(as.formula(paste(
                    'norm.data',
                    paste0(lm.formua, collapse = '+') ,
                    sep = '~'
                )))
                norm.data <- t(norm.data$residuals)
                colnames(norm.data) <- colnames(norm.data)
                row.names(norm.data) <- row.names(norm.data)

            }
            ## Regressing out variables ####
            if (is.null(normalization) & !is.null(regress.out.variables)){
                if (isTRUE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- Applying log transformation and then regressing out the ',
                            paste0(regress.out.variables, collapse = '&'),
                            ' variable(s) on the samples ',
                            x,
                            '" group from the data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    if (!is.null(pseudo.count)){
                        norm.data <- log2(assay(se.obj[, selected.samples], assay.name) + pseudo.count)
                    } else {
                        norm.data <- log2(assay(se.obj[, selected.samples], i = assay.name))
                    }

                } else if (isFALSE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- Regressing out ',
                            paste0(regress.out.variables, collapse = '&'),
                            x,
                            '" group from the data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    norm.data <- assay(se.obj[, selected.samples], assay)
                }
                ### regression ####
                sample.info <- as.data.frame(colData(se.obj[, selected.samples]))
                norm.data <- t(norm.data)
                lm.formua <- paste('sample.info', regress.out.variables, sep = '$')
                norm.data <- lm(as.formula(paste(
                    'norm.data',
                    paste0(lm.formua, collapse = '+') ,
                    sep = '~'
                )))
                norm.data <- t(norm.data$residuals)
                colnames(norm.data) <- colnames(norm.data)
                row.names(norm.data) <- row.names(norm.data)
            }
            ## Applying log transformation ####
            if (is.null(normalization) & is.null(regress.out.variables)) {
                if (isTRUE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- Applying the log2 within the samples from "',
                            x,
                            '" group data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    if(!is.null(pseudo.count)){
                        norm.data <- log2(assay(x = se.obj[, selected.samples], i = assay.name) + pseudo.count)
                    } else {
                        norm.data <- log2(assay(se.obj[, selected.samples], i = assay.name))
                    }

                } else if (isFALSE(apply.log)){
                    printColoredMessage(
                        message = paste0(
                            '- No library size normalization and transformation is applied on data from ',
                            x,
                            '" group data.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    norm.data <- assay(se.obj[, selected.samples], assay)
                }
            }
            return(norm.data)
        })
    names(all.norm.data) <- sub.group.sample.size

    # Selecting input data for KNN analysis ####
    printColoredMessage(
        message = '-- Selecting the data input for the knn analysis',
        color = 'magenta',
        verbose = verbose
        )
    all.data.input <- lapply(
        sub.group.sample.size,
        function(x){
            norm.data <- all.norm.data[[x]]
            ## Using expression data with hvg #####
            if (data.input == 'expr' & !is.null(hvg)) {
                printColoredMessage(
                    message = '- Selecting the gene expression matrix with the highly variable genes as the data input.',
                    color = 'blue',
                    verbose = verbose
                    )
                norm.data <- t(norm.data[hvg,])
            }
            ## Using gene expression data with all genes #####
            if (data.input == 'expr' & is.null(hvg)) {
                printColoredMessage(
                    message = paste0(
                        '- Selecting the gene expression matrix with all genes as the data input for',
                        ' the sub-group: ',
                        x,
                        '.'
                        ),
                    color = 'blue',
                    verbose = verbose
                    )
                norm.data <- t(norm.data)
            }
            ## Using PCA with hvg #####
            if (data.input == 'pcs' & !is.null(hvg)) {
                printColoredMessage(
                    message = paste0(
                        '- Performing PCA on the gene expression matrix using highly',
                        ' variable genes and using PCs as the data input.'),
                    color = 'blue',
                    verbose = verbose
                    )
                sv.dec <- BiocSingular::runSVD(
                    x = t(norm.data[hvg,]),
                    k = nb.pcs,
                    BSPARAM = svd.bsparam,
                    center = center,
                    scale = scale
                    )
                norm.data <- sv.dec$u
            }
            ## Using PCA with all genes #####
            if (data.input == 'pcs' & is.null(hvg)) {
                printColoredMessage(
                    message = '- Performing PCA on the gene expression matrix and selecting PCs as the data input.',
                    color = 'blue',
                    verbose = verbose
                    )
                sv.dec <- BiocSingular::runSVD(
                    x = t(norm.data),
                    k = nb.pcs,
                    BSPARAM = svd.bsparam,
                    center = center,
                    scale = scale
                    )
                norm.data <- sv.dec$u
            }
            return(norm.data)
        })
    names(all.data.input) <- sub.group.sample.size

    # Finding k nearest neighbors ####
    printColoredMessage(
        message = '-- Finding k nearest neighbors within each sub-group of the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    printColoredMessage(
        message = paste0(
            '- For individual samples within the selected sub-group(s) of the variable "',
            uv.variable,
            '", k = ',
             nb.knn,
            ' nearest neighbours will be found.'),
        color = 'orange',
        verbose = verbose
        )
    pb <- utils::txtProgressBar(
        min = 0,
        max = length(sub.group.sample.size),
        style = 3
        )
    all.knn <- lapply(
        1:length(sub.group.sample.size),
        function(x) {
            norm.data <- all.data.input[[sub.group.sample.size[x]]]
            ## Finding knn ####
            printColoredMessage(
                message = paste0(
                    '- Finding k nearest neighbors within the "',
                    sub.group.sample.size[x],
                    '" sub-group.'),
                color = 'blue',
                verbose = verbose
                )
            knn.samples <- RANN::nn2(
                data = norm.data,
                query = norm.data,
                k = c (nb.knn + 1),
                treetype = 'bd'
                )
            knn.index <- as.data.frame(knn.samples$nn.idx)
            colnames(knn.index) <- paste0('dataset.index', c(1:c(nb.knn + 1)))
            selected.samples <- se.obj[[uv.variable]] == sub.group.sample.size[x]
            ovral.cell.no <- sapply(
                1:ncol(knn.index),
                function(x) {
                    all.samples.index.temp <- all.samples.index[selected.samples]
                    all.samples.index.temp[knn.index[ , x]]
                    }
                )
            ## Applying a sanity check ####
            if ( length(unique(se.obj[[uv.variable]][as.vector(ovral.cell.no)])) != 1 ){
                stop('There are something wrong with the sample annotation.')
            }
            if (!unique(se.obj[[uv.variable]][as.vector(ovral.cell.no)]) == sub.group.sample.size[x] ){
                stop('There are something wrong with the sample annotation.')
            }

            colnames(ovral.cell.no) <- paste0('overal.index', 1:c( nb.knn + 1))
            knn.index <- as.data.frame(cbind(ovral.cell.no , knn.index))
            # Computing distance between all knn ####
            printColoredMessage(
                message = paste0(
                    '- Calculating all pairwise distances between the k nearest neighbors within the "',
                    sub.group.sample.size[x],
                    '" subgroup.'),
                color = 'blue',
                verbose = verbose
                )
            knn.dis <- round(as.data.frame(knn.samples$nn.dists), digits = 3)
            knn.dis <- knn.dis[,-1, drop = FALSE]
            colnames(knn.dis) <- paste0('distance1_', 2:c(nb.knn + 1))
            knn.index.dist <- cbind(knn.index, knn.dis)
            if (isTRUE(nb.knn > 1)) {
                all.comb <- combn(
                    x = paste0('dataset.index', 2:c(nb.knn + 1)),
                    m = 2
                    )
                all.comb.names <- combn(x = 2:c(nb.knn + 1), m = 2)
                for (z in 1:ncol(all.comb)) {
                    pair.dist <- unlist(lapply(
                        1:nrow(knn.index.dist),
                        function(y) {
                            col1 <- knn.index.dist[, all.comb[, z][1]][y]
                            col2 <- knn.index.dist[, all.comb[, z][2]][y]
                            stats::dist(
                                x = norm.data[c(col1, col2) ,],
                                method = "euclidean",
                                diag = FALSE,
                                upper = FALSE
                            )
                        }))
                    name <- paste0('dist',
                               all.comb.names[, z][1],
                               '_',
                               all.comb.names[, z][2])
                    knn.index.dist[[name]] <- pair.dist
                }
            }
            cols.index <- grep('distance', colnames(knn.index.dist))
            knn.index.dist$aver.dist <- rowMeans(knn.index.dist[, cols.index, drop = FALSE])
            set.seed(4589)
            knn.index.dist$rank.aver.dist <- rank(
                x = knn.index.dist$aver.dist,
                ties.method = 'random'
                )
            knn.index.dist$group <- sub.group.sample.size[x]
            setTxtProgressBar(pb, x)
            message(' ')
            return(knn.index.dist)
        })
    all.knn <- do.call(rbind, all.knn)
    all.knn <- all.knn[order(all.knn$overal.index1) , ]
    row.names(all.knn) <- c(1:ncol(se.obj))
    sample.ids <- unlist(lapply(colnames(se.obj), function(x) sub("_[^_]*$", "", x)))
    all.knn$sample.ids.1 <- sample.ids[all.knn$overal.index1]
    samples.ids <- sapply(
        1:nb.knn,
        function(x){
            sample.ids[all.knn[ , x+1]]
        })
    colnames(samples.ids) <- paste0('sample.ids.', c(c(1:nb.knn) + 1))
    all.knn <- cbind(all.knn , samples.ids)
    se.obj[[uv.variable]] <- initial.variable
    colnames(se.obj) <- sample.ids

    # Saving the results ####
    ## Selecting knn group name ####
    if(is.null(knn.group.name)){
        knn.group.name <- uv.variable
    }
    ## Selecting knn sets name ####
    if (is.null(knn.sets.name)){
        if (is.numeric(se.obj[[uv.variable]])){
            knn.sets.name <- paste0(assay.name, '|', nb.clusters ,'groups|', nb.knn, 'knn')
        } else knn.sets.name <- paste0(assay.name, '|', length(unique(se.obj[[uv.variable]])) ,'groups|', nb.knn, 'knn')
    }
    ## Saving the results in the SummarizedExperiment object ####
    message(' ')
    printColoredMessage(
        message = '-- Saving the results:',
        color = 'magenta',
        verbose = verbose
        )
    if (isTRUE(save.se.obj)) {
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

        printColoredMessage(
            message = '- All the knn results are saved in the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The findKnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## Outputing the results as matrix  ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '- All the knn results are outputed as matrix.',
                            color = 'blue',
                            verbose = verbose)
        printColoredMessage(message = '------------The findKnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(all.knn)
    }
}
