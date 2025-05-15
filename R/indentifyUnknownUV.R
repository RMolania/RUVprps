#' Identify potential unknown sources of unwanted variation in RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function identifies unknown sources of unwanted variation in RNA-seq data using different robust statistical
#' approaches.
#'
#' @details
#' Identification of sources of unwanted variation is essential for creating PRPS data for RUV-III normalization.
#' There may be instances where certain sources of unwanted variation remain unrecorded; these are referred to as
#' "unknown" sources. This function uses three different approaches—`rle`, `pca`, and `sample.scoring`—to identify
#' potential unknown sources of unwanted variation when none are known.
#'
#' - In the `rle` approach, a clustering method specified by `clustering.methods` is applied to the RLE medians,
#' IQRs, or both. In the absence of unwanted variation, no distinguishable clusters should form.
#' - In the `pca` approach, principal component analysis is applied to either a set of negative control genes or
#' all genes. The first PCs are then clustered to detect unknown sources of variation.
#' - In the `sample.scoring` approach, samples are scored against gene sets (e.g., housekeeping genes) whose variation
#' may indicate unwanted variation. A clustering method is then applied to the scores.
#'
#' @references
#' 1. Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.
#' 2. Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character string. Specifies the assay name to be used for identifying potential sources of unwanted
#'  variation.
#' @param approach Character string. One of `rle`, `pca`, or `sample.scoring`, specifying the method to detect unknown variation.
#' @param rle.comp Character string. One of `median`, `iqr`, or `both`. Specifies which RLE summary statistic to use.
#'The default is set to `median`.
#' @param regress.out.bio.variables Character string or vector. Column name(s) of biological variables in the sample
#' annotation.The default is set to `NULL`.
#' @param regress.out.bio.gene.sets List. Biological gene signatures to regress out before identifying unwanted variation.
#' The default is set to `NULL`.
#' @param uv.gene.sets List. Gene sets related to unwanted variation for use in `sample.scoring`.The default is set to `NULL`.
#' @param ncg Vector. Negative control genes. If not `NULL`, analysis will be restricted to these genes.The default is set to `NULL`.
#' @param clustering.methods Character string. Clustering method: one of `kmeans`, `cut`, `quantile`, or `nbClust`.
#'The default is set to `nbClust`.
#' @param nbClust.diss Dissimilarity matrix. If provided, `nbClust.distance` must be `NULL`.The default is set to `NULL`.
#' @param nbClust.distance Character string. Distance measure: `euclidean`, `maximum`, `manhattan`, `canberra`, `binary`,
#'  `minkowski`, or `NULL`.The default is set to `euclidean`.
#' @param nbClust.min.nc Numeric. Minimum number of clusters. Must be between 1 and (n - 1).
#' @param nbClust.max.nc Numeric. Maximum number of clusters. Must be ≥ `nbClust.min.nc`.The default is set to 15.
#' @param nbClust.method Character string. Clustering method: `ward.D`, `ward.D2`, `single`, `complete`, `average`,
#' `mcquitty`, `median`, `centroid`, or `kmeans`.
#' @param nbClust.index Character string. Clustering index to evaluate: e.g., `silhouette`, `gap`, `ch`, `db`, `all`, etc.
#' @param nbClust.alphaBeale Numeric. Significance threshold for Beale's index.
#' @param max.samples.per.batch Numeric. Max proportion of samples per cluster (only applies to `nbClust`).The default is set to 0.1.
#' @param nb.clusters Numeric. Number of clusters when using `kmeans`, `cut`, or `quantile`.The default is set to 3.
#' @param apply.log Logical. Apply log-transformation before analysis.The default is set to `TRUE`.
#' @param pseudo.count Numeric. Pseudo count added before log-transforming data.The default is set to 1.
#' @param nb.pcs Numeric. Number of principal components to use in `pca` approach.The default is set to 2.
#' @param center Logical. Whether to center data before PCA.The default is set to `TRUE`.
#' @param scale Logical. Whether to scale data before PCA.The default is set to `FALSE`.
#' @param svd.bsparam A `BiocParallelParam` object. Controls parallelization for SVD computation.The default is set to `bsparam()`.
#' @param remove.current.estimates Character string. Whether to remove current estimates of unknown batches.The default is set to
#'  `TRUE`.
#' @param output.name Character string. Output file name. If `NULL`, a name is automatically generated.
#' @param assess.se.obj Logical. If `TRUE`, the `checkSeobj` function is applied to validate the object.The default is set to
#' `TRUE`.
#' @param remove.na Character string. Where to remove `NA` values: `assays`, `sample.annotation`, `both`, or `none`.
#'The default is set to `both`.
#' @param save.se.obj Logical. If `TRUE`, results are saved in the metadata under `metadata$UV$Unknown`.The default is set to `TRUE`.
#' @param plot.output Logical. If `TRUE`, generates clustering input plots colored by identified groups.The default is set to `TRUE`.
#' @param order.batches Logical. If `TRUE`, orders estimated batches in the plot.The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, displays messages throughout function execution.


#' @importFrom SummarizedExperiment assay colData
#' @importFrom singscore rankGenes simpleScore
#' @importFrom BiocSingular bsparam runSVD
#' @importFrom GGally ggpairs wrap
#' @importFrom stats as.formula
#' @importFrom NbClust NbClust
#' @importFrom dplyr arrange
#' @import RColorBrewer
#' @import ggplot2
#' @export

identifyUnknownUV <- function(
        se.obj,
        assay.name,
        approach = 'rle',
        rle.comp = 'median',
        regress.out.bio.variables = NULL,
        regress.out.bio.gene.sets = NULL,
        uv.gene.sets = NULL,
        ncg = NULL,
        clustering.methods = 'nbClust',
        nbClust.diss = NULL,
        nbClust.distance = "euclidean",
        nbClust.min.nc = 2,
        nbClust.max.nc = 5,
        nbClust.method = 'kmeans',
        nbClust.index = 'silhouette',
        nbClust.alphaBeale = 0.1,
        max.samples.per.batch = NULL,
        nb.clusters = 3,
        apply.log = TRUE,
        pseudo.count = 1,
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        remove.current.estimates = FALSE,
        output.name = NULL,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        plot.output = TRUE,
        order.batches = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The indentifyUnknownUV function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking function inputs ####
    if (is.null(assay.name)){
        stop('The "assay.name" cannot be empty.')
    }
    if (length(assay.name) > 1) {
        stop('The "assay.name" must be a single assay name.')
    }
    if (is.null(approach)){
        stop('The "approach" cannot be empty.')
    }
    if (length(approach) > 1){
        stop('The approach must be one of the "rle", "pca" or "sample.scoring".')
    }
    if (!approach %in% c('rle', 'pca', 'sample.scoring') ) {
        stop('The approach must be one of the "rle", "pca" or "sample.scoring".')
    }
    if (is.logical(regress.out.bio.variables)){
        stop('The "regress.out.bio.variables" should be either NULL or a vector of variable names.')
    }
    if (!is.null(regress.out.bio.variables)) {
        if (sum(regress.out.bio.variables %in% colnames(colData(se.obj))) != length(regress.out.bio.variables) )
            stop('Some or all "regress.out.bio.variables" variables are not found in the SummarizedExperiment object.')
    }
    if (is.logical(regress.out.bio.gene.sets)){
        stop('The "regress.out.bio.gene.sets" should be either NULL or a list of genes.')
    }
    if (!is.null(regress.out.bio.gene.sets)){
        lapply(
            regress.out.bio.gene.sets,
            function(x){
                if (sum(regress.out.bio.gene.sets %in% row.names(se.obj)) == 0)
                    stop('The "regress.out.bio.gene.sets" are not found in the SummarizedExperiment object.')
            })
    }
    if (is.logical(uv.gene.sets)){
        stop('The "uv.gene.sets" should be either NULL or a list of genes.')
    }
    if (!is.null(uv.gene.sets)){
        lapply(
            names(uv.gene.sets),
            function(x){
                if (sum(uv.gene.sets[[x]] %in% row.names(se.obj)) < 2)
                    stop('Some or all genes of "uv.gene.sets" are not found in the SummarizedExperiment object.')
            })
    }
    if (is.logical(ncg)){
        stop('The "ncg" should be either NULL or vector of genes ids.')
    }
    if (!is.null(ncg)) {
        if (sum(ncg %in% row.names(se.obj)) == 0)
            stop('Some or "ncg" genes are found in the SummarizedExperiment object.')
    }
    if (length(clustering.methods) > 1) {
        stop('A single method should be provided for the "clustering.methods".')
    }
    if (clustering.methods == 'nbClust'){
        if (is.null(nbClust.min.nc)){
            stop('The "nbClust.min.nc" must be specified, when the clustering.methods is nbClust.')
        } else if (nbClust.min.nc < 0 | nbClust.min.nc == 1){
            stop('The "nbClust.min.nc" must be equal or more than 2, when the clustering.methods is nbClust.')
        } else if (is.null(nbClust.max.nc)){
            stop('The "nbClust.max.nc" must be specified, when the clustering.methods is nbClust.')
        } else if (nbClust.max.nc < 0 | nbClust.max.nc < 2){
            stop('The "nbClust.max.nc" must be equal or more than 2, when the clustering.methods is nbClust.')
        } else if (!nbClust.method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans")){
            stop('The "nbClust.method" must be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans."')
        } else if (!nbClust.index %in% c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db",
                                         "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
                                         "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw",
                                         "all", "alllong")){
            stop('The nbClust.index must of one of: "kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db",
                                         "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
                                         "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw",
                                         "all", "alllong"')
        } else if (!is.null(max.samples.per.batch)){
            if (max.samples.per.batch == 0 | max.samples.per.batch < 0  | max.samples.per.batch >= 1){
                stop('The value of max.samples.per.batch must be between 0<max.samples.per.batch<1.')
            }
        }
    }

    if (!clustering.methods %in% c('kmeans', 'cut', 'quantile', 'nbClust')) {
        stop('The clustering.methods should be one of "kmeans", "cut", "quantile" or "nbClust".')
    }
    if (clustering.methods %in% c('kmeans', 'cut', 'quantile')){
        if (is.null(nb.clusters))
            stop('The "nb.clusters" must be specified when the "clustering.methods" is kmeans, cut or quantile.')
    }
    if (approach == 'rle'){
        if (!rle.comp %in% c('median', 'iqr', 'both'))
            stop('The "rle.comp" should be one of "median", "iqr" or "both".')
    }
    if (approach == 'pca'){
        if (nb.pcs == 0 | is.null(max.samples.per.batch))
            stop('The value of "nb.pcs" should be more than 0 when the arroach is equal to pca.')
    }
    if (isTRUE(apply.log)){
        if (pseudo.count < 0)
            stop('The valuse of pseudo.count cannot be negative.')
    }
    if (approach == 'pca' & nb.pcs > 1 & clustering.methods %in% c('cut', 'quantile')){
        stop(paste0('The nb.pcs should be 1 to use the ', clustering.methods, ' method for clustering.'))
    }
    if (is.null(regress.out.bio.variables) & remove.na == 'both'){
        stop('The "remove.na" cannot be set to "both" when the "regress.out.bio.variables = NULL".')
    } else if (is.null(regress.out.bio.variables) & remove.na == 'sample.annotation'){
        stop('The "remove.na" cannot be set to "sample.annotation" when the "regress.out.bio.variables = NULL".')
    }
    if (is.logical(output.name)){
        stop('The "output.name" should be eitehr NULL or a character.')
    }

    # Removing current unwanted variation estimates for the assay ####
    if (isTRUE(remove.current.estimates)){
        printColoredMessage(
            message = paste0('The current estimated unknown batches:'),
            color = 'magenta',
            verbose = verbose
            )
        if (!'UnKnownUV' %in%  names(se.obj@metadata)) {
            printColoredMessage(
                message = paste0('- There is not any estimated unknown batches in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
                )
        } else  if (assay.name %in% names(se.obj@metadata[['UnKnownUV']])) {
            printColoredMessage(
                message = paste0('- The current estimated unknown batches for the  ', assay.name, ' data is removed.'),
                color = 'blue',
                verbose = verbose
                )
            se.obj@metadata[['UnKnownUV']][[assay.name]] <- list()
        } else {
            printColoredMessage(
                message = paste0('- There is not any estimated unknown batches for the  ', assay.name, ' data.'),
                color = 'blue',
                verbose = verbose
                )
        }
    }
    # Checking the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(regress.out.bio.variables),
            remove.na = remove.na,
            verbose = verbose
            )
    }
    # Data log transformation ####
    if (isTRUE(apply.log)){
        printColoredMessage(
            message = '-- Applying log transformation on all the specified assay(s):',
            color = 'magenta',
            verbose = verbose
            )
        expr.data <- applyLog(
            se.obj = se.obj,
            assay.names = assay.name,
            pseudo.count = pseudo.count,
            )[[assay.name]]
    }
    if (isFALSE(apply.log)){
        printColoredMessage(
            message = '-- The specified assay will be used without applying log transformation.',
            color = 'blue',
            verbose = verbose
            )
        expr.data <- assay(x = se.obj, i = assay.name)

    }
    # Regressing out variables and biological gene sets ####
    ## regressing out biological variables ####
    if (!is.null(regress.out.bio.variables) & is.null(regress.out.bio.gene.sets)){
        printColoredMessage(
            message = paste0(
                '- The ',
                paste0(regress.out.bio.variables, collapse = ' & '),
                ' variables will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- t(expr.data)
        lm.formula <- paste('se.obj', regress.out.bio.variables, sep = '$')
        adjusted.data <- lm(as.formula(paste('expr.data', paste0(lm.formula, collapse = '+') , sep = '~')))
        expr.data <- t(adjusted.data$residuals)
        colnames(expr.data) <- colnames(se.obj)
        row.names(expr.data) <- row.names(se.obj)
        }
    ## regressing out biological gene sets ####
    if (is.null(regress.out.bio.variables) & !is.null(regress.out.bio.gene.sets)){
        printColoredMessage(
            message = '- Calculating sample scores for individual gene sets of the "regress.out.bio.gene.sets" list.',
            color = 'blue',
            verbose = verbose
        )
        ranked.data <- rankGenes(expr.data)
        regress.out.bio.gene.sets <- sapply(
            regress.out.bio.gene.sets,
            function(x) singscore::simpleScore(
                rankData = ranked.data,
                upSet = x)$TotalScore
        )
        rm(ranked.data)
        printColoredMessage(
            message = paste0(
                '- The sample scores of individual gene list of ',
                'regress.out.bio.gene.sets',
                ' will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose)
        expr.data <- t(expr.data)
        adjusted.data <- lm(expr.data~ regress.out.bio.gene.sets)
        expr.data <- t(adjusted.data$residuals)
        colnames(expr.data) <- colnames(se.obj)
        row.names(expr.data) <- row.names(se.obj)
        }
    ## regressing out biological variable and biological gene sets ####
    if (!is.null(regress.out.bio.variables) & !is.null(regress.out.bio.gene.sets)){
        printColoredMessage(
            message = paste0(
                '- The sample scores of individual gene list of ',
                ' "regress.out.bio.gene.sets" and the ',
                paste0(regress.out.bio.variables, collapse = ' & '),
                ' variables will be regressed out from the data,',
                ' please make sure your data is log transformed.'),
            color = 'blue',
            verbose = verbose
            )
        ranked.data <- rankGenes(expr.data)
        regress.out.bio.gene.sets <- sapply(
            regress.out.bio.gene.sets,
            function(x) singscore::simpleScore(
                rankData = ranked.data,
                upSet = x)$TotalScore
            )
        all.variables <- as.data.frame(cbind(
            regress.out.bio.gene.sets,
            as.data.frame(colData(se.obj)[, regress.out.bio.variables, drop = FALSE]))
            )
        lm.formula <- paste('all.variables', colnames(all.variables), sep = '$')
        expr.data <- t(expr.data)
        adjusted.data <- lm(as.formula(paste('expr.data', paste0(lm.formula, collapse = '+') , sep = '~')))
        expr.data <- t(adjusted.data$residuals)
        colnames(expr.data) <- colnames(se.obj)
        row.names(expr.data) <- row.names(se.obj)
        }

    # Selecting data input for clustering ####
    printColoredMessage(
        message = '-- Selecting input data for clustering:',
        color = 'magenta',
        verbose = verbose
        )
    ## PCA approach ####
    if (approach == 'pca'){
        if (is.null(ncg)){
            ### PCA on all genes ####
            printColoredMessage(
                message = paste0(
                    '- Applying PCA on the data and use the first ',
                    nb.pcs,
                    ' PCs as an input for clustering.'),
                color = 'blue',
                verbose = verbose
            )
            set.seed(2233)
            sv.dec <- runSVD(
                x = t(expr.data),
                k = nb.pcs,
                BSPARAM = svd.bsparam,
                center = center,
                scale = scale
            )
            input.data <- sv.dec$u
            colnames(input.data) <- c(paste0('PC', 1:ncol(input.data)))
            if (clustering.methods == 'nbClust'){
                input.data.name <- paste0(
                    approach,
                    '|AllGenes_nbClust.',
                    nbClust.method,
                    'Clustering')
            } else input.data.name <- paste0(
                approach,
                '|AllGenes|',
                clustering.methods,
                'Clustering'
            )
        }
        if (!is.null(ncg)){
            ### PCA on all NCG ####
            printColoredMessage(
                message = paste0(
                    '- Applying PCA on the data using the "ncg" gene only, and use the first ',
                    nb.pcs,
                    ' PCs as an input for clustering.'),
                color = 'blue',
                verbose = verbose
            )
            set.seed(2233)
            sv.dec <- runSVD(
                x = t(expr.data[ncg , ]),
                k = nb.pcs,
                BSPARAM = svd.bsparam,
                center = center,
                scale = scale)
            input.data = sv.dec$u
            if (clustering.methods == 'nbClust'){
                input.data.name <- paste0(
                    approach,
                    '|NCG|nbClust.',
                    nbClust.method,
                    'Clustering'
                )
            } else input.data.name <- paste0(
                approach,
                '|NCG|',
                clustering.methods,
                'Clustering'
            )
        }
    }
    ## RLE approach ####
    if (approach == 'rle'){
        if (is.null(ncg)){
            ### RLE on all the genes ####
            printColoredMessage(
                message = paste0('- Applying tge RLE on the data.'),
                color = 'blue',
                verbose = verbose
                )
            rle.data <- expr.data - rowMedians(expr.data)
            if (clustering.methods == 'nbClust'){
                input.data.name <- paste0(
                    approach,
                    '.',
                    rle.comp,
                    '|AllGenes|nbClust.',
                    nbClust.method,
                    'Clustering'
                    )
            } else input.data.name <- paste0(approach, '.', rle.comp,'|AllGenes|', clustering.methods, 'Clustering')
        } else if (!is.null(ncg)){
            ### RLE on all the NCG ####
            printColoredMessage(
                message = '- Applying the RLE on the data using only "ncg" genes.',
                color = 'blue',
                verbose = verbose
                )
            rle.data <- expr.data[ncg , ] - rowMedians(expr.data[ncg , ])
            if (clustering.methods == 'nbClust'){
                input.data.name <- paste0(
                    approach,
                    '.',
                    rle.comp,
                    '|NCG|nbClust.',
                    nbClust.method,
                    'Clustering'
                    )
            } else input.data.name <- paste0(
                approach,
                '.',
                rle.comp,
                '|NCG|',
                clustering.methods,
                'Clustering'
                )
        }
        if (rle.comp == 'median'){
            printColoredMessage(
                message = '- Useing the RLE medians as input data.',
                color = 'blue',
                verbose = verbose
            )
            input.data <- colMedians(rle.data)
        } else if (rle.comp == 'iqr'){
            printColoredMessage(
                message = '- Use the RLE IQRs as input data for clustering.',
                color = 'blue',
                verbose = verbose
            )
            input.data <- colIQRs(rle.data)
        }
    }
    ## Sample scoring approach ####
    if (approach == 'sample.scoring'){
        if (is.null(uv.gene.sets)){
            all.uv.gene.sets <- list(ncg = ncg)
        } else if (!is.null(uv.gene.sets) & !is.null(ncg)){
            all.uv.gene.sets <- uv.gene.sets[['ncg']] <- ncg
        } else if (!is.null(uv.gene.sets) & is.null(ncg)){
            all.uv.gene.sets <- uv.gene.sets
        }
        printColoredMessage(
            message = paste0(
                '- Calculating sample scores for individual gene set(s) of the ',
                '"uv.gene.sets" and using the scores as input for clustering.'),
            color = 'blue',
            verbose = verbose
            )
        ranked.data <- rankGenes(expr.data)
        input.data <- sapply(
            names(all.uv.gene.sets),
            function(x) singscore::simpleScore(
                rankData = ranked.data,
                upSet = all.uv.gene.sets[[x]])$TotalScore
            )
        # names(input.data) <- names(all.uv.gene.sets)
        rm(ranked.data)
        if (clustering.methods == 'nbClust'){
            input.data.name <- paste0(
                approach,
                '|nbClust.',
                nbClust.method,
                'Clustering'
                )
        } else input.data.name <- paste0(
            approach,
            '|',
            clustering.methods,
            'Clustering'
            )
    }
    # Clustering ####
    printColoredMessage(
        message = '- Clustering the inpute data',
        color = 'magenta',
        verbose = verbose
        )
    ## kmeans ####
    if (clustering.methods == 'kmeans'){
        printColoredMessage(
            message = paste0(
                '- Applying kmeans with centers = ',
                nb.clusters,
                ' on the data'),
            color = 'blue',
            verbose = verbose
            )
        set.seed(3344)
        groups <- kmeans(
            x = input.data,
            centers = nb.clusters,
            iter.max = 10000)$cluster
        uv.sources <- paste0('Batch' , groups)
    }
    ## cut ####
    if (clustering.methods == 'cut'){
        printColoredMessage(
            message = paste0(
                '- Applying the cut method with breaks = ',
                nb.clusters,
                ' on the data'),
            color = 'blue',
            verbose = verbose
            )
        groups <- as.numeric(cut(
                x = input.data,
                breaks = nb.clusters,
                include.lowest = TRUE
            ))
        uv.sources <- paste0('Batch' , groups)
    }
    ## quantile ####
    if (clustering.methods == 'quantile'){
        printColoredMessage(
            message = paste0(
                '- Apply the quantile method with probs = ',
                paste0(round(seq(0, 1, 1/nb.clusters), digits = 2)),
                ' on the data'),
            color = 'blue',
            verbose = verbose
            )
        quantiles <- quantile(x = input.data, probs = seq(0, 1, 1 / nb.clusters))
        groups <- as.numeric(cut(
                x = input.data,
                breaks = quantiles,
                include.lowest = TRUE
            ))
        uv.sources <- paste0('Batch' , groups)
    }
    ## nbClust ####
    if (clustering.methods == 'nbClust'){
        if (is.numeric(max.samples.per.batch)){
            printColoredMessage(
                message = '- Applying the nbClust method on the summarized data.',
                color = 'blue',
                verbose = verbose
            )
            initial.clusters <- NbClust(
                data = input.data,
                diss = nbClust.diss,
                distance = nbClust.distance,
                min.nc = nbClust.min.nc,
                max.nc = nbClust.max.nc,
                method = nbClust.method,
                index = nbClust.index,
                alphaBeale = nbClust.alphaBeale
            )
            batch.samples <- data.frame(
                id = colnames(se.obj),
                batch = initial.clusters$Best.partition
            )
            selected.clusters <- findRepeatingPatterns(
                vec = batch.samples$batch,
                n.repeat = round(max.samples.per.batch * ncol(se.obj), digits = 0)
            )
            while (length(selected.clusters) > 0) {
                more.clusters <- lapply(
                    selected.clusters,
                    function(x) {
                        index <- batch.samples$batch == x
                        if (is.matrix(input.data)) {
                            sub.input.data <- input.data[index , ]
                        } else
                            sub.input.data <- input.data[index]
                        sub.clusters <- NbClust(
                            data = sub.input.data,
                            diss = nbClust.diss,
                            distance = nbClust.distance,
                            min.nc = nbClust.min.nc,
                            max.nc = nbClust.max.nc,
                            method = nbClust.method,
                            index = nbClust.index,
                            alphaBeale = nbClust.alphaBeale
                        )
                        data.frame(
                            id = batch.samples$id[index],
                            batch = paste0(x, sub.clusters$Best.partition)
                        )
                    })
                more.clusters <- do.call(rbind, more.clusters)
                batch.samples$batch[match(more.clusters$id, batch.samples$id)] <-
                    more.clusters$batch
                selected.clusters <- findRepeatingPatterns(
                    vec = batch.samples$batch,
                    n.repeat = round(max.samples.per.batch * ncol(se.obj), digits = 0)
                )
            }
            uv.sources <- paste0('Batch', as.numeric(as.factor(batch.samples$batch)))
        }
        if (is.null(max.samples.per.batch)){
            initial.clusters <- NbClust(
                data = input.data,
                diss = nbClust.diss,
                distance = nbClust.distance,
                min.nc = nbClust.min.nc,
                max.nc = nbClust.max.nc,
                method = nbClust.method,
                index = nbClust.index,
                alphaBeale = nbClust.alphaBeale
                )
            uv.sources <- paste0('Batch', initial.clusters$Best.partition)
        }

    }

    # Reporting the number of the possible batches ####
    printColoredMessage(
        message = paste0(
            length(unique(uv.sources)),
            ' potential batches are found in the ',
            assay.name,
            ' data.'),
        color = 'blue',
        verbose = verbose
        )

    # Plotting the outputs ####
    currentCols <-  c(
        RColorBrewer::brewer.pal(8, "Dark2")[-5],
        RColorBrewer::brewer.pal(10, "Paired"),
        RColorBrewer::brewer.pal(12, "Set3"),
        RColorBrewer::brewer.pal(9, "Blues")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Oranges")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Greens")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Purples")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Reds")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "Greys")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "BuGn")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "PuRd")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "BuPu")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(9, "YlGn")[c(8, 3, 7, 4, 6, 9, 5)],
        RColorBrewer::brewer.pal(10, "Paired")
    )
    colors.selected <- currentCols[1:length(unique(uv.sources))]
    names(colors.selected) <- sort(unique(uv.sources))
    if (!is.matrix(input.data)) {
        data.to.plot <- data.frame(
            input.data = input.data,
            batches = factor(x = paste0('Batch', as.numeric(as.factor(uv.sources))),
                             levels = paste0('Batch', sort(unique(
                                 as.numeric(as.factor(uv.sources))
                             )))))
        if (isTRUE(order.batches))
            data.to.plot <- arrange(data.to.plot, batches)
        data.to.plot$samples <- c(1:ncol(se.obj))
        p <- ggplot(data = data.to.plot, aes(x = samples, y = input.data, color = batches)) +
            geom_point() +
            xlab('Samples') +
            ylab(paste0('Input data (', approach, ')')) +
            scale_color_manual(values = colors.selected, name = 'Batches') +
            theme(
                panel.background = element_blank(),
                legend.key = element_blank(),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12)
            ) +
            guides(colour = guide_legend(override.aes = list(size = 5)))
        if (isTRUE(plot.output))
            print(p)

    } else{
        if (ncol(input.data) == 1) {
            data.to.plot <- as.data.frame(input.data)
            input <- samples <- batches <- NULL
            colnames(data.to.plot) <- 'input'
            data.to.plot$batches <- factor(
                x = uv.sources,
                levels = sort(unique(uv.sources))
                )
            if (isTRUE(order.batches))
                data.to.plot <- arrange(data.to.plot, batches)
            data.to.plot$samples <- c(1:ncol(se.obj))
            p <- ggplot(data = data.to.plot, aes(
                    x = samples,
                    y = input,
                    color = batches
                )) +
                geom_point() +
                xlab('Samples') +
                ylab(paste0('Input data (', approach, ')')) +
                scale_color_manual(values = colors.selected, name = 'Batch') +
                theme(
                    panel.background = element_blank(),
                    legend.key = element_blank(),
                    legend.text = element_text(size = 12),
                    legend.title = element_text(size = 14),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12),
                    axis.text.x = element_text(size = 9),
                    axis.text.y = element_text(size = 9)
                ) +
                guides(colour = guide_legend(override.aes = list(size = 5)))
            if (isTRUE(plot.output))
                print(p)

        } else {
            data.to.plot <- as.data.frame(input.data)
            data.to.plot$batches <- factor(x = paste0('Batch', as.numeric(as.factor(uv.sources))),
                                           levels = paste0('Batch', sort(unique(
                                               as.numeric(as.factor(uv.sources))
                                           ))))
            p <- GGally::ggpairs(
                data = data.to.plot[, 1:(ncol(data.to.plot) - 1)],
                mapping = ggplot2::aes(colour = data.to.plot[, ncol(data.to.plot)]),
                showStrips = FALSE,
                switch = 'y',
                labeller = NULL,
                diag = list(continuous = wrap("diagAxis", labelSize = 8, diagAxis = 0)),
                upper = "blank") +
                theme(
                    panel.background = element_blank(),
                    legend.key = element_blank(),
                    legend.text = element_text(size = 12),
                    legend.title = element_text(size = 14),
                    panel.grid.major = element_blank(),
                    axis.ticks = element_blank(),
                    strip.background = element_blank(),
                    strip.text.x.bottom = element_text(size = 0),
                    strip.text = element_text(size = 0),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 2),
                    axis.text.x = element_text(size = 0),
                    axis.text.y = element_text(size = 0)
                ) +
                scale_color_manual(values = colors.selected)
            if (isTRUE(plot.output))
                print(p)
        }
    }

    # Saving the results ####
    printColoredMessage(message = '- Save the the results:',
                        color = 'magenta',
                        verbose = verbose)
    # out put name ####
    if (is.null(output.name)){
        input.data.name <- paste0(length(unique(uv.sources)), 'batches|', input.data.name)
    } else input.data.name <- output.name

    if (isTRUE(save.se.obj)){
        if (!'UnKnownUV' %in%  names(se.obj@metadata)) {
            se.obj@metadata[['UnKnownUV']] <- list()
        }
        if (!assay.name %in%  names(se.obj@metadata[['UnKnownUV']])){
            se.obj@metadata[['UnKnownUV']][[assay.name]] <- list()
        }
        if (!input.data.name %in%  names(se.obj@metadata[['UnKnownUV']][[assay.name]])){
            se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]] <- list()
        }
        se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['batches']] <-
            paste0('Batch', as.numeric(as.factor(uv.sources)))
        se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['input.data']] <- input.data
        se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['plot']] <- p
        printColoredMessage(
            message = 'The potentail unknow sources of variation are saved to the metadata of the SummarizedExperiment object',
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The indentifyUnknownUV function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj)
    }
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = 'The results are outputed as list.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The indentifyUnknownUV function finished.',
            color = 'white',
            verbose = verbose
            )
        return(list(
            batches = paste0('Batch', as.numeric(as.factor(uv.sources))),
            input.data = input.data,
            plot = p )
            )
    }
}
