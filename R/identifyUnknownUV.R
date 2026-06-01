#' Identify unknown sources of unwanted variation in RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @references
#' - Gandolfo L.C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. *PLoS ONE*. \url{https://doi.org/10.1371/journal.pone.0191629}
#' - Molania R., et al. Removing unwanted variation from large-scale RNA sequencing data with PRPS. *Nature Biotechnology*. \url{https://doi.org/10.1038/s41587-022-01440-w}
#'
#' @description
#' This function estimates unknown sources of unwanted variation in RNA-seq data using different robust statistical
#' approaches.
#'
#' @details
#' Identification of sources of unwanted variation is essential for creating PRPS data for RUV-III normalization.
#' There may be instances where certain sources of unwanted variation remain unrecorded; these are referred to as
#' "unknown" sources. This function uses three different approaches, `rle`, `pca`, and `sample.scoring`, to identify
#' potential unknown sources of unwanted variation when none are known.
#'
#' - In the `rle` approach, a clustering method specified by the `clustering.methods` argument is applied to the RLE medians,
#' IQRs, or both. In the absence of unwanted variation, no distinguishable clusters should form. Any clusters will be considered
#' as soruces of unwanted variation e.g., batch effects.
#' - In the `pca` approach, principal component analysis is applied to either a set of negative control genes or
#' all genes. The first PCs are then clustered to detect unknown sources of variation.
#' - In the `sample.scoring` approach, samples are scored against gene sets (e.g., housekeeping genes) whose variation
#' may indicate unwanted variation. A clustering method is then applied to the scores.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character that specifies a data (assay) name in the SummarizedExperiment object to be
#' used for identifying potential sources of unwanted variation.
#' @param approach Character. A character string that specifies which approach should be used to estimate unwanted variation
#' in the data. The options are `rle`, `pca`, or `sample.scoring`. The default is set to `rle`.
#' @param rle.comp Character.  A character string that specifies which RLE summary statistic to be used when the  `approach = rle`.
#' The options are `median`, `iqr`, or `both`. The default is set to `median`.
#' @param uv.gene.sets List. A list of genes sets that capture unwanted variation and to be used when  `approach = sample.scoring`
#' .The default is set to `NULL`. This must be specified when `approach = sample.scoring`.
#' @param chronological.detection Logical. If samples are ordered based on a chronological e.g. time or sequencing, ...,
#' then a chronological analysis can find any patterns in the data. If this information is not available, this argument
#' mus be set to `FALSE`. The default is set to `FALSE`. Refer to the details for more information.
#' @param changepoint.type Character. A character that specifies which chronological test should be used. The options are
#' `mean` and `meanvar`. If `mean` selected, the function will apply the `cpt.mean` function, and if `meanvar` selected
#' the function will apply `cpt.meanvar` function. Refer to the `changepoint` package for more information.
#' @param clustering.methods Character. A character that indicates which clustering method should be used to estimate
#' unwanted variation. The options are  `kmeans`, `cut`, `quantile`, or `nbClust`. The default is set to `nbClust`.
#' @param ncg Vector. A vector of negative control genes that are assumed to be affected by unwanted variation. If it is
#' specified, the `rle` and `pca` approaches will be performed using only these genes. he default is et to `NULL`.
#' @param regress.out.bio.variables Character. A character or a vector of characters that indifcates name(s) of the column(s)
#' in the SummarizedExperiment that contains biological variable(s). These variable (s) will be regressed out before applying
#' any data stigmatization process. The default is set to `NULL`.
#' @param regress.out.bio.gene.sets List. A list of biological gene signatures to be used to score all samples. Then the
#' score(s) wuill be regress out from the data before any data stigmatization process. The default is set to `NULL`.
#' @param nbClust.diss Dissimilarity matrix. If provided, `nbClust.distance` must be `NULL`.The default is set to `NULL`.
#' @param nbClust.distance Character. Distance measure: `euclidean`, `maximum`, `manhattan`, `canberra`, `binary`,
#'  `minkowski`, or `NULL`.The default is set to `euclidean`.
#' @param nbClust.min.nc Numeric. Minimum number of clusters. Must be between 1 and (n - 1).
#' @param nbClust.max.nc Numeric. Maximum number of clusters. Must be ≥ `nbClust.min.nc`.The default is set to 15.
#' @param nbClust.method Character. Clustering method: `ward.D`, `ward.D2`, `single`, `complete`, `average`,
#' `mcquitty`, `median`, `centroid`, or `kmeans`.
#' @param nbClust.index Character. Clustering index to evaluate: e.g., `silhouette`, `gap`, `ch`, `db`, `all`, etc.
#' @param nbClust.alphaBeale Numeric. Significance threshold for Beale's index.
#' @param max.samples.per.batch Numeric. Max proportion of samples per cluster (only applies to `nbClust`).The default
#' is set to `NULL`.
#' @param mclust.max.clusters Numeric. A numeric value that specified the maximum number of mixture components (clusters)
#' for which the BIC is to be calculated by the `Mclust` function. The default is set to 20. Refer to the `mclust` package
#' for more details.
#' @param nb.clusters Numeric. A numeric value that species the number of clusters when using the `approach` is set to one
#' of the `kmeans`, `cut`, or `quantile`.The default is set to 3.
#' @param cpt.penalty Character. A character that specifies the penalty approach for the `cpt.mean` and `cpt.meanvar`
#' functions. The options are `None`, `SIC`, `BIC`, `MBIC`, `AIC`, `Hannan-Quinn`, `Asymptotic`, `Manual` and `CROPS`. The
#' default is set to `MBIC`. Refer to the `changepoint` package for more information.
#' @param cpt.pen.value  Numeric. from the the `changepoint` package: The theoretical type I error e.g.0.05 when using
#' the Asymptotic penalty. A vector of length 2 (min,max) if using the CROPS penalty. The value of the penalty when using
#' the Manual penalty option - this can be a numeric value or text giving the formula to use. Available variables are,
#' n=length of original data, null=null likelihood, alt=alternative likelihood, tau=proposed changepoint,
#' diffparam=difference in number of alternatve and null parameters. The defualt is set to 0.
#' @param cpt.method Character. The options are `AMOC`, `PELT`, `SegNeigh`, or `BinSeg`. Default is set to `PELT`.
#' @param cpt.q Numeric. From the `changepoint` package: the maximum number of changepoints to search for using the
#' "BinSeg" method. The maximum number of segments (number of changepoints + 1) to search for using the "SegNeigh" method.
#' @param cpt.test.stat Character. The assumed test statistic / distribution of the data. The options are `Normal` and
#' `CUSUM`. The default is set to `Normal`.
#' @param cpt.minseglen Numeric. Positive integer giving the minimum segment length (no. of observations between changes),
#'  default is the minimum allowed by theory.
#' @param apply.log Logical. Indicating whether to apply log-transformation before analysis.The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value that indicates a pseudo count to be added to all the genes expression
#' before log-transforming data.The default is set to 1.
#' @param nb.pcs Numeric. A uumberical value that indicates the number of principal components to be calculated when the
#' approach is set to `pca`. The default is set to 2.
#' @param center Logical. Indicates whether to center the data before applying SVD. If center is `TRUE`, centering is
#' performed by subtracting the column means of the data from their corresponding columns. The default is set to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data before applying SVD. If scale is set to `TRUE`, scaling is
#' done by dividing the (centered) columns of the assays by their standard deviations,  if center is `TRUE`, and the root
#' mean square otherwise. The default is set to `FALSE`.
#' @param svd.bsparam A `BiocParallelParam` object specifying how parallelization should be performed. The default is set
#' to `bsparam()`. See the `runSVD()` function from the BiocSingular R package for more details.
#' @param remove.current.estimates Logical. Indicating whether to remove current estimates of unknown batches in the
#' SummarizedExperiment object.The default is set to `TRUE`.
#' @param assess.bio.association Logical. Indicating whether to assess the association of the estimated batch effects with
#' any known biological variables. The default is set to `FALSE`.
#' @param bio.variables Character. A Character or a vector of characters specifying the column names of biological
#' variables in the sample annotation of the 'SummarizedExperiment  object. These 'bio.variables' can be either categorical
#' or continuous variables. The default is set to `NULL`.
#' @param bio.clustering.method Character. A Character specifying the clustering method to be applied for grouping each
#' continuous biological variable. Options include `kmeans`, `cut`, and `quantile`. The default is set to `kmeans` clustering.
#' @param nb.bio.clusters Numeric. A value indicating the number of groups for continuous sources of biological variation.
#' The default is 3. This implies that each continuous variable will be split into 3 groups using the specified
#' `clustering.method`.
#' @param assess.uv.association Logical. Indicating whether to assess the association of the estimated batch effects with
#' any known unwanted variables. The default is set to `FALSE`.
#' @param uv.variables Character. A character or a vector of characters specifying the column names of unwanted variation
#' variables in the sample annotation of the 'SummarizedExperiment  object. These 'uv.variables' can be either categorical
#' or continuous variables.
#' @param uv.clustering.method Character. A character specifying the clustering method to be applied for grouping each
#' continuous unwanted variation variable. Options include `kmeans`, `cut`, and `quantile`. The default is set to `kmeans`
#' clustering.
#' @param nb.uv.clusters Numeric. A numeric value value indicating the number of groups for continuous sources of unwanted
#' variation variables. The default is set to 3. This implies that each continuous variable will be split into 3 groups
#' using the specified `uv.clustering.method`.
#' @param generate.association.plot Logical. Indicating whether to generate a heatmap of the estimated unwanted variation
#' across any specified `bio.variables` and `uv.variables`. The default is set to `FALSE`
#' @param plot.output Logical. Indicating whether to plot the outputs. The default is set to `TRUE`. Individual plots are
#' stored in the SummarizedExperiment by default.
#' @param color.palette Character. A character string indicating which color palette should be used. The options are
#' `nrc`, `pan.selection.a`, `pan.selection.b` and `pan.selection.c`, The default is set to `nrc`.
#' @param order.batches Logical. If `TRUE`, orders estimated batches in the plot. This can better visualize the patterns
#'  of the estimated unwanted variation. The default is set to `TRUE`.
#' @param add.to.sample.annotation Logical. If `TRUE`, the estimated unknown unwanted variation
#' will be stored in the sample annotation of the SummarizedExperiment object. If `FALSE`, the estimated unknown unwanted
#' variation will be stored in the metadata or returned as a list.
#' @param col.name Character. A character string specifying the name of the new column  in the sample annotation where
#' the estimated unknown unwanted variation will be stored. The default is `Estimated.batches`.
#' @param output.name Character. A character that specifies the name of the out files in the SummarizedExperiment object.
#' If `NULL`, a name is automatically generated base on: `paste0(length(unique(uv.sources)),'batches|',input.data.name)`.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. The default is set to `TRUE`.
#' See the `checkSeObj()` function for more details.
#' @param remove.na Character. A character that specifies whether to remove NA or missing values from the data sets (assays).
#' The options are `assays`, `sample.annotation` and `none`. The default is set to `assays`.
#' @param save.se.obj Logical. If `TRUE`, results are saved in the metadata under `metadata$UV$Unknown`.The default is set
#' to `TRUE`.
#' @param verbose Logical. If `TRUE`, displays messages for the different steps of the function.
#'
#' @importFrom changepoint cpt.mean cpt.meanvar cpts
#' @importFrom SummarizedExperiment assay colData
#' @importFrom singscore rankGenes simpleScore
#' @importFrom BiocSingular bsparam runSVD
#' @importFrom GGally ggpairs wrap
#' @importFrom stats as.formula
#' @importFrom NbClust NbClust
#' @importFrom dplyr arrange
#' @importFrom grid gpar
#' @import RColorBrewer
#' @import ggplot2
#'
#' @export

identifyUnknownUV <- function(
        se.obj,
        assay.name,
        approach                  = 'rle',
        rle.comp                  = 'median',
        uv.gene.sets              = NULL,
        chronological.detection   = FALSE,
        changepoint.type          = 'meanvar',
        clustering.methods        = 'nbClust',
        ncg                       = NULL,
        regress.out.bio.variables = NULL,
        regress.out.bio.gene.sets = NULL,
        nbClust.diss              = NULL,
        nbClust.distance          = "euclidean",
        nbClust.min.nc            = 3,
        nbClust.max.nc            = 30,
        nbClust.method            = 'kmeans',
        nbClust.index             = 'silhouette',
        nbClust.alphaBeale        = 0.1,
        max.samples.per.batch     = NULL,
        mclust.max.clusters       = 20,
        nb.clusters               = 3,
        cpt.penalty               = 'MBIC',
        cpt.pen.value             = 0,
        cpt.method                = 'PELT',
        cpt.q                     = 5,
        cpt.test.stat             = 'Normal',
        cpt.minseglen             = 1,
        apply.log                 = TRUE,
        pseudo.count              = 1,
        nb.pcs                    = 2,
        center                    = TRUE,
        scale                     = FALSE,
        svd.bsparam               = bsparam(),
        remove.current.estimates  = FALSE,
        assess.bio.association    = FALSE,
        bio.variables             = NULL,
        nb.bio.clusters           = 3,
        bio.clustering.method     = 'kmeans',
        assess.uv.association     = FALSE,
        uv.variables              = NULL,
        nb.uv.clusters            = 3,
        uv.clustering.method      = 'kmeans',
        generate.association.plot = FALSE,
        plot.output               = TRUE,
        color.palette             = 'pan.selection.a',
        order.batches             = FALSE,
        add.to.sample.annotation  = TRUE,
        col.name                  = 'Estimated.batches',
        output.name               = NULL,
        check.se.obj              = TRUE,
        remove.na                 = 'none',
        save.se.obj               = TRUE,
        verbose                   = TRUE
        ){
    printColoredMessage(
        message = '------------The indentifyUnknownUV function starts:',
        color   = 'white',
        verbose = verbose
        )
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
    if (!is.logical(chronological.detection)){
        stop('The "chronological.detection" must be logical')
    }
    if (!changepoint.type %in% c('mean', 'meanvar') ) {
        stop('The "changepoint.type" must be one of the "mean" or "meanvar".')
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

    if (!clustering.methods %in% c('kmeans', 'cut', 'quantile', 'nbClust', 'mclust')) {
        stop('The clustering.methods should be one of "kmeans", "cut", "quantile", "nbClust" and "mclust".')
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
        if (nb.pcs == 0 | is.null(nb.pcs))
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
    }
    if (is.null(regress.out.bio.variables) & remove.na == 'sample.annotation'){
        stop('The "remove.na" cannot be set to "sample.annotation" when the "regress.out.bio.variables = NULL".')
    }
    if (is.logical(output.name)){
        stop('The "output.name" should be eitehr NULL or a character.')
    }

    # Removing current unwanted variation estimates for the assay ####
    if (isTRUE(remove.current.estimates)){
        printColoredMessage(
            message = paste0('- The current estimated unknown batches:'),
            color   = 'magenta',
            verbose = verbose
            )
        if (!'UnKnownUV' %in%  names(se.obj@metadata)) {
            printColoredMessage(
                message = paste0('- There is not any estimated unknown batches in the SummarizedExperiment object.'),
                color   = 'blue',
                verbose = verbose
                )
        } else if (assay.name %in% names(se.obj@metadata[['UnKnownUV']])) {
            printColoredMessage(
                message = paste0(
                    '- The current estimated unknown batches for the  ',
                    assay.name,
                    ' data is removed.'),
                color   = 'blue',
                verbose = verbose
                )
            se.obj@metadata[['UnKnownUV']][[assay.name]] <- list()
        } else {
            printColoredMessage(
                message = paste0(
                    '- There is not any estimated unknown batches for the  ',
                    assay.name,
                    ' data.'),
                color   = 'blue',
                verbose = verbose
                )
        }
    }

    # Checking the SummarizedExperiment object ####
    if (isTRUE(check.se.obj)) {
        se.obj <- checkSeObj(
            se.obj      = se.obj,
            assay.names = assay.name,
            variables   = c(regress.out.bio.variables),
            remove.na   = remove.na,
            verbose     = verbose
            )
    }
    # Applying data log transformation ####
    if (isTRUE(apply.log)){
        printColoredMessage(
            message = '- Applying log transformation on all the specified assay(s):',
            color = 'magenta',
            verbose = verbose
            )
        expr.data <- applyLog(
            se.obj       = se.obj,
            assay.names  = assay.name,
            check.se.obj = FALSE,
            pseudo.count = pseudo.count,
            verbose      = verbose
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
        expr.data            <- t(expr.data)
        lm.formula           <- paste('se.obj', regress.out.bio.variables, sep = '$')
        adjusted.data        <- lm(as.formula(paste('expr.data', paste0(lm.formula, collapse = '+') , sep = '~')))
        expr.data            <- t(adjusted.data$residuals)
        colnames(expr.data)  <- colnames(se.obj)
        row.names(expr.data) <- row.names(se.obj)
        }
    ## Regressing out biological gene sets ####
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
                upSet    = x)$TotalScore
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
        expr.data            <- t(expr.data)
        adjusted.data        <- lm(expr.datab~ regress.out.bio.gene.sets)
        expr.data            <- t(adjusted.data$residuals)
        colnames(expr.data)  <- colnames(se.obj)
        row.names(expr.data) <- row.names(se.obj)
        }
    ## Regressing out biological variable and biological gene sets ####
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
                upSet    = x)$TotalScore
            )
        all.variables <- as.data.frame(cbind(
            regress.out.bio.gene.sets,
            as.data.frame(colData(se.obj)[, regress.out.bio.variables, drop = FALSE]))
            )
        lm.formula           <- paste('all.variables', colnames(all.variables), sep = '$')
        expr.data            <- t(expr.data)
        adjusted.data        <- lm(as.formula(paste('expr.data', paste0(lm.formula, collapse = '+') , sep = '~')))
        expr.data            <- t(adjusted.data$residuals)
        colnames(expr.data)  <- colnames(se.obj)
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
                x       = t(expr.data),
                k       = nb.pcs,
                BSPARAM = svd.bsparam,
                center  = center,
                scale   = scale
            )
            input.data <- sv.dec$u
            colnames(input.data) <- c(paste0('PC', 1:ncol(input.data)))
            if (clustering.methods == 'nbClust'){
                input.data.name <- paste0(
                    approach,
                    '|AllGenes_nbClust.',
                    nbClust.method,
                    'Clustering'
                    )
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
                x       = t(expr.data[ncg , ]),
                k       = nb.pcs,
                BSPARAM = svd.bsparam,
                center  = center,
                scale   = scale
                )
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
                message = paste0('- Applying the RLE on the data.'),
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
            } else input.data.name <- paste0(
                approach,
                '.',
                rle.comp,
                '|AllGenes|',
                clustering.methods,
                'Clustering'
                )
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
                message = '- Using the RLE medians as input data.',
                color = 'blue',
                verbose = verbose
            )
            input.data <- colMedians(rle.data)
        } else if (rle.comp == 'iqr'){
            printColoredMessage(
                message = '- Using the RLE IQRs as input data for clustering.',
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
            color   = 'blue',
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
    # Clustering analysis ####
    if (isFALSE(chronological.detection)){
        printColoredMessage(
            message = '- Clustering the inpute data',
            color   = 'magenta',
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
            quantiles <- quantile(
                x     = input.data,
                probs = seq(0, 1, 1 / nb.clusters)
                )
            groups <- as.numeric(cut(
                x              = input.data,
                breaks         = quantiles,
                include.lowest = TRUE
            ))
            uv.sources <- paste0('Batch' , groups)
        }
        ## mclust ####
        if (clustering.methods == 'mclust'){
            printColoredMessage(
                message = '- Applying the mclust function :',
                color = 'blue',
                verbose = verbose
            )
            set.seed(3344)
            groups <- Mclust(
                data    = input.data,
                G       = 2:mclust.max.clusters,
                verbose = verbose)
            uv.sources <- paste0('Batch' , groups$classification)
        }
        ## nbClust ####
        if (clustering.methods == 'nbClust'){
            ### Considering maximum samples per clusters ####
            if (is.numeric(max.samples.per.batch)){
                printColoredMessage(
                    message = '- Applying the nbClust method on the summarized data.',
                    color = 'blue',
                    verbose = verbose
                )
                initial.clusters <- NbClust(
                    data       = input.data,
                    diss       = nbClust.diss,
                    distance   = nbClust.distance,
                    min.nc     = nbClust.min.nc,
                    max.nc     = nbClust.max.nc,
                    method     = nbClust.method,
                    index      = nbClust.index,
                    alphaBeale = nbClust.alphaBeale
                )
                batch.samples <- data.frame(
                    id    = colnames(se.obj),
                    batch = initial.clusters$Best.partition
                )
                selected.clusters <- findRepeatingPatterns(
                    vec      = batch.samples$batch,
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
                                data       = sub.input.data,
                                diss       = nbClust.diss,
                                distance   = nbClust.distance,
                                min.nc     = nbClust.min.nc,
                                max.nc     = nbClust.max.nc,
                                method     = nbClust.method,
                                index      = nbClust.index,
                                alphaBeale = nbClust.alphaBeale
                            )
                            data.frame(
                                id    = batch.samples$id[index],
                                batch = paste0(x, sub.clusters$Best.partition)
                            )
                        })
                    more.clusters <- do.call(rbind, more.clusters)
                    batch.samples$batch[match(more.clusters$id, batch.samples$id)] <-
                        more.clusters$batch
                    selected.clusters <- findRepeatingPatterns(
                        vec      = batch.samples$batch,
                        n.repeat = round(max.samples.per.batch * ncol(se.obj), digits = 0)
                    )
                }
                uv.sources <- paste0('Batch', as.numeric(as.factor(batch.samples$batch)))
            }
            ## Without considering maximum samples per clusters ####
            if (is.null(max.samples.per.batch)){
                nb.clusters <- NbClust(
                    data       = input.data,
                    diss       = nbClust.diss,
                    distance   = nbClust.distance,
                    min.nc     = nbClust.min.nc,
                    max.nc     = nbClust.max.nc,
                    method     = nbClust.method,
                    index      = nbClust.index,
                    alphaBeale = nbClust.alphaBeale
                )
                uv.sources <- paste0('Batch', nb.clusters$Best.partition)
            }
        }
    }
    ## Chronological order analysis ####
    if (isTRUE(chronological.detection)){
        printColoredMessage(
            message = '- Finiding chronological patterns in the inpute data',
            color = 'magenta',
            verbose = verbose
        )
        ### Chronological order analysis using mean ####
        if (changepoint.type == 'mean'){
            cpt <- cpt.mean(
                data            = input.data,
                penalty         = cpt.penalty,
                pen.value       = cpt.pen.value,
                method          = cpt.method,
                Q               = cpt.q,
                test.stat       = cpt.test.stat,
                minseglen       = cpt.minseglen,
                class           = TRUE,
                param.estimates = TRUE
                )
            cpts <- cpts(cpt)
            if (length(cpts) == 0){
                stop('The "cpt.mean" function cannot find any chronological patterns in the data.')
            }
            uv.sources <- rep(
                x = 1:(length(cpts) + 1),
                times = diff(c(0, cpts, length(input.data)))
            )
            uv.sources <- paste0('Batch', uv.sources)
        }
        ### Chronological order analysis using mean and variance ####
        if (changepoint.type == 'meanvar'){
            cpt <- cpt.meanvar(
                data            = input.data,
                penalty         = cpt.penalty,
                pen.value       = cpt.pen.value,
                method          = cpt.method,
                Q               = cpt.q,
                test.stat       = cpt.test.stat,
                minseglen       = cpt.minseglen,
                class           = TRUE,
                param.estimates = TRUE
            )
            cpts <- cpts(cpt)
            if (length(cpts) == 0){
                stop('The "cpt.meanvar" function cannot find any chronological patterns in the data.')
            }
            uv.sources <- rep(
                x = 1:(length(cpts) + 1),
                times = diff(c(0, cpts, length(input.data)))
            )
            uv.sources <- paste0('Batch', uv.sources)
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
    colors.selected <- selectColors(
        nb.color = 1:length(unique(uv.sources)) ,
        group = color.palette
        )
    names(colors.selected) <- sort(unique(uv.sources))
    if (!is.matrix(input.data)) {
        data.to.plot <- data.frame(
            input.data = input.data,
            batches = factor(
                x = paste0('Batch', as.numeric(as.factor(uv.sources))),
                levels = paste0('Batch', sort(unique(as.numeric(as.factor(uv.sources)))))
                ))
        if (isTRUE(order.batches))
            data.to.plot <- arrange(data.to.plot, batches)
        data.to.plot$samples <- c(1:ncol(se.obj))
        p.batches <- ggplot(data = data.to.plot, aes(x = samples, y = input.data, color = batches)) +
            geom_point() +
            geom_smooth(aes(group = 1), method = "loess", se = FALSE, span = 0.2, color = "grey") +
            scale_color_manual(values = colors.selected, name = 'Batches') +
            xlab('Samples') +
            ylab(paste0('Input data (', approach, ')')) +
            theme(
                panel.background = element_blank(),
                legend.key       = element_blank(),
                legend.text      = element_text(size = 12),
                legend.title     = element_text(size = 14),
                axis.line        = element_line(colour = 'black', linewidth = 1),
                axis.title.x     = element_text(size = 14),
                axis.title.y     = element_text(size = 14),
                axis.text.x      = element_text(size = 12),
                axis.text.y      = element_text(size = 12)) +
            guides(colour        = guide_legend(override.aes = list(size = 5)))
        if (isTRUE(plot.output)) print(p.batches)
    } else {
        if (ncol(input.data) == 1) {
            data.to.plot <- as.data.frame(input.data)
            input <- samples <- batches <- NULL
            colnames(data.to.plot) <- 'input'
            data.to.plot$batches <- factor(
                x = uv.sources,
                levels = sort(unique(uv.sources))
                )
            if (isTRUE(order.batches)) data.to.plot <- arrange(data.to.plot, batches)
            data.to.plot$samples <- c(1:ncol(se.obj))
            p.batches <- ggplot(data = data.to.plot, aes(
                    x = samples,
                    y = input,
                    color = batches )) +
                geom_point() +
                xlab('Samples') +
                ylab(paste0('Input data (', approach, ')')) +
                scale_color_manual(values = colors.selected, name = 'Batch') +
                theme(
                    panel.background = element_blank(),
                    legend.key       = element_blank(),
                    legend.text      = element_text(size = 12),
                    legend.title     = element_text(size = 14),
                    axis.line        = element_line(colour = 'black', linewidth = 1),
                    axis.title.x     = element_text(size = 12),
                    axis.title.y     = element_text(size = 12),
                    axis.text.x      = element_text(size = 9),
                    axis.text.y      = element_text(size = 9)) +
                guides(colour = guide_legend(override.aes = list(size = 5)))
            if (isTRUE(plot.output)) print(p.batches)
        } else {
            data.to.plot <- as.data.frame(input.data)
            data.to.plot$batches <- factor(
                x      = paste0('Batch', as.numeric(as.factor(uv.sources))),
                levels = paste0('Batch', sort(unique(as.numeric(as.factor(uv.sources)))))
                )
            p.batches <- GGally::ggpairs(
                data       = data.to.plot[, 1:(ncol(data.to.plot) - 1)],
                mapping    = ggplot2::aes(colour = data.to.plot[, ncol(data.to.plot)]),
                showStrips = FALSE,
                switch     = 'y',
                labeller   = NULL,
                diag       = list(continuous = wrap("diagAxis", labelSize = 8, diagAxis = 0)),
                upper      = "blank") +
                theme(
                    panel.background    = element_blank(),
                    legend.key          = element_blank(),
                    legend.text         = element_text(size = 12),
                    legend.title        = element_text(size = 14),
                    panel.grid.major    = element_blank(),
                    axis.ticks          = element_blank(),
                    strip.background    = element_blank(),
                    strip.text.x.bottom = element_text(size = 0),
                    strip.text          = element_text(size = 0),
                    axis.line           = element_line(colour = 'black', linewidth = 1),
                    axis.title.x        = element_text(size = 12),
                    axis.title.y        = element_text(size = 2),
                    axis.text.x         = element_text(size = 0),
                    axis.text.y         = element_text(size = 0)) +
                scale_color_manual(values = colors.selected)
            if (isTRUE(plot.output)) print(p.batches)
        }
    }

    # Assessing association between the estimates batches and variable ####
    ## Biological variables ####
    if (isTRUE(assess.bio.association)){
        homo.bio.groups <- createHomogeneousBioGroups(
            se.obj            = se.obj,
            bio.variables     = bio.variables,
            nb.clusters       = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            check.se.obj      = FALSE,
            save.se.obj       = FALSE,
            remove.na         = 'none',
            verbose           = verbose
            )
        bio.association <- DescTools::CramerV(
            x = uv.sources,
            y = homo.bio.groups
            )
        if (isTRUE(generate.association.plot)){
            batches.bio <- table(
                uv.sources,
                homo.bio.groups
            )
            h.bio <- Heatmap(
                matrix = batches.bio,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                col = c('grey90', 'darkgreen'),
                name = 'Frequency',
                row_names_gp = gpar(fontsize = 18),
                column_names_gp = gpar(fontsize = 16),
                heatmap_legend_param = list(
                    labels_gp = gpar(fontsize = 22),
                    title_gp = gpar(fontsize = 22),
                    legend_direction = 'horizontal', legend_width = unit(7, "cm"))
            )
            if (isTRUE(plot.output)) print(h.bio)
        } else h.bio <- NULL
    }

    ## Unwanted variables ####
    if (isTRUE(assess.uv.association)){
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj            = se.obj,
            uv.variables      = uv.variables,
            nb.clusters       = nb.uv.clusters,
            clustering.method = uv.clustering.method,
            check.se.obj      = FALSE,
            save.se.obj       = FALSE,
            remove.na         = 'none',
            verbose           = verbose
            )
        uv.association <- DescTools::CramerV(
            x = uv.sources,
            y = homo.uv.groups
            )
        if (isTRUE(generate.association.plot)){
            batches <- table(
                uv.sources,
                homo.uv.groups
                )
            h.uv <- Heatmap(
                matrix               = batches,
                cluster_rows         = FALSE,
                cluster_columns      = FALSE,
                col                  = c('grey90', 'darkgreen'),
                name                 = 'Frequency',
                row_names_gp         = gpar(fontsize = 18),
                column_names_gp      = gpar(fontsize = 16),
                heatmap_legend_param = list(
                    labels_gp        = gpar(fontsize = 22),
                    title_gp         = gpar(fontsize = 22),
                    legend_direction = 'horizontal',
                    legend_width     = unit(7, "cm"))
                )
            if (isTRUE(plot.output)) print(h.uv)
        } else h.uv <- NULL
    }
    if (isTRUE(add.to.sample.annotation)){
        se.obj[[col.name]] <- paste0('Batch', as.numeric(as.factor(uv.sources)))
    }
    # Saving the results ####
    printColoredMessage(
        message = '- Saving the the results:',
        color = 'magenta',
        verbose = verbose
        )
    ## Selecting the out put name ####
    if (is.null(output.name)){
        input.data.name <- paste0(
            length(unique(uv.sources)),
            'batches|',
            input.data.name
            )
    } else input.data.name <- output.name
    ## Saving in the SummarizedExperiment object ####
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

        if (!'plots' %in%  names(se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]])){
            se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['plots']] <- list()
        }
        se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['plots']][['batch.plot']] <- p.batches
        if (isTRUE(assess.bio.association)){
            se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['plots']][['bio.plot']] <- h.bio
        }
        if (isTRUE(assess.uv.association)){
            se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['plots']][['uv.plot']] <- h.uv
        }

        if (!'variable.association' %in%  names(se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]])){
            se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['variable.association']] <- list()
        }
        if (isTRUE(assess.bio.association)){
            if (!'bio.association' %in% names(se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['variable.association']])){
                se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['variable.association']][['bio.association']] <- list()
            }
            se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['variable.association']][['bio.association']] <- bio.association
        }
        if (isTRUE(assess.uv.association)){
            if (!'uv.association' %in% names(se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['variable.association']])){
                se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['variable.association']][['uv.association']] <- list()
            }
            se.obj@metadata[['UnKnownUV']][[assay.name]][[input.data.name]][['variable.association']][['uv.association']] <- uv.association
        }
        printColoredMessage(
            message = 'The potentail unknow sources of variation are saved to the metadata of the SummarizedExperiment object',
            color   = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The indentifyUnknownUV function finished.',
            color   = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ## Saving as a list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = 'The results are outputed as list.',
            color   = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The indentifyUnknownUV function finished.',
            color   = 'white',
            verbose = verbose
            )
        return(list(
            batches    = uv.sources,
            input.data = input.data,
            plot       = p.batches)
            )
    }
}


