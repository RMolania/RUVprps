#' Finds NCGs using an unsupervised approach.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function identifies a set of genes to be used as negative control genes (NCGs) when no biological variation is
#' known. This function applies gene-level correlation, ANOVA, mean-variance relationship, and MAD analyses across and
#' between sample groups to select NCGs for different purposes, including RUV-III normalization.
#'
#' @details
#' The function initially uses gene-level correlation and ANOVA to identify genes significantly influenced by continuous
#' and categorical sources of variation, respectively. Then, it performs a mean-variance relationship or median absolute
#' deviation (MAD) analysis on each gene within homogeneous sample groups, considering unwanted variables, to detect genes
#' that are highly variable due to biological factors. Finally, various methods are applied to consolidate the statistical
#' findings and determine an appropriate set of genes as NCGs for different purposes,
#' including RUV-III normalization.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character that indicates the name of an data (assay) in the SummarizedExperiment object.
#' The selected assay should be the one that will be used for the RUV-III-PRPS normalization.
#' @param uv.variables Character. A character or a vector of characters indicating the name of the columns in the
#' SummarizedExperiment object that contain unwanted variables. These variables can be both categorical or continouse or
#' a combination of these. If all unwanted variation is unknown, use the  `identifyUnknownUV()` function to estimate them
#' first, then run `findNcgsUnSupervised()`.
#' @param clustering.method Character. The clustering method that should be used to group any continuous unwanted variable.
#' The options are: `kmeans`, `cut`, and `quantile`. The default is set to `kmeans`.
#' @param nb.clusters Numeric. A numeric value that indicates the number of clusters for grouping continuous unwanted
#' variation. The default is set to 3.
#' @param ncg.selection.method Character. A character that indicates how to summarize different statistics and select a
#' set of genes as negative control genes. The options are: `prod`, `average`, `sum`, `non.overlap`, `auto`, and `quantile`.
#' The default is set to `non.overlap`. For more information, refer to the details of the function.
#' @param nb.ncg Numeric. A numeric value that specifies the number of genes to be chosen as negative control genes (NCG)
#' when the `ncg.selection.method` parameter is set to `auto`. This value, `nb.ncg`, corresponds to a fraction of the total
#' genes in the SummarizedExperiment object. The default is set to 0.1.
#' @param hvg.method Character. A character indicating the method to select highly variable genes. The options are: `var`
#' and `mad`. `var` models gene variance based on a mean-variance trend. `mad` uses the median absolute deviation. The
#' default is set to `var`.
#' @param top.rank.bio.genes Numeric. A numeric value that indicates the percentage of top-ranked genes that are highly
#' affected by biological variation. This is required to be specified when the `ncg.selection.method` is either `auto`
#' or `non.overlap`. The default is set to 0.2.
#' @param top.rank.uv.genes Numeric. A numeric value that indicates the percentage of top-ranked genes that are highly
#' affected by unwanted variables. This is required to be specified when the `ncg.selection.method` is either `auto` or
#' `non.overlap`. The default is set to 0.2.
#' @param bio.percentile Numeric. A numeric value that specifies the percentile cut-off to select genes that are highly
#' affected by biological variation. This is required to be specified when the `ncg.selection.method` is set to `quantile`.
#' The default is set to 0.8.
#' @param uv.percentile Numeric. A numeric value that specifies the percentile cut-off to select genes that are highly
#' affected by unwanted variation. This is required to be specified when the `ncg.selection.method` is set to `quantile`.
#' The default is set to 0.8.
#' @param grid.group Character. A character that indicates whether the grid search should be performed on biological,
#' unwanted, or both factors when the `ncg.selection.method` is set to `auto`. The options are `bio`, `uv`, or `both`.
#' The default is set to `uv`. For more details, refer to the function documentation.
#' @param grid.direction Character. A character that indicates whether the grid search should be performed in decreasing
#' or increasing order when the `ncg.selection.method` is set to `auto`. The options are: `increase` and `decrease`. The
#' default is set to `decrease`.
#' @param grid.nb Numeric. A numeric value that indicates the number of genes for grid search when the `ncg.selection.method`
#' is set to `auto`. In the `auto` approach, the grid search starts with the initial `top.rank.uv.genes` and
#' `top.rank.bio.genes` values and adds or drops the `grid.nb` in each loop to find `nb.ncg` of genes as negative control
#' genes. The default is set to 20.
#' @param min.sample.for.mad Numeric. A numeric value that indicates minimum number of samples required per group for MAD
#' analysis. The default is set to 3.
#' @param min.sample.for.var Numeric. A numeric value that indicates minimum number of samples required per group for
#' gene-variance relationship analysis. The default is set to 15.
#' @param min.sample.for.aov Numeric. A numeric value that indicates minimum number of samples required per group  for
#' ANOVA. The default is set to 3.
#' @param min.sample.for.correlation Numeric. A numeric value that indicates minimum number of samples required for
#' correlation analysis. The default is set to 10
#' @param corr.method Character. A character that indicates which correlation methods should be used for the correlation
#' analyses. The options are `pearson` or `spearman`. The default is set to `spearman`.
#' @param a Numeric. The significance level used for the confidence intervals in the correlation; by default, it is set
#' to 0.05. Refer to the function `correls` from the **Rfast** R package for more details.
#' @param rho Numeric. The value of the hypothesized correlation to be used in the hypothesis testing. The default is
#' set to 0. Refer to the function `correls` from the **Rfast** R package for more details.
#' @param anova.method Character. A character that  indicates which ANOVA method to use. The options are `aov` or `welch`.
#' The default is se to `aov`. Refer to the function `row_oneway_equalvar()` or `row_oneway_welch()` from the
#' **matrixTests** R package for more details.
#' @param create.ncg.rank.plot Logical. Indicates whether to generate a heatmap that shows the rank of the all genes
#' with respect to their biological and unwanted variation effects. The default is set to `FALSE`.
#' @param plot.ncg.rank Logical. Indicates whether to plot a heatmap that shows the rank of the all genes
#' with respect to their biological and unwanted variation effects, while function is running. The default is set to `FALSE`.
#' @param filter.ncgs Logical. Whether to filter selected NCGs based on public human housekeeping gene sets. The default
#' is set to `FALSE`.
#' @param common.hk Character. Specifies group of housekeeping genes to use: `cancer` or `non.cancer`. The default is set
#' to `cancer`.
#' @param nb.stable.genes Numeric. A numeric value that specifies the number of top stable genes to be obtained from the
#' `getStableGenes()` function in the **singescore** R package. The default is set to 2000.
#' @param hk.group Character. Column name in the gene annotation containing non-cancer housekeeping genes. Options include:
#' `bulk.rnaseq.hk.genes.v1`, `bulk.rnaseq.hk.genes.v2`, `micorarray.hk.genes`, `nanostring.pan.cancer.hk.genes`,
#' `singscore.pan.cancer.hk.genes`. The default is set to `micorarray.hk.genes`.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected genes as negative control or not.
#' This analysis involves principal component analysis on the selected genes, followed by exploration of the R^2 or vector
#' correlation between the first `nb.pcs` principal components and the biological and unwanted variables. The default is
#' set to `TRUE`.
#' @param variables.to.assess.ncg Character. A character string or vector of strings indicating the column names in sample
#' annotation of of the  SummarizedExperiment object that contain variables whose association with the selected genes as
#' NCG needs to be evaluated. The default is set to `NULL`. This means all the variables specified in the `bio.variables`
#' and `uv.variables` will be assessed.
#' @param nb.pcs Numeric. A numeric value that indicates the number of the first principal components of selected negative
#' control genes to be used to assess their performance. The default is set to 10.
#' @param center Logical. Indicates whether to scale the data before applying SVD. If `TRUE`, centering is done by subtracting
#' the column means of the assay from their corresponding columns. The default is set to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data before applying SVD. If `TRUE`, scaling is done by dividing the
#' (centered) columns of the assays by their standard deviations if centering is `TRUE`, and by the root mean square otherwise.
#' The default is set to `FALSE`.
#' @param plot.ncg.assessment Logical. Indicates whether to plot the output of the NCG assessment while function is running
#' . The default is set to `TRUE`.
#' @param regress.out.variables Character. A character or a vector of character indicating the names of the columns in
#' the SummarizedExperiment object that contain variables to be regressed out from the data before identifying biologically
#' variable genes. The default is set to `NULL`.
#' @param normalization Character. A character that indicates which normalization method should be use to mitigate the
#' variation in library size before finding genes that are highly affected by biological variation. The options are :
#' `CPM`, `TMM`, `VST`, `upper`, `full` and `medium`. The default is set to  `CPM`. Refer to the `applyOtherNormalization()`
#' function for more details.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before any statistical analyses.
#' The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value to be added as a pseudo count to all measurements before applying log
#' transformation. The default is set to 1.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object before any analysis. If `TRUE`,
#'  the function `checkSeObj()` will be used. The default is set to `TRUE`.
#' @param remove.na Character set. Indicates whether to remove NA or missing values from the SummarizedExperiment object
#' The options are: `assays`, the `sample.annotation`, `both`, or `none`. If `assays` is selected, genes containing NA or
#' missing values will be excluded. If `sample.annotation` is selected, the samples containing NA or missing values for
#' any `bio.variables` or `uv.variables` will be excluded. The default is set to `none`.
#' @param ncg.group.name Character. A character to be used as name of the group of NCG. The default is set to `NULL`, then
#' the function create a names as following: `paste0('ncg|unsupervised')`. We refer to the details of the function for
#' more details.
#' @param ncg.set.name Character. A character to be used as name of the NCG set based on current variables and parameters
#' The default is set to `NULL`, then the function create a names as following:
#' `paste0(sum(ncg.selected),'|',paste0(bio.variables, collapse = '&'),'|',paste0(uv.variables, collapse = '&'),'|AnoCorrAs:',
#' ncg.selection.method,'|',assay.name)`.We refer to the details of the function for more details.
#' @param save.imf Logical. Indicates whether to save the intermediate file. If `TRUE`, the function saves the results
#' of the statistical analyses in the metadata of the SummarizedExperiment object. If users want to change the parameters
#' including `nb.ncg`, `ncg.selection.method`, `top.rank.bio.genes`, and `top.rank.uv.genes`, the analyses will not be
#' re-calculated. The default is set to `FALSE`.
#' @param use.imf Logical. Indicates whether to use the intermediate file. The default is set to `FALSE`.
#' @param imf.name Character string. A name to save the intermediate file. If `NULL`, the function generates a name.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the
#' SummarizedExperiment object or output the result. The default is `TRUE`.
#' @param verbose Logical. If `TRUE`, shows messages of different steps of the function.
#' @param use.rank TTTT
#' @param samples.to.use TTTT
#' @param regress.out.rle.med TTT
#'
#' @return A `SummarizedExperiment` object containing the selected negative control genes and optional assessment plots,
#' or a list of the results.
#'
#' @importFrom matrixTests row_oneway_equalvar row_oneway_welch
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom BiocSingular bsparam runSVD
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr progress_estimated
#' @importFrom fastDummies dummy_cols
#' @importFrom Rfast correls rowcvs
#' @importFrom matrixStats rowMads
#' @importFrom tidyr pivot_longer
#' @importFrom scran modelGeneVar
#' @importFrom ruv design.matrix
#' @importFrom Biobase rowMax
#' @importFrom stats aov
#' @import ggplot2
#' @export

findNcgsUnSupervisedByAnovaCorr <- function(
        se.obj,
        assay.name,
        uv.variables,
        ncg.selection.method = 'quantile',
        use.rank = FALSE,
        samples.to.use = 'all',
        nb.ncg = 0.1,
        hvg.method = 'mad',
        top.rank.bio.genes = 0.8,
        top.rank.uv.genes = 0.2,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'decrease',
        grid.nb = 40,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        min.sample.for.mad = 3,
        min.sample.for.var = 15,
        min.sample.for.aov = 3,
        min.sample.for.correlation = 10,
        corr.method = "spearman",
        a = 0.05,
        rho = 0,
        anova.method = 'aov',
        create.ncg.rank.plot = FALSE,
        plot.ncg.rank = FALSE,
        filter.ncgs = FALSE,
        common.hk = 'cancer',
        nb.stable.genes = 2000,
        hk.group = 'micorarray.hk.genes',
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        plot.ncg.assessment = TRUE,
        regress.out.variables = NULL,
        regress.out.rle.med = FALSE,
        normalization = 'CPM',
        apply.log = TRUE,
        pseudo.count = 1,
        check.se.obj = TRUE,
        remove.na = 'both',
        ncg.group.name = NULL,
        ncg.set.name = NULL,
        save.imf = FALSE,
        use.imf = FALSE,
        imf.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findNcgsUnSupervisedByAnovaCorr function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the  functions inputs ####
    if (length(assay.name) > 1 | is.logical(assay.name)){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    }
    if (sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables)){
        stop('Some or all the "uv.variables" cannot be found in the SummarizedExperiment object.')
    }
    if (nb.ncg >= 1 | nb.ncg <= 0){
        stop('The "nb.ncg" should be a positve value  0 < nb.ncg < 1.')
    }
    if (!ncg.selection.method %in% c('auto', 'non.overlap', 'quantile', 'prod', 'sum', 'average')){
        stop('The "ncg.selection.method" must be one of "auto", "quantile", "prod", "sum", "average", "non.overlap".')
    }
    if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
        stop('The "top.rank.bio.genes" msut be a positve value  0 < top.rank.bio.genes < 1.')
    }
    if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
        stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes < 1.')
    }
    if (grid.nb < 1 | grid.nb > nrow(se.obj)){
        stop(paste0('The "grid.nb" must be a positve value  0 < grid.nb < ', nrow(se.obj), '.'))
    }
    if (!grid.group %in% c('bio', 'uv', 'both')){
        stop('The "grid.group" must be one of "bio", "uv" or "non.overlap".')
    }
    if (!grid.direction %in% c('increase', 'decrease', 'auto')){
        stop('The "grid.direction" must be one of "increase", "decrease" or "auto".')
    }
    if (is.null(min.sample.for.aov)){
        stop('The "min.sample.for.aov" cannot be empty.')
    }
    if (min.sample.for.aov <= 2){
        stop('The "min.sample.for.aov" should be at least 3.')
    }
    if (is.null(min.sample.for.correlation)){
        stop('The min.sample.for.correlation cannot be empty.')
    }
    if (min.sample.for.correlation >= ncol(se.obj) | min.sample.for.correlation < 3){
        stop('The "min.sample.for.correlation" msut be more than 2 and less than the total number of samples in the data.')
    }
    if (!anova.method %in% c('aov', 'welch')){
        stop('The anova.method must be one of the "aov" or "welch".')
    }
    if (isFALSE(is.logical(assess.ncg))){
        stop('The "assess.ncg" must be "TRUE" or "FALSE.')
    }
    if (length(nb.pcs) > 1){
        stop('The "nb.pcs" must be a postive integer value.')
    }
    if (nb.pcs < 0){
        stop('The "nb.pcs" must be a postive integer value.')
    }
    if (isFALSE(is.logical(scale))) {
        stop('The "scale" must be "TRUE" or "FALSE.')
    }
    if (isFALSE(is.logical(center))) {
        stop('The "center" must be "TRUE" or "FALSE.')
    }
    if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE.')
    }
    if (length(pseudo.count) > 1){
        stop('The "pseudo.count" must be 0 or a postive integer value.')
    }
    if (pseudo.count < 0){
        stop('The "pseudo.count" must be 0 or a postive integer value.')
    }
    if (isFALSE(is.logical(check.se.obj))) {
        stop('The "check.se.obj" must be "TRUE" or "FALSE.')
    }
    if (is.null(check.se.obj)) {
        if (isTRUE(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables))) {
            stop('All or some of "uv.variables" cannot be found in the SummarizedExperiment object.')
        } else if (!is.null(variables.to.assess.ncg)) {
            if (isTRUE(sum(variables.to.assess.ncg %in% colnames(colData(se.obj))) != length(variables.to.assess.ncg))) {
                stop('All or some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
            }
        }
    }
    if (!is.null(regress.out.variables)){
        if (isTRUE(sum(regress.out.variables %in% colnames(colData(se.obj))) != length(regress.out.variables))) {
            stop('All or some of "regress.out.variables" cannot be found in the SummarizedExperiment object.')
        }
    }

    if (isTRUE(filter.ncgs)){
        if (!is.character(common.hk) | length(common.hk) > 1){
            stop('The "common.hk" must be one of the "cancer" or "non.cancer".')
        }
        if (!common.hk %in% c('cancer', 'non.cancer')){
            stop('The "common.hk" must be one of the "cancer" or "non.cancer".')
        }
        if (!is.numeric(nb.stable.genes) | length(nb.stable.genes) > 1 ){
            stop('The "nb.stable.genes" must be a numeric postive value.')
        }
        if (nb.stable.genes < 0){
            stop('The "nb.stable.genes" must be a numeric postive value.')
        }
        if (!is.character(hk.group) | length(hk.group) > 1){
            stop('The "hk.group" must be one of the "bulk.rnaseq.hk.genes.v1", "bulk.rnaseq.hk.genes.v2", "micorarray.hk.genes", "nanostring.pan.cancer.hk.genes", "singscore.pan.cancer.hk.genes".')
        }
        if (!hk.group %in% c("bulk.rnaseq.hk.genes.v1", "bulk.rnaseq.hk.genes.v2", "micorarray.hk.genes", "nanostring.pan.cancer.hk.genes", "singscore.pan.cancer.hk.genes")){
            stop('The "hk.group" must be one of the "bulk.rnaseq.hk.genes.v1", "bulk.rnaseq.hk.genes.v2", "micorarray.hk.genes", "nanostring.pan.cancer.hk.genes", "singscore.pan.cancer.hk.genes".')
        }
    }
    if (isTRUE(use.rank)){
        if (isTRUE(ncg.selection.method == 'quantile')){
            if (is.null(bio.percentile) | is.null(uv.percentile))
                stop('The "bio.percentile" or "uv.percentile" cannot be NULL.')
            if (bio.percentile > 1 | bio.percentile < 0)
                stop('The "bio.percentile" must be a postive value between 0 and 1.')
            if (uv.percentile > 1 | uv.percentile < 0)
                stop('The "uv.percentile" must be a postive value between 0 and 1.')
        }
    }

    if (is.logical(samples.to.use)){
        if (length(samples.to.use) != ncol(se.obj)){
            stop('The length of the "samples.to.use" must be the same as the number of columns in the SummarizedExperiment object.')
        }
        se.obj.all <- se.obj
        se.obj <- se.obj[ , samples.to.use]
    }

    # Checking the SummarizedExperiment object ####
    if (isTRUE(check.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = unique(c(uv.variables, regress.out.variables, variables.to.assess.ncg)),
            remove.na = remove.na,
            verbose = verbose
            )
    }
    # Applying data normalization and transformation ####
    printColoredMessage(
        message = '-- Applying data normalization and transformation:',
        color = 'magenta',
        verbose = verbose
        )
    ## Applying library size normalization ####
    if (!is.null(normalization)){
        expr.data.nor <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            pseudo.count = pseudo.count,
            apply.log = apply.log,
            check.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
    }
    ## Applying log transformation ####
    if (isTRUE(apply.log) & !is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                '- Applying log2 + ',
                pseudo.count,
                ' (pseudo.count) on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose
            )
        expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
    } else if (isTRUE(apply.log) & is.null(pseudo.count)){
        printColoredMessage(
            message = paste0(
                '- Applying log2 on the ',
                assay.name,
                ' data.'),
            color = 'blue',
            verbose = verbose
            )
        expr.data <- log2(assay(x = se.obj, i = assay.name))
    } else if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                'The ',
                assay.name,
                ' data will be used without any transformation.'),
            color = 'blue',
            verbose = verbose
            )
        expr.data <- assay(x = se.obj, i = assay.name)
    }

    # Finding negative control genes ####
    if (isFALSE(use.imf)){
        printColoredMessage(
            message = '-- Finding a subset of genes as negative control genes (NCGs):',
            color = 'magenta',
            verbose = verbose
            )
        ## Identifying genes that are highly affected by unwanted variation ####
        printColoredMessage(
            message = '-- Finding genes that are highly affected by each specified source(s) of unwanted variation:',
            color = 'orange',
            verbose = verbose
            )
        ### Finding the classes of variables ####
        uv.var.class <- unlist(lapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]]))
            )
        categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
        continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]

        ### Performing ANOVA between genes and categorical sources of unwanted variation ####
        if (isTRUE(length(categorical.uv) > 0)){
            printColoredMessage(
                message = paste0(
                    '- Performing ANOVA between individual gene-level ',
                    'expression and each categorical source of unwanted variation: ',
                    paste0(categorical.uv, collapse = ' & '),
                    '.'),
                color = 'blue',
                verbose = verbose
                )
            anova.genes.uv <- lapply(
                categorical.uv,
                function(x) {
                    keep.samples <- findRepeatingPatterns(
                        vec = colData(se.obj)[[x]],
                        n.repeat = min.sample.for.aov
                        )
                    if (isTRUE(length(keep.samples) == 0)){
                        stop(paste0(
                            'There are not enough samples to perform ANOVA between individual gene expression and the ',
                            x,
                            ' variable. Possible solutions is to lower min.sample.for.aov or remove the "',
                            x,
                            '" from the uv.variables and re-run the function.')
                            )
                    } else if (isTRUE(length(keep.samples) == 1)){
                        stop(paste0(
                            'There is only a single batch from in the ',
                            x,
                            ' variable that have enough samples ',
                            min.sample.for.aov,
                            ' (min.sample.for.aov). Possible solutions is to lower min.sample.for.aov or remove the "',
                            x,
                            '" from the uv.variables and re-run the function.')
                            )
                    } else if (isTRUE(length(keep.samples) != length(unique(colData(se.obj)[[x]]))) ){
                        not.coverd <- unique(colData(se.obj)[[x]])[!unique(colData(se.obj)[[x]]) %in% keep.samples]
                        printColoredMessage(
                            message = paste0(
                                '- Note, the ',
                                paste0(not.coverd, collapse = '&'),
                                ' batches do not have enough samples for the ANOVA analysis.'),
                            color = 'red',
                            verbose = verbose)
                    }
                    keep.samples <- colData(se.obj)[[x]] %in% keep.samples
                    if (anova.method == 'aov'){
                        anova.gene.batch <- as.data.frame(row_oneway_equalvar(
                            x = expr.data[ , keep.samples],
                            g = se.obj@colData[[x]][keep.samples])
                            )
                    } else if (anova.method == 'welch'){
                        anova.gene.batch <- as.data.frame(row_oneway_welch(
                            x = expr.data[ , keep.samples],
                            g = se.obj@colData[[x]][keep.samples])
                            )
                    }
                    #### Calculating adjusted p-value and effect size ####
                    anova.gene.batch$p.adjusted <- p.adjust(anova.gene.batch$pvalue, method = "BH")
                    anova.gene.batch$eta.squared <- (anova.gene.batch$statistic * anova.gene.batch$df.between) /
                        (anova.gene.batch$statistic * anova.gene.batch$df.between + anova.gene.batch$df.within)
                    anova.gene.batch$g.statistic <- anova.gene.batch$eta.squared

                    #### Ranking the F-statistics ####
                    set.seed(2233)
                    anova.gene.batch$ranked.genes <- rank(
                        x = -anova.gene.batch[ , 'statistic'],
                        ties.method = 'random'
                        )
                    anova.gene.batch
                })
            names(anova.genes.uv) <- categorical.uv
        } else anova.genes.uv <- NULL

        ### Applying correlation between genes and continuous sources of unwanted variation ####
        if (isTRUE(length(continuous.uv) > 0)){
            printColoredMessage(
                message = paste0(
                    '- Performing ',
                    corr.method,
                    ' correlation between individual gene-level ',
                    'expression and each continuous source of unwanted variations: ',
                    paste0(continuous.uv, collapse = '&'),
                    '.'),
                color = 'blue',
                verbose = verbose
                )
            if (isTRUE(ncol(se.obj) <= min.sample.for.correlation)){
                stop(paste0(
                    'There are not enough samples (min.sample.for.correlation:',
                    min.sample.for.correlation,
                    ') to perform correlation analysis.',
                    ' A possible solution is to lower the "min.sample.for.correlation" value.'))
            }
            corr.genes.uv <- lapply(
                continuous.uv,
                function(x) {
                    corr.genes.var <- as.data.frame(correls(
                        y = se.obj@colData[, x],
                        x = t(expr.data),
                        type = corr.method,
                        a = a ,
                        rho = rho)
                        )
                    corr.genes.var <- cbind(
                        round(x = corr.genes.var[, 1:4], digits = 3),
                        corr.genes.var[, 5, drop = FALSE]
                        )
                    #### Calculating adjusted p-value and making obs of correlation ####
                    corr.genes.var$p.adjusted <- p.adjust(corr.genes.var[['p-value']], method = "BH")
                    corr.genes.var$g.statistic <- abs(corr.genes.var$correlation)

                    #### Ranking the correlation coefficent ####
                    set.seed(2233)
                    colnames(corr.genes.var)[colnames(corr.genes.var) == 'correlation' ] <- 'statistic'
                    corr.genes.var$ranked.genes <- rank(
                        x = -abs(corr.genes.var[, 'statistic']),
                        ties.method = 'random'
                        )
                    row.names(corr.genes.var) <- row.names(expr.data)
                    corr.genes.var
                })
            names(corr.genes.uv) <- continuous.uv
        } else corr.genes.uv <- NULL

        ## Identifying genes that are highly affected by possible biological variation ####
        printColoredMessage(
            message = '-- Finding genes that are potentially highly affected by biological variation:',
            color = 'orange',
            verbose = verbose
            )
        #### Selecting the data input ####
        if (!is.null(normalization)) {
            data.to.use <- expr.data.nor
        } else if (is.null(normalization)) data.to.use <- expr.data

        #### Regressing out variables ####
        if (!is.null(regress.out.variables)){
            printColoredMessage(
                message = paste0(
                    '- The ',
                    paste0(regress.out.variables, collapse = ' & '),
                    ' variable(s) will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose
                )
            data.to.use <- t(data.to.use)
            lm.formula <- paste('se.obj', regress.out.variables, sep = '$')
            adjusted.data <- lm(as.formula(
                paste('data.to.use', paste0(lm.formula, collapse = '+') , sep = '~'))
                )
            data.to.use <- t(adjusted.data$residuals)
            colnames(data.to.use) <- colnames(se.obj)
            row.names(data.to.use) <- row.names(se.obj)
        }
        ## Applying different statistical test within each homogeneous sample groups with respect to the unwanted variable ####
        ##### Finding all possible sample groups with respect to the unwanted variables ####
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = uv.variables,
            clustering.method =  clustering.method,
            nb.clusters = nb.clusters,
            check.se.obj = FALSE,
            save.se.obj = FALSE,
            verbose = verbose
            )
        ##### Applying MAD ####
        if (hvg.method == 'mad'){
            printColoredMessage(
                message = paste0(
                    '- Performing MAD on individual gene expression',
                    ' within each homogeneous sample groups with respect to the unwanted variables.'),
                color = 'blue',
                verbose = verbose
                )
            selected.homo.uv.groups <- findRepeatingPatterns(
                vec = homo.uv.groups,
                n.repeat = min.sample.for.mad
                )
            if (isTRUE(length(selected.homo.uv.groups) > 0)){
                bio.genes <- sapply(
                    selected.homo.uv.groups,
                    function(x){
                        index.samples <- homo.uv.groups == x
                        if (isTRUE(regress.out.rle.med)){
                            temp.data <- data.to.use[ , index.samples, drop = FALSE]
                            rle.med <- colMedians(temp.data - rowMedians(temp.data))
                            regress.med <- lm(t(temp.data) ~ rle.med)
                            temp.data <- t(regress.med$residuals)
                            colnames(temp.data) <- colnames(se.obj)[index.samples]
                            row.names(temp.data) <- row.names(se.obj)
                            temp.data <- temp.data - min(temp.data)
                            rowMads(temp.data)
                        } else rowMads(x = data.to.use[ , index.samples, drop = FALSE])
                    })
                bio.genes <- matrixStats::rowMedians(bio.genes)
                bio.genes <- data.frame(bio.stat = bio.genes)
                bio.genes$bio.stat.scaled <- (bio.genes$bio.stat - min(bio.genes$bio.stat, na.rm = TRUE)) /
                    (max(bio.genes$bio.stat, na.rm = TRUE) - min(bio.genes$bio.stat, na.rm = TRUE))
                set.seed(3322)
                bio.genes$bio.ranks <- rank(x = bio.genes$bio.stat, ties.method = 'random')
                row.names(bio.genes) <- row.names(se.obj)
            } else{
                stop(paste0(
                    'There is no any homogenous sample groups with at least ',
                    min.sample.for.mad,
                    ' samples to perform MAD.')
                )
            }
        }
        ##### Applying  MAD and CV ####
        if (hvg.method == 'mad.cv'){
            printColoredMessage(
                message = paste0(
                    '- Performing MAD and CV on individual gene expression',
                    ' within each homogeneous sample groups with respect to the unwanted variables.'),
                color = 'blue',
                verbose = verbose
                )
            selected.homo.uv.groups <- findRepeatingPatterns(
                vec = homo.uv.groups,
                n.repeat = min.sample.for.mad
                )
            if (isTRUE(length(selected.homo.uv.groups) > 0)){
                bio.genes <- lapply(
                    selected.homo.uv.groups,
                    function(x){
                        index.samples <- homo.uv.groups == x
                        if (isTRUE(regress.out.rle.med)){
                            temp.data <- data.to.use[ , index.samples, drop = FALSE]
                            rle.med <- colMedians(temp.data - rowMedians(temp.data))
                            regress.med <- lm(t(temp.data) ~ rle.med)
                            temp.data <- t(regress.med$residuals)
                            colnames(temp.data) <- colnames(se.obj)[index.samples]
                            row.names(temp.data) <- row.names(se.obj)
                            temp.data <- temp.data - min(temp.data)
                            data.frame(cv = rowcvs(temp.data), mad = rowMads(temp.data))
                        } else {
                            data.frame(
                                cv = rowcvs(x = data.to.use[ , index.samples, drop = FALSE]),
                                mad = rowMads(x = data.to.use[ , index.samples, drop = FALSE])
                                )}
                    })
                bio.genes <- as.matrix(do.call(cbind, bio.genes))
                bio.genes.cv <- rowMedians(bio.genes[ , grep('cv', colnames(bio.genes))])
                bio.genes.mad <- rowMedians(bio.genes[ , grep('mad', colnames(bio.genes))])
                set.seed(3322)
                bio.genes.cv.ranked <- rank(
                    x = rowMedians(bio.genes[ , grep('cv', colnames(bio.genes))]),
                    ties.method = 'random'
                    )
                bio.genes.mad.ranked <- rank(
                    x = rowMedians(bio.genes[ , grep('mad', colnames(bio.genes))]),
                    ties.method = 'random'
                    )
                bio.genes <- data.frame(
                    bio.genes.cv = bio.genes.cv,
                    bio.genes.cv.scaled = (bio.genes.cv - min(bio.genes.cv, na.rm = TRUE)) /
                        (max(bio.genes.cv, na.rm = TRUE) - min(bio.genes.cv, na.rm = TRUE)),
                    bio.genes.mad = bio.genes.mad,
                    bio.genes.mad.scaled = (bio.genes.mad - min(bio.genes.mad, na.rm = TRUE)) /
                        (max(bio.genes.mad, na.rm = TRUE) - min(bio.genes.mad, na.rm = TRUE)),
                    bio.genes.cv.ranked = bio.genes.cv.ranked,
                    bio.genes.mad.ranked = bio.genes.mad.ranked
                    )
                bio.genes$bio.stat = pmax(bio.genes$bio.genes.mad.scaled, bio.genes$bio.genes.cv.scaled)
                bio.genes$gene.pro <- bio.genes$bio.genes.cv.ranked * bio.genes$bio.genes.mad.ranked
                set.seed(3322)
                bio.genes$bio.ranks <- rank(
                    x = bio.genes$gene.pro,
                    ties.method = 'random'
                    )
                row.names(bio.genes) <- row.names(se.obj)
            } else{
                stop(paste0(
                    'There is no any homogenous sample groups with at least ',
                    min.sample.for.mad,
                    ' samples to perform MAD.')
                )
            }
        }
        # ##### Applying gene-variance relationship analysis ####
        # if (hvg.method == 'var'){
        #     printColoredMessage(
        #         message = paste0(
        #             '- Performing gene-variance relationship analysis on individual gene expression',
        #             ' within each homogeneous sample groups with respect to the unwanted variables.'),
        #         color = 'blue',
        #         verbose = verbose
        #         )
        #     selected.homo.uv.groups <- findRepeatingPatterns(
        #         vec = homo.uv.groups,
        #         n.repeat = min.sample.for.var
        #         )
        #     if (isTRUE(length(selected.homo.uv.groups) > 0)){
        #         selected.samples <- homo.uv.groups %in% selected.homo.uv.groups
        #         batch.design <- design.matrix(a = homo.uv.groups[selected.samples])
        #         bio.genes <- modelGeneVar(x = data.to.use[ , selected.samples], design = batch.design)
        #         bio.genes$bio.ranks <- rank(x = bio.genes$bio, ties.method = 'random')
        #     } else {
        #         stop(paste0(
        #             'There is no any homogenous sample groups with at least ',
        #             min.sample.for.var,
        #             ' samples to perform mean-variacne relationship analysis.')
        #         )
        #     }
        # }
    }
    # Checking and reading the intermediate file ####
    ## read intermediate file ####
    if (isTRUE(use.imf)){
        if (is.null(imf.name)){
            imf.name <- paste0(assay.name, '|un.supervised|', ncg.selection.method)
        }
        if (is.null(se.obj@metadata$IMF$NCG[[imf.name]]))
            stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
        all.tests <- se.obj@metadata$IMF$NCG[[imf.name]]
        bio.genes <- all.tests$bio.genes
        anova.genes.uv <- all.tests$anova.genes.uv
        corr.genes.uv <- all.tests$corr.genes.uv
    }

    ## save intermediate file ####
    if (isTRUE(save.imf)){
        if (length(se.obj@metadata$IMF) == 0 ) {
            se.obj@metadata[['IMF']] <- list()
        }
        if (!'NCG' %in% names(se.obj@metadata[['IMF']])){
            se.obj@metadata[['IMF']][['NCG']] <- list()
        }
        if (is.null(imf.name)){
            imf.name <- paste0(assay.name, '|un.supervised|', ncg.selection.method)
        }
        if (!imf.name %in% names(se.obj@metadata[['IMF']][['NCG']])){
            se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list()
        }
        se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list(
            bio.genes = bio.genes,
            anova.genes.uv = anova.genes.uv,
            corr.genes.uv = corr.genes.uv)
        }
    # Summarizing the statistical results to select NCGs ####
    printColoredMessage(
        message = '-- Summarizing the statistical results to select a set of genes as NCGs:',
        color = 'magenta',
        verbose = verbose
        )
    ## Ratio ####
    if (isFALSE(use.rank)){
        if (!is.null(anova.genes.uv) & is.null(corr.genes.uv)){
            g.statistic <- sapply(
                names(anova.genes.uv),
                function(x){
                    anova.genes.uv[[x]]$g.statistic
                })
            g.statistic <- rowMax(g.statistic)
            all.stats <- data.frame(
                genes = row.names(se.obj),
                bio = bio.genes$bio.stat,
                uv = g.statistic
                )
            all.stats$uv.bio <- all.stats$uv /c(bio + 10e-6)
        }
        if (is.null(anova.genes.uv) & !is.null(corr.genes.uv)){
            g.statistic <- sapply(
                names(corr.genes.uv),
                function(x){
                    corr.genes.uv[[x]]$g.statistic
                })
            g.statistic <- rowMax(g.statistic)
            all.stats <- data.frame(
                genes = row.names(se.obj),
                bio = bio.genes$bio.stat,
                uv = corr.genes.uv$g.statistic
                )
            all.stats$uv.bio <- all.stats$uv /c(bio + 10e-6)
        }
        if (!is.null(anova.genes.uv) & !is.null(corr.genes.uv)){
            g.statistic.corr <- sapply(
                names(corr.genes.uv),
                function(x){
                    corr.genes.uv[[x]]$g.statistic
                })
            g.statistic.corr <- rowMax(g.statistic.corr)
            g.statistic.aov <- sapply(
                names(anova.genes.uv),
                function(x){
                    anova.genes.uv[[x]]$g.statistic
                })
            g.statistic.aov <- rowMax(g.statistic.aov)

            all.stats <- data.frame(
                gene = row.names(se.obj),
                bio = bio.genes$bio.stat,
                uv = pmax(g.statistic.aov, g.statistic.corr)
                )
            all.stats$tech.var <- all.stats$uv /c(all.stats$bio + 10e-6)
        }
        if (ncg.selection.method == 'quantile'){
            if (!is.null(uv.percentile) & !is.null(bio.percentile)){
                all.stats <- all.stats[all.stats$tech.var > quantile(x = all.stats$tech.var, probs = uv.percentile) , ]
                if (quantile(x = all.stats$bio, probs = bio.percentile) == 0){
                    all.stats <- all.stats
                } else {
                    all.stats <- all.stats[all.stats$bio < quantile(x = all.stats$bio, probs = bio.percentile) , ]
                }
                nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                ncg.selected <- row.names(se.obj) %in% all.stats$gene[1:nb.ncg]
            } else if (!is.null(uv.percentile) & is.null(bio.percentile)){
                all.stats <- all.stats[all.stats$tech.var > quantile(x = all.stats$tech.var, probs = uv.percentile) , ]
                nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                ncg.selected <- row.names(se.obj) %in% all.stats$gene[1:nb.ncg]
            } else if (is.null(uv.percentile) & !is.null(bio.percentile)){
                all.stats <- all.stats[all.stats$bio < quantile(x = all.stats$bio, probs = bio.percentile) , ]
                nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                ncg.selected <- row.names(se.obj) %in% all.stats$gene[1:nb.ncg]
            } else if (is.null(uv.percentile) & is.null(bio.percentile)){
                nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                ncg.selected <- row.names(se.obj) %in% all.stats$gene[1:nb.ncg]
            }
        }
    }
    ## Product, sum or average of ranks ####
    if (ncg.selection.method %in% c('prod', 'sum', 'average')) {
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        all.uv.ranks <- lapply(
            all.uv.tests,
            function(x){
                temp <- get(x)
                if (length(names(temp))!=0){
                    ranks.data <- lapply(
                        names(temp),
                        function(y) temp[[y]]$ranked.genes)
                    ranks.data <- do.call(cbind, ranks.data)
                    colnames(ranks.data) <- names(temp)
                    ranks.data}
            })
        all.uv.ranks <- do.call(cbind, all.uv.ranks)
        row.names(all.uv.ranks) <- row.names(se.obj)
        all.ranks <- cbind(all.uv.ranks, bio.genes[ , 'bio.ranks', drop = FALSE])
        ### Applying the product of ranks ####
        if (ncg.selection.method == 'prod'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the product of ranks.',
                color = 'blue',
                verbose = verbose
                )
            stat.summary <- rowProds(as.matrix(all.ranks))
            if (sum(is.infinite(stat.summary)) > 0)
                stop('The product of ranks results in infinity values.')
        }
        ## Applying the average of ranks ####
        if (ncg.selection.method == 'sum'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the sum of ranks.',
                color = 'blue',
                verbose = verbose
                )
            stat.summary <- rowSums(as.matrix(all.ranks))
        }
        ## Applying the sum of ranks ####
        if (ncg.selection.method == 'average'){
            printColoredMessage(
                message = '- A set of NCG will be selected based on the average of ranks.',
                color = 'blue',
                verbose = verbose
                )
            stat.summary <- rowMeans(as.matrix(all.ranks))
        }
        ## Selecting top genes as NCGS ####
        all.ranks$stat.summary <- stat.summary
        set.seed(112233)
        all.ranks$rank.stat.summary <- rank(x = all.ranks$stat.summary, ties.method = 'random')
        all.ranks <- all.ranks[order(all.ranks$rank.stat.summary, decreasing = FALSE) , ]
        ncg.selected <- row.names(all.ranks[1:round(nb.ncg* nrow(se.obj), digits = 0) , ])
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    ## Non-overlap approach ####
    if (ncg.selection.method == 'non.overlap' & isTRUE(use.imf)){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs using the "non.overlap" approach.',
            color = 'orange',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '-1: selecting top ',
                top.rank.bio.genes *100,
                '% of highly affected genes by possible bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
        )
        ### Selecting genes affected by biological variation ####
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
        printColoredMessage(
            message = paste0(
                '-- ',
                length(top.bio.genes),
                ' genes are selected.'),
            color = 'blue',
            verbose = verbose
        )

        ## Selecting genes affected by unwanted variation ####
        printColoredMessage(
            message = paste0(
                '-2: selecting top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by each unwanted variation.'),
            color = 'blue',
            verbose = verbose
        )
        top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if (!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                            row.names(temp.data[[y]])[index] })
                        ))
                }
            })))
        printColoredMessage(
            message = paste0(
                '-- ',
                length(top.uv.genes),
                ' genes are selected.'),
            color = 'blue',
            verbose = verbose
        )
        ## Selecting of NCGS ####
        printColoredMessage(
            message = '- all genes found in 1 will be excluded from ones found in 2.',
            color = 'blue',
            verbose = verbose
        )
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if (isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    ## Quantile approach ####
    if (isTRUE(ncg.selection.method == 'quantile') & isTRUE(use.imf)){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "quantile" approach..',
            color = 'blue',
            verbose = verbose)
        ### find biological percentile ####
        bio.quan <- quantile(x = bio.genes$bio.mad, probs = bio.percentile)
        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.mad > bio.quan]

        ## find UV percentile ####
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if (!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$statistic < quantile(x = temp.data[[y]]$statistic , probs = uv.percentile)
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        printColoredMessage(
            message = paste0(
                '- Selecting ',
                length(top.uv.genes),
                ' genes with the UV F-statistics higher than ',
                ' (' , uv.percentile* 100,
                '% percentile), and exclude any genes presents in ',
                length(top.bio.genes),
                ' genes with the biological F-statistics higher than ',
                ' (' ,
                bio.percentile* 100,
                '% percentile).'),
            color = 'blue',
            verbose = verbose)
        top.uv.genes <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if (isTRUE(length(top.uv.genes) == 0)) stop('No NCGs can be found based on the current parameters.')
        ncg.selected <- row.names(se.obj) %in% top.uv.genes
    }

    ## Auto approach ####
    if (ncg.selection.method == 'auto' & isTRUE(use.imf)){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs using the "auto" approach.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- Selecting top ',
                top.rank.uv.genes * 100,
                '% of the most affected genes by unwanted variation, and then exclude any genes that are also in the top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
            )
        ### select genes affected by biological variation ####
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]

        ## select genes affected by unwanted variation ####
        top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        top.uv.genes <- unique(unlist(lapply(
            all.uv.tests,
            function(x){
                if (!is.null(x)){
                    temp.data <- get(x)
                    ranks.data <- unique(unlist(lapply(
                        names(temp.data),
                        function(y){
                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                            row.names(temp.data[[y]])[index] })))
                }
            })))
        ## select NCG ####
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        # if (isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
        printColoredMessage(
            message = paste0('- ', length(ncg.selected), ' genes are found.'),
            color = 'blue',
            verbose = verbose
            )
        ## assess the need for grid search ####
        nb.ncg <- round(c(nb.ncg * nrow(se.obj)), digits = 0)
        ncg.ranges <- round(x = 0.01 *nb.ncg, digits = 0)
        if (length(ncg.selected) > c(nb.ncg + ncg.ranges) | length(ncg.selected) < c(nb.ncg - ncg.ranges)) {
            if (isTRUE(nb.ncg > length(ncg.selected))){
                con <- parse(text = paste0("nb.ncg", ">", "length(ncg.selected)"))
                printColoredMessage(
                    message = paste0(
                        '- The number of selected genes ',
                        length(ncg.selected),
                        ' is less than the number (',
                        nb.ncg ,
                        ') of specified genes ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose)
            }
            if (isTRUE(nb.ncg < length(ncg.selected))){
                con <- parse(text = paste0("length(ncg.selected)", ">", "nb.ncg"))
                printColoredMessage(
                    message = paste0(
                        '- The number of selected genes ', length(ncg.selected),
                        ' is larger than the number (', nb.ncg ,') of specified genes ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose)
            }
            ## grid search ####
            ### grid group: both bio and uv variable ####
            if (grid.group == 'both'){
                printColoredMessage(
                    message = '- The grid search will be applied on both biological and unwanted factors. ',
                    color = 'blue',
                    verbose = verbose
                    )
                #### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the values of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- min(
                        nrow(se.obj) - top.rank.uv.genes.nb,
                        top.rank.bio.genes.nb
                        )
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if (top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x){
                                if (!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if (top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- min(top.rank.uv.genes.nb, c(nrow(se.obj) - top.rank.bio.genes.nb))
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb > 1 & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if (top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x){
                                if (!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if (top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ### check selection ####
                if (length(ncg.selected) == 0)
                    stop('No NCGs can be found based on the current parameters.')
                ### update numbers ####
                # bio
                top.rank.bio.genes.nb <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes.nb/nrow(se.obj) * 100, digits = 2)
                if (top.rank.bio.genes >= 100) top.rank.bio.genes = 100
                # uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if (top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Updating the selection. Select top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
                ncg.selected <- row.names(se.obj) %in% ncg.selected
            }
            ##### grid group: bio ####
            if (grid.group == 'bio'){
                printColoredMessage(
                    message = '- The grid search will be applied on biological factor. ',
                    color = 'blue',
                    verbose = verbose
                    )
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the value of "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose)
                    lo <- top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if (top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- nrow(se.obj) - top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if (top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(bio.genes)[bio.genes$bio.ranks > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ##### check selection ####
                if (length(ncg.selected) == 0)
                    stop('No NCGs can be found based on the current parameters.')
                # gene selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # bio
                top.rank.bio.genes.nb <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes.nb/nrow(se.obj) * 100, digits = 0)
                if (top.rank.bio.genes >= 100) top.rank.bio.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Update the selection. Select top ',
                        top.rank.uv.genes * 100,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose
                    )
            }
            ##### grid group: uv ####
            if (grid.group == 'uv'){
                printColoredMessage(
                    message = '- The grid search will be applied on genes that are affected by unwanted variation. ',
                    color = 'blue',
                    verbose = verbose)
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the value of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- nrow(se.obj) - top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if (top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x){
                                if (!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the value of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if (top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                        top.uv.genes <- unique(unlist(lapply(
                            all.uv.tests,
                            function(x){
                                if (!is.null(x)){
                                    temp.data <- get(x)
                                    ranks.data <- unique(unlist(lapply(
                                        names(temp.data),
                                        function(y){
                                            index <- temp.data[[y]]$ranked.genes < top.rank.uv.genes.nb
                                            row.names(temp.data[[y]])[index] })))
                                }
                            })))
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                }
                ##### check selection ####
                if (length(ncg.selected) == 0)
                    stop('No NCGs can be found based on the current parameters.')
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if (top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Update the selection. Select top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes * 100,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
            }else {
                printColoredMessage(
                    message = paste0('- ', length(ncg.selected), ' genes are selected as NCGs.'),
                    color = 'blue',
                    verbose = verbose)
            }
        } else ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    # Filtering selected negative control genes ######
    if (isTRUE(filter.ncgs)){
        printColoredMessage(
            message = '- Filtering the selected NCGs based on publicly available stable genes:',
            color = 'blue',
            verbose = verbose
        )
        if (common.hk == 'cancer'){
            common.hks <- singscore::getStableGenes(n_stable = nb.stable.genes)
            common.hks <- intersect(common.hks, row.names(se.obj)[ncg.selected])
            ncg.selected <- row.names(se.obj) %in% common.hks
        }
        if (common.hk == 'non.cancer'){
            common.hks <- row.names(se.obj)[rowData(se.obj)[[hk.group]]]
            common.hks <- intersect(common.hks, row.names(se.obj)[ncg.selected])
            ncg.selected <- row.names(se.obj) %in% common.hks
        }
    }

    printColoredMessage(
        message = paste0(
            '- A set of ',
            sum(ncg.selected),
            ' genes are selected for NCG.'),
        color = 'blue',
        verbose = verbose
        )

    # Plotting ####
    if (isTRUE(create.ncg.rank.plot)){
        all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
        all.test.res <- lapply(
            all.uv.tests,
            function(x){
                if (!is.null(x)){
                    temp.data <- get(x)
                    all.test.res <- lapply(
                        names(temp.data),
                        function(y){
                            test.res <- temp.data[[y]][ , 'ranked.genes', drop = FALSE]
                            test.res$group <- rep(y , nrow(test.res))
                            test.res
                        })
                    all.test.res <- do.call(cbind, all.test.res)
                }
                all.test.res
            })
        all.test.res <- Filter(Negate(is.null), all.test.res)
        all.test.res <- do.call(cbind, all.test.res)
        temp.data <- lapply(
            seq(1, ncol(all.test.res), 2),
            function(x){
                temp.data <- all.test.res[ , x, drop = FALSE]
                colnames(temp.data) <- all.test.res[ , x+1][1]
                temp.data
            })
        temp.data <- do.call(cbind, temp.data)
        temp.data$Biology <- bio.genes$bio.ranks
        temp.data$NCG <- ncg.selected
        ha <- ComplexHeatmap::rowAnnotation(
            NCG = temp.data$NCG,
            col = list(NCG = c('TRUE' = 'gray10', 'FALSE' = 'gray'))
        )
        ncg.rank.plot <- ComplexHeatmap::Heatmap(
            temp.data[ , seq_len(ncol(temp.data) - 1)],
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_row_names = FALSE,
            right_annotation = ha,
            column_names_rot = 45,
            col = viridis::magma(n = 20),
            heatmap_legend_param = list(
                title = 'Ranks',
                title_gp = grid::gpar(fontsize = 14),
                by_row = TRUE,
                ncol = 1)
        )
        if (isTRUE(plot.ncg.rank)) print(ncg.rank.plot)
    }


    # Performance assessment of the selected NCG ####
    ## pca ####
    if (isTRUE(assess.ncg)){
        printColoredMessage(
            message = '-- Assessing the performance of the selected NCG set:',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '- Performing PCA on only the selected genes as NCG.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- Exploring the association of the first ',
                nb.pcs,
                ' with the ',
                paste0(uv.variables, collapse = ' & '),
                ' variables.'),
            color = 'blue',
            verbose = verbose
            )
        pca.data <- BiocSingular::runSVD(
            x = t(expr.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM = bsparam(),
            center = center,
            scale = scale
            )
        pc.var <- (pca.data$d^2) / (ncol(expr.data) - 1)
        centered.data <- scale(
            t(expr.data[ncg.selected , ]),
            center = center,
            scale = scale
            )
        total.var <- sum(colVars(centered.data))
        percentage <- round(x = c(pc.var / total.var) * 100, digits = 2)
        if (is.null(variables.to.assess.ncg))
            variables.to.assess.ncg <- uv.variables
        ## regression and vector correlations ####
        all.corr <- lapply(
            variables.to.assess.ncg,
            function(x){
                if (class(se.obj[[x]]) %in% c('numeric', 'integer')){
                    rSquared <- sapply(
                        1:nb.pcs,
                        function(y) summary(lm(se.obj[[x]] ~ pca.data$u[, 1:y]))$r.squared)
                } else if (class(se.obj[[x]]) %in% c('factor', 'character')){
                    catvar.dummies <- dummy_cols(se.obj[[x]])
                    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca.pcs <- sapply(
                        1:nb.pcs,
                        function(y){ cca <- cancor(
                            x = pca.data$u[, 1:y, drop = FALSE],
                            y = catvar.dummies)
                        1 - prod(1 - cca$cor^2)})
                }
            })
        names(all.corr) <- variables.to.assess.ncg
        pca.ncg <- as.data.frame(do.call(cbind, all.corr)) %>%
            dplyr::mutate(pcs = c(1:nb.pcs)) %>%
            tidyr::pivot_longer(
                -pcs,
                names_to = 'Groups',
                values_to = 'ls'
                )
        ncg.assessment.plot <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
            geom_line(aes(color = Groups), linewidth = 1) +
            geom_point(aes(color = Groups), size = 2) +
            xlab('PCs') +
            ylab ('Correlations') +
            ggtitle('') +
            scale_x_continuous(breaks = (1:nb.pcs), labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 16),
                strip.text.x = element_text(size = 10),
                plot.title = element_text(size = 16)
                )
        p.pca.percentage <- data.frame(var = percentage, no = 1:nb.pcs) %>%
            ggplot(., aes(x = no, y = percentage)) +
            geom_point(size = 3) +
            ylab('Variation(%)') +
            ylim(c(0,100)) +
            geom_line() +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.line.x  = element_blank(),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 10),
                legend.title = element_text(size = 14),
                plot.title = element_text(size = 16),
                plot.margin = unit(c(0, 0, 3, 0), "pt")
                )
        ncg.assessment.plot <- ncg.assessment.plot / p.pca.percentage + plot_layout(heights = c(3, 1))
        if (isTRUE(plot.ncg.assessment)) print(ncg.assessment.plot)
    }
    # Saving the NCGs ####
    ## Adding results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Saving the selected NCGs:',
        color = 'magenta',
        verbose = verbose
        )
    if (is.null(ncg.group.name)){
        ncg.group.name <- paste0('ncg|unsupervised')
    }
    if (is.null(ncg.set.name)){
        ncg.set.name <- paste0(
            sum(ncg.selected),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|AnoCorrMad:',
            ncg.selection.method,
            '|',
            assay.name
            )
    }
    if (is.logical(samples.to.use)){
        se.obj <- se.obj.all
    }
    if (isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- Saving the selected set of NCG to the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        ## Check
        if (!'NCG' %in% names(se.obj@metadata)) {
            se.obj@metadata[['NCG']] <- list()
        }
        ## check
        if (!'un.supervised' %in% names(se.obj@metadata[['NCG']])) {
            se.obj@metadata[['NCG']][['un.supervised']] <- list()
        }
        ## check
        if (!ncg.group.name %in% names(se.obj@metadata[['NCG']][['un.supervised']])) {
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]] <- list()
        }
        ## check
        if (!ncg.set.name %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]])) {
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]] <- list()
        }
        ## check
        if (!'ncg.set' %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]])) {
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- list()
        }
        se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- ncg.selected

        if (isTRUE(assess.ncg)){
            ## check
            if (!'assessment.plot' %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]])) {
                se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- list()
            }
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <-
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][['assessment.plot']][[ncg.set.name]] <- ncg.assessment.plot
        }
        if (isTRUE(create.ncg.rank.plot)){
            ## check
            if (!'rank.plot' %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]])) {
                se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['rank.plot']] <- list()
            }
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['rank.plot']] <- ncg.rank.plot
        }
        printColoredMessage(
            message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The findNcgsUnSupervised function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ## export output as vector ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = '- The NCGs are outpputed as a logical vector.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The findNcgsUnSupervised function finished.',
            color = 'white',
            verbose = verbose
            )
        return(list(ncg.selected = ncg.selected))
    }
}

