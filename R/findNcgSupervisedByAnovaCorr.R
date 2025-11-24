#' Finds NCGs using ANOVA and correlation.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function employs gene level correlation and ANOVA analyses across all samples to identify a set of suitable genes
#' as negative control genes (NCGs) for RUV-III normalization. Both biological and unwanted variation sources is necessary
#' and should be specified.
#'
#' @details
#' The function utilizes correlation analysis to identify genes that are highly affected by continuous sources of variation,
#' while it uses ANOVA to identify genes that are highly affected by categorical sources of variation. The function selects
#' genes as negative control genes (NCG) based on high correlation coefficients and F-statistics for unwanted sources of
#' variation, and low correlation coefficients and F-statistics for biological sources of variation. Various approaches
#' are employed for the final gene selection.
#' The function uses 5 ways to summarize two gene-level F-statistics obtained for the biological and unwanted variation.
#' The function uses either the values or the ranks of F-statistics for NCGs selection. The function ranks the negative
#' of F-statistics values for unwanted variation.
#' The lower the ranks, the greater the impact of unwanted variation on genes. The function ranks the F-statistics for
#' biological variation. The higher the ranks, the greater the impact of biological variation on genes. The options are
#' `prod`, `sum`, `average`, `auto` or `non.overlap` and `quantile`.
#'
#' If `prod`, `sum` and `average` is set:
#'
#' The product, sum or average of ranks of F-statistics is calculated. Then, the function selects `nb.ncg` `numbers of
#' genes as negative control genes that have the lowest ranks.
#'
#' If `non.overlap` is selected:
#' \enumerate{
#'    \item The function selects the top `top.rank.bio.genes` genes that have the highest ranks of F-statistics
#'    for biological variation.
#'    \item The function selects the top `top.rank.uv.genes` genes that have the lowest ranks of F-statistics for
#'    unwanted variation.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1. This will be a set of genes as
#'    negative control genes.
#' }
#'
#' If `auto` is selected:
#' \enumerate{
#'    \item The function selects the top `top.rank.bio.genes` genes that have the highest ranks of F-statistics for
#'    biological variation.
#'    \item  The function selects the top `top.rank.uv.genes` genes that have the lowest ranks of F-statistics  for
#'    unwanted variation.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1.
#'    \item If the number of selected genes is larger or smaller than the specified `nb.ncg`, the function applies an
#'    auto search to find approximate `nb.ncg` of genes as negative control genes as follow. The auto search will either
#'    decrease or increase the values of either `top.rank.bio.genes` or `top.rank.uv.genes` or both till to find
#'    approximate `nb.ncg` of genes as negative control genes.
#' }
#' If `quantile` is selected:
#' \enumerate{
#'    \item The function selects the `bio.percentile` percentile of F-statistics for biological variation. Then, selects
#'    all the genes that have F-statistics larger the calculated percentile.
#'    \item The function selects the `uv.percentile` percentile of F-statistics for unwanted variation. Then, selects
#'    all the genes that have F-statistics larger the calculated percentile.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1.
#' }
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character that indicates the name of an data (assay) in the `SummarizedExperiment` object.
#' The selected assay should be the one that will be used for the RUV-III-PRPS normalization.
#' @param bio.variables Character. A character string or vector of strings indicating the column name(s) of the biological
#' variable(s) in the SummarizedExperiment object. These variable can be categorical or continuous or a combination. This
#' argument cannot be `NULL`.
#' @param uv.variables Character. A character string or vector of strings indicating the column name(s) of the unwanted
#' variable(s) in the SummarizedExperiment object. These variable can be categorical or continuous or a combination.This
#'  argument cannot be `NULL`.
#' @param approach TTT
#' @param use.rank TTT
#' @param ncg.selection.method Character. A character that indicates how to summarize different statistics and select a
#' set of genes as negative control genes. The options are: `prod`, `average`, `sum`, `non.overlap`, `auto`, and `quantile`.
#' The default is set to `non.overlap`. For more information, refer to the details of the function.
#' @param nb.ncg Numeric. A numeric value that specifies the number of genes to be chosen as negative control genes (NCG)
#' when the `ncg.selection.method` parameter is set to `auto`. This value, `nb.ncg`, corresponds to a fraction of the total
#' genes in the SummarizedExperiment object. The default is set to 0.1.
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
#' @param  rank.continuous.var TTTT
#' @param  rank.categorical.var TTTT
#' @param filter.ncgs Logical. Whether to filter selected NCGs based on public human housekeeping gene sets. The default
#' is set to `FALSE`.
#' @param common.hk Character. Specifies group of housekeeping genes to use: `cancer` or `non.cancer`. The default is set
#' to `cancer`.
#' @param nb.stable.genes Numeric. A numeric value that specifies the number of top stable genes to be obtained from the
#' `getStableGenes()` function in the **singescore** R package. The default is set to 2000.
#' @param hk.group Character. Column name in the gene annotation containing non-cancer housekeeping genes. Options include:
#' `bulk.rnaseq.hk.genes.v1`, `bulk.rnaseq.hk.genes.v2`, `micorarray.hk.genes`, `nanostring.pan.cancer.hk.genes`,
#' `singscore.pan.cancer.hk.genes`. The default is set to `micorarray.hk.genes`.
#' @param grid.direction Character. A character that indicates whether the grid search should be performed in decreasing
#' or increasing order when the `ncg.selection.method` is set to `auto`. The options are: `increase` and `decrease`. The
#' default is set to `decrease`.
#' @param grid.nb Numeric. A numeric value that indicates the number of genes for grid search when the `ncg.selection.method`
#' is set to `auto`. In the `auto` approach, the grid search starts with the initial `top.rank.uv.genes` and
#' `top.rank.bio.genes` values and adds or drops the `grid.nb` in each loop to find `nb.ncg` of genes as negative control
#' genes. The default is set to 20.
#' @param create.ncg.rank.plot Logical. Indicates whether to generate a heatmap that shows the rank of the all genes
#' with respect to their biological and unwanted variation effects. The default is set to `FALSE`.
#' @param plot.ncg.rank Logical. Indicates whether to plot a heatmap that shows the rank of the all genes
#' with respect to their biological and unwanted variation effects, while function is running. The default is set to `FALSE`.
#' @param min.sample.for.aov Numeric. A numeric value that indicates the minimum number of samples that are required to
#' perform ANOVA analyses between continuous sources of variation (biological and unwanted variation) with individual
#' gene expression. The default is set to 3. The minimum value is 3.
#' @param min.sample.for.correlation Numeric. A numeric value that indicates the minimum number of samples that are required
#' to perform correlation analyses between continuous sources of variation (biological and unwanted variation) with
#'individual gene expression. The default is set to 10. The minimum value can be 3.
#' @param regress.out.bio.variables Character. A character string or vector of strings that indicate the column names of
#' biological variables in the SummarizedExperiment object that will be regressed out from the data before performing
#' correlation and ANOVA. Regressing out biological variables might help better identify genes that are highly affected
#' by unwanted variation. The default is set to `NULL`.
#' @param regress.out.uv.variables Character. A character string or vector of strings that indicate the column names of
#' unwanted variables in the SummarizedExperiment object that will be regressed out from the data before performing
#' correlation and ANOVA. Regressing out unwanted variables might help better identify genes that are highly affected by
#' biological variation. The default is set to `NULL`.
#' @param normalization Character. A character that indicates which normalization method should be use to mitigate the
#' variation in library size before finding genes that are highly affected by biological variation. The options are :
#' `CPM`, `TMM`, `VST`, `upper`, `full` and `medium`. The default is set to  `CPM`. Refer to the `applyOtherNormalization()`
#' function for more details.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before any statistical analyses.
#' The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value to be added as a pseudo count to all measurements before applying log
#' transformation. The default is set to 1.
#' @param corr.method Character. A character that indicates which correlation methods should be used for the correlation
#' analyses. The options are `pearson` or `spearman`. The default is set to `spearman`.
#' @param a Numeric. The significance level used for the confidence intervals in the correlation; by default, it is set
#' to 0.05. Refer to the function `correls` from the **Rfast** R package for more details.
#' @param rho Numeric. The value of the hypothesized correlation to be used in the hypothesis testing. The default is
#' set to 0. Refer to the function `correls` from the **Rfast** R package for more details.
#' @param anova.method Character. A character that  indicates which ANOVA method to use. The options are `aov` or `welch`.
#' The default is se to `aov`. Refer to the function `row_oneway_equalvar()` or `row_oneway_welch()` from the
#' **matrixTests** R package for more details.
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
#' @param imf.name Character string. A name to save the intermediate file. If `NULL`, the function generates a name.
#' @param use.imf Logical. Indicates whether to use the intermediate file. The default is set to `FALSE`.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the SummarizedExperiment
#' object or output the result. The default is `TRUE`.
#' @param verbose Logical. If `TRUE`, shows messages of different steps of the function.
#' @param samples.to.use TTTT
#' @param bio.clustering.method TTTT
#' @param bio.groups TTTT
#' @param nb.bio.clusters TTTT
#' @param uv.groups TTTT
#' @param uv.clustering.method TTTT
#' @param nb.uv.clusters TTTT
#' @param svd.bsparam TTTT
#'
#' @return Either the SummarizedExperiment object containing a set of negative control genes in the metadata or a
#' logical vector of the selected negative control genes.
#'
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @importFrom matrixTests row_oneway_welch row_oneway_equalvar
#' @importFrom BiocSingular runSVD bsparam
#' @importFrom fastDummies dummy_cols
#' @importFrom matrixStats rowProds
#' @importFrom tidyr pivot_longer
#' @importFrom stats quantile
#' @importFrom Rfast correls
#' @importFrom dplyr mutate
#' @import ggplot2
#' @export

findNcgSupervisedByAnovaCorr <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        approach = 'AnovaCorr.AcrossAllSamples',
        nb.ncg = 0.1,
        samples.to.use = 'all',
        use.rank = FALSE,
        rank.continuous.var = 'correlation',
        rank.categorical.var = 'fvalue',
        ncg.selection.method = 'non.overlap',
        top.rank.bio.genes = 0.5,
        top.rank.uv.genes = 0.5,
        bio.percentile = 0.2,
        uv.percentile = 0.8,
        grid.group = 'uv',
        grid.direction = 'decrease',
        grid.nb = 20,
        filter.ncgs = FALSE,
        common.hk = 'cancer',
        nb.stable.genes = 2000,
        hk.group = 'micorarray.hk.genes',
        create.ncg.rank.plot = FALSE,
        plot.ncg.rank = FALSE,
        min.sample.for.aov = 3,
        min.sample.for.correlation = 10,
        regress.out.bio.variables = NULL,
        regress.out.uv.variables = NULL,
        bio.groups = NULL,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.groups = NULL,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        normalization = 'CPM',
        apply.log = TRUE,
        pseudo.count = 1,
        corr.method = "spearman",
        a = 0.05,
        rho = 0,
        anova.method = 'aov',
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 10,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        plot.ncg.assessment = TRUE,
        ncg.group.name = NULL,
        ncg.set.name = NULL,
        save.imf = FALSE,
        imf.name = NULL,
        use.imf = FALSE,
        check.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    # Applying across all samples ####
    if(approach == 'AnovaCorr.AcrossAllSamples' | approach == 'both'){
        printColoredMessage(message = '------------The findNcgByAnovaCorr function starts with "AcrossAllSamples" mode:',
                            color = 'white',
                            verbose = verbose)
        # Checking the function inputs ####
        if (!is.vector(assay.name) | length(assay.name) > 1 | is.logical(assay.name)){
            stop('The "assay.name" must be a single data(assay) name in the SummarizedExperiment object.')
        }
        if (is.null(bio.variables)){
            stop('The "bio.variables" cannot be empty or "NULL".')
        }
        if (is.null(uv.variables)){
            stop('The "uv.variables" cannot be empty or "NULL".')
        }
        if (!is.vector(bio.variables) | !is.vector(uv.variables) ){
            stop('The "uv.variables" and "bio.variables" must be a vector of variables names in the SummarizedExperiment object.')
        }
        if (length(intersect(bio.variables, uv.variables)) > 0){
            stop('Individual specified variable must be either in "bio.variables" or "uv.variables".')
        }
        if (nb.ncg >= 1 | nb.ncg <= 0){
            stop('The "nb.ncg" must be a positve value 0 < nb.ncg < 1.')
        }
        if (!ncg.selection.method %in% c('prod', 'sum', 'average', 'auto', 'non.overlap', 'quantile')){
            stop('The "ncg.selection.method" muat be one of "prod", "sum", "average", "auto", "non.overlap" or "quantile".')
        }
        if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
            stop('The "top.rank.bio.genes" must be a positve value  0 < top.rank.bio.genes < 1.')
        }
        if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
            stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes < 1.')
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
            stop('The "anova.method" must be one of the "aov" or "welch".')
        }
        if (isFALSE(is.logical(assess.ncg))){
            stop('The "assess.ncg" must be "TRUE" or "FALSE.')
        }
        if (nb.pcs < 0){
            stop('The "nb.pcs" must be a postive integer value.')
        }
        if (isFALSE(is.logical(scale))) {
            stop('The "scale" must be "TRUE" or "FALSE".')
        }
        if (isFALSE(is.logical(center))) {
            stop('The "center" must be "TRUE" or "FALSE".')
        }
        if (isFALSE(is.logical(check.se.obj))) {
            stop('The "check.se.obj" must be "TRUE" or "FALSE".')
        }
        if (isFALSE(is.logical(apply.log))) {
            stop('The "apply.log" must be "TRUE" or "FALSE".')
        }
        if (isFALSE(is.logical(save.se.obj))) {
            stop('The "save.se.obj" must be "TRUE" or "FALSE".')
        }
        if (isFALSE(is.logical(use.imf))) {
            stop('The "use.imf" must be "TRUE" or "FALSE".')
        }
        if (isFALSE(is.logical(save.imf))) {
            stop('The "save.imf" must be "TRUE" or "FALSE".')
        }
        if (isFALSE(is.logical(verbose))) {
            stop('The "verbose" must be "TRUE" or "FALSE".')
        }
        if (!is.null(regress.out.bio.variables) | !is.null(regress.out.uv.variables)){
            if (is.logical(regress.out.bio.variables) | is.logical(regress.out.uv.variables))
                stop(paste0('The "regress.out.bio.variables" or "regress.out.bio.variables" ',
                            'must names of columns in the the SummarizedExperiment object.'))
        }
        if (isTRUE(ncg.selection.method == 'auto')){
            if (isFALSE(is.numeric(grid.nb))){
                stop('The "grid.nb" must be a postive integer value.')
            } else if (grid.nb < 0 | length(grid.nb) > 1 ){
                stop('The "grid.nb" must be a postive integer value.')
            } else if (isTRUE(is.logical(grid.group))){
                stop('The "grid.group" must be on of the "both", "uv" or "bio".')
            } else if (isTRUE(length(grid.group) > 1)){
                stop('The "grid.group" must be on of the "both", "uv" or "bio".')
            } else if (isTRUE(!grid.group %in% c('both', 'uv', 'bio'))){
                stop('The "grid.group" must be on of the "both", "uv" or "bio".')
            } else if (isTRUE(is.logical(grid.direction))){
                stop('The "grid.direction" must be on of the "decrease",or "increase".')
            } else if (isTRUE(length(grid.direction) > 1)){
                stop('The "grid.direction" must be on of the "decrease" or "increase".')
            } else if (isTRUE(!grid.direction %in% c('decrease', 'increase'))){
                stop('The "grid.direction" must be on of the "decrease" or "increase".')
            }
        }
        if (!is.null(normalization)){
            if (!is.null(regress.out.uv.variables))
                printColoredMessage(
                    message = paste0('Both normalization and regress.out.uv.variables are selected.',
                                     'The function will perfom normalization first and the regression the UV variables.'),
                    color = 'magenta',
                    verbose = verbose)
        }
        if (isTRUE(apply.log)){
            if (length(pseudo.count) > 1 | pseudo.count < 0 | is.null(pseudo.count))
                stop('The "pseudo.count" must be 0 or a postive integer value.')
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
        if (isFALSE(check.se.obj)) {
            if (isTRUE(sum(bio.variables %in% colnames(colData(se.obj))) != length(bio.variables))) {
                stop('All or some of "bio.variables" cannot be found in the SummarizedExperiment object.')
            }
            if (isTRUE(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables))) {
                stop('All or some of "uv.variables" cannot be found in the SummarizedExperiment object.')
            }
            if (!is.null(variables.to.assess.ncg)) {
                if (isTRUE(sum(variables.to.assess.ncg %in% colnames(colData(se.obj))) != length(variables.to.assess.ncg)))
                    stop('All or some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
            }
        }

        # Checking the SummarizedExperiment object ####
        if (isTRUE(check.se.obj)) {
            se.obj <- checkSeObj(
                se.obj = se.obj,
                assay.names = assay.name,
                variables = unique(c(bio.variables, uv.variables, variables.to.assess.ncg)),
                remove.na = remove.na,
                verbose = verbose
            )
        }
        if (remove.na == 'none'){
            if (is.null(variables.to.assess.ncg))
                variables.to.assess.ncg <- c(bio.variables, uv.variables)
            mout <- lapply(
                variables.to.assess.ncg,
                function(x){
                    if (sum(is.na(se.obj[[x]])) > 0)
                        stop('There are NA or missing values in the specified variables.')
                })
        }

        # Selecting a subset of samples for analysis ####
        if (is.logical(samples.to.use)){
            se.obj.all <- se.obj
            se.obj <- se.obj[ , samples.to.use]
        }
        # Applying data transformation and normalization ####
        printColoredMessage(
            message = '-- Applying data transformation and normalization:',
            color = 'magenta',
            verbose = verbose
        )
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
                    '- The ',
                    assay.name,
                    ' data will be used without any log transformation.'),
                color = 'blue',
                verbose = verbose)
            expr.data <- assay(x = se.obj, i = assay.name)
        }
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
        ## Regressing out unwanted variables ####
        if (!is.null(regress.out.uv.variables) & !is.null(normalization)){
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(regress.out.uv.variables, collapse = ' & '),
                    ' will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    'Note, we do not recommend regressing out the ',
                    paste0(regress.out.uv.variables, collapse = ' & '),
                    ' if they are largely associated with the ',
                    paste0(bio.variables, collapse = ' & '),
                    '.'),
                color = 'red',
                verbose = verbose
            )
            expr.data.reg.uv <- t(expr.data.nor)
            uv.variables.all <- paste('se.obj', regress.out.uv.variables, sep = '$')
            expr.data.reg.uv <- lm(as.formula(paste(
                'expr.data.reg.uv',
                paste0(uv.variables.all, collapse = '+') ,
                sep = '~')))
            expr.data.reg.uv <- t(expr.data.reg.uv$residuals)
            colnames(expr.data.reg.uv) <- colnames(se.obj)
            row.names(expr.data.reg.uv) <- row.names(se.obj)
        }
        if (!is.null(regress.out.uv.variables) & is.null(normalization)){
            printColoredMessage(
                message = paste0(
                    'The',
                    paste0(regress.out.uv.variables, collapse = ' & '),
                    ' will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    'Note: we do not recommend regressing out ',
                    paste0(regress.out.uv.variables, collapse = ' & '),
                    'if they are largely associated with the ',
                    paste0(bio.variables, collapse = ' & '), '.'),
                color = 'red',
                verbose = verbose
            )
            expr.data.reg.uv <- t(expr.data)
            uv.variables.all <- paste('se.obj', regress.out.uv.variables, sep = '$')
            expr.data.reg.uv <- lm(as.formula(paste(
                'expr.data.reg.uv',
                paste0(uv.variables.all, collapse = '+') ,
                sep = '~'
            )))
            expr.data.reg.uv <- t(expr.data.reg.uv$residuals)
            colnames(expr.data.reg.uv) <- colnames(se.obj)
            row.names(expr.data.reg.uv) <- row.names(se.obj)
        }

        ## Regressing out biological variables ####
        if (!is.null(regress.out.bio.variables)){
            printColoredMessage(
                message = paste0(
                    paste0(regress.out.bio.variables, collapse = ' & '),
                    ' will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    'We do not recommend regressing out the ',
                    paste0(regress.out.bio.variables, collapse = ' & '),
                    'if they are largely associated with the ',
                    paste0(uv.variables, collapse = ' & '), '.'),
                color = 'red',
                verbose = verbose
            )
            expr.data.reg.bio <- t(expr.data)
            bio.variables.all <- paste('se.obj', regress.out.bio.variables, sep = '$')
            expr.data.reg.bio <- lm(as.formula(paste(
                'expr.data.reg.bio',
                paste0(bio.variables.all, collapse = '+') ,
                sep = '~')))
            expr.data.reg.bio <- t(expr.data.reg.bio$residuals)
            colnames(expr.data.reg.bio) <- colnames(se.obj)
            row.names(expr.data.reg.bio) <- row.names(se.obj)
        }

        # Applying gene-level ANOVA and correlation analyses  ####
        if (isFALSE(use.imf)){
            ## Selecting genes that are highly affected by unwanted variation ####
            printColoredMessage(
                message = '-- Finding genes that are highly affected by each specified source(s) of unwnated variation:',
                color = 'magenta',
                verbose = verbose
            )
            uv.var.class <- unlist(lapply(
                uv.variables,
                function(x) class(colData(se.obj)[[x]]))
            )
            categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
            continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
            ### anova between genes and categorical sources of unwanted variation ####
            if (length(categorical.uv) > 0 ){
                if (!is.null(regress.out.bio.variables)){
                    data.to.use <- expr.data.reg.bio
                } else data.to.use <- expr.data
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
                        if (length(keep.samples) == 0){
                            stop(paste0(
                                'There are not enough samples to perform ANOVA between individual gene expression and the ',
                                x,
                                ' variable. A possible solution is to lower "min.sample.for.aov" or remove',
                                x,
                                'from the "uv.variables" and re-run the function.')
                            )
                        } else if (length(keep.samples) == 1 ){
                            stop(paste0(
                                'There is only a single batch from ',
                                x,
                                ' that have enough samples ',
                                min.sample.for.aov,
                                '(min.sample.for.aov). A possible solution is to lower min.sample.for.aov or remove')
                            )
                        } else if (length(keep.samples) != length(unique(colData(se.obj)[[x]])) ){
                            not.coverd <- unique(colData(se.obj)[[x]])[!unique(colData(se.obj)[[x]]) %in% keep.samples]
                            printColoredMessage(
                                message = paste0(
                                    '- Note, the ',
                                    paste0(not.coverd, collapse = '&'),
                                    ' batches do not have enough samples for the ANOVA analysis.'),
                                color = 'red',
                                verbose = verbose
                                )
                        }
                        keep.samples <- colData(se.obj)[[x]] %in% keep.samples
                        if (anova.method == 'aov'){
                            anova.gene.batch <- as.data.frame(row_oneway_equalvar(
                                x = data.to.use[ , keep.samples],
                                g = se.obj@colData[, x][keep.samples])
                            )
                        }
                        if (anova.method == 'welch.correction'){
                            anova.gene.batch <- as.data.frame(row_oneway_welch(
                                x = data.to.use[ , keep.samples],
                                g = se.obj@colData[, x][keep.samples])
                            )
                        }
                        # effect size
                        anova.gene.batch$p.adjusted <- p.adjust(anova.gene.batch$pvalue, method = "BH")
                        anova.gene.batch$eta.squared <- (anova.gene.batch$statistic * anova.gene.batch$df.between) /
                            (anova.gene.batch$statistic * anova.gene.batch$df.between + anova.gene.batch$df.within)
                        anova.gene.batch$g.statistic <- anova.gene.batch$eta.squared
                        set.seed(2233)
                        anova.gene.batch$ranked.genes <- rank(
                            x = -anova.gene.batch[, 'statistic'],
                            ties.method = 'random'
                        )
                        anova.gene.batch
                    })
                names(anova.genes.uv) <- categorical.uv
                # rm(data.to.use)
            } else anova.genes.uv <- NULL

            ### correlation between genes and categorical sources of unwanted variation ####
            if (length(continuous.uv) > 0 ){
                if (!is.null(regress.out.bio.variables)){
                    data.to.use <- expr.data.reg.bio
                } else data.to.use <- expr.data
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
                if (ncol(se.obj) <= min.sample.for.correlation){
                    stop(paste0(
                        'There are not enough samples (min.sample.for.correlation:',
                        min.sample.for.correlation,
                        ') to perform correlation analysis.',
                        ' A possible soultion in to lower min.sample.for.correlation.')
                    )
                }
                corr.genes.uv <- lapply(
                    continuous.uv,
                    function(x) {
                        corr.genes.var <- as.data.frame(correls(
                            y = se.obj@colData[, x],
                            x = t(data.to.use),
                            type = corr.method,
                            a = a ,
                            rho = rho)
                            )
                        corr.genes.var <- cbind(
                            round(x = corr.genes.var[, 1:4], digits = 3),
                            corr.genes.var[, 5, drop = FALSE]
                            )
                        corr.genes.var$p.adjusted <- p.adjust(corr.genes.var[['p-value']], method = "BH")
                        corr.genes.var$g.statistic <- abs(corr.genes.var$correlation)
                        set.seed(2233)
                        colnames(corr.genes.var)[colnames(corr.genes.var) == 'correlation' ] <- 'statistic'
                        corr.genes.var$ranked.genes <- rank(
                            -abs(corr.genes.var[, 'statistic']),
                            ties.method = 'random'
                            )
                        row.names(corr.genes.var) <- row.names(data.to.use)
                        corr.genes.var
                    })
                names(corr.genes.uv) <- continuous.uv
                rm(data.to.use)
            } else corr.genes.uv <- NULL

            ## Selecting genes that are  highly affected by biology ####
            printColoredMessage(
                message = '-- Finding genes that are highly affected by each specified source(s) of biological variation:',
                color = 'magenta',
                verbose = verbose
                )
            bio.var.class <- unlist(lapply(
                bio.variables,
                function(x) class(colData(se.obj)[[x]]))
                )
            continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
            categorical.bio <- bio.variables[bio.var.class %in% c('factor', 'character')]
            ### anova between genes and categorical sources of biological variation ####
            if (length(categorical.bio) > 0 ){
                if (!is.null(normalization) & is.null(regress.out.uv.variables)){
                    data.to.use <- expr.data.nor
                } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)){
                    data.to.use <- expr.data.reg.uv
                } else if (!is.null(regress.out.uv.variables) & is.null(normalization)){
                    data.to.use <- expr.data.reg.uv
                } else if (is.null(regress.out.uv.variables) & is.null(normalization)){
                    data.to.use <- expr.data
                }
                printColoredMessage(
                    message = paste0(
                        '- Performing ANOVA between individual gene-level ',
                        'expression and each categorical source of biological variation: ',
                        paste0(categorical.bio, collapse = ' & '),
                        '.'),
                    color = 'blue',
                    verbose = verbose
                )
                anova.genes.bio <- lapply(
                    categorical.bio,
                    function(x) {
                        keep.samples <- findRepeatingPatterns(
                            vec = colData(se.obj)[[x]],
                            n.repeat = min.sample.for.aov
                            )
                        if ( length(keep.samples) == 0){
                            stop(paste0(
                                'There are not enough samples to perfrom ANOVA between individual genes expression and the ',
                                x,
                                ' variable. Possible solutions is to lower min.sample.for.aov or remove',
                                x,
                                'from the bio.variables and re-run the function.'))
                        } else if (length(keep.samples) == 1 ){
                            stop(paste0(
                                'There is only a single batch from ',
                                x,
                                ' that have enough samples ',
                                min.sample.for.aov,
                                '(min.sample.for.aov). Possible solutions is to lower min.sample.for.aov or remove',
                                x,
                                'from the bio.variables and re-run the function'))
                        } else if (length(keep.samples) != length(unique(colData(se.obj)[[x]])) ){
                            not.coverd <- unique(colData(se.obj)[[x]])[!unique(colData(se.obj)[[x]]) %in% keep.samples]
                            printColoredMessage(
                                message = paste0(
                                    'Note, the',
                                    paste0(not.coverd, collapse = '&'),
                                    ' groups do not have enough samples for the ANOVA analysis.'),
                                color = 'red',
                                verbose = verbose
                                )
                        }
                        keep.samples <- colData(se.obj)[[x]] %in% keep.samples
                        if (anova.method == 'aov'){
                            anova.genes <- as.data.frame(row_oneway_equalvar(
                                x = data.to.use[ , keep.samples],
                                g = se.obj@colData[, x][keep.samples])
                                )
                        } else if (anova.method == 'welch.correction'){
                            anova.genes <- as.data.frame(row_oneway_equalvar(
                                x = data.to.use[ , keep.samples],
                                g = se.obj@colData[, x][keep.samples])
                                )
                        }
                        anova.genes$p.adjusted <- p.adjust(anova.genes$pvalue, method = "BH")
                        anova.genes$eta.squared <- (anova.genes$statistic * anova.genes$df.between) /
                            (anova.genes$statistic * anova.genes$df.between + anova.genes$df.within)
                        anova.genes$g.statistic <- anova.genes$eta.squared
                        set.seed(2233)
                        anova.genes$ranked.genes <- rank(anova.genes[, 'statistic'], ties.method = 'random')
                        anova.genes
                    })
                names(anova.genes.bio) <- categorical.bio
                rm(data.to.use)
            } else anova.genes.bio <- NULL

            ### correlation between genes and continuous sources of biological variation ####
            if (length(continuous.bio) > 0 ){
                if (!is.null(normalization) & is.null(regress.out.uv.variables)){
                    data.to.use <- expr.data.nor
                } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)){
                    data.to.use <- expr.data.reg.uv
                } else if (!is.null(regress.out.uv.variables) & is.null(normalization)){
                    data.to.use <- expr.data.reg.uv
                } else if (is.null(regress.out.uv.variables) & is.null(normalization)){
                    data.to.use <- expr.data
                }
                ### gene-batch anova
                printColoredMessage(
                    message = paste0(
                        '- Performing ',
                        corr.method,
                        ' correlation between individual gene-level ',
                        'expression and each continuous sources of biological variation: ',
                        paste0(continuous.bio, collapse = '&'),
                        '.'),
                    color = 'blue',
                    verbose = verbose
                )
                if (ncol(se.obj) <= min.sample.for.correlation){
                    stop(paste0(
                        'There are not enough samples (min.sample.for.correlation:',
                        min.sample.for.correlation,
                        ') to perform correlation analysis.',
                        ' A possible soultion in to lower min.sample.for.correlation.')
                    )
                }
                corr.genes.bio <- lapply(
                    continuous.bio,
                    function(x) {
                        corr.genes.bio.var <- as.data.frame(correls(
                            y = se.obj@colData[, x],
                            x = t(data.to.use),
                            type = corr.method,
                            a = a ,
                            rho = rho))
                        corr.genes.bio.var <- cbind(
                            round(x = corr.genes.bio.var[, 1:4], digits = 3),
                            corr.genes.bio.var[, 5, drop = FALSE]
                            )
                        row.names(corr.genes.bio.var) <- row.names(data.to.use)
                        corr.genes.bio.var$p.adjusted <- p.adjust(corr.genes.bio.var[['p-value']], method = "BH")
                        corr.genes.bio.var$g.statistic <- abs(corr.genes.bio.var$correlation)
                        set.seed(2233)
                        colnames(corr.genes.bio.var)[colnames(corr.genes.bio.var) == 'correlation' ] <- 'statistic'
                        corr.genes.bio.var$ranked.genes <- rank(
                            x = abs(corr.genes.bio.var[, 'statistic']),
                            ties.method = 'random')
                        corr.genes.bio.var
                    })
                names(corr.genes.bio) <- continuous.bio
            } else corr.genes.bio <- NULL
        }
        # Reading the intermediate file ####
        if (isTRUE(use.imf)){
            printColoredMessage(
                message = '- Retrieving the results of ANOVA and correlations from the the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            if (is.null(imf.name)){
                imf.name <- paste0(assay.name, '|AcrossSamples|', ncg.selection.method)
            }
            if (is.null(se.obj@metadata$IMF$NCG[[imf.name]]))
                stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
            all.tests <- se.obj@metadata$IMF$NCG[[imf.name]]
            anova.genes.bio <- all.tests$anova.genes.bio
            corr.genes.bio <- all.tests$corr.genes.bio
            anova.genes.uv <- all.tests$anova.genes.uv
            corr.genes.uv <- all.tests$corr.genes.uv
        }
        # Saving the intermediate file ####
        if (isTRUE(save.imf)){
            printColoredMessage(
                message = '-- Save a intermediate file:',
                color = 'magenta',
                verbose = verbose)
            printColoredMessage(
                message = '- The results of ANOVA and correlations are saved in the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            if (length(se.obj@metadata$IMF) == 0 ) {
                se.obj@metadata[['IMF']] <- list()
            }
            if (!'NCG' %in% names(se.obj@metadata[['IMF']])){
                se.obj@metadata[['IMF']][['NCG']] <- list()
            }
            if (is.null(imf.name)){
                imf.name <- paste0(assay.name, '|AcrossSamples|', ncg.selection.method)
            }
            if (!imf.name %in% names(se.obj@metadata[['IMF']][['NCG']])){
                se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list()
            }
            se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list(
                anova.genes.bio = anova.genes.bio,
                corr.genes.bio = corr.genes.bio,
                anova.genes.uv = anova.genes.uv,
                corr.genes.uv = corr.genes.uv)
        }

        # Selecting a set of genes as NCG ####
        printColoredMessage(
            message = '-- Selecting a set of genes as NCG:',
            color = 'magenta',
            verbose = verbose
            )
        ## Ratio ####
        if (isFALSE(use.rank)){
            all.tests <- c(
                'anova.genes.bio',
                'corr.genes.bio',
                'anova.genes.uv',
                'corr.genes.uv'
                )
            var.partition <- lapply(
                all.tests,
                function(x){
                    temp <- get(x)
                    if (length(names(temp))!=0){
                        ranks.data <- lapply(
                            names(temp),
                            function(y) temp[[y]]$g.statistic)
                        ranks.data <- do.call(cbind, ranks.data)
                        colnames(ranks.data) <- names(temp)
                        ranks.data }
                })
            var.partition <- as.data.frame(do.call(cbind, var.partition))
            row.names(var.partition) <- row.names(se.obj)
            var.partition$bio <- rowSums(var.partition[ , bio.variables, drop = FALSE])
            var.partition$uv <- rowSums(var.partition[ , uv.variables, drop = FALSE])
            var.partition$ratio <- var.partition$uv / c(var.partition$bio + 1e-6)
            var.partition <- var.partition[order(-var.partition$ratio), ]
            if (ncg.selection.method == 'quantile'){
                if (!is.null(uv.percentile) & !is.null(bio.percentile)){
                    var.partition <- var.partition[var.partition$uv > quantile(x = var.partition$uv, probs = uv.percentile) , ]
                    if (quantile(x = var.partition$bio, probs = bio.percentile) == 0){
                        var.partition <- var.partition
                    } else {
                        var.partition <- var.partition[var.partition$bio < quantile(x = var.partition$bio, probs = bio.percentile) , ]
                    }
                    nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                    ncg.selected <- row.names(se.obj) %in% row.names(var.partition)[1:nb.ncg]
                } else if (!is.null(uv.percentile) & is.null(bio.percentile)){
                    var.partition <- var.partition[var.partition$uv > quantile(x = var.partition$uv, probs = uv.percentile) , ]
                    nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                    ncg.selected <- row.names(se.obj) %in% row.names(var.partition)[1:nb.ncg]
                } else if (is.null(uv.percentile) & !is.null(bio.percentile)){
                    var.partition <- var.partition[var.partition$bio < quantile(x = var.partition$bio, probs = bio.percentile) , ]
                    nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                    ncg.selected <- row.names(se.obj) %in% row.names(var.partition)[1:nb.ncg]
                } else if (is.null(uv.percentile) & is.null(bio.percentile)){
                    nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                    ncg.selected <- row.names(se.obj) %in% row.names(var.partition)[1:nb.ncg]
                }
            }
        }
        ## prod, sum or average of ranks ####
        if (ncg.selection.method %in% c('prod', 'sum', 'average') & isTRUE(use.rank)){
            all.tests <- c(
                'anova.genes.bio',
                'corr.genes.bio',
                'anova.genes.uv',
                'corr.genes.uv'
                )
            all.stats <- lapply(
                all.tests,
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
            all.stats <- do.call(cbind, all.stats)
            row.names(all.stats) <- row.names(se.obj)
            ### product of ranks ####
            if (ncg.selection.method == 'prod'){
                printColoredMessage(
                    message = '- A set of NCG will be selected based on the product of ranks.',
                    color = 'blue',
                    verbose = verbose
                )
                stat.summary <- rowProds(all.stats)
                if (sum(is.infinite(stat.summary)) > 0)
                    stop('The product of ranks results in infinity values.')
            }
            ## average of ranks ####
            if (ncg.selection.method == 'sum'){
                printColoredMessage(
                    message = '- A set of NCG will be selected based on the sum of ranks.',
                    color = 'blue',
                    verbose = verbose
                )
                stat.summary <- rowSums(all.stats)
            }
            ## sum of ranks ####
            if (ncg.selection.method == 'average'){
                printColoredMessage(
                    message = '- A set of NCG will be selected based on the average of ranks.',
                    color = 'blue',
                    verbose = verbose
                )
                stat.summary <- rowMeans(all.stats)
            }
            ## select top genes as NCGS ####
            all.stats <- as.data.frame(all.stats)
            row.names(all.stats) <- row.names(se.obj)
            all.stats$stat.summary <- stat.summary
            set.seed(112233)
            all.stats$rank.stat.summary <- rank(x = all.stats$stat.summary, ties.method = 'random')
            all.stats <- all.stats[order(all.stats$rank.stat.summary, decreasing = FALSE) , ]
            ncg.selected <- row.names(all.stats[1:round(nb.ncg* nrow(se.obj), digits = 0) , ])
            ncg.selected <- row.names(se.obj) %in% ncg.selected
        }

        ## non.overlap approach ####
        if (ncg.selection.method == 'non.overlap' & isTRUE(use.rank)){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the "non.overlap" approach.',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    '- Selecting top ',
                    top.rank.uv.genes * 100,
                    '% of highly affected genes by the unwanted variation, and then exclude all top ',
                    top.rank.bio.genes *100,
                    '% of highly affected genes by the bioloigcal variation.'),
                color = 'blue',
                verbose = verbose
            )
            ### select genes affected by biological variation ####
            top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
            top.bio.genes <- unique(unlist(lapply(
                all.bio.tests,
                function(x){
                    if (!is.null(x)){
                        temp.data <- get(x)
                        ranks.data <- unique(unlist(lapply(
                            names(temp.data),
                            function(y){
                                index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                row.names(temp.data[[y]])[index] })))
                    }
                })))

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
            ## select of NCGS ####
            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
            if (isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
            ncg.selected <- row.names(se.obj) %in% ncg.selected
        }

        ## quantile approach ####
        if (isTRUE(ncg.selection.method == 'quantile') & isTRUE(use.rank)){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the "quantile" approach..',
                color = 'blue',
                verbose = verbose
            )
            ### find biological percentile ####
            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
            top.bio.genes <- unique(unlist(lapply(
                all.bio.tests,
                function(x){
                    if (!is.null(x)){
                        temp.data <- get(x)
                        ranks.data <- unique(unlist(lapply(
                            names(temp.data),
                            function(y){
                                index <- temp.data[[y]]$statistic > quantile(x = temp.data[[y]]$statistic , probs = bio.percentile)
                                row.names(temp.data[[y]])[index] })))
                    }
                })))

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
                    ' (' ,
                    uv.percentile* 100,
                    '% percentile), and exclude any genes presents in ',
                    length(top.bio.genes),
                    ' genes with the biological F-statistics higher than ',
                    ' (' ,
                    bio.percentile* 100,
                    '% percentile).'),
                color = 'blue',
                verbose = verbose
            )
            top.uv.genes <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
            if (isTRUE(length(top.uv.genes) == 0)) stop('No NCGs can be found based on the current parameters.')
            ncg.selected <- row.names(se.obj) %in% top.uv.genes
        }

        ## auto approach ####
        if (ncg.selection.method == 'auto' & isTRUE(use.rank)){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the "auto" approach.',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    '- Selecting top ',
                    top.rank.uv.genes * 100,
                    '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                    top.rank.bio.genes * 100,
                    '% of highly affected genes by the bioloigcal variation.'),
                color = 'blue',
                verbose = verbose
            )
            ### select genes affected by biological variation ####
            top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
            top.bio.genes <- unique(unlist(lapply(
                all.bio.tests,
                function(x){
                    if (!is.null(x)){
                        temp.data <- get(x)
                        ranks.data <- unique(unlist(lapply(
                            names(temp.data),
                            function(y){
                                index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                row.names(temp.data[[y]])[index] })))
                    }
                })))
            ## select genes affected by unwanted variation ####
            top.rank.uv.genes.nb <- round(c(top.rank.uv.genes * nrow(se.obj)), digits = 0)
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
                message = paste0(
                    '- ',
                    length(ncg.selected),
                    ' genes are found.'),
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
                        verbose = verbose
                    )
                }
                if (isTRUE(nb.ncg < length(ncg.selected))){
                    con <- parse(text = paste0("length(ncg.selected)", ">", "nb.ncg"))
                    printColoredMessage(
                        message = paste0(
                            '- The number of selected genes ',
                            length(ncg.selected),
                            ' is larger than the number (',
                            nb.ncg ,
                            ') of specified genes ',
                            'by "nb.ncg". A grid search will be performed.'),
                        color = 'blue',
                        verbose = verbose
                    )
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
                            verbose = verbose
                        )
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
                            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                            top.bio.genes <- unique(unlist(lapply(
                                all.bio.tests,
                                function(x){
                                    if (!is.null(x)){
                                        temp.data <- get(x)
                                        ranks.data <- unique(unlist(lapply(
                                            names(temp.data),
                                            function(y){
                                                index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                                row.names(temp.data[[y]])[index] })))
                                    }
                                })))
                            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                        }
                    }
                    ### decreasing order ####
                    if (grid.direction == 'decrease'){
                        printColoredMessage(
                            message = '- The grid search will decrease the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                            color = 'blue',
                            verbose = verbose
                        )
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
                            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                            top.bio.genes <- unique(unlist(lapply(
                                all.bio.tests,
                                function(x){
                                    if (!is.null(x)){
                                        temp.data <- get(x)
                                        ranks.data <- unique(unlist(lapply(
                                            names(temp.data),
                                            function(y){
                                                index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                                row.names(temp.data[[y]])[index] })))
                                    }
                                })))
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
                            '- Update the selection. Select top ',
                            top.rank.uv.genes,
                            '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                            top.rank.bio.genes,
                            '% of highly affected genes by the bioloigcal variation.'),
                        color = 'blue',
                        verbose = verbose
                    )
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
                            verbose = verbose
                        )
                        lo <- top.rank.bio.genes.nb
                        pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                        while(eval(con) & top.rank.bio.genes.nb > 1){
                            pro.bar$pause(0.1)$tick()$print()
                            # bio genes
                            top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                            if (top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                            top.bio.genes <- unique(unlist(lapply(
                                all.bio.tests,
                                function(x){
                                    if (!is.null(x)){
                                        temp.data <- get(x)
                                        ranks.data <- unique(unlist(lapply(
                                            names(temp.data),
                                            function(y){
                                                index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                                row.names(temp.data[[y]])[index] })))
                                    }
                                })))
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
                            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                            top.bio.genes <- unique(unlist(lapply(
                                all.bio.tests,
                                function(x){
                                    if (!is.null(x)){
                                        temp.data <- get(x)
                                        ranks.data <- unique(unlist(lapply(
                                            names(temp.data),
                                            function(y){
                                                index <- temp.data[[y]]$ranked.genes > top.rank.bio.genes.nb
                                                row.names(temp.data[[y]])[index] })))
                                    }
                                })))
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
                            '- Updating the selection. Select top ',
                            top.rank.uv.genes * 100,
                            '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                            top.rank.bio.genes,
                            '% of highly affected genes by the bioloigcal variation.'),
                        color = 'blue',
                        verbose = verbose)
                }
                ##### grid group: uv ####
                if (grid.group == 'uv'){
                    printColoredMessage(
                        message = '- The grid search will be applied on unwanted factor. ',
                        color = 'blue',
                        verbose = verbose
                    )
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
                            '- Updating the selection. Select top ',
                            top.rank.uv.genes,
                            '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                            top.rank.bio.genes * 100,
                            '% of highly affected genes by the bioloigcal variation.'),
                        color = 'blue',
                        verbose = verbose
                    )
                }
            } else {
                printColoredMessage(
                    message = paste0(
                        '- ',
                        sum(ncg.selected),
                        ' genes are selected as NCGs.'),
                    color = 'blue',
                    verbose = verbose)
            }
        }

        # Filtering selected negative control genes ######
        if (isTRUE(filter.ncgs)){
            printColoredMessage(
                message = '-- Filtering the selected NCGs based on publicly available stable genes:',
                color = 'magenta',
                verbose = verbose
            )
            if (common.hk == 'cancer'){
                printColoredMessage(
                    message = '- Using stable genes proposed for pan-cancer solid tissues:',
                    color = 'blue',
                    verbose = verbose
                )
                common.hks <- singscore::getStableGenes(n_stable = nb.stable.genes)
                common.hks <- intersect(common.hks, row.names(se.obj)[ncg.selected])
                ncg.selected <- row.names(se.obj) %in% common.hks
            }
            if (common.hk == 'non.cancer'){
                printColoredMessage(
                    message = '- Using stable genes proposed for human normal tissues:',
                    color = 'blue',
                    verbose = verbose
                )
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
        # Plotting tha rank of all genes ####
        if (isTRUE(create.ncg.rank.plot)){
            printColoredMessage(
                message = '- Generating a heatmap plot of all the ranks of all the NCGs across all the variables',
                color = 'magenta',
                verbose = verbose
            )
            all.uv.bio.tests <- lapply(
                c(all.bio.tests, all.uv.tests),
                function(x){
                    if (!is.null(x)){
                        temp.data <- get(x)
                        temp.data <- lapply(
                            names(temp.data),
                            function(y) {
                                tm <- temp.data[[y]][ , c('statistic', 'ranked.genes')]
                                tm$group <- rep(y, nrow(tm))
                                tm
                            })
                        temp.data <- do.call(cbind, temp.data)
                        temp.data
                    }
                })
            all.uv.bio.tests <- Filter(Negate(is.null), all.uv.bio.tests)
            all.uv.bio.tests <- do.call(cbind, all.uv.bio.tests)
            temp.data <- lapply(
                seq(3, ncol(all.uv.bio.tests), 3),
                function(x){
                    temp.data <- all.uv.bio.tests[ , x-1, drop = FALSE]
                    colnames(temp.data) <- all.uv.bio.tests[ , x][1]
                    temp.data
                })
            temp.data <- do.call(cbind, temp.data)
            temp.data$ncg <- ncg.selected
            ha <- ComplexHeatmap::rowAnnotation(
                NCG = temp.data$ncg,
                col = list(ncg = c('TRUE' = 'gray10', 'FALSE' = 'gray'))
            )
            ncg.rak.plot <- ComplexHeatmap::Heatmap(
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
            if (isTRUE(plot.ncg.rank)) print(ncg.rak.plot)
        }
        # Assessing of the selected set of NCG ####
        ### Applying PCA using only the NCGs ####
        if (isTRUE(assess.ncg)){
            printColoredMessage(
                message = '-- Assessing the performance of the selected NCG set:',
                color = 'magenta',
                verbose = verbose
            )
            if (is.null(variables.to.assess.ncg)){
                variables.to.assess.ncg <- c(bio.variables, uv.variables)
            }
            assess.ncg.plot <- assessNCGs(
                se.obj = se.obj,
                assay.name = assay.name,
                variables.to.assess.ncg = variables.to.assess.ncg,
                ncg = ncg.selected,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                nb.pcs = nb.pcs,
                svd.bsparam = svd.bsparam,
                center = center,
                scale = scale,
                plot.output = plot.ncg.assessment,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose
            )
        }
        ### Saving the results ####
        printColoredMessage(
            message = '-- Saving the selected NCG to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose
        )
        if (is.null(ncg.set.name)){
            ncg.set.name <- paste0(
                sum(ncg.selected),
                '|',
                paste0(bio.variables, collapse = '&'),
                '|',
                paste0(uv.variables, collapse = '&'),
                '|AnoCorrAs:',
                ncg.selection.method,
                '|',
                assay.name)
        }
        if (is.null(ncg.group.name)){
            ncg.group.name <- 'NcgAcrossAllSamples'
        }
        ### Adding the results to the SummarizedExperiment object ####
        if (is.logical(samples.to.use)){
            se.obj <- se.obj.all
        }
        if (isTRUE(save.se.obj)){
            ## Check if metadata NCG already exists
            if (length(se.obj@metadata$NCG) == 0 ) {
                se.obj@metadata[['NCG']] <- list()
            }
            if (!'supervised' %in% names(se.obj@metadata[['NCG']])){
                se.obj@metadata[['NCG']][['supervised']] <- list()
            }
            if (!ncg.group.name %in% names(se.obj@metadata[['NCG']][['supervised']])){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]] <- list()
            }
            if (!ncg.set.name %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]] )){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]] <- list()
            }
            if (!'ncg.set' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]])){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- list()
            }
            se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- ncg.selected

            if (isTRUE(create.ncg.rank.plot)){
                if (!'rank.plot' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]])){
                    se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['rank.plot']] <- list()
                }
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['rank.plot']] <- ncg.rak.plot
            }
            if (isTRUE(assess.ncg)){
                if (!'rank.plot' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]])){
                    se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- list()
                }
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- assess.ncg.plot
            }
            printColoredMessage(
                message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = '------------The findNcgAcrossSamples function finished.',
                color = 'white',
                verbose = verbose
            )
            return(se.obj)
        }
        ### Exporting the results as logical vector ####
        if (isFALSE(save.se.obj)){
            printColoredMessage(
                message = '-- The NCGs and assessment plot are outputed as a list.',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = '------------The findNcgAcrossSamples function finished.',
                color = 'white',
                verbose = verbose
            )
            if (isTRUE(assess.ncg)){
                return(list(ncg.selected = ncg.selected, assess.ncg.plot = assess.ncg.plot))
            }
            if (isFALSE(assess.ncg)){
                return(list(ncg.selected = ncg.selected))
            }
        }
    }

    # Applying per sample groups ####
    if(approach == 'AnovaCorr.PerBatchPerBiology'){
        printColoredMessage(message = '------------The findNcgByAnovaCorr function starts with "PerBatchPerBiology" mode:',
                            color = 'white',
                            verbose = verbose)

        ### Checking the functions inputs ####
        if (length(assay.name) > 1 | is.logical(assay.name)){
            stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
        }
        if (nb.ncg >= 1 | nb.ncg <= 0){
            stop('The "nb.ncg" should be a positve value  0 < nb.ncg < 1.')
        }
        if (!ncg.selection.method %in% c('prod', 'sum', 'avergae', 'auto', 'non.overlap', 'quantile')){
            stop('The "ncg.selection.method" must be one of "prod", "sum", "avergae", "auto" or "non.overlap".')
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
        if (isFALSE(check.se.obj)) {
            if (isTRUE(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables))) {
                stop('All or some of "uv.variables" cannot be found in the SummarizedExperiment object.')
            } else if (!is.null(variables.to.assess.ncg)) {
                if (isTRUE(sum(variables.to.assess.ncg %in% colnames(colData(se.obj))) != length(variables.to.assess.ncg))) {
                    stop('All or some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
                }
            }
        }
        if (!is.null(regress.out.uv.variables)){
            if (isTRUE(sum(regress.out.uv.variables %in% colnames(colData(se.obj))) != length(regress.out.uv.variables))) {
                stop('All or some of "regress.out.uv.variables" cannot be found in the SummarizedExperiment object.')
            }
        }
        if (!is.null(regress.out.bio.variables)){
            if (isTRUE(sum(regress.out.bio.variables %in% colnames(colData(se.obj))) != length(regress.out.bio.variables))) {
                stop('All or some of "regress.out.bio.variables" cannot be found in the SummarizedExperiment object.')
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

        # Checking the SummarizedExperiment object ####
        if (isTRUE(check.se.obj)) {
            se.obj <- checkSeObj(
                se.obj = se.obj,
                assay.names = assay.name,
                variables = unique(c(bio.variables, uv.variables, bio.groups, uv.groups, variables.to.assess.ncg)),
                remove.na = remove.na,
                verbose = verbose)
        }

        if (remove.na == 'none'){
            if (is.null(variables.to.assess.ncg))
                variables.to.assess.ncg <- c(bio.variables, uv.variables)
            mout <- lapply(
                variables.to.assess.ncg,
                function(x){
                    if (sum(is.na(se.obj[[x]])) > 0)
                        stop('There are NA or missing values in the specified variables.')
                })
        }
        if (is.logical(samples.to.use)){
            se.obj.all <- se.obj
            se.obj <- se.obj[ , samples.to.use]
        }

        # Data transformation and normalization ####
        printColoredMessage(message = '-- Applying data transformation and normalization:',
                            color = 'magenta',
                            verbose = verbose
        )
        ## apply log ####
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
                verbose = verbose)
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log)) {
            printColoredMessage(
                message = paste0(
                    'The ',
                    assay.name,
                    ' data will be used without any log transformation.'),
                color = 'blue',
                verbose = verbose)
            expr.data <- assay(x = se.obj, i = assay.name)
        }

        ## normalization ####
        if (!is.null(normalization)) {
            printColoredMessage(
                message = '-- Applying data normalization:',
                color = 'magenta',
                verbose = verbose
            )
            expr.data.nor <- applyOtherNormalizations(
                se.obj = se.obj,
                assay.name = assay.name,
                method = normalization,
                pseudo.count = pseudo.count,
                apply.log = apply.log,
                check.se.obj = FALSE,
                save.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose)
        }
        # Regressing out variables ####
        if (!is.null(regress.out.uv.variables) | !is.null(regress.out.bio.variables)){
            printColoredMessage(
                message = '-- Regressing out unwanted or biological variables:',
                color = 'magenta',
                verbose = verbose
            )
        }
        ## regress out unwanted variables ####
        if (!is.null(regress.out.uv.variables)){
            printColoredMessage(
                message = '- Regressing out the specified unwanted variables:',
                color = 'blue',
                verbose = verbose
            )
            if (!is.null(normalization)){
                expr.data.reg.uv <- expr.data.nor
            } else expr.data.reg.uv <- expr.data
            printColoredMessage(
                message = paste0(
                    'The ',
                    paste0(regress.out.uv.variables, collapse = ' & '),
                    ' will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    'We do not recommend regressing out ',
                    paste0(regress.out.uv.variables, collapse = ' & '),
                    ' if they are largely associated with the ',
                    paste0(bio.variables, collapse = ' & '),
                    ' variables.'),
                color = 'red',
                verbose = verbose
            )
            expr.data.reg.uv <- t(expr.data.reg.uv)
            uv.variables.all <- paste('se.obj', regress.out.uv.variables, sep = '$')
            expr.data.reg.uv <- lm(as.formula(paste(
                'expr.data.reg.uv',
                paste0(uv.variables.all, collapse = '+') ,
                sep = '~')))
            expr.data.reg.uv <- t(expr.data.reg.uv$residuals)
            colnames(expr.data.reg.uv) <- colnames(se.obj)
            row.names(expr.data.reg.uv) <- row.names(se.obj)
        }

        ## regressing out biological variables ####
        if (!is.null(regress.out.bio.variables)){
            printColoredMessage(
                message = '- Regressing the specified biological variables:',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    paste0(regress.out.bio.variables, collapse = ' & '),
                    ' will be regressed out from the data,',
                    ' please make sure your data is log transformed.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    'We do not recommend regressing out ',
                    paste0(regress.out.bio.variables, collapse = ' & '),
                    ' if they are largely associated with the ',
                    paste0(uv.variables, collapse = ' & '),
                    ' variable(s).'),
                color = 'red',
                verbose = verbose
            )
            expr.data.reg.bio <- t(expr.data)
            bio.variables.all <- paste('se.obj', regress.out.bio.variables, sep = '$')
            expr.data.reg.bio <- lm(as.formula(paste(
                'expr.data.reg.bio',
                paste0(bio.variables.all, collapse = '+') ,
                sep = '~')))
            expr.data.reg.bio <- t(expr.data.reg.bio$residuals)
            colnames(expr.data.reg.bio) <- colnames(se.obj)
            row.names(expr.data.reg.bio) <- row.names(se.obj)
        }

        # Statistical analyses ####
        if (isFALSE(use.imf)){
            printColoredMessage(
                message = '-- Finding a subset of genes as negative control genes:',
                color = 'magenta',
                verbose = verbose
            )
            ## select genes that are highly affected by sources of unwanted variation ####
            printColoredMessage(
                message = '-- Selecting genes that are highly affected by each source(s) of unwanted variation:',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = '- Note, this step will be performed within each homogeneous sample groups with respect to the biological variation.',
                color = 'blue',
                verbose = verbose
            )

            ### create all possible major homogeneous biological groups ####
            printColoredMessage(
                message = '- Creating all possible major homogeneous sample groups with respect to biological variables:',
                color = 'blue',
                verbose = verbose
            )
            if (is.null(bio.groups)){
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(bio.variables, collapse = ' & '),
                        ' variable(s) will be used to create all possible major homogeneous biological groups.'),
                    color = 'blue',
                    verbose = verbose
                )
                all.bio.groups <- createHomogeneousBioGroups(
                    se.obj = se.obj,
                    bio.variables = bio.variables,
                    nb.clusters = nb.bio.clusters,
                    clustering.method = bio.clustering.method,
                    check.se.obj = FALSE,
                    save.se.obj = FALSE,
                    remove.na = 'none',
                    verbose = verbose
                )
            } else if (!is.null(bio.groups)) {
                printColoredMessage(
                    message = paste0(
                        '- The',
                        paste0(bio.groups, collapse = ' & '),
                        ' variables will be used to create all possible major homogeneous biological groups.'),
                    color = 'blue',
                    verbose = verbose
                )
                all.bio.groups <- createHomogeneousBioGroups(
                    se.obj = se.obj,
                    bio.variables = bio.groups,
                    nb.clusters = nb.bio.clusters,
                    clustering.method = bio.clustering.method,
                    check.se.obj = FALSE,
                    save.se.obj = FALSE,
                    remove.na = 'none',
                    verbose = verbose)
            }

            ### correlation between gene expression and all continuous source of unwanted variation with each biological groups ####
            uv.var.class <- unlist(lapply(
                uv.variables,
                function(x) class(colData(se.obj)[[x]]))
            )
            continuous.uv <- uv.variables[uv.var.class %in% c('numeric', 'integer')]
            if (isTRUE(length(continuous.uv) > 0)) {
                printColoredMessage(
                    message = paste0(
                        '-- Performing correlation analysis between gene expression and ',
                        'each specified continuous sources of unwanted variation:'),
                    color = 'blue',
                    verbose = verbose
                )
                selected.bio.groups <- findRepeatingPatterns(
                    vec = all.bio.groups,
                    n.repeat = min.sample.for.correlation
                )
                if (isTRUE(length(selected.bio.groups) > 0)) {
                    if (is.null(regress.out.bio.variables)) {
                        data.to.use <- expr.data
                    } else data.to.use <- expr.data.reg.bio
                    corr.genes.uv <- lapply(
                        continuous.uv,
                        function(x) {
                            all.corr <- lapply(
                                selected.bio.groups,
                                function(y) {
                                    selected.samples <- all.bio.groups == y
                                    corr.genes <- as.data.frame(correls(
                                        x = t(data.to.use[, selected.samples]),
                                        y = se.obj@colData[[x]][selected.samples],
                                        type = corr.method,
                                        a = a ,
                                        rho = rho
                                    ))
                                    corr.genes <- cbind(
                                        round(x = corr.genes[, 1:4], digits = 3),
                                        corr.genes[, 5, drop = FALSE]
                                    )
                                    corr.genes$g.statistic <- abs(corr.genes$correlation)
                                    set.seed(2233)
                                    corr.genes$ranked.genes <- rank(
                                        -abs(corr.genes[, 'correlation']),
                                        ties.method = 'random'
                                    )
                                    row.names(corr.genes) <- row.names(data.to.use)
                                    corr.genes
                                })
                            names(all.corr) <- selected.bio.groups
                            all.corr
                        })
                    names(corr.genes.uv) <- continuous.uv
                } else if (isTRUE(length(selected.bio.groups) == 0))
                    stop(paste0(
                        'There are not homogeneous biological groups that have at least ',
                        min.sample.for.correlation,
                        ' (min.sample.for.correlation) samples for correlation analysis.')
                    )
            } else corr.genes.uv <- NULL

            ### anova between gene expression and all categorical source of variation within each biological groups ####
            categorical.uv <- uv.variables[uv.var.class %in% c('factor', 'character')]
            if (isTRUE(length(categorical.uv) > 0)) {
                printColoredMessage(
                    message = paste0(
                        '-- Performing ANOVA between individual gene expression and each ',
                        'specified categorical sources of unwanted variation:'),
                    color = 'blue',
                    verbose = verbose
                )
                anova.genes.uv <- lapply(
                    categorical.uv,
                    function(x) {
                        bio.batch <- table(all.bio.groups, colData(se.obj)[[x]])
                        cover.sample.groups <- rowSums(bio.batch >= min.sample.for.aov) == length(unique(se.obj[[x]]))
                        if (isTRUE(sum(cover.sample.groups) > 0)) {
                            printColoredMessage(
                                message = paste0(
                                    sum(cover.sample.groups),
                                    ' homogeneous biological group(s) have at least ',
                                    min.sample.for.aov,
                                    ' (min.sample.for.aov) samples within individual batches of the ',
                                    x,
                                    ' variable.'),
                                color = 'blue',
                                verbose = verbose)
                        }
                        if (isTRUE(sum(cover.sample.groups) == 0)){
                            printColoredMessage(
                                message = paste0(
                                    'There are not homogeneous biological groups that have at least ',
                                    min.sample.for.aov ,
                                    ' (min.sample.for.aov) samples within each batches of the ',
                                    x,
                                    ' variable. This may result in unsatisfactory NCG selection.'),
                                color = 'red',
                                verbose = verbose)
                        }
                        selected.bio.groups <- names(which(rowSums(bio.batch >= min.sample.for.aov) > 1))
                        if (isTRUE(length(selected.bio.groups) == 0)){
                            stop('There is not enough groups to perform ANOVA.')
                        }
                        if (is.null(regress.out.bio.variables)) {
                            data.to.use <- expr.data
                        } else data.to.use <- expr.data.reg.bio
                        all.anova <- lapply(
                            selected.bio.groups,
                            function(i) {
                                selected.samples <- all.bio.groups == i
                                if (anova.method == 'aov') {
                                    anova.genes.batch <- as.data.frame(row_oneway_equalvar(
                                        x = data.to.use[, selected.samples],
                                        g = se.obj@colData[, x][selected.samples])
                                        )
                                } else if (anova.method == 'welch.correction') {
                                    anova.genes.batch <- as.data.frame(row_oneway_welch(
                                        x = data.to.use[, selected.samples],
                                        g = se.obj@colData[, x][selected.samples])
                                        )
                                }
                                anova.genes.batch$p.adjusted <- p.adjust(anova.genes.batch$pvalue, method = "BH")
                                anova.genes.batch$eta.squared <- (anova.genes.batch$statistic * anova.genes.batch$df.between) /
                                    (anova.genes.batch$statistic * anova.genes.batch$df.between + anova.genes.batch$df.within)
                                anova.genes.batch$g.statistic <- anova.genes.batch$eta.squared
                                set.seed(2233)
                                anova.genes.batch$ranked.genes <- rank(
                                    -anova.genes.batch[, 'statistic'],
                                    ties.method = 'random'
                                    )
                                anova.genes.batch
                            })
                        names(all.anova) <- selected.bio.groups
                        all.anova
                    })
                names(anova.genes.uv) <- categorical.uv
            } else anova.genes.uv <- NULL

            ## select genes that are highly affected by biological variation  ####
            printColoredMessage(
                message = '-- Selecting genes that are highly affected by each specified source(s) of biological variation:',
                color = 'blue',
                verbose = verbose
                )
            printColoredMessage(
                message = '- This step will be performed within each possible homogeneous unwanted groups.',
                color = 'blue',
                verbose = verbose
                )
            ### create all possible homogeneous uv groups ####
            printColoredMessage(
                message = '- Creating all possible homogeneous sample groups with respect to the specified unwanted variables:',
                color = 'blue',
                verbose = verbose
                )
            if (is.null(uv.groups)){
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(uv.variables, collapse = ' & '),
                        ' variables will be used as a major sources of unwanted variation',
                        ' to find all possible groups.'),
                    color = 'blue',
                    verbose = verbose
                    )
                all.uv.groups <- createHomogeneousUVGroups(
                    se.obj = se.obj,
                    uv.variables = uv.variables,
                    nb.clusters = nb.uv.clusters,
                    clustering.method = uv.clustering.method,
                    check.se.obj = FALSE,
                    save.se.obj = FALSE,
                    verbose = verbose
                )
            } else if (!is.null(uv.groups)){
                printColoredMessage(
                    message = paste0(
                        'The ',
                        paste0(uv.groups, collapse = ' & '),
                        ' variables will be used as a major sources of unwanted variation',
                        ' to find all possible groups.'),
                    color = 'blue',
                    verbose = verbose
                )
                all.uv.groups <- createHomogeneousUVGroups(
                    se.obj = se.obj,
                    uv.variables = uv.groups,
                    nb.clusters = nb.uv.clusters,
                    clustering.method = uv.clustering.method,
                    check.se.obj = FALSE,
                    save.se.obj = FALSE,
                    verbose = verbose)
            }
            ### correlation between gene expression and all continuous source of biological variation with each uv groups ####
            bio.var.class <- unlist(lapply(
                bio.variables,
                function(x) class(colData(se.obj)[[x]]))
                )
            continuous.bio <- bio.variables[bio.var.class %in% c('numeric', 'integer')]
            if (length(continuous.bio) > 0) {
                printColoredMessage(
                    message = '-- Correlation analyses:',
                    color = 'magenta',
                    verbose = verbose
                )
                selected.uv.groups <- findRepeatingPatterns(
                    vec = all.uv.groups,
                    n.repeat = min.sample.for.correlation
                    )
                if (length(selected.uv.groups) > 0) {
                    if (length(selected.uv.groups) == 1) {
                        group = 'group has'
                    } else group = 'groups have'
                    printColoredMessage(
                        message = paste0(
                            length(selected.uv.groups),
                            ' homogeneous groups with respect to the sources of unwanted variation ',
                            group,
                            ' at least ',
                            min.sample.for.correlation,
                            ' (min.sample.for.correlation) samples to pefrom correlation between gene-level',
                            'expression and all the continuous sources of bioloical variation.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    if (is.null(regress.out.uv.variables) & is.null(normalization)) {
                        data.to.use <- expr.data
                    } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)) {
                        data.to.use <- expr.data.reg.uv
                    } else if (is.null(regress.out.uv.variables) & !is.null(normalization)) {
                        data.to.use <- expr.data.nor
                    }
                    corr.genes.bio <- lapply(
                        continuous.bio,
                        function(x) {
                            all.corr <- lapply(
                                selected.uv.groups,
                                function(y) {
                                    selected.samples <- all.uv.groups == y
                                    corr.genes <- as.data.frame(correls(
                                        y = se.obj@colData[, x][selected.samples],
                                        x = t(data.to.use[, selected.samples]),
                                        type = corr.method,
                                        a = a ,
                                        rho = rho
                                    ))
                                    corr.genes <- cbind(
                                        round(x = corr.genes[, 1:4], digits = 3),
                                        corr.genes[, 5, drop = FALSE]
                                        )
                                    corr.genes$g.statistic <- abs(corr.genes$correlation)
                                    set.seed(2233)
                                    corr.genes$ranked.genes <- rank(
                                        abs(corr.genes[, 'correlation']),
                                        ties.method = 'random'
                                    )
                                    row.names(corr.genes) <- row.names(data.to.use)
                                    corr.genes
                                })
                            names(all.corr) <- selected.uv.groups
                            all.corr
                        })
                    names(corr.genes.bio) <- continuous.bio
                } else {
                    stop(
                        paste0(
                            'There are not homogeneous groups with respect to sources of unwanted variation that have at least ',
                            min.sample.for.correlation,
                            ' (min.sample.for.correlation) samples for correlation analysis between gene-level expression and all',
                            ' the continuous sources of bioloical variation.'))
                }
            } else corr.genes.bio <- NULL

            ### anova between gene expression and all categorical source of biological variation with each uv groups ####
            categorical.bio <- bio.variables[bio.var.class %in% c('factor', 'character')]
            if (length(categorical.bio) > 0) {
                printColoredMessage(
                    message = '-- ANOV analyses:',
                    color = 'magenta',
                    verbose = verbose
                )
                anova.genes.bio <- lapply(
                    categorical.bio,
                    function(x) {
                        bio.batch <- table(all.uv.groups, colData(se.obj)[[x]])
                        cover.sample.groups <- rowSums(bio.batch >= min.sample.for.aov) == length(unique(se.obj[[x]]))
                        if (isTRUE(sum(cover.sample.groups) > 0)) {
                            printColoredMessage(
                                message = paste0(
                                    sum(cover.sample.groups),
                                    ' homogeneous unwanted group(s) have at least ',
                                    min.sample.for.aov,
                                    ' (min.sample.for.aov) samples within individual groups of the ', x, ' variable.'),
                                color = 'blue',
                                verbose = verbose)
                        }
                        if (isTRUE(sum(cover.sample.groups) == 0)){
                            printColoredMessage(
                                message = paste0(
                                    'There are not homogeneous unwanted groups that have at least ',
                                    min.sample.for.aov , ' (min.sample.for.aov) samples within each batches of the ',
                                    x,
                                    ' variable. This may result in unsatisfactory NCG selection.'),
                                color = 'red',
                                verbose = verbose)
                        }
                        selected.uv.groups <- names(which(rowSums(bio.batch >= min.sample.for.aov) > 1))
                        if (length(selected.uv.groups) == 0) {
                            stop(paste0(
                                'It seems there is complete association between ',x,
                                ' homogeneous groups with respect to unwanted variation.'))
                        }
                        if (is.null(regress.out.uv.variables) &is.null(normalization)) {
                            data.to.use <- expr.data
                        } else if (!is.null(regress.out.uv.variables) & !is.null(normalization)) {
                            data.to.use <- expr.data.reg.uv
                        } else if (is.null(regress.out.uv.variables) & !is.null(normalization)) {
                            data.to.use <- expr.data.nor
                        }
                        all.anova <- lapply(
                            selected.uv.groups,
                            function(i) {
                                selected.samples <- all.uv.groups == i
                                if (anova.method == 'aov') {
                                    anova.gene.bio <- as.data.frame(row_oneway_equalvar(
                                        x = data.to.use[, selected.samples],
                                        g = se.obj@colData[, x][selected.samples])
                                        )
                                }
                                if (anova.method == 'welch.correction') {
                                    anova.gene.bio <- as.data.frame(row_oneway_welch(
                                        x = data.to.use[, selected.samples],
                                        g = se.obj@colData[, x][selected.samples])
                                        )
                                }
                                anova.gene.bio$p.adjusted <- p.adjust(anova.gene.bio$pvalue, method = "BH")
                                anova.gene.bio$eta.squared <- (anova.gene.bio$statistic * anova.gene.bio$df.between) /
                                    (anova.gene.bio$statistic * anova.gene.bio$df.between + anova.gene.bio$df.within)
                                anova.gene.bio$g.statistic <- anova.gene.bio$eta.squared
                                set.seed(2233)
                                anova.gene.bio$ranked.genes <- rank(
                                    anova.gene.bio[, 'statistic'],
                                    ties.method = 'random'
                                    )
                                anova.gene.bio
                            })
                        names(all.anova) <- selected.uv.groups
                        all.anova
                    })
                names(anova.genes.bio) <- categorical.bio
            } else anova.genes.bio <- NULL
        }

        # Intermediate file ####
        ## read intermediate file ####
        if (isTRUE(use.imf)){
            if (is.null(imf.name)){
                imf.name <- paste0(assay.name, '|PerBiologyPerBatch|', ncg.selection.method)
            }
            if (is.null(se.obj@metadata$IMF$NCG[[imf.name]]))
                stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
            all.tests <- se.obj@metadata$IMF$NCG[[imf.name]]
            anova.genes.bio <- all.tests$anova.genes.bio
            corr.genes.bio <- all.tests$corr.genes.bio
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
                imf.name <- paste0(assay.name, '|PerBiologyPerBatch|', ncg.selection.method)
            }
            if (!imf.name %in% names(se.obj@metadata[['IMF']][['NCG']])){
                se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list()
            }
            se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list(
                anova.genes.bio = anova.genes.bio,
                corr.genes.bio = corr.genes.bio,
                anova.genes.uv = anova.genes.uv,
                corr.genes.uv = corr.genes.uv)
        }

        # Selection of NCG ####
        printColoredMessage(message = '-- Selection a set of genes as NCG:',
                            color = 'magenta',
                            verbose = verbose)

        ## product, sum or average of ranks ####
        if (isFALSE(use.rank)){
            all.tests <- c(
                'anova.genes.bio',
                'corr.genes.bio',
                'anova.genes.uv',
                'corr.genes.uv'
            )
            var.partition <- lapply(
                all.tests,
                function(x) {
                    if (!is.null(x)) {
                        temp.data <- get(x)
                        ranks.data <- lapply(
                            names(temp.data),
                            function(y) {
                                all.ranks <- lapply(
                                    names(temp.data[[y]]),
                                    function(i) temp.data[[y]][[i]]$g.statistic
                                    )
                                names(all.ranks) <- paste0(y , '||', names(temp.data[[y]]))
                                all.ranks <- do.call(cbind, all.ranks)
                                all.ranks})
                        ranks.data <- do.call(cbind, ranks.data)
                        names(ranks.data) <- names(temp.data)
                        ranks.data
                    }
                })
            var.partition <- as.data.frame(do.call(cbind, var.partition))
            row.names(var.partition) <- row.names(se.obj)
            colnames(var.partition) <- sub("\\|\\|.*", "", colnames(var.partition))
            matched.cols.bio <- grep(paste(bio.variables, collapse = "|"), colnames(var.partition))
            matched.cols.uv <- grep(paste(uv.variables, collapse = "|"), colnames(var.partition))

            var.partition$bio <- rowSums(var.partition[ , matched.cols.bio, drop = FALSE])
            var.partition$uv <- rowSums(var.partition[ , matched.cols.uv, drop = FALSE])
            var.partition$ratio <- var.partition$uv / c(var.partition$bio + 1e-6)
            var.partition <- var.partition[order(-var.partition$ratio), ]

            if (ncg.selection.method == 'quantile'){
                if (!is.null(uv.percentile) & !is.null(bio.percentile)){
                    var.partition <- var.partition[var.partition$uv > quantile(x = var.partition$uv, probs = uv.percentile, na.rm = TRUE) , ]
                    if (quantile(x = var.partition$bio, probs = bio.percentile) == 0){
                        var.partition <- var.partition
                    } else {
                        var.partition <- var.partition[var.partition$bio < quantile(x = var.partition$bio, probs = bio.percentile, na.rm = TRUE) , ]
                    }
                    nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                    ncg.selected <- row.names(se.obj) %in% row.names(var.partition)[1:nb.ncg]
                } else if (!is.null(uv.percentile) & is.null(bio.percentile)){
                    var.partition <- var.partition[var.partition$uv > quantile(x = var.partition$uv, probs = uv.percentile, na.rm = TRUE) , ]
                    nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                    ncg.selected <- row.names(se.obj) %in% row.names(var.partition)[1:nb.ncg]
                } else if (is.null(uv.percentile) & !is.null(bio.percentile)){
                    var.partition <- var.partition[var.partition$bio < quantile(x = var.partition$bio, probs = bio.percentile, na.rm = TRUE) , ]
                    nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                    ncg.selected <- row.names(se.obj) %in% row.names(var.partition)[1:nb.ncg]
                } else if (is.null(uv.percentile) & is.null(bio.percentile)){
                    nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
                    ncg.selected <- row.names(se.obj) %in% row.names(var.partition)[1:nb.ncg]
                }
            }

        }
        if (ncg.selection.method %in% c('prod', 'sum', 'average') & isTRUE(use.rank)) {
            all.tests <-c(
                'anova.genes.bio',
                'corr.genes.bio',
                'anova.genes.uv',
                'corr.genes.uv'
                )
            all.stats <- lapply(
                all.tests,
                function(x) {
                    if (!is.null(x)) {
                        temp.data <- get(x)
                        ranks.data <-
                            lapply(
                                names(temp.data),
                                function(y) {
                                    all.ranks <- lapply(
                                        names(temp.data[[y]]),
                                        function(i) temp.data[[y]][[i]]$ranked.genes)
                                    names(all.ranks) <- names(temp.data[[y]])
                                    all.ranks <- do.call(cbind, all.ranks)
                                    all.ranks
                                })
                        ranks.data <- do.call(cbind, ranks.data)
                        names(ranks.data) <- names(temp.data)
                        ranks.data
                    }
                })
            all.stats <- as.data.frame(do.call(cbind, all.stats))
            row.names(all.stats) <- row.names(se.obj)
            ### product of ranks ####
            if (ncg.selection.method == 'prod') {
                printColoredMessage(
                    message = 'A set of NCG will be selected based on the product of ranks.',
                    color = 'blue',
                    verbose = verbose
                )
                all.stats$all.rank <- 10 ^ (rowSums(log(all.stats, base = 10)))
                if (sum(is.infinite(all.stats$all.rank)) > 0) {
                    stop('The product of ranks results in infinity values.')
                }
            }
            ## sum of ranks ####
            if (ncg.selection.method == 'sum') {
                printColoredMessage(
                    message = 'A set of NCG will be selected based on the sum of ranks.',
                    color = 'blue',
                    verbose = verbose)
                all.stats$all.rank <- rowSums(x = all.stats, na.rm = TRUE)
            }
            ## average of ranks ####
            if (ncg.selection.method == 'average') {
                printColoredMessage(
                    message = 'A set of NCG will be selected based on the average of ranks.',
                    color = 'blue',
                    verbose = verbose
                )
                all.stats$all.rank <- rowMeans(x = all.stats, na.rm = TRUE)
            }

            all.stats <- all.stats[order(all.stats$all.rank, decreasing = FALSE), ]
            ncg.selected <- row.names(all.stats)[1:round(nb.ncg * nrow(se.obj), digits = 0)]
            ncg.selected <- row.names(se.obj) %in% ncg.selected
        }

        ## non.overlap approach ####
        if (ncg.selection.method == 'non.overlap' & isTRUE(use.rank) ) {
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the "non.overlap" approach.',
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    '- Selecting top ',
                    top.rank.uv.genes * 100,
                    '% of highly affected genes by the unwanted variation, and then exclude all top ',
                    top.rank.bio.genes *100,
                    '% of highly affected genes by the bioloigcal variation.'),
                color = 'blue',
                verbose = verbose
            )
            ### select genes affected by biological variation ####
            top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
            top.bio.genes <- unique(unlist(lapply(
                all.bio.tests,
                function(x) {
                    if (isTRUE(!is.null(x))) {
                        temp.data <- get(x)
                        ranks.data <- unique(unlist(lapply(
                            names(temp.data),
                            function(y) {
                                all.ranks <- sapply(
                                    names(temp.data[[y]]),
                                    function(z) temp.data[[y]][[z]]$ranked.genes)
                                set.seed(2233)
                                all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                            })))
                    }
                })))
            ## select genes affected by unwanted variation ####
            top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
            all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
            top.uv.genes <- unique(unlist(lapply(
                all.uv.tests,
                function(x) {
                    temp <- get(x)
                    if (length(names(temp)) != 0) {
                        ranks.data <- lapply(
                            names(temp),
                            function(y) {
                                unlist(lapply(names(temp[[y]]),
                                              function(z) {
                                                  index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                  row.names(temp[[y]][[z]])[index]
                                              }))
                            })
                    }
                })))
            ## select of NCGS ####
            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
            if (isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
            ncg.selected <- row.names(se.obj) %in% ncg.selected
        }

        ## auto approach ####
        if (isTRUE(ncg.selection.method == 'auto') & isTRUE(use.rank)){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the "auto" approach.',
                color = 'blue',
                verbose = verbose)
            printColoredMessage(
                message = paste0(
                    '- Selecting top ',
                    top.rank.uv.genes * 100,
                    '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                    top.rank.bio.genes * 100,
                    '% of highly affected genes by the bioloigcal variation.'),
                color = 'blue',
                verbose = verbose
            )
            ### find highly affected genes by biology ####
            top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
            top.bio.genes <- unique(unlist(lapply(
                all.bio.tests,
                function(x) {
                    if (isTRUE(!is.null(x))) {
                        temp.data <- get(x)
                        ranks.data <- unique(unlist(lapply(
                            names(temp.data),
                            function(y) {
                                all.ranks <- sapply(
                                    names(temp.data[[y]]),
                                    function(z) temp.data[[y]][[z]]$ranked.genes)
                                set.seed(2233)
                                all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                            })))
                    }
                })))
            ## select genes affected by unwanted variation ####
            top.rank.uv.genes.nb <- round(c(top.rank.uv.genes * nrow(se.obj)), digits = 0)
            all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
            top.uv.genes <- unique(unlist(lapply(
                all.uv.tests,
                function(x) {
                    temp <- get(x)
                    if (length(names(temp)) != 0) {
                        ranks.data <- lapply(
                            names(temp),
                            function(y) {
                                unlist(lapply(names(temp[[y]]),
                                              function(z) {
                                                  index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                  row.names(temp[[y]][[z]])[index]
                                              }))
                            })
                    }
                })))
            ## select NCG ####
            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
            # if (isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
            printColoredMessage(
                message = paste0(
                    '- ',
                    length(ncg.selected),
                    ' genes are found.'),
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
                        verbose = verbose
                    )
                }
                if (isTRUE(nb.ncg < length(ncg.selected))){
                    con <- parse(text = paste0("length(ncg.selected)", ">", "nb.ncg"))
                    printColoredMessage(
                        message = paste0(
                            '- The number of selected genes ',
                            length(ncg.selected),
                            ' is larger than the number (',
                            nb.ncg ,
                            ') of specified genes ',
                            'by "nb.ncg". A grid search will be performed.'),
                        color = 'blue',
                        verbose = verbose)
                }
                ## grid search ####
                ### grid group: both bio and uv variable ####
                if (grid.group == 'both'){
                    #### increasing order ####
                    if (grid.direction == 'increase'){
                        lo <- min(
                            nrow(se.obj) - top.rank.uv.genes.nb,
                            top.rank.bio.genes.nb
                        )
                        pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                        while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj) & top.rank.bio.genes.nb > 1){
                            pro.bar$pause(0.1)$tick()$print()
                            # uv
                            all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                            top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                            if (top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                            top.uv.genes <- unique(unlist(lapply(
                                all.uv.tests,
                                function(x) {
                                    temp <- get(x)
                                    if (length(names(temp)) != 0) {
                                        ranks.data <- lapply(
                                            names(temp),
                                            function(y) {
                                                unlist(lapply(names(temp[[y]]),
                                                              function(z) {
                                                                  index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                                  row.names(temp[[y]][[z]])[index]
                                                              }))
                                            })
                                    }
                                })))
                            # bio
                            top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                            if (top.rank.bio.genes.nb < 0) top.rank.bio.genes.nb = 1
                            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                            top.bio.genes <- unique(unlist(lapply(
                                all.bio.tests,
                                function(x) {
                                    if (isTRUE(!is.null(x))) {
                                        temp.data <- get(x)
                                        ranks.data <- unique(unlist(lapply(
                                            names(temp.data),
                                            function(y) {
                                                all.ranks <- sapply(
                                                    names(temp.data[[y]]),
                                                    function(z) temp.data[[y]][[z]]$ranked.genes)
                                                set.seed(2233)
                                                all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                                row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                                            })))
                                    }
                                })))
                            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                        }
                        if (length(ncg.selected) == 0)
                            stop('No NCGs can be found based on the current parameters.')
                    }
                    ### decreasing order ####
                    if (grid.direction == 'decrease'){
                        lo <- min(top.rank.uv.genes.nb, c(nrow(se.obj) - top.rank.bio.genes.nb))
                        pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                        while(eval(con) & top.rank.uv.genes.nb > 1 & top.rank.bio.genes.nb < nrow(se.obj)){
                            pro.bar$pause(0.1)$tick()$print()
                            # uv
                            all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                            top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                            if (top.rank.uv.genes.nb < 0 ) top.rank.uv.genes.nb = 1
                            top.uv.genes <- unique(unlist(lapply(
                                all.uv.tests,
                                function(x) {
                                    temp <- get(x)
                                    if (length(names(temp)) != 0) {
                                        ranks.data <- lapply(
                                            names(temp),
                                            function(y) {
                                                unlist(lapply(names(temp[[y]]),
                                                              function(z) {
                                                                  index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                                  row.names(temp[[y]][[z]])[index]
                                                              }))
                                            })
                                    }
                                })))
                            # bio
                            top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                            if (top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                            top.bio.genes <- unique(unlist(lapply(
                                all.bio.tests,
                                function(x) {
                                    if (isTRUE(!is.null(x))) {
                                        temp.data <- get(x)
                                        ranks.data <- unique(unlist(lapply(
                                            names(temp.data),
                                            function(y) {
                                                all.ranks <- sapply(
                                                    names(temp.data[[y]]),
                                                    function(z) temp.data[[y]][[z]]$ranked.genes)
                                                set.seed(2233)
                                                all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                                row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                                            })))
                                    }
                                })))
                            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                        }
                        if (length(ncg.selected) == 0)
                            stop('No NCGs can be found based on the current parameters.')
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
                    ###### increasing order ####
                    if (grid.direction == 'increase'){
                        lo <- top.rank.bio.genes.nb
                        pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                        while(eval(con) & top.rank.bio.genes.nb > 1){
                            pro.bar$pause(0.1)$tick()$print()
                            top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                            if (top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                            top.bio.genes <- unique(unlist(lapply(
                                all.bio.tests,
                                function(x) {
                                    if (isTRUE(!is.null(x))) {
                                        temp.data <- get(x)
                                        ranks.data <- unique(unlist(lapply(
                                            names(temp.data),
                                            function(y) {
                                                all.ranks <- sapply(
                                                    names(temp.data[[y]]),
                                                    function(z) temp.data[[y]][[z]]$ranked.genes)
                                                set.seed(2233)
                                                all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                                row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                                            })))
                                    }
                                })))
                            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                            ncg.selected
                        }
                    }
                    ##### decreasing order ####
                    if (grid.direction == 'decrease'){
                        lo <- nrow(se.obj) - top.rank.bio.genes.nb
                        pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                        while(eval(con) & top.rank.bio.genes.nb < nrow(se.obj)){
                            pro.bar$pause(0.1)$tick()$print()
                            top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                            if (top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                            all.bio.tests <- c('anova.genes.bio', 'corr.genes.bio')
                            top.bio.genes <- unique(unlist(lapply(
                                all.bio.tests,
                                function(x) {
                                    if (isTRUE(!is.null(x))) {
                                        temp.data <- get(x)
                                        ranks.data <- unique(unlist(lapply(
                                            names(temp.data),
                                            function(y) {
                                                all.ranks <- sapply(
                                                    names(temp.data[[y]]),
                                                    function(z) temp.data[[y]][[z]]$ranked.genes)
                                                set.seed(2233)
                                                all.ranks <- rank(x = rowMeans(all.ranks), ties.method = 'random')
                                                row.names(se.obj)[all.ranks > top.rank.bio.genes.nb]
                                            })))
                                    }
                                })))
                        }
                        ##### check selection ####
                        if (length(ncg.selected) == 0) stop('NCGs cannot be found based on the current parameters.')
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
                            verbose = verbose)
                    }
                }
                ##### grid group: bio ####
                if (grid.group == 'uv'){
                    printColoredMessage(
                        message = '- The grid search will be applied on unwanted factor. ',
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
                            all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                            top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                            if (top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                            top.uv.genes <- unique(unlist(lapply(
                                all.uv.tests,
                                function(x) {
                                    temp <- get(x)
                                    if (length(names(temp)) != 0) {
                                        ranks.data <- lapply(
                                            names(temp),
                                            function(y) {
                                                unlist(lapply(names(temp[[y]]),
                                                              function(z) {
                                                                  index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                                  row.names(temp[[y]][[z]])[index]
                                                              }))
                                            })
                                    }
                                })))
                            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                            ncg.selected
                        }
                    }
                    ##### decreasing order ####
                    if (grid.direction == 'decrease'){
                        printColoredMessage(
                            message = '- The grid search will decrease the value of "top.rank.uv.genes". ',
                            color = 'blue',
                            verbose = verbose)
                        lo <- top.rank.uv.genes.nb
                        pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                        while(eval(con) & top.rank.uv.genes.nb > 1){
                            pro.bar$pause(0.1)$tick()$print()
                            all.uv.tests <- c('anova.genes.uv', 'corr.genes.uv')
                            top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                            if (top.rank.uv.genes.nb < 0 ) top.rank.uv.genes.nb = 1
                            top.uv.genes <- unique(unlist(lapply(
                                all.uv.tests,
                                function(x) {
                                    temp <- get(x)
                                    if (length(names(temp)) != 0) {
                                        ranks.data <- lapply(
                                            names(temp),
                                            function(y) {
                                                unlist(lapply(names(temp[[y]]),
                                                              function(z) {
                                                                  index <- temp[[y]][[z]]$ranked.genes < top.rank.uv.genes.nb
                                                                  row.names(temp[[y]][[z]])[index]
                                                              }))
                                            })
                                    }
                                })))
                            ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                            ncg.selected
                        }
                    }
                    ##### check selection ####
                    if (length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
                    ncg.selected <- row.names(se.obj) %in% ncg.selected
                    ##### update numbers ####
                    # uv
                    top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 0)
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
                }
            }
        }
        # Filtering selected negative control genes ######
        if (isTRUE(filter.ncgs)){
            printColoredMessage(
                message = '-- Filtering the selected NCGs based on publicly available stable genes:',
                color = 'magenta',
                verbose = verbose
                )
            if (common.hk == 'cancer'){
                printColoredMessage(
                    message = '- Using stable genes proposed for pan-cancer solid tissues:',
                    color = 'blue',
                    verbose = verbose
                    )
                common.hks <- singscore::getStableGenes(n_stable = nb.stable.genes)
                common.hks <- intersect(common.hks, row.names(se.obj)[ncg.selected])
                ncg.selected <- row.names(se.obj) %in% common.hks
            }
            if (common.hk == 'non.cancer'){
                printColoredMessage(
                    message = '- Using stable genes proposed for human normal tissues:',
                    color = 'blue',
                    verbose = verbose
                )
                common.hks <- row.names(se.obj)[rowData(se.obj)[[hk.group]]]
                common.hks <- intersect(common.hks, row.names(se.obj)[ncg.selected])
                ncg.selected <- row.names(se.obj) %in% common.hks
            }
        }

        printColoredMessage(
            message = paste0(
                '- A set of ',
                sum(ncg.selected),
                ' genes are selected as NCG.'),
            color = 'blue',
            verbose = verbose
        )
        # Plotting ####
        if (isTRUE(create.ncg.rank.plot)){
            printColoredMessage(
                message = '- Generating a heatmap plot of all the ranks of all the NCGs across all the variables',
                color = 'magenta',
                verbose = verbose
            )
            all.uv.bio.tests <- lapply(
                c(all.bio.tests, all.uv.tests),
                function(x){
                    if (!is.null(x)){
                        temp.data <- get(x)
                        temp.data <- lapply(
                            names(temp.data),
                            function(y) {
                                tests.data <- sapply(
                                    names(temp.data[[y]]),
                                    function(z){
                                        tm <- temp.data[[y]][[z]][ , 'ranked.genes']
                                    })
                                tests.data <- rowMeans(tests.data)
                                tests.data <- data.frame(
                                    ranks = tests.data,
                                    group = rep(y, length(tests.data))
                                )
                                tests.data
                            })
                        temp.data <- do.call(cbind , temp.data)
                        temp.data
                    }
                })
            all.uv.bio.tests <- Filter(Negate(is.null), all.uv.bio.tests)
            all.uv.bio.tests <- do.call(cbind, all.uv.bio.tests)
            temp.data <- lapply(
                seq(1, ncol(all.uv.bio.tests), 2),
                function(x){
                    temp.data <- all.uv.bio.tests[ , x, drop = FALSE]
                    colnames(temp.data) <- all.uv.bio.tests[ , x+1][1]
                    temp.data
                })
            temp.data <- do.call(cbind, temp.data)
            temp.data$ncg <- ncg.selected
            ha <- ComplexHeatmap::rowAnnotation(
                NCG = temp.data$ncg,
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

        # Assessing the performance of selected set of NCG  ####
        ## Performing PCA on only selected genes as NCG ####
        if (isTRUE(assess.ncg)) {
            printColoredMessage(
                message = '-- Assessing the performance of selected NCG set:',
                color = 'magenta',
                verbose = verbose
            )
            if (is.null(variables.to.assess.ncg)){
                variables.to.assess.ncg <- c(bio.variables, uv.variables)
            }
            assess.ncg.plot <- assessNCGs(
                se.obj = se.obj,
                assay.name = assay.name,
                variables.to.assess.ncg = variables.to.assess.ncg,
                ncg = ncg.selected,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                nb.pcs = nb.pcs,
                svd.bsparam = svd.bsparam,
                center = center,
                scale = scale,
                plot.output = plot.ncg.assessment,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose
            )
        }
        # Saving results ####
        ## Adding results to the SummarizedExperiment object ####
        if (is.null(ncg.set.name)){
            ncg.set.name <- paste0(
                sum(ncg.selected),
                '|',
                paste0(bio.variables, collapse = '&'),
                '|',
                paste0(uv.variables, collapse = '&'),
                '|PbPbio:',
                ncg.selection.method,
                '|',
                assay.name)
        }
        if (is.null(ncg.group.name)){
            ncg.group.name <- 'NcgPerBioPerBatch'
        }
        ### Adding the results to the SummarizedExperiment object ####
        if (is.logical(samples.to.use)){
            se.obj <- se.obj.all
        }
        if (isTRUE(save.se.obj)){
            ## Check if metadata NCG already exists
            if (length(se.obj@metadata$NCG) == 0 ) {
                se.obj@metadata[['NCG']] <- list()
            }
            if (!'supervised' %in% names(se.obj@metadata[['NCG']])){
                se.obj@metadata[['NCG']][['supervised']] <- list()
            }
            if (!ncg.group.name %in% names(se.obj@metadata[['NCG']][['supervised']])){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]] <- list()
            }
            if (!ncg.set.name %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]] )){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]] <- list()
            }
            if (!'ncg.set' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]])){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- list()
            }
            se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- ncg.selected

            if (isTRUE(create.ncg.rank.plot)){
                if (!'ranl.plot' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]])){
                    se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['rank.plot']] <- list()
                }
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['rank.plot']] <- ncg.rank.plot
            }
            if (isTRUE(assess.ncg)){
                if (!'ranl.plot' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]])){
                    se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- list()
                }
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- assess.ncg.plot
            }
            printColoredMessage(
                message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
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
        ## add results to the SummarizedExperiment object ####
        if (isFALSE(save.se.obj)){
            printColoredMessage(message = '------------The findNcgPerBiologyPerBatch function finished.',
                                color = 'white',
                                verbose = verbose)
            return(list(ncg.selected = ncg.selected))
        }
    }
}
