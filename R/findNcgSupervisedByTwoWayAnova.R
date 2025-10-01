#' Finds NCGs using two way ANOVA.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function utilizes two-way ANOVA to identify a set of suitable genes as negative control genes (NCG) for RUV-III
#' normalization. Both biological and unwanted variation sources is necessary and should be specified.
#'
#' @details
#' The function begins by creating all possible sample groups based on biological and unwanted variation separately.
#' Subsequently, these groups are used as factors in two-way ANOVA to identify genes highly influenced by biological and
#' unwanted variation. Finally, the function selects genes with the possible highest F-statistics for unwanted variation and
#' lowest F-statistics for biological variation. Various approaches are employed for the final gene selection; please refer
#' to the details for more information.
#' The function uses 5 ways to summarize two gene-level F-statistics obtained for the biological and unwanted variation.
#' The function uses either the values or the ranks of F-statistics for NCGs selection. The function ranks the
#' negative of F-statistics values for unwanted variation. The lower the ranks, the greater the impact of unwanted
#' variation on genes. The function ranks the F-statistics for biological variation. The higher the ranks, the greater
#' the impact of biological variation on genes. The options are `prod`, `sum`, `average`, `auto`, `non.overlap` and
#' `quantile`.
#'
#' If `prod`, `sum` and `average` is set:
#'
#' * The product, sum or average of ranks of F-statistics is calculated. Then, the function selects `nb.ncg` numbers of
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
#'    \item The function selects the top `top.rank.uv.genes` genes that have the lowest ranks of F-statistics for
#'    unwanted variation.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1.
#'    \item If the number of selected genes is larger or smaller than the specified `nb.ncg`, the function applies an
#'    auto search to find approximate `nb.ncg` of genes as negative control genes as follow. The auto search will either
#'    decrease or increase the values of either `top.rank.bio.genes` or `top.rank.uv.genes` or both till to find
#'    approximate `nb.ncg` of genes as negative control genes.
#' }
#'
#' If `quantile` is selected:
#' \enumerate{
#'    \item The function selects the `bio.percentile` percentile of F-statistics for biological variation. Then, selects
#'    all the genes that have F-statistics larger the calculated percentile.
#'    \item The function selects the `uv.percentile` percentile of F-statistics for unwanted variation. Then, selects
#'    all the genes that have F-statistics larger the calculated percentile.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1.
#' }
#'
#' Assess the performance of NCGS:
#' * The function can assess the initial performance of selected NCGs. This analysis involves principal component analysis
#' on only the selected NCG and then explore the R^2 or vector correlation between the `nb.pcs` first principal components
#' and with the specified variables. Ideal NCGS, should show high and low R^2 or vector correlation for unwanted and
#' biological variation respectively.
#'
#' @references
#' * Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.
#' * Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character string indicating the name of the data (assay) in the `SummarizedExperiment`
#' object. This data should be the one that will be used as input data for the RUV-III normalization.
#' @param bio.variables Character. A character string or vector of strings indicating the column name(s) of the biological
#' variable(s) in the SummarizedExperiment object. These variable can be categorical or continuous or a combination. This
#' argument cannot be `NULL`.
#' @param uv.variables Character. A character string or vector of strings indicating the column name(s) of the unwanted
#' variable(s) in the SummarizedExperiment object. These variable can be categorical or continuous or a combination.This
#' argument cannot be `NULL`.
#' @param nb.ncg Numeric. A numeric value that specifies the number of genes to be chosen as negative control genes (NCG)
#' when the `ncg.selection.method` parameter is set to `auto`. This value, `nb.ncg`, corresponds to a fraction of the total
#' genes in the SummarizedExperiment object. The default is set to 0.1.
#' @param samples.to.use TTT
#' @param use.rank TTT
#' @param ratio.variable TTT
#' @param ncg.selection.method Character. A character that indicates how to summarize different statistics and select a
#' set of genes as negative control genes. The options are: `prod`, `average`, `sum`, `non.overlap`, `auto`, and `quantile`.
#' The default is set to `non.overlap`. For more information, refer to the details of the function.
#' @param rank.variable TTTT

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
#' @param filter.ncgs Logical. Whether to filter selected NCGs based on public human housekeeping gene sets. The default
#' is set to `FALSE`.
#' @param common.hk Character. Specifies group of housekeeping genes to use: `cancer` or `non.cancer`. The default is set
#' to `cancer`.
#' @param nb.stable.genes Numeric. A numeric value that specifies the number of top stable genes to be obtained from the
#' `getStableGenes()` function in the **singescore** R package. The default is set to 2000.
#' @param hk.group Character. Column name in the gene annotation containing non-cancer housekeeping genes. Options include:
#' `bulk.rnaseq.hk.genes.v1`, `bulk.rnaseq.hk.genes.v2`, `micorarray.hk.genes`, `nanostring.pan.cancer.hk.genes`,
#' `singscore.pan.cancer.hk.genes`. The default is set to `micorarray.hk.genes`.
#' @param create.ncg.rank.plot Logical. Indicates whether to generate a heatmap that shows the rank of the all genes
#' with respect to their biological and unwanted variation effects. The default is set to `FALSE`.
#' @param plot.ncg.rank Logical. Indicates whether to plot a heatmap that shows the rank of the all genes
#' with respect to their biological and unwanted variation effects, while function is running. The default is set to `FALSE`.
#' @param bio.clustering.method Character. A character Indicating which clustering methods should be used to group continuous
#' sources of biological variation if any is provided. The options are: `kmeans`, `cut`, and `quantile`. The default is
#' set to `kmeans` clustering. Refer to the `createHomogeneousBioGroups()` function for more details.
#' @param nb.bio.clusters Numeric. A numeric value indicating the number of clusters for each continuous source of biological
#' variation. The default is set to 2.
#' @param uv.clustering.method Character. Indicates which clustering methods should be used to group continuous sources
#' of unwanted variation if any is provided. The options are: `kmeans`, `cut`, and `quantile`. The default is
#' set to `kmeans` clustering. Refer to the `createHomogeneousUvGroups()` function for more details.
#' @param nb.uv.clusters Numeric. A numeric that indicates the number of clusters for each continuous source of unwanted
#' variation. The default is set to 2.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing any statistical
#' analysis. The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value to be added as a pseudo count to all measurements before log transformation
#' .The default is se to 1.
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
#' @param scale Logical. Indicates whether to scale the data before applying SVD. If `TRUE`, scaling is done by dividing
#' the(centered) columns of the assays by their standard deviations if centering is `TRUE`, and by the root mean square
#' otherwise. The default is set to `FALSE`.
#' @param svd.bsparam TTTT
#' @param plot.ncg.assessment Logical. Indicates whether to plot the output of the NCG assessment while function is running
#' . The default is set to `TRUE`.
#' @param nb.cores Numeric. A numeric value to specify number of cores for palatalization. The default is set to 1.
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
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object before any analysis. If `TRUE`,
#'  the function `checkSeObj()` will be used. The default is set to `TRUE`.
#' @param remove.na Character set. Indicates whether to remove NA or missing values from the SummarizedExperiment object
#' The options are: `assays`, the `sample.annotation`, `both`, or `none`. If `assays` is selected, genes containing NA or
#' missing values will be excluded. If `sample.annotation` is selected, the samples containing NA or missing values for
#' any `bio.variables` or `uv.variables` will be excluded. The default is set to `none`.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the SummarizedExperiment
#' object or output the result. The default is `TRUE`.
#' @param verbose Logical. If `TRUE`, shows messages of different steps of the function.
#' @param use.fvalues TTT
#'
#' @importFrom dplyr mutate progress_estimated
#' @importFrom BiocSingular runSVD bsparam
#' @importFrom SummarizedExperiment assay
#' @importFrom fastDummies dummy_cols
#' @importFrom tidyr pivot_longer
#' @importFrom parallel mclapply
#' @importFrom ggpubr ggarrange
#' @importFrom stats aov
#' @importFrom car Anova
#' @import ggplot2
#' @export

findNcgSupervisedByTwoWayAnova <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        nb.ncg = 0.05,
        samples.to.use = 'all',
        use.rank = FALSE,
        ratio.variable = 'eta.squared',
        ncg.selection.method = 'non.overlap',
        rank.variable = 'fvalue',
        top.rank.bio.genes = 0.8,
        top.rank.uv.genes = 0.2,
        bio.percentile = 0.2,
        uv.percentile = 0.8,
        grid.group = 'uv',
        grid.direction = 'increase',
        grid.nb = 20,
        filter.ncgs = FALSE,
        common.hk = 'cancer',
        nb.stable.genes = 2000,
        hk.group = 'micorarray.hk.genes',
        create.ncg.rank.plot = FALSE,
        use.fvalues = FALSE,
        plot.ncg.rank = FALSE,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 10,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        plot.ncg.assessment = TRUE,
        nb.cores = 1,
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
    printColoredMessage(
        message = '------------The findNcgByTwoWayAnova function starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking inputs ####
    if (!is.vector(assay.name) | length(assay.name) > 1 | is.logical(assay.name) | assay.name == 'all'){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    }
    if (is.null(bio.variables)){
        stop('The "bio.variables" cannot be empty or "NULL".')
    }
    if (is.null(uv.variables)){
        stop('The "uv.variables" cannot be empty or "NULL".')
    }
    if (!is.vector(bio.variables) | !is.vector(uv.variables) ){
        stop('The "uv.variables" and "bio.variables" must be a vector of variables name(s) in the SummarizedExperiment object.')
    }
    if (length(intersect(bio.variables, uv.variables)) > 0){
        stop('Individual specified variable must be either in the "bio.variables" or "uv.variables".')
    }
    if (!is.numeric(nb.ncg)){
        stop('The "nb.ncg" must be a positve numeric value 0 < nb.ncg < 1.')
    }
    if (nb.ncg >= 1 | nb.ncg <= 0){
        stop('The "nb.ncg" must be a positve value 0 < nb.ncg < 1.')
    }
    if (length(ncg.selection.method) > 1 | is.logical(ncg.selection.method)){
        stop('The "ncg.selection.method" muat be one of "ratio", prod", "sum", "average", "auto", "non.overlap" or "quantile".')
    }
    if (!ncg.selection.method %in% c('ratio', 'prod', 'sum', 'average', 'auto', 'non.overlap', 'quantile')){
        stop('The "ncg.selection.method" muat be one of "ratio", "prod", "sum", "average", "auto", "non.overlap" or "quantile".')
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
    if (ncg.selection.method %in% c('non.overlap' , 'auto')){
        if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
            stop('The "top.rank.bio.genes" must be a positve value  0 < top.rank.bio.genes =< 1.')
        } else if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
            stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes =< 1.')
        }
    }
    if (ncg.selection.method == 'ratio'){
        if (!ratio.variable %in% c('fvalue', 'variance')){
            stop('The "ratio.variable" must be one of the "fvalue" or "variance".')
        }
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

    # Performing gene level two_way ANOVA ####
    if (isFALSE(use.imf)){
        ## Data transformation ####
        printColoredMessage(
            message = '-- Applying data transformation:',
            color = 'magenta',
            verbose = verbose
            )
        ### Applying log2 + pseudo count transformation ####
        if (isTRUE(apply.log)){
            expr.data <- applyLog(
                se.obj = se.obj,
                assay.names = assay.name,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose
                )[[assay.name]]
        }
        if (isFALSE(apply.log)) {
            printColoredMessage(
                message = paste0(
                    '- The ',
                    assay.name,
                    ' data will be used without any log transformation.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- assay(x = se.obj, i = assay.name)
        }

        ## Creating all possible homogeneous sample groups ####
        ### Biological groups ####
        printColoredMessage(
            message = '-- Creating all possible groups with respect to the specified sources of biological variation:',
            color = 'magenta',
            verbose = verbose
            )
        all.bio.groups <- createHomogeneousBioGroups(
            se.obj = se.obj,
            bio.variables = bio.variables,
            nb.clusters = nb.bio.clusters,
            clustering.method = bio.clustering.method,
            check.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = remove.na,
            verbose = verbose
            )
        ### Unwanted groups ####
        printColoredMessage(
            message = '-- Creating all possible groups with respect to the specified sources of unwanted variation:',
            color = 'magenta',
            verbose = verbose
            )
        all.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = uv.variables,
            nb.clusters = nb.uv.clusters,
            clustering.method = uv.clustering.method,
            check.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
        ## Applying Two way ANOVA ####
        printColoredMessage(
            message = '-- Performing two way ANOVA:',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = paste0(
                '- This is between all individual gene-level expression',
                ' and considering both biological and unwanted variables created above as factors.'),
            color = 'blue',
            verbose = verbose
            )
        all.aov <- mclapply(
            1:nrow(expr.data),
            function(x){
                sub.data <- data.frame(batch = all.uv.groups, bio = all.bio.groups, gene = expr.data[x, ])
                lm.fit <- lm(gene ~ batch + bio , sub.data)
                result.twa <- car::Anova(lm.fit, type = "II")
                ## Calculating partial eta squared ####
                f.values <- result.twa$`F value`
                df.effect <- result.twa$Df
                df.error <- result.twa["Residuals", "Df"]
                partial.eta.sq <- (f.values * df.effect) / (f.values * df.effect + df.error)
                total.sum.sq <- sum(result.twa$`Sum Sq`)
                data.frame(
                    uv.pct = result.twa$`Sum Sq`[1] / total.sum.sq,
                    bio.pct = result.twa$`Sum Sq`[2] / total.sum.sq,
                    residual.pct = result.twa$`Sum Sq`[3] / total.sum.sq,
                    uv.eta.p = partial.eta.sq[1],
                    bio.eta.p = partial.eta.sq[2],
                    uv.fvalue = result.twa$`F value`[1],
                    bio.fvalue = result.twa$`F value`[2]
                )
            }, mc.cores = nb.cores
            )
        all.aov <- do.call(rbind, all.aov)
        all.aov$uv.fvalue.scaled <- (all.aov$uv.fvalue - min(all.aov$uv.fvalue)) /
            (max(all.aov$uv.fvalue) - min(all.aov$uv.fvalue))
        all.aov$bio.fvalue.scaled <- (all.aov$bio.fvalue - min(all.aov$bio.fvalue)) /
            (max(all.aov$bio.fvalue) - min(all.aov$bio.fvalue))
        all.aov <- round(x = all.aov, digits = 3)
        row.names(all.aov) <- row.names(se.obj)

        ### Ranking of the F-statistics ####
        if (isTRUE(use.rank)){
            if (rank.variable == 'fvalue'){
                set.seed(2190)
                all.aov$uv.rank <- rank(x = -all.aov$uv.fvalue, ties.method = 'random')
                set.seed(2190)
                all.aov$bio.rank <- rank(x = all.aov$bio.fvalue, ties.method = 'random')
                all.aov$uv.rank.plot <- rank(x = all.aov$uv.fvalue, ties.method = 'random')
            }
            ### Ranking of the percentage variation ####
            if (rank.variable == 'variance'){
                set.seed(2190)
                all.aov$uv.rank <- rank(x = -all.aov$uv.pct, ties.method = 'random')
                set.seed(2190)
                all.aov$bio.rank <- rank(x = all.aov$bio.pct, ties.method = 'random')
                all.aov$uv.rank.plot <- rank(x = all.aov$uv.pct, ties.method = 'random')
            }
            ### Ranking of the partial eta sqaured ####
            if (rank.variable == 'eta.squared'){
                set.seed(2190)
                all.aov$uv.rank <- rank(x = -all.aov$uv.eta.p, ties.method = 'random')
                set.seed(2190)
                all.aov$bio.rank <- rank(x = all.aov$bio.eta.p, ties.method = 'random')
                all.aov$uv.rank.plot <- rank(x = all.aov$uv.eta.p, ties.method = 'random')
            }
        }
    }
    # Reading the intermediate file ####
    if (isTRUE(use.imf)){
        printColoredMessage(
            message = '- Retrieving the results of two-way ANOVA from the the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        if (is.null(imf.name)){
            imf.name <- paste0(
                assay.name,
                '|TwoWayAnova|',
                ncg.selection.method
                )
        }
        if (is.null(se.obj@metadata$IMF$NCG[[imf.name]]))
            stop('The intermediate file cannot be found in the metadata of the SummarizedExperiment object.')
        all.aov <- se.obj@metadata$IMF$NCG[[imf.name]]
    }
    # Saving the intermediate file ####
    if (isTRUE(save.imf)){
        printColoredMessage(
            message = '-- Saving a intermediate file:',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '- The results of two-way ANOVA is saved in the SummarizedExperiment object.',
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
            imf.name <- paste0(
                assay.name,
                '|TwoWayAnova|',
                ncg.selection.method)
        }
        if (!imf.name %in% names(se.obj@metadata[['IMF']][['NCG']])){
            se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- list()
        }
        se.obj@metadata[['IMF']][['NCG']][[imf.name]] <- all.aov
    }

    # Selecting a set of genes as NCG ####
    printColoredMessage(
        message = '-- Selecting a set of genes as NCG:',
        color = 'magenta',
        verbose = verbose
        )

    ## Applying ratio approach ####
    if (isFALSE(use.rank)) {
        ### Ratio of F value ####
        if (ratio.variable == 'fvalue'){
            all.aov$ratio.fvalue <- all.aov$uv.fvalue.scaled / c(all.aov$bio.fvalue.scaled + 1e-6)
            all.aov <- all.aov[order(-all.aov$ratio.fvalue), , drop = FALSE]
            all.aov <- all.aov[all.aov$uv.fvalue.scaled > quantile(x = all.aov$uv.fvalue.scaled, probs = uv.percentile) , ]
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% row.names(all.aov)[1:nb.ncg]
        }
        ## Ratio of variance ####
        if (ratio.variable == 'variance'){
            all.aov$ratio.var <- all.aov$uv.pct / c(all.aov$bio.pct + 1e-6)
            all.aov <- all.aov[order(-all.aov$ratio.var), ]
            all.aov <- all.aov[all.aov$uv.pct > quantile(x = all.aov$uv.pct, probs = uv.percentile) , ]
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% row.names(all.aov)[1:nb.ncg]
        }
        ## Ratio of partial eta ####
        if (ratio.variable == 'eta.squared'){
            all.aov$ratio.var <- all.aov$uv.eta.p / c(all.aov$bio.eta.p + 1e-6)
            all.aov <- all.aov[order(-all.aov$ratio.var), ]
            all.aov <- all.aov[all.aov$uv.eta.p > quantile(x = all.aov$uv.eta.p, probs = uv.percentile) , ]
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% row.names(all.aov)[1:nb.ncg]
        }
    }
    ## Applying product, average and sum of ranks ####
    if (ncg.selection.method %in% c('prod', 'average', 'sum') & isTRUE(use.rank)){
        ### Product of ranks ####
        if (isTRUE(ncg.selection.method == 'prod')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the product of ranks.',
                color = 'blue',
                verbose = verbose)
            all.aov$all.rank <- all.aov$bio.rank * all.aov$uv.rank
            if (sum(is.infinite(all.aov$all.rank)) > 0)
                stop('The product of ranks results in infinity values.')
        }
        ## Average of ranks ####
        if (isTRUE(ncg.selection.method == 'average')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the average of ranks.',
                color = 'blue',
                verbose = verbose)
            all.aov$all.rank <- rowMeans(all.aov[ , c('bio.rank', 'uv.rank')])
        }
        ## Sum of ranks ####
        if (isTRUE(ncg.selection.method == 'sum')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the sum of ranks.',
                color = 'blue',
                verbose = verbose
                )
            all.aov$all.rank <- rowSums(all.aov[ , c('bio.rank', 'uv.rank')])
        }
        ## Selecting top genes as NCGS ####
        nb.ncg <- round(x = nb.ncg * nrow(se.obj), digits = 0)
        printColoredMessage(
            message = paste0(
                '- Selecting ',
                nb.ncg ,
                ' genes as NCGS.'),
            color = 'blue',
            verbose = verbose
            )
        all.aov <- all.aov[order(all.aov$all.rank, decreasing = FALSE), ]
        ncg.selected <- row.names(all.aov)[1:nb.ncg]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }
    ## Applying quantile approach ####
    if (ncg.selection.method == 'quantile' & isTRUE(use.rank)){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "quantile" approach.',
            color = 'blue',
            verbose = verbose
            )
        ### Selecting biological percentile ####
        bio.quan <- quantile(x = all.aov$bio.rank , probs = bio.percentile)
        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > bio.quan]

        ## Selecting UV percentile ####
        uv.quan <- quantile(x = all.aov$uv.rank , probs = uv.percentile)
        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank < uv.quan]
        printColoredMessage(
            message = paste0(
                '- Selecting ',
                length(top.uv.genes),
                ' genes with the uv ranked F-statistics lower than ',
                uv.quan,
                ' (' ,
                uv.percentile* 100,
                '% percentile), and exclude any genes presents in ',
                length(top.bio.genes),
                ' genes with the biological ranked F-statistics higher than ',
                bio.quan,
                ' (' ,
                bio.percentile* 100,
                '% percentile).'),
            color = 'blue',
            verbose = verbose)
        ## Selecting top genes as NCGS ####
        top.uv.genes <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if (isTRUE(length(top.uv.genes) == 0)) stop('No NCGs can be found based on the current parameters.')
        ncg.selected <- row.names(se.obj) %in% top.uv.genes
    }
    ## Applying non overlap approach ####
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
                '% of highly affected genes by the unwanted variation, and then exclude top ',
                top.rank.bio.genes *100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
            )
        ### Selecting genes highly affected by biological variation ####
        top.rank.bio.genes.nb <- round(c(1-top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
        ## Selecting genes highly affected by unwanted variation ####
        top.rank.uv.genes <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes]
        ## Selecting top genes as NCGS ####
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
        if (isTRUE(sum(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
    }
    ## Applying auto approach ####
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
                '% of highly affected genes by the unwanted variation, and then exclude all genes in top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
            )
        ### Selecting genes affected by biological variation ####
        nb.ncg <- round(nb.ncg * nrow(se.obj), digits = 0)
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
        ## Selecting genes affected by unwanted variation ####
        top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank < top.rank.uv.genes.nb]
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
        ncg.ranges <- round(x = 0.01 *nb.ncg, digits = 0)
        if (length(ncg.selected) > c(nb.ncg + ncg.ranges) | length(ncg.selected) < c(nb.ncg - ncg.ranges) ){
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
                        '- The number of selected genes ',
                        length(ncg.selected),
                        ' is larger than ',
                        nb.ncg ,
                        ', which is the specified number of NCG ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose)
            }
            #### Applying grid search ####
            ##### grid group: both bio and uv variable ####
            if (grid.group == 'both'){
                printColoredMessage(
                    message = '- The grid search will be applied on both biological and unwanted factors. ',
                    color = 'blue',
                    verbose = verbose
                    )
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- min(
                        nrow(se.obj) - top.rank.uv.genes.nb,
                        top.rank.bio.genes.nb
                        )
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while (eval(con) & top.rank.uv.genes.nb < nrow(se.obj) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if (top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes.nb]
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if (top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- min(top.rank.uv.genes.nb, c(nrow(se.obj) - top.rank.bio.genes.nb))
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while (eval(con) & top.rank.uv.genes.nb > 1 & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if (top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes.nb]
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if (top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                # genes selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                ## bio
                top.rank.bio.genes <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes/nrow(se.obj) * 100, digits = 2)
                if (top.rank.bio.genes >= 100) top.rank.bio.genes = 100
                ## uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if (top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Updating the selection. Selecting top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then excluding any genes in top ',
                        top.rank.bio.genes,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
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
                        message = '- The grid search will increase the number of "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                        )
                    lo <- top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if (top.rank.bio.genes.nb < 1 ) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(all.aov)[all.aov$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
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
                        if (top.rank.bio.genes.nb > nrow(se.obj) ) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(all.aov)[ all.aov$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
                }
                # gene selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # bio
                top.rank.bio.genes.nb <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes.nb/nrow(se.obj) * 100, digits = 2)
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
                        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
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
                        top.uv.genes <- row.names(all.aov)[all.aov$uv.rank <  top.rank.uv.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
                }
                # gene selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if (top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Updating the selection. Selecting top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes * 100,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
                }
            } else {
            printColoredMessage(
                message = paste0(
                    length(ncg.selected),
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
            '- ',
            sum(ncg.selected),
            ' genes are selected as negative control genes.'),
        color = 'blue',
        verbose = verbose
        )

    # Plotting the F-statistics ####
    if (isTRUE(create.ncg.rank.plot)){
        if (isTRUE(use.fvalues)){
            all.aov$ncg <- factor(x = ncg.selected, levels = c('TRUE', 'FALSE'))
            p.fvalues <- ggplot(all.aov, aes(x = log2(bio.fvalue), y = log2(uv.fvalue), color =  ncg)) +
                geom_point(alpha = .1) +
                scale_color_manual(values = c('darkgreen', 'grey10')) +
                xlab('Biology (log2 of F-values)') +
                ylab('Unwanted variation (log2 of F-values)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            p.fvalues.rank <- ggplot(all.aov, aes(x = bio.rank, y = uv.rank.plot , color =  ncg)) +
                geom_point(alpha = .1) +
                scale_color_manual(values = c('darkgreen', 'grey10')) +
                xlab('Biology (rank of F-values)') +
                ylab('Unwanted variation (rank of F-values)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            all.plots <- ggarrange(p.fvalues + p.fvalues.rank)
            print(all.plots)
            rm(p.fvalues, p.fvalues.rank)
        }
        if (isFALSE(use.fvalues)){
            all.aov$ncg <- factor(x = ncg.selected, levels = c('TRUE', 'FALSE'))
            p.fvalues <- ggplot(all.aov, aes(x = bio.pct, y = uv.pct, color =  ncg)) +
                geom_point(alpha = .1) +
                scale_color_manual(values = c('darkgreen', 'grey10')) +
                xlab('Biology (percentage of variation)') +
                ylab('Unwanted variation (percentage of variation)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            p.fvalues.rank <- ggplot(all.aov, aes(x = bio.rank, y = uv.rank.plot, color =  ncg)) +
                geom_point(alpha = .1) +
                scale_color_manual(values = c('darkgreen', 'grey10')) +
                xlab('Biology (rank of percentage of variation') +
                ylab('Unwanted variation (rank of percentage of variation)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            all.plots <- ggarrange(p.fvalues + p.fvalues.rank)
            print(all.plots)
            rm(p.fvalues, p.fvalues.rank)
        }
    } else all.plots = NULL
    # Assessing the performance of selected NCG ####
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
    } else assess.ncg.plot = NULL

    # Saving the results ####
    ### Adding the results to the SummarizedExperiment object ####
    if (is.logical(samples.to.use)){
        se.obj <- se.obj.all
    }
    #### Selecting output name ####
    if (is.null(ncg.group.name)){
        ncg.group.name <- 'NcgTwoWayAnova'
    }
    if (is.null(ncg.set.name)){
        ncg.set.name <- paste0(
            sum(ncg.selected),
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|TWAnova:',
            ncg.selection.method,
            '|',
            assay.name)
    }
    #### Saving all the results ####
    if (isTRUE(save.se.obj)){
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
            se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['rank.plot']] <- all.plots
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
    ### Export results as logical vector ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = '-- The set of NCGs is outpputed as a logical vector.',
            color = 'magenta',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The findNcgByTwoWayAnova function finished.',
            color = 'white',
            verbose = verbose
            )
        return(list(
            ncg = ncg.selected,
            assess.ncg.plot = assess.ncg.plot,
            rank.plots = all.plots))
    }
}




