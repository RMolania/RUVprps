#' Finds a set of negative control genes using supervised approaches.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function includes three different methods — `findNcgByTwoWayAnova`, `findNcgAcrossSamples`, and
#' `findNcgPerBiologyPerBatch` — to identify a set of negative control genes (NCG) for `RUV-III-PRPS` normalization.
#' See each function's documentation for additional details.
#'
#' @param se.obj A `SummarizedExperiment` object.
#' @param assay.name Character. The name of the assay in the `SummarizedExperiment` object to be used for `RUV-III-PRPS` normalization.
#' Raw (unnormalized) data is recommended.
#' @param bio.variables Character. Column names in the `SummarizedExperiment` object that contain biological variables.
#' These may be categorical or continuous. Continuous variables will be grouped into `nb.bio.clusters` using the method specified
#' in `bio.clustering.method`. Cannot be `NULL`.
#' @param uv.variables Character. Column names representing unwanted variables in the `SummarizedExperiment` object.
#' Continuous variables will be grouped into `nb.uv.clusters` using `uv.clustering.method`. Cannot be `NULL`.
#' @param approach Character. Method for selecting NCGs. Options: `AnovaCorr.PerBatchPerBiology`, `AnovaCorr.AcrossAllSamples`,
#' and `TwoWayAnova`. The default is set to `TwoWayAnova`.
#' @param ncg.selection.method Character. Strategy to summarize F-statistics from two-way ANOVA and select NCGs.
#' Options: `prod`, `average`, `sum`, `non.overlap`, `auto`, and `quantile`. The default is set to `non.overlap`.
#' @param nb.ncg Numeric. Proportion of total genes to select as NCGs. The default is set to `0.1`.
#' @param top.rank.bio.genes Numeric. Fraction of genes highly influenced by biological variation, required for
#' `non.overlap` or `auto` methods. The default is set to `0.5`.
#' @param top.rank.uv.genes Numeric. Fraction of genes highly influenced by unwanted variation, required for
#' `non.overlap` or `auto`. The default is set to `0.5`.
#' @param bio.percentile Numeric. F-statistic percentile threshold for identifying biologically variable genes. The default is set to `0.8`.
#' @param uv.percentile Numeric. F-statistic percentile threshold for identifying genes affected by unwanted variation. The default is set to `0.8`.
#' @param grid.group Character. Which factor to use in grid search for `auto` method: `bio`, `uv`, or `both`. The default is set to `uv`.
#' @param grid.direction Character. Order of grid search: `increase` or `decrease`. The default is set to `decrease`.
#' @param grid.nb Numeric. Number of genes to test during grid search. The default is set to `20`.
#' @param bio.groups Character. Column name(s) used to group samples by biological variables. If `NULL`, `bio.variables` is used.
#' @param bio.clustering.method Character. Clustering method for grouping continuous biological variation.
#' See `createHomogeneousBioGroups` for options. The default is set to `kmeans`.
#' @param nb.bio.clusters Numeric. Number of clusters per continuous biological variable. The default is set to `2`.
#' @param uv.groups Character. Column name(s) for grouping samples by unwanted variables. If `NULL`, `uv.variables` is used.
#' @param uv.clustering.method Character. Clustering method for grouping continuous unwanted variation.
#' See `createHomogeneousUvGroups` for options. The default is set to `kmeans`.
#' @param nb.uv.clusters Numeric. Number of clusters per continuous unwanted variable. The default is set to `2`.
#' @param normalization Character. Normalization method before assessing biological variation. The default is set to `CPM`.
#' If `NULL`, no normalization is applied. See `applyOtherNormalizations` for details.
#' @param regress.out.bio.variables Character. Column names of biological variables to regress out before identifying
#' unwanted-variable-associated genes. The default is set to `NULL`.
#' @param regress.out.uv.variables Character. Column names of unwanted variables to regress out before identifying
#' biologically associated genes. The default is set to `NULL`.
#' @param apply.log Logical. Whether to log-transform the data prior to analysis. The default is set to `TRUE`.
#' @param pseudo.count Numeric. Pseudo-count added before log transformation. The default is set to `1`.
#' @param anova.method Character. Method for ANOVA. See specific method documentation.
#' @param min.sample.for.aov Numeric. Minimum samples required per group for ANOVA. The default is set to `3`.
#' @param corr.method Character. Correlation method for association analysis. Options: `pearson`, `spearman`. The default is set to `spearman`.
#' @param a Numeric. Significance level (alpha) for correlation confidence intervals. The default is set to `0.05`.
#' @param rho Numeric. Hypothesized correlation value. The default is set to `0`.
#' @param min.sample.for.correlation Numeric. Minimum number of samples per group for correlation analysis. The default is set to `10`.
#' @param assess.ncg Logical. Whether to evaluate selected NCGs using PCA and correlation with variables. The default is set to `TRUE`.
#' @param variables.to.assess.ncg Character. Variables used to assess selected NCGs. If `NULL`, both `bio.variables` and
#' `uv.variables` are used. The default is set to `NULL`.
#' @param nb.pcs Numeric. Number of principal components used in NCG performance evaluation. The default is set to `5`.
#' @param center Logical. Whether to center the data before PCA. See `computePCA` for details. The default is set to `TRUE`.
#' @param scale Logical. Whether to scale the data before PCA. The default is set to `FALSE`.
#' @param assess.se.obj Logical. Whether to validate the `SummarizedExperiment` object using `checkSeObj()`. The default is set to `TRUE`.
#' @param remove.na Character. Indicates whether to remove `NA` values from `assays`, `sample.annotation`, `both`, or `none`.
#' The default is set to `none`.
#' @param save.se.obj Logical. Whether to save results in `se.obj@metadata$NCG$supervised$output.name`. The default is set to `TRUE`.
#' @param output.name Character. Name to store results under. If `NULL`, auto-generated using:
#' `paste0(sum(ncg.selected), '|', paste0(bio.variables, collapse = '&'), '|', paste0(uv.variables, collapse = '&'), '|TWAnova:', ncg.selection.method, '|', assay.name)`.
#' @param ncg.group Character. Label for the group of selected NCGs.
#' @param plot.output Character. Whether and what type of plot to produce.
#' @param use.imf Logical. Whether to use an intermediate file. The default is set to `FALSE`.
#' @param save.imf Logical. Whether to save the intermediate file (results from two-way ANOVA). Speeds up tuning when reusing results.
#' The default is set to `FALSE`.
#' @param imf.name Character. Name for the intermediate file. If `NULL`, auto-generated as:
#' `paste0(assay.name, '|TwoWayAnova|', ncg.selection.method)`.
#' @param verbose Logical. If `TRUE`, display messages during execution.
#'
#' @return Either the updated `SummarizedExperiment` object containing the selected negative control genes
#' or a logical vector identifying the NCGs.



#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @export

findNcgSupervised <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        approach = 'TwoWayAnova',
        ncg.selection.method = 'non.overlap',
        nb.ncg = 0.1,
        top.rank.bio.genes = 0.5,
        top.rank.uv.genes = 0.5,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'increase',
        grid.nb = 20,
        bio.groups = NULL,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.groups = NULL,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        normalization = 'CPM',
        regress.out.uv.variables = NULL,
        regress.out.bio.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        anova.method = 'aov',
        min.sample.for.aov = 3,
        corr.method = "spearman",
        a = 0.05,
        rho = 0,
        min.sample.for.correlation = 10,
        assess.ncg = TRUE,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        assess.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        output.name = NULL,
        ncg.group = NULL,
        plot.output = TRUE,
        use.imf = FALSE,
        save.imf = FALSE,
        imf.name = NULL,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The findNcgSupervised function starts:',
                        color = 'white',
                        verbose = verbose)
    # check inputs ####
    if (!approach %in% c('AnovaCorr.PerBatchPerBiology', 'AnovaCorr.AcrossAllSamples', 'TwoWayAnova')){
        stop('The approach must be one of the "AnovaCorr.PerBatchPerBiology", "AnovaCorr.AcrossAllSamples" or "TwoWayAnova".')
    }

    # find NCGs ####
    ## find NCGs using AnovaCorr.PerBatchPerBio approach ####
    if (approach == 'AnovaCorr.PerBatchPerBiology'){
        se.obj <- findNcgPerBiologyPerBatch(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            ncg.selection.method = ncg.selection.method,
            nb.ncg = nb.ncg,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.percentile = bio.percentile,
            uv.percentile = uv.percentile,
            grid.group = grid.group,
            grid.direction = grid.direction,
            grid.nb = grid.nb,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            regress.out.bio.variables = regress.out.bio.variables,
            bio.groups = bio.groups,
            bio.clustering.method = bio.clustering.method,
            nb.bio.clusters = nb.bio.clusters,
            regress.out.uv.variables = regress.out.uv.variables,
            uv.groups = uv.groups,
            uv.clustering.method = uv.clustering.method,
            nb.uv.clusters = nb.uv.clusters,
            normalization = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            corr.method = corr.method,
            a = a,
            rho = rho,
            anova.method = anova.method,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            output.name = output.name,
            ncg.group = ncg.group,
            plot.output = plot.output,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            verbose = verbose
            )
    }
    ## find NCGs using AnovaCorr.AcrossAllSamples approach ####
    if (approach == 'AnovaCorr.AcrossAllSamples'){
        se.obj <- findNcgAcrossSamples(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            ncg.selection.method = ncg.selection.method,
            nb.ncg = nb.ncg,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.percentile = bio.percentile,
            uv.percentile = uv.percentile,
            grid.group = grid.group,
            grid.direction = grid.direction,
            grid.nb = grid.nb,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            regress.out.bio.variables = regress.out.bio.variables,
            regress.out.uv.variables = regress.out.uv.variables,
            normalization = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            corr.method = corr.method,
            a = a,
            rho = rho,
            anova.method = anova.method,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            output.name = output.name,
            ncg.group = ncg.group,
            plot.output = plot.output,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            verbose = verbose
            )
    }
    ## find NCGs using TwoWayAnova approach ####
    if (approach == 'TwoWayAnova'){
        se.obj <- findNcgByTwoWayAnova(
            se.obj = se.obj,
            assay.name = assay.name,
            bio.variables = bio.variables,
            uv.variables = uv.variables,
            ncg.selection.method = ncg.selection.method,
            nb.ncg = nb.ncg,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.percentile = bio.percentile,
            uv.percentile = uv.percentile,
            grid.group = grid.group,
            grid.direction = grid.direction,
            grid.nb = grid.nb,
            bio.clustering.method = bio.clustering.method,
            nb.bio.clusters = nb.bio.clusters,
            uv.clustering.method = uv.clustering.method,
            nb.uv.clusters = nb.uv.clusters,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            save.se.obj = save.se.obj,
            output.name = output.name,
            ncg.group = ncg.group,
            plot.output = plot.output,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            verbose = verbose
            )
    }
    printColoredMessage(message = '------------The findNcgSupervised function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
