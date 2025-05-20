#' Finds NCGs using supervised approaches.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function includes three different functions incuding `findNcgByTwoWayAnova()`, `findNcgAcrossSamples()`, and
#' `findNcgPerBiologyPerBatch()`, to identify a set of negative control genes (NCG) for `RUV-III-PRPS` normalization.
#' See each function's documentation for additional details.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character that indicates the name of an data (assay) in the SummarizedExperiment object.
#' The selected assay should be the one that will be used for the RUV-III-PRPS normalization.
#' @param bio.variables Character. A character string or vector of strings indicating the column name(s) of the biological
#' variable(s) in the SummarizedExperiment object. These variable can be categorical or continuous or a combination. This
#' argument cannot be `NULL`.
#' @param uv.variables Character. A character string or vector of strings indicating the column name(s) of the unwanted
#' variable(s) in the SummarizedExperiment object. These variable can be categorical or continuous or a combination.This
#'  argument cannot be `NULL`.
#' @param approach Character. A character that indicates which NCG approaches should be used. The options are :
#' `AnovaCorr.PerBatchPerBiology`, `AnovaCorr.AcrossAllSamples` and  `TwoWayAnova`. The default is set to `AnovaCorr.AcrossAllSamples`.
#' See the details for mroe information.
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
#' @param bio.groups A character vector or a symbol. Indicates the column names that contain biological variables in the
#' SummarizedExperiment object. If not NULL, bio.groups' will be used to group samples into different homogeneous biological groups.
#' @param bio.clustering.method A character string. Indicates which clustering method should be used to group continuous
#' sources of biological variation if any are provided. The default is 'kmeans' clustering.
#' @param nb.bio.clusters Numeric. Indicates the number of clusters for each continuous source of biological variation.
#' The default is set to 3.
#' @param uv.groups A character vector or a symbol. Indicates the column names that contain unwanted variation variables in the
#' SummarizedExperiment object. If not NULL, `uv.groups` will be used to group samples into different homogeneous unwanted
#' variation groups.
#' @param uv.clustering.method A character string. Indicates which clustering method should be used to group continuous
#' sources of unwanted variation. The default is set to 'kmeans' clustering.
#' @param nb.uv.clusters Numeric. Indicates the number of clusters for each continuous source of unwanted variation (UV).
#' By default, it is set to 2.
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
#' @param anova.method Character. A character that  indicates which ANOVA method to use. The options are `aov` or `welch`.
#' The default is set to `aov`. Refer to the function `row_oneway_equalvar()` or `row_oneway_welch()` from the
#' **matrixTests** R package for more details.
#' @param min.sample.for.aov Numeric. A numeric value that indicates the minimum number of samples that are required to
#' perform ANOVA analyses between continuous sources of variation (biological and unwanted variation) with individual
#' gene expression. The default is set to 3. The minimum value is 3.
#' @param corr.method Character. A character that indicates which correlation methods should be used for the correlation
#' analyses. The options are `pearson` or `spearman`. The default is set to `spearman`.
#' @param a Numeric. The significance level used for the confidence intervals in the correlation; by default, it is set
#' to 0.05. Refer to the function `correls` from the **Rfast** R package for more details.
#' @param rho Numeric. The value of the hypothesized correlation to be used in the hypothesis testing. The default is
#' set to 0. Refer to the function `correls` from the **Rfast** R package for more details.
#' @param min.sample.for.correlation Numeric. A numeric value that indicates the minimum number of samples that are required
#' to perform correlation analyses between continuous sources of variation (biological and unwanted variation) with
#'individual gene expression. The default is set to 10. The minimum value can be 3.
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
#' @param use.imf Logical. Indicates whether to use the intermediate file. The default is set to `FALSE`.
#' @param imf.name Character string. A name to save the intermediate file. If `NULL`, the function generates a name.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the
#' SummarizedExperiment object or output the result. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, shows messages of different steps of the function.
#'
#' @return Either the SummarizedExperiment object containing a set of negative control genes in the metadata or a
#' logical vector of the selected negative control genes.
#'
#' @importFrom SummarizedExperiment assay SummarizedExperiment
#' @export

findNcgSupervised <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        approach = 'AnovaCorr.AcrossAllSamples',
        ncg.selection.method = 'non.overlap',
        nb.ncg = 0.1,
        top.rank.bio.genes = 0.5,
        top.rank.uv.genes = 0.5,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'increase',
        grid.nb = 20,
        create.ncg.rank.plot = FALSE,
        plot.ncg.rank = FALSE,
        bio.groups = NULL,
        bio.clustering.method = 'kmeans',
        nb.bio.clusters = 3,
        uv.groups = NULL,
        uv.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        regress.out.uv.variables = NULL,
        regress.out.bio.variables = NULL,
        normalization = 'CPM',
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
        nb.pcs = 10,
        center = TRUE,
        scale = FALSE,
        plot.ncg.assessment = TRUE,
        check.se.obj = TRUE,
        remove.na = 'none',
        ncg.group.name = NULL,
        ncg.set.name = NULL,
        save.imf = FALSE,
        use.imf = FALSE,
        imf.name = NULL,
        save.se.obj = TRUE,
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
            create.ncg.rank.plot = create.ncg.rank.plot,
            plot.ncg.rank = plot.ncg.rank,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            regress.out.bio.variables = regress.out.bio.variables,
            regress.out.uv.variables = regress.out.uv.variables,
            bio.groups = bio.groups,
            bio.clustering.method = bio.clustering.method,
            nb.bio.clusters = nb.bio.clusters,
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
            plot.ncg.assessment = plot.ncg.assessment,
            check.se.obj = check.se.obj,
            remove.na = remove.na,
            ncg.group.name = ncg.group.name,
            ncg.set.name = ncg.set.name,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            save.se.obj = save.se.obj,
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
            create.ncg.rank.plot = create.ncg.rank.plot,
            plot.ncg.rank = plot.ncg.rank,
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
            plot.ncg.assessment = plot.ncg.assessment,
            check.se.obj = check.se.obj,
            remove.na = remove.na,
            ncg.group.name = ncg.group.name,
            ncg.set.name = ncg.set.name,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            save.se.obj = save.se.obj,
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
            create.ncg.rank.plot = create.ncg.rank.plot,
            plot.ncg.rank = plot.ncg.rank,
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
            plot.ncg.assessment = plot.ncg.assessment,
            check.se.obj = check.se.obj,
            remove.na = remove.na,
            ncg.group.name = ncg.group.name,
            ncg.set.name = ncg.set.name,
            save.imf = save.imf,
            imf.name = imf.name,
            use.imf = use.imf,
            save.se.obj = save.se.obj,
            verbose = verbose
            )
    }
    printColoredMessage(message = '------------The findNcgSupervised function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}
