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
#' @importFrom tidyr pivot_longer
#' @importFrom scran modelGeneVar
#' @importFrom ruv design.matrix
#' @importFrom Rfast correls
#' @importFrom stats aov
#' @import ggplot2
#' @export

findNcgUnSupervised <- function(
        se.obj,
        assay.name,
        uv.variables,
        approach = 'LinearMixedModel',
        clustering.method = 'kmeans',
        nb.clusters = 3,
        ncg.idenfitication.approach = 'LMM.BioUvAdjustment',
        form,
        nb.bio.pcs = 3,
        nb.uv.pcs = 3,
        use.rank = FALSE,
        ncg.selection.method = 'non.overlap',
        nb.ncg = 0.1,
        samples.to.use = NULL,
        hvg.method = 'mad',
        top.rank.bio.genes = 0.8,
        top.rank.uv.genes = 0.2,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'decrease',
        grid.nb = 20,
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
        nb.pcs = 10,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        plot.ncg.assessment = TRUE,
        regress.out.variables = NULL,
        regress.out.rle.med = TRUE,
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
    if (approach == 'LinearMixedModel'){
        se.obj <- findNcgUnSupervisedByLinearMixedModel(
            se.obj = se.obj,
            assay.name = assay.name,
            uv.variables = uv.variables,
            form = form,
            ncg.idenfitication.approach = ncg.idenfitication.approach,
            nb.bio.pcs = nb.bio.pcs,
            nb.uv.pcs = nb.uv.pcs,
            ncg.selection.method = ncg.selection.method,
            use.rank = use.rank,
            regress.out.rle.med = regress.out.rle.med,
            samples.to.use = samples.to.use,
            nb.ncg = nb.ncg,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.percentile = bio.percentile,
            uv.percentile = uv.percentile,
            grid.group = grid.group,
            grid.direction = grid.direction,
            grid.nb = grid.nb,
            filter.ncgs = filter.ncgs,
            common.hk = common.hk,
            nb.stable.genes = nb.stable.genes,
            hk.group = hk.group,
            create.ncg.rank.plot = create.ncg.rank.plot,
            assess.ncg = assess.ncg,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            svd.bsparam = svd.bsparam,
            plot.ncg.assessment = plot.ncg.assessment,
            nb.cores = nb.cores,
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
    if (approach == 'AnovaCorr'){
        se.obj <- findNcgsUnSupervisedByAnovaCorr(
            se.obj = se.obj,
            assay.name = assay.name,
            uv.variables = uv.variables,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            ncg.selection.method = ncg.selection.method,
            nb.ncg = nb.ncg,
            use.rank = use.rank,
            samples.to.use = samples.to.use,
            hvg.method = hvg.method,
            top.rank.bio.genes = top.rank.bio.genes,
            top.rank.uv.genes = top.rank.uv.genes,
            bio.percentile = bio.percentile,
            uv.percentile = uv.percentile,
            grid.group = grid.group,
            grid.direction = grid.direction,
            grid.nb = grid.nb,
            min.sample.for.mad = min.sample.for.mad,
            min.sample.for.var = min.sample.for.var,
            min.sample.for.aov = min.sample.for.aov,
            min.sample.for.correlation = min.sample.for.correlation,
            corr.method = corr.method,
            a = a ,
            rho = rho,
            anova.method = anova.method,
            create.ncg.rank.plot = create.ncg.rank.plot,
            plot.ncg.rank = plot.ncg.rank,
            filter.ncgs = filter.ncgs,
            common.hk = common.hk,
            nb.stable.genes = nb.stable.genes,
            hk.group = hk.group,
            assess.ncg = assess.ncg,
            variables.to.assess.ncg = variables.to.assess.ncg,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            plot.ncg.assessment = plot.ncg.assessment,
            regress.out.variables = regress.out.variables,
            regress.out.rle.med = regress.out.rle.med,
            normalization = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            check.se.obj = check.se.obj,
            remove.na = remove.na,
            ncg.group.name = ncg.group.name,
            ncg.set.name = ncg.set.name,
            save.imf = save.imf,
            use.imf = use.imf,
            imf.name = imf.name,
            save.se.obj = save.se.obj,
            verbose = verbose
        )
    }
    return(se.obj)
}
