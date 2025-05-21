#' Assesses the variation in variables.
#'
#' @author Ramyar Molania
#'
#'#' @references
#' Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research, 2019.
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS, Nature Biotechnology, 2023
#'
#' @description
#' This function applies a range of global and gene-level metrics to assess the variation of specified
#' variables in the assay(s) of a SummarizedExperiment object.
#'
#' @details
#' Several assessments will be performed:
#' For each categorical variable:
#' - PCA plot of the categorical variable.
#' - Silhouette and ARI computed on the categorical variable.
#' - Differential analysis based on ANOVA between the gene expression and the categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and the categorical variable.
#'
#' For each continuous variable:
#' - Linear regression between the first cumulative PC and the continuous variable.
#' - Correlation between gene expression and the continuous variable.
#'
#' It will output the following plots:
#' - PCA plot of each categorical variable.
#' - Boxplot of the F-test distribution from ANOVA between gene expression and each categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and each categorical variable.
#' - Combined silhouette plot of all categorical variable pairs.
#' - Linear regression between the first cumulative PC and each continuous variable.
#' - Boxplot of the correlation between gene expression and each continuous variable.
#' - RLE plot distribution.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or vector for selecting assay name(s) in the SummarizedExperiment object. The default is set to `all`.
#' @param variables Character. Column names containing variables to assess for variation. Can be categorical or continuous.
#' @param plots.to.exclude Character. A character string or vector indicating which metrics to exclude. Use `getAssessmentMetrics()` to view options. The default is set to `NULL`.
#' @param apply.log Logical. Whether to apply a log-transformation. The default is set to `FALSE`.
#' @param pseudo.count Numeric. Pseudo count added to all measurements before log-transformation. The default is set to 1.
#' @param general.points.size Numeric. Size of points in scatter plots. The default is set to 1.
#' @param rle.iqr.width Numeric. Width of IQR in RLE plots. The default is set to 1.
#' @param rle.median.points.size Numeric. Size of median points in RLE plots. The default is set to 1.
#' @param rle.median.points.color Character. Color of median points in RLE plots. The default is set to `red`.
#' @param rle.geom.hline.color Character. Color of the horizontal line in RLE plots. The default is set to `cyan`.
#' @param rle.plot.ncol Numeric. Number of columns in the RLE boxplot grid. The default is set to 1.
#' @param rle.plot.nrow Numeric. Number of rows in the RLE boxplot grid. The default is set to 3.
#' @param rle.var.plot.ncol Numeric. Number of columns in scatter plots of RLE medians against variables. The default is set to 3.
#' @param rle.var.plot.nrow Numeric. Number of rows in scatter plots of RLE medians against variables. The default is set to 3.
#' @param rle.colors List. List of colors for variables in RLE plots. The default is set to `NULL`.
#' @param fast.pca Logical. Use fewer PCs for faster PCA. The default is set to `TRUE`.
#' @param compute.nb.pcs Numeric. Number of PCs to compute for fast PCA. The default is set to 4.
#' @param nb.pcs.toplot.pca Numeric. Number of PCs to plot against each other. The default is set to 3. Must be ≤ `compute.nb.pcs`.
#' @param center Logical. Whether to center data before SVD. The default is set to `TRUE`.
#' @param scale Logical. Whether to scale data before SVD. The default is set to `FALSE`.
#' @param svd.bsparam A `BiocParallelParam` object for parallelization. The default is set to `bsparam()`. See `runSVD()` in BiocSingular.
#' @param pca.variables.colors List. List of colors for categorical variables in PCA plots. The default is set to `NULL`.
#' @param pca.plot.ncol Numeric. Number of columns in PCA plot grid. The default is set to 3.
#' @param pca.plot.nrow Numeric. Number of rows in PCA plot grid. The default is set to 3.
#' @param pca.stroke.size Numeric. Stroke size in PCA plots. The default is set to 0.1.
#' @param pca.stroke.color Character. Color of stroke in PCA plots. The default is set to `gray`.
#' @param pca.points.alpha Numeric. Transparency of PCA points. The default is set to 0.5.
#' @param pca.densities.alpha Numeric. Transparency of densities in PCA plots. The default is set to 0.5.
#' @param sil.dist.measure Character. Distance measure for silhouette: `euclidean`, `maximum`, `manhattan`, `canberra`,
#' `binary`, or `minkowski`. The default is set to `euclidean`.
#' @param sli.nb.pcs Numeric. Number of PCs to use for silhouette. The default is set to 3. Must be ≤ `compute.nb.pcs`.
#' @param ari.clustering.method Character. Clustering method for ARI: `mclust` or `hclust`. The default is set to `hclust`.
#' @param ari.hclust.method Character. Agglomeration method for `hclust`: `ward.D`, `ward.D2`, `single`, `complete`,
#' `average`, `mcquitty`, `median`, `centroid`. The default is set to `complete`.
#' @param ari.hclust.dist.measure Character. Distance measure for `hclust`: `euclidean`, `maximum`, `manhattan`, `canberra`,
#'  `binary`, `minkowski`. The default is set to `euclidean`.
#' @param ari.nb.pcs Numeric. Number of PCs for ARI. The default is set to 3. Must be ≤ `compute.nb.pcs`.
#' @param vca.nb.pcs Numeric. Number of PCs for vector correlation. The default is set to 3.
#' @param lra.nb.pcs Numeric. Number of PCs for linear regression. The default is set to 3.
#' @param corr.method Character. Correlation method: `pearson`, `kendall`, `spearman`. The default is set to `spearman`.
#' @param a Numeric. Significance level for confidence intervals in correlation. The default is set to 0.05.
#' @param rho Numeric. Hypothesized correlation for testing. The default is set to 0.
#' @param correlation.plot.ncol Numeric. Number of columns in correlation boxplot grid. The default is set to 1.
#' @param correlation.plot.nrow Numeric. Number of rows in correlation boxplot grid. The default is set to 3.
#' @param anova.method Character. ANOVA method: `aov` or `welch`. The default is set to `aov`.
#' @param anova.plot.ncol Numeric. Columns in ANOVA F-value boxplot grid. The default is set to 1.
#' @param anova.plot.nrow Numeric. Rows in ANOVA F-value boxplot grid. The default is set to 3.
#' @param pcorr.method Character. Method for gene-gene partial correlation: `pearson`, `kendall`, `spearman`. The default is set to `spearman`.
#' @param pcorr.genes Vector. Genes for pairwise correlation. Logical, numeric, or gene names. The default is set to `NULL`.
#' @param pcorr.select.genes Logical. Whether to pre-filter genes by correlation cutoff. The default is set to `TRUE`.
#' @param pcorr.reference.data Character. Assay to use for gene selection. The default is set to `NULL`.
#' @param pcorr.corr.cutoff Numeric. Correlation cutoff for selecting genes. The default is set to 0.7.
#' @param pcorr.filter.genes Logical. Whether to filter pairwise correlations before plotting. The default is set to `TRUE`.
#' @param pcorr.corr.dif.cutoff Numeric. Cutoff for difference in correlation. The default is set to 0.3.
#' @param pcorr.plot.ncol Numeric. Columns in correlation plot grid. Used when multiple assays selected.
#' @param pcorr.plot.nrow Numeric. Rows in correlation plot grid. Used when multiple assays selected.
#' @param deg.method Character. Differential expression method. The default is set to `limma`.
#' @param deg.plot.ncol Numeric. Columns in DEG p-value plot grid. The default is set to 1.
#' @param deg.plot.nrow Numeric. Rows in DEG p-value plot grid. The default is set to 3.
#' @param gene.set.score.reference.data Character. Column name or assay name used as reference for gene set scoring.
#' @param gene.set.score.regress.out.variables Character. Column(s) to regress out before scoring. The default is set to `NULL`.
#' @param gene.set.score.list List. Gene sets used for enrichment analysis.
#' @param gene.set.score.normalization Character. Normalization method for gene set scoring. The default is set to `NULL`.
#' @param gene.set.score.assessment TTT
#' @param gene.set.score.variables.to.assess TTT
#' @param check.se.obj Logical. Whether to check the SummarizedExperiment object.
#' @param remove.na Character. Whether to remove `NA` values from assays. Options: `assays`, `none`. The default is set to `assays`.
#' @param override.check Logical. Skip recalculating metrics if already present. The default is set to `FALSE`.
#' @param verbose Logical. Print progress messages. The default is set to `FALSE`.
#'
#' @return A SummarizedExperiment object containing assessment metrics and plots. Optionally saves a PDF of results.
#'
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom SummarizedExperiment assays colData
#' @importFrom gridExtra grid.arrange grid.table
#' @importFrom ggforestplot geom_stripes
#' @importFrom graphics plot.new text
#' @export

assessVariation <- function(
        se.obj,
        assay.names = 'all',
        variables,
        plots.to.exclude = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        general.points.size = 1.5,
        rle.iqr.width = 2,
        rle.median.points.size = 1,
        rle.median.points.color = 'red',
        rle.geom.hline.color = "cyan",
        rle.plot.ncol = 2,
        rle.plot.nrow = 3,
        rle.var.plot.ncol = 3,
        rle.var.plot.nrow = 2,
        rle.colors = NULL,
        fast.pca = TRUE,
        compute.nb.pcs = 10,
        nb.pcs.toplot.pca = 3,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        pca.variables.colors = NULL,
        pca.plot.nrow = 2,
        pca.plot.ncol = 3,
        pca.stroke.size = 0.05,
        pca.stroke.color = 'grey',
        pca.points.alpha = 0.5,
        pca.densities.alpha = 0.5,
        sil.dist.measure = 'euclidian',
        sli.nb.pcs = 3,
        ari.clustering.method = "hclust",
        ari.hclust.method = "complete",
        ari.hclust.dist.measure = "euclidian",
        ari.nb.pcs = 3,
        vca.nb.pcs = 10,
        lra.nb.pcs = 10,
        corr.method = 'spearman',
        a = 0.05,
        rho = 0,
        correlation.plot.ncol = 1,
        correlation.plot.nrow = 3,
        anova.method = 'aov',
        anova.plot.ncol = 1,
        anova.plot.nrow = 3,
        pcorr.method = 'spearman',
        pcorr.genes = NULL,
        pcorr.select.genes = FALSE,
        pcorr.reference.data = NULL,
        pcorr.corr.cutoff = 0.6,
        pcorr.filter.genes = TRUE,
        pcorr.corr.dif.cutoff = 0.1,
        pcorr.plot.ncol = 2,
        pcorr.plot.nrow = 2,
        deg.method = 'limma',
        deg.plot.ncol = 1,
        deg.plot.nrow = 1,
        gene.set.score.reference.data = NULL,
        gene.set.score.regress.out.variables = NULL,
        gene.set.score.list,
        gene.set.score.normalization = NULL,
        gene.set.score.assessment = FALSE,
        gene.set.score.variables.to.assess = NULL,
        check.se.obj = TRUE,
        remove.na = 'none',
        override.check = FALSE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The assessVariation function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the inputs of function ####
    if(!is.vector(assay.names)  | is.logical(assay.names) ){
        stop('The "assay.names" must be a single assay name or assay.names = "all" in the SummarizedExperiment object.')
    }
    if (length(assay.names) == 1 && assay.names != 'all') {
        if (!assay.names %in% names(assays(se.obj)))
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    if (length(assay.names) > 1) {
        if (length(setdiff(assay.names, names(assays(se.obj)))) > 0)
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    if (is.null(variables)){
        stop('The "variables" cannot be empty or NULL.')
    }
    if (is.logical(plots.to.exclude)){
        stop('The "to.exclude" must be a vector or NULL.')
    }
    if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE".')
    }
    if (isTRUE(apply.log)){
        if(length(pseudo.count) > 1 | pseudo.count < 0 | is.null(pseudo.count))
            stop('The "pseudo.count" must be 0 or a postive numeric value.')
    }
    if (isFALSE(is.logical(fast.pca))) {
        stop('The "fast.pca" must be "TRUE" or "FALSE".')
    }
    if (compute.nb.pcs < 0 | is.logical(compute.nb.pcs)){
        stop('The "compute.nb.pcs" must be a postive numeric value.')
    }
    if (isFALSE(is.logical(scale))) {
        stop('The "scale" must be "TRUE" or "FALSE".')
    } else if (isFALSE(is.logical(center))) {
        stop('The "center" must be "TRUE" or "FALSE".')
    }

    if (isTRUE(sli.nb.pcs) > compute.nb.pcs){
        stop('The "sli.nb.pcs" cannot be larger than "compute.nb.pcs".')
    }
    if (isTRUE(ari.nb.pcs) > compute.nb.pcs){
        stop('The "ari.nb.pcs" cannot be larger than "compute.nb.pcs".')
    }
    if (isTRUE(nb.pcs.toplot.pca) > compute.nb.pcs){
        stop('The "nb.pcs.toplot.pca" cannot be larger than "compute.nb.pcs".')
    }
    if (isTRUE(vca.nb.pcs) > compute.nb.pcs){
        stop('The "vca.nb.pcs" cannot be larger than "compute.nb.pcs".')
    }
    if (isTRUE(lra.nb.pcs) > compute.nb.pcs){
        stop('The "lra.nb.pcs" cannot be larger than "compute.nb.pcs".')
    }
    if (fast.pca & is.null(compute.nb.pcs)) {
        stop('To perform fast PCA, the number of PCs (compute.nb.pcs) must specified.')
    } else if (fast.pca & compute.nb.pcs == 0) {
        stop('To perform fast PCA, the number of PCs (compute.nb.pcs) must specified.')
    }
    if (fast.pca & is.null(compute.nb.pcs)) {
        stop('To perform fast PCA, the number of PCs (compute.nb.pcs) must specified.')
    } else if (fast.pca & compute.nb.pcs == 0) {
        stop('To perform fast PCA, the number of PCs (compute.nb.pcs) must specified.')
    }

    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Assessing the SummarizedExperiment object ####
    if (isTRUE(check.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            variables = variables,
            remove.na = remove.na,
            verbose = verbose)
    }
    # Getting all possible metrics for each variable #####
    printColoredMessage(
        message = '-- Get all possible assessment metrics and plot:',
        color = 'magenta',
        verbose = verbose)
    se.obj <- getAssessmentMetrics(
        se.obj = se.obj,
        variables = variables,
        plot.output = FALSE,
        save.se.obj = TRUE
        )
    ## Metrics and plots to generate #####
    metrics.table <- se.obj@metadata$AssessmentMetrics$metrics.table
    if(!is.null(plots.to.exclude)){
        printColoredMessage(
            message = paste0('- Exclude all the specified codes'),
            color = 'blue',
            verbose = verbose
            )
        metrics.table <- metrics.table[!metrics.table$Code %in% plots.to.exclude, ]
    }
    n.plots <- length(unique(paste0(metrics.table$Metrics, metrics.table$PlotTypes)))
    printColoredMessage(
        message = paste0('- In total, ',n.plots, ' assessment plots will be generated.'),
        color = 'blue',
        verbose = verbose)

    # RLE #####
    ## compute rle #####
    if('RLE' %in% metrics.table$Metrics){
        if('rlePlot' %in% metrics.table$PlotTypes ){
            rle.outputs.to.return <- 'all'
        } else rle.outputs.to.return <- 'rle.med.iqr'
        se.obj <- computeRLE(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            outputs.to.return = rle.outputs.to.return,
            check.se.obj = FALSE,
            remove.na = 'none',
            save.se.obj = TRUE,
            override.check = override.check,
            verbose = verbose
            )
    }

    ## plot general rle #####
    if('rlePlot' %in% metrics.table$PlotTypes){
        se.obj <- plotRLE(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            variable = NULL,
            variable.colors = NULL,
            ylim.rle.plot = NULL,
            iqr.width = rle.iqr.width,
            median.points.size = rle.median.points.size,
            median.points.color = rle.median.points.color,
            geom.hline.color = rle.geom.hline.color,
            plot.ncol = rle.plot.ncol,
            plot.nrow = rle.plot.nrow,
            plot.output = FALSE,
            save.se.obj = TRUE,
            verbose = verbose
            )
    }

    ## plot colored rle #####
    if('coloredRLEplot' %in% metrics.table$PlotTypes){
        coloredRLEplot.vars <- metrics.table$PlotTypes == 'coloredRLEplot'
        coloredRLEplot.vars <- metrics.table$Variables[coloredRLEplot.vars]
        for(i in coloredRLEplot.vars){
            se.obj <- plotRLE(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                variable.colors = rle.colors,
                ylim.rle.plot = NULL,
                iqr.width = rle.iqr.width,
                median.points.size = rle.median.points.size,
                median.points.color = "grey",
                geom.hline.color = rle.geom.hline.color,
                plot.ncol = rle.plot.ncol,
                plot.nrow = rle.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose
                )
        }
    }
    ## plot rle medians with variable #####
    if('rleMedians' %in% metrics.table$Factors){
        rleMedplot.vars <- metrics.table$Factors == 'rleMedians'
        rleMedplot.vars <- metrics.table$Variables[rleMedplot.vars]
        for(i in rleMedplot.vars){
            se.obj <- plotRleVariable(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                rle.data.type = 'rle.medians',
                ylim.rle.med.plot = NULL,
                ylim.rle.iqr.plot = NULL,
                points.size = general.points.size,
                plot.ncol = rle.var.plot.ncol,
                plot.nrow = rle.var.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose
                )
        }
    }

    ## plot rle iqr with variable #####
    if('rleIqr' %in% metrics.table$Factors){
        rleIqrplot.vars <- metrics.table$Factors == 'rleIqr'
        rleIqrplot.vars <- metrics.table$Variables[rleIqrplot.vars]
        for(i in rleIqrplot.vars){
            se.obj <- plotRleVariable(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                rle.data.type = 'rle.iqrs',
                ylim.rle.med.plot = NULL,
                ylim.rle.iqr.plot = NULL,
                points.size = general.points.size,
                plot.ncol = rle.var.plot.ncol,
                plot.nrow = rle.var.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    # PCA ####
    ## compute pca ####
    if(sum(c('PCA', 'LRA', 'VCA', 'ARI', 'Silhouette')  %in% metrics.table$Metrics) > 0 ) {
        se.obj <- RUVIIIPRPS::computePCA(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            fast.pca = fast.pca,
            nb.pcs = compute.nb.pcs,
            center = center,
            scale = scale,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            svd.bsparam = svd.bsparam,
            check.se.obj = FALSE,
            remove.na = 'none',
            save.se.obj = TRUE,
            override.check = override.check,
            verbose = verbose)
    }
    ## scatter plot pca ####
    if('PCA' %in% metrics.table$Metrics & 'scatterPlot' %in% metrics.table$PlotTypes){
        pca.scatter.vars <- metrics.table$Metrics == 'PCA' & metrics.table$PlotTypes == 'scatterPlot'
        pca.scatter.vars <- metrics.table$Variables[pca.scatter.vars]
        for(i in pca.scatter.vars){
            se.obj <- RUVIIIPRPS::plotPCA(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = nb.pcs.toplot.pca,
                plot.type = "scatter",
                variable.colors = pca.variables.colors,
                points.size = general.points.size,
                stroke.color = pca.stroke.color,
                stroke.size = pca.stroke.size,
                points.alpha = pca.points.alpha,
                densities.alpha = pca.densities.alpha,
                plot.ncol = pca.plot.nrow,
                plot.nrow = pca.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = TRUE)
        }
    }
    ## boxplot pca ####
    if('PCA' %in% metrics.table$Metrics & 'boxPlot' %in% metrics.table$PlotTypes){
        pca.boxplot.vars <- metrics.table$Metrics == 'PCA' & metrics.table$PlotTypes == 'boxPlot'
        pca.boxplot.vars <- metrics.table$Variables[pca.boxplot.vars]
        for(i in pca.boxplot.vars){
            se.obj <- RUVIIIPRPS::plotPCA(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = nb.pcs.toplot.pca,
                plot.type = "boxplot",
                variable.colors = pca.variables.colors,
                points.size = general.points.size,
                stroke.color = pca.stroke.color,
                stroke.size = pca.stroke.size,
                points.alpha = pca.points.alpha,
                densities.alpha = pca.densities.alpha,
                plot.ncol = pca.plot.nrow,
                plot.nrow = pca.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = TRUE)
        }
    }
    # Vector correlation ####
    ## compute vector correlation ####
    if('VCA' %in% metrics.table$Metrics){
        pc.vec.corr.vars <- metrics.table$Metrics == 'VCA'
        pc.vec.corr.vars <- metrics.table$Variables[pc.vec.corr.vars]
        for(i in pc.vec.corr.vars){
            se.obj <- computePCVariableCorrelation(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = vca.nb.pcs,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    ## plot vector correlation ####
    if('VCA' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'VCA'
        pc.vec.corr.vars <- metrics.table$Variables[index]
        for(i in pc.vec.corr.vars){
            se.obj <- plotPCVariableCorrelation(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = vca.nb.pcs,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    # Linear regression ####
    ## compute linear regression ####
    if('LRA' %in% metrics.table$Metrics){
        pc.reg.vars <- metrics.table$Metrics == 'LRA'
        pc.reg.vars <- metrics.table$Variables[pc.reg.vars]
        for(i in pc.reg.vars){
            se.obj <- computePCVariableRegression(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = lra.nb.pcs,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    ## plot linear regression ####
    if('LRA' %in% metrics.table$Metrics){
        pc.reg.vars <- metrics.table$Metrics == 'LRA'
        pc.reg.vars <- metrics.table$Variables[pc.reg.vars]
        for(i in pc.reg.vars){
            se.obj <- plotPCVariableRegression(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = lra.nb.pcs,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }

    # Silhouette coefficient ####
    ## compute silhouette coefficients ####
    if('Silhouette' %in% metrics.table$Metrics){
        index.single <- metrics.table$Metrics == 'Silhouette' &
            metrics.table$PlotTypes == 'barPlot'
        sil.single.vars <- metrics.table$Variables[index.single]
        index.combined <- metrics.table$Metrics == 'Silhouette' &
            metrics.table$PlotTypes == 'combinedPlot'
        sil.combined.vars <- metrics.table$Variables[index.combined]
        sil.combined.vars <- unlist(strsplit(x = sil.combined.vars, split = '&'))
        all.sil.vars <- unique(c(sil.single.vars, sil.combined.vars))
        for(i in all.sil.vars){
            se.obj <- computeSilhouette(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                dist.measure = sil.dist.measure,
                fast.pca = fast.pca,
                nb.pcs = sli.nb.pcs,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    ## barplot silhouette coefficients  ####
    if('Silhouette' %in% metrics.table$Metrics & 'barPlot' %in% metrics.table$PlotTypes){
        index <- metrics.table$Metrics == 'Silhouette' & metrics.table$PlotTypes == 'barPlot'
        sil.single.vars <- metrics.table$Variables[index]
        for(i in sil.single.vars){
            se.obj <- plotSilhouette(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = i,
                plot.type = 'single.plot',
                silhouette.method = paste0('sil.', sil.dist.measure),
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    ## combined plot silhouette coefficients  ####
    if('Silhouette' %in% metrics.table$Metrics & 'combinedPlot' %in% metrics.table$PlotTypes){
        index <- metrics.table$Metrics == 'Silhouette' & metrics.table$PlotTypes == 'combinedPlot'
        sil.combined.vars <- metrics.table$Variables[index]
        for(i in sil.combined.vars){
            se.obj <- plotSilhouette(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = strsplit(x = i, split = '&')[[1]],
                plot.type = 'combined.plot',
                silhouette.method = paste0('sil.', sil.dist.measure),
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }

    # ARI ####
    ## compute adjusted rand index ####
    if('ARI' %in% metrics.table$Metrics){
        index.single <- metrics.table$Metrics == 'ARI' &
            metrics.table$PlotTypes == 'barPlot'
        ari.single.vars <- metrics.table$Variables[index.single]
        index.combined <- metrics.table$Metrics == 'ARI' &
            metrics.table$PlotTypes == 'combinedPlot'
        ari.combined.vars <- metrics.table$Variables[index.combined]
        ari.combined.vars <- unlist(strsplit(x = ari.combined.vars, split = '&'))
        all.ari.vars <- unique(c(ari.single.vars, ari.combined.vars))
        for(i in all.ari.vars){
            se.obj <- computeARI(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                clustering.method = ari.clustering.method,
                hclust.method = ari.hclust.method,
                hclust.dist.measure = ari.hclust.dist.measure,
                fast.pca = fast.pca,
                nb.pcs = ari.nb.pcs,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    ## barplot adjusted rand index  ####
    if('ARI' %in% metrics.table$Metrics & 'barPlot' %in% metrics.table$PlotTypes){
        index <- metrics.table$Metrics == 'ARI' & metrics.table$PlotTypes == 'barPlot'
        ari.single.vars <- metrics.table$Variables[index]
        if(ari.clustering.method == 'mclust'){
            ari.method <- 'mclust'
        } else ari.method <- paste0('hclust.', ari.hclust.method, '.', ari.hclust.dist.measure)
        for(i in ari.single.vars){
            se.obj <- plotARI(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = i,
                plot.type = 'single.plot',
                ari.method = ari.method,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    ## combined adjusted rand index  ####
    if('ARI' %in% metrics.table$Metrics & 'combinedPlot' %in% metrics.table$PlotTypes){
        index <- metrics.table$Metrics == 'ARI' & metrics.table$PlotTypes == 'combinedPlot'
        ari.combined.vars <- metrics.table$Variables[index]
        if(ari.clustering.method == 'mclust'){
            ari.method <- 'mclust'
        } else ari.method <- paste0('hclust.', ari.hclust.method, '.', ari.hclust.dist.measure)
        for(i in ari.combined.vars){
            se.obj <- plotARI(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = strsplit(x = i, split = '&')[[1]],
                plot.type = 'combined.plot',
                ari.method = ari.method,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    # Gene variable correlation ####
    ## compute gene variable correlations ####
    if('Correlation' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'Correlation'
        gene.var.corr.vars <- unique(metrics.table$Variables[index])
        for(i in gene.var.corr.vars){
            se.obj <- computeGenesVariableCorrelation(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                method = corr.method,
                a = a,
                rho = rho,
                plot.top.genes = FALSE,
                nb.top.genes = NULL,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                apply.round = TRUE,
                check.se.obj = FALSE,
                override.check = override.check,
                remove.na = 'none',
                save.se.obj = TRUE)
        }
    }
    ## plot gene variable correlations ####
    if('Correlation' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'Correlation'
        gene.var.corr.vars <- unique(metrics.table$Variables[index])
        for(i in gene.var.corr.vars){
            se.obj <- plotGenesVariableCorrelation(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                correlation.method = corr.method,
                plot.ncol = correlation.plot.ncol,
                plot.nrow = correlation.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
            }
    }

    # Gene variable anova ####
    ## compute gene variable anova ####
    if('ANOVA' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'ANOVA'
        gene.var.anova.vars <- unique(metrics.table$Variables[index])
        for(i in gene.var.anova.vars){
            se.obj <- computeGenesVariableAnova(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                method = anova.method,
                plot.top.genes = FALSE,
                nb.top.genes = NULL,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                apply.round = TRUE,
                check.se.obj = FALSE,
                remove.na = 'none',
                override.check = override.check,
                save.se.obj = TRUE)
        }
    }
    ## plot gene variable anova ####
    if('ANOVA' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'ANOVA'
        gene.var.corr.vars <- unique(metrics.table$Variables[index])
        for(i in gene.var.corr.vars){
            se.obj <- plotGenesVariableAnova(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                anova.method = anova.method,
                plot.ncol = anova.plot.ncol,
                plot.nrow = anova.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = TRUE)
        }
    }

    # DGE ####
    ## compute dge ####
    if('DGE' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'DGE'
        dge.vars <- unique(metrics.table$Variables[index])
        for(i in dge.vars){
            se.obj <- computeDGE(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                method = deg.method,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                save.se.obj = TRUE,
                override.check = override.check,
                verbose = verbose
                )
        }
    }
    ## plot p-value hist ####
    if('DGE' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'DGE'
        dge.vars <- unique(metrics.table$Variables[index])
        for(i in dge.vars){
            se.obj <- plotDGE(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                method = deg.method,
                plot.ncol = deg.plot.ncol,
                plot.nrow = deg.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose
                )
        }
    }

    # Partial correlation ####
    ## compute partial correlation ####
    if('PartialCorrelation' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'PartialCorrelation'
        pcorr.vars <- unique(metrics.table$Variables[index])
        for(i in pcorr.vars){
            se.obj <- computeGenesPartialCorrelation(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                method = pcorr.method,
                genes = pcorr.genes[[i]],
                select.genes = pcorr.select.genes,
                reference.data = pcorr.reference.data,
                corr.coff.cutoff = pcorr.corr.cutoff,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                apply.round = TRUE,
                check.se.obj = FALSE,
                remove.na = 'none',
                save.se.obj = TRUE,
                override.check = override.check,
                verbose = verbose
            )
        }
    }
    ## plot partial correlation ####
    if('PartialCorrelation' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'PartialCorrelation'
        pcorr.vars <- unique(metrics.table$Variables[index])
        for(i in pcorr.vars){
            index <- metrics.table$Metrics == 'PartialCorrelation' & metrics.table$Variables == i
            plot.types <- metrics.table$PlotTypes[index]
                for(j in plot.types){
                    if(j == 'scatterPlot'){
                        plot.type = 'scatter.plot'
                    } else if(j == 'barPlot'){
                        plot.type = 'barplot'
                    } else if (j == 'histogram') plot.type = 'histogram'
                    se.obj <- plotGenesPartialCorrelation(
                        se.obj = se.obj,
                        assay.names = levels(assay.names),
                        variable = i,
                        method = pcorr.method,
                        plot.type = plot.type,
                        filter.genes = pcorr.filter.genes,
                        corr.dif.cutoff = pcorr.corr.dif.cutoff,
                        plot.ncol = pcorr.plot.ncol,
                        plot.nrow = pcorr.plot.nrow,
                        plot.output = FALSE,
                        save.se.obj = TRUE,
                        verbose = verbose)
                }
            }
    }
    # Gene set scoring ####
    ## compute gene set scoring ####
    if('GeneSetScore' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'GeneSetScore'
        gene.set.vars <- unique(metrics.table$Variables[index])
        gene.set.vars <- gene.set.vars[gene.set.vars %in% names(gene.set.score.list)]
        if(length(gene.set.vars) > 0){
            for(i in gene.set.vars){
                se.obj <- computeGeneSetScore(
                    se.obj = se.obj,
                    assay.names = levels(assay.names),
                    upset.genes = gene.set.score.list[[i]]$upset.genes,
                    downset.genes = gene.set.score.list[[i]]$downset.genes,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    normalization = gene.set.score.normalization,
                    regress.out.variables = gene.set.score.regress.out.variables,
                    assess.score = gene.set.score.assessment ,
                    variables.to.assess = gene.set.score.variables.to.assess,
                    corr.method = 'spearman',
                    gene.set.name = i,
                    plot.output = FALSE,
                    check.se.obj = FALSE,
                    save.se.obj = TRUE,
                    verbose = verbose
                )
            }
        }
    }

    ## plot partial correlation ####
    if('GeneSetScore' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'GeneSetScore'
        gene.set.vars <- unique(metrics.table$Variables[index])
        gene.set.vars <- gene.set.vars[gene.set.vars %in% names(gene.set.score.list)]
        for(i in gene.set.vars){
            se.obj <- plotGeneSetScore(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                reference.score = gene.set.score.reference.data,
                gene.set.name = i,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose
                )
        }
    }
    printColoredMessage(message = '------------The assessVariation function finished.',
                        color = 'white',
                        verbose = verbose)
    return(se.obj)
}

