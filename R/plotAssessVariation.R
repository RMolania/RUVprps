#' Plot the results of the assess variation .
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character or character vector. Specifies the name(s) of the assay(s) in the
#' SummarizedExperiment object. The default is set to `all`, so all the assays in the SummarizedExperiment object will
#' be selected.
#' @param variables Character or character vector. Specifies the name(s) of the column(s) in the SummarizedExperiment
#' object for which the specified metrics below have been calculated.
#' @param fast.pca Logical. Indicates whether the fast SVD approach was used in the `computePCA` function. The default
#' is set to `TRUE`.
#' @param anova.method Character. Indicates which ANOVA method was used to compute the association between
#' gene-level expression and a categorical variable. Options are: `aov` and `welch`. Default is `aov`.
#' @param corr.method Character. Specifies which correlation method to use to compute correlation between
#' gene-level expression and a continuous variable. Options are: `pearson`, `kendall`, or `spearman`. Default is `spearman`.
#' @param pcorr.method Character. Indicates which correlation method to use to compute gene-gene partial
#' correlation. Options are: `pearson`, `kendall`, or `spearman`. Default is `spearman`.
#' @param sil.dist.measure Character. Indicates which distance measure to apply on the PCs to calculate silhouette scores.
#' Options are: `euclidean`, `maximum`, `manhattan`, `canberra`, `binary`, or `minkowski`. Default is `euclidean`.
#' Refer to the function `computeSilhouette` for more details.
#' @param ari.clustering.method Character. Indicates which clustering method to apply on the PCs to
#' calculate ARI. Options are: `mclust` or `hclust`. Default is `hclust`. Refer to the `computeARI` function for more details.
#' @param ari.hclust.method Character. Specifies the agglomeration method to use when `clustering.method`
#' is set to `hclust` for computing ARI. Options are: `ward.D`, `ward.D2`, `single`, `complete`, `average` (UPGMA),
#' `mcquitty` (WPGMA), `median` (WPGMC), or `centroid` (UPGMC). Default is `complete`. Refer to the `computeARI`
#' function for more details.
#' @param ari.hclust.dist.measure Character. Specifies the distance measure to use in the `dist` function when
#' `clustering.method` is set to `hclust` for computing ARI. Options are: `euclidean`, `maximum`, `manhattan`,
#' `canberra`, `binary`, or `minkowski`. Default is `euclidean`. Refer to the `computeARI` function for more details.
#' @param output.file.name Character. The path and name of the output file where the assessment plots will be saved in
#' PDF format.
#' @param pdf.width Numeric. Specifies the width (in inches) of the output PDF file. Default is typically around 8.
#' @param pdf.height Numeric. Specifies the height (in inches) of the output PDF file. Default is typically around 6.
#' @param verbose Logical. If `TRUE`, displays messages describing different steps of the function.
#' @importFrom SummarizedExperiment assays
#' @export


plotAssessVariation <- function(
        se.obj,
        assay.names = 'all',
        variables,
        fast.pca = TRUE,
        anova.method = 'aov',
        corr.method = 'spearman',
        pcorr.method = 'spearman',
        sil.dist.measure = 'euclidian',
        ari.clustering.method = "hclust",
        ari.hclust.method = "complete",
        ari.hclust.dist.measure = "euclidian",
        output.file.name = NULL,
        pdf.width = 12,
        pdf.height = 12,
        verbose = TRUE
        ){
    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    metrics.table <- se.obj@metadata$AssessmentMetrics$metrics.table

    # Putting all plots together ####
    ## find classes of different variables ####
    if (!is.null(variables)) {
        vars.class <- sapply(
            variables,
            function(x) class(colData(se.obj)[[x]]))
        categorical.vars <- names(vars.class[vars.class %in% c('character', 'factor')])
        continuous.vars <- names(vars.class[vars.class %in% c('numeric', 'integer')])
    }
    ## selecting output file names ####
    if (is.null(output.file.name)){
        output.file.name <- 'RUVIIIPRPS_AssessVariation'
    }
    pdf(paste0(output.file.name, '.pdf'),
        width = pdf.width,
        height = pdf.height
        )
    plot.new()
    text(.5, .7, "Assess variation", font = 2, cex = 2.5)
    text(.5, .6, "variables:", font = 2, cex = 2.5)
    text(.5, .4, paste0(variables, collapse = '\n'), font = 2, cex = 2)
    print(se.obj@metadata$AssessmentMetrics$plot)

    ## general RLE plot ####
    if ('rlePlot' %in% metrics.table$PlotTypes) {
        if (length(assay.names) > 1){
            print(
                se.obj@metadata$Plots$global.level$RLE$gene.median.center$general$un.colored
            )
        }
    }
    ## continuous variables ####
    for(i in continuous.vars){
        plot.new()
        text(.5, .7, paste0("Assess variation \n in the variable: \n ", i ), font = 2, cex = 2.5)
        metrics.table.var <- metrics.table[metrics.table$Variables == i, ]
        ### scatter plot between RLE medians and variable ####
        if ('rleMedians' %in% metrics.table.var$Factors){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$corr.medians.variable[[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$corr.medians.variable[[i]]$scatter.plot
            )
        }
        ## scatter plot between RLE iqrs and variable ####
        if ('rleIqr' %in% metrics.table.var$Factors ){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$corr.iqrs.variable[[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$corr.iqr.variable[[i]]$scatter.plot
            )
        }
        ## scatter plot between PCs and variable ####
        if (isTRUE(fast.pca)){
            svd.method <- 'fast.svd'
        } else svd.method <- 'ordinary.svd'
        if ('PCA' %in% metrics.table.var$Metrics){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$PCA[[svd.method]][[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$PCA[[svd.method]][[i]]$scatter.plot
            )
        }
        ## line-dot plot for linear regression ####
        if ('LRA' %in% metrics.table.var$Metrics){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$LRA[[svd.method]][[i]]$line.dotplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$LRA[[svd.method]][[i]]$plot
            )
        }
        ## boxplot of gene variable correlation coefficients ####
        if ('CorrelationboxPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$Correlation[[corr.method]][[i]]$cor.coef.boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$Correlation[[corr.method]][[i]]$cor.coef.boxplot
            )
        }
        ## histograms of gene variable correlation coefficients ####
        if ('CorrelationpvalHist' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$Correlation[[corr.method]][[i]]$cor.coef.boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$Correlation[[corr.method]][[i]]$cor.coef.boxplot
            )
        }
        ## scatter plot of gene variable partial correlation coefficients ####
        if ('PartialCorrelationscatterPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$PPcorr[[pcorr.method]][[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$PPcorr[[pcorr.method]][[i]]$scatter.plot
            )
        }
        ## barplot of gene variable partial correlation coefficients ####
        if ('PartialCorrelationbarPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$PPcorr[[pcorr.method]][[i]]$barplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$PPcorr[[pcorr.method]][[i]]$barplot
            )
        }
        ## histogram of gene variable partial correlation coefficients ####
        if ('PartialCorrelationhistogram' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$PPcorr[[pcorr.method]][[i]]$histogram
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$PPcorr[[pcorr.method]][[i]]$histogram
            )
        }
        ## histogram of gene variable partial correlation coefficients ####
        if ('GeneSetScore' %in% metrics.table.var$Metrics){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$GeneSetSocore$singscore[[i]]$general
                )
            }
        }

        if ('LISI' %in% metrics.table.var$Metrics){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$LISI[[svd.method]][[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$LISI[[svd.method]][[i]]$boxplot
            )
        }
    }

    ## categorical variables ####
    for(i in categorical.vars){
        plot.new()
        text(.5, .7, paste0("Assess variation \n in the variable: \n ", i ), font = 2, cex = 2.5)
        metrics.table.var <- metrics.table[metrics.table$Variables == i, ]
        ### rle plots colored by the variable ####
        if ('coloredRLEplot' %in% metrics.table.var$PlotTypes){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$gene.median.center[[i]]$colored
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$gene.median.center[[i]]$colored
            )
        }
        ## boxplot between RLE medians and variable ####
        if ('rleMedians' %in% metrics.table.var$Factors){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$corr.medians.variable[[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$corr.medians.variable[[i]]$boxplot
            )
        }
        ## boxplot between RLE iqrs and variable ####
        if ('rleIqr' %in% metrics.table.var$Factors){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$RLE$corr.iqrs.variable[[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$RLE$corr.iqr.variable[[i]]$boxplot
            )
        }
        ## scatter plot of PCs colored by the variable ####
        if (isTRUE(fast.pca)){
            svd.method <- 'fast.svd'
        } else svd.method <- 'ordinary.svd'
        if ('pcsscatterPlot' %in% paste0(metrics.table.var$Factors, metrics.table.var$PlotTypes)){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$PCA[[svd.method]][[i]]$boxplot.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$PCA[[svd.method]][[i]]$boxplot.plot
            )
        }
        ## boxplot plot of PCs colored by the variable ####
        if ('pcsboxPlot' %in% paste0(metrics.table.var$Factors, metrics.table.var$PlotTypes) ){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$PCA[[svd.method]][[i]]$scatter.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$PCA[[svd.method]][[i]]$scatter.plot
            )
        }
        ## line-dot plot for vector correlation ####
        if ('VCA' %in% metrics.table.var$Metrics ){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$VCA[[svd.method]][[i]]$line.dotplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$VCA[[svd.method]][[i]]$plot
            )
        }
        ## barplot of ari ####
        if ('ARIbarPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)) {
            if (ari.clustering.method == 'mclust'){
                ari.method <- 'mclust'
            } else ari.method <- paste0('hclust.', ari.hclust.method, '.', ari.hclust.dist.measure)
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$ARI[[ari.method]][[i]]$single.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$ARI[[ari.method]][[i]]$single.plot
            )
        }
        ## barplot of silhouette ####
        if ('SilhouettebarPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes)){
            silhouette.method <- paste0('sil.', sil.dist.measure)
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$Silhouette[[silhouette.method]][[i]]$single.plot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$Silhouette[[silhouette.method]][[i]]$single.plot
            )
        }
        ## boxplot of gene variable ANOVA F-stat ####
        if ('ANOVAboxPlot' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes) ){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$ANOVA[[anova.method]][[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$ANOVA[[anova.method]][[i]]$boxplot
            )
        }
        ## histograms of gene variable ANOVA p-values ####
        if ('ANOVApvalHist' %in% paste0(metrics.table.var$Metrics, metrics.table.var$PlotTypes) ){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$ANOVA[[anova.method]][[i]]$histogram
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$ANOVA[[anova.method]][[i]]$histogram
            )
        }
        ## histograms of differential gene expression p-values  ####
        if ('DGE' %in% metrics.table.var$Metrics ){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$gene.level$DGE$Wilcoxon[[i]]$histogram
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$gene.level$DGE$Wilcoxon$time.interval$plot
            )
        }
        if ('LISI' %in% metrics.table.var$Metrics){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$LISI[[svd.method]][[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$LISI[[svd.method]][[i]]$boxplot
            )
        }
        if ('KBET' %in% metrics.table.var$Metrics){
            if (length(assay.names) > 1){
                print(
                    se.obj@metadata$Plots$global.level$KBET[[svd.method]][[i]]$boxplot
                )
            } else print(
                se.obj@metadata$Metrics[[assay.names]]$global.level$KBET[[svd.method]][[i]]$boxplot
            )
        }
    }
    dev.off()
}


