#' Assess the variation in variables.
#'
#' @author Ramyar Molania
#'
#' @references
#' Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids Research, 2019.
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS, Nature Biotechnology, 2023
#'
#' @description
#' This function applies a range of global and gene-level metrics to assess the impact variables on gene expression data
#' in a SummarizedExperiment object.
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
#' @param assay.names Character. A character string or vector for selecting dataset name(s) in the SummarizedExperiment
#' object. The default is set to `all`. All available dataset(s) wiil be selected.
#' @param bio.variables Character. A character string or vector of column names in the sample annotation that containing
#' biological variables to assess their impact on the dataset(s). The variables can be categorical or continuous.
#' @param uv.variables Character. A character string or vector of column names in the sample annotation that containing
#' biological variables to assess their impact on the dataset(s). The variables can be categorical or continuous.
#' @param assessment.level Character. A character string specifying the level of assessment. The options are `L1` amd `L2`.
#' The default is set to `L2`. Refer to the details of the function for more information.
#' @param plots.to.exclude Character. A character string or vector indicating which metrics to exclude. Use the
#' `getAssessmentMetrics()` function to view all the possible metrics. The default is set to `NULL`.
#' @param apply.log Logical. Whether to apply a log-transformation on the dataset(s) before applying the metrics. The
#' default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value as pseudo count added to all measurements before log-transformation. The
#' default is set to 1.
#' @param general.points.size Numeric. Size of points in scatter plots. The default is set to 1.5.
#' @param rle.iqr.width Numeric. Width of IQR in RLE plots. The default is set to 2.
#' @param rle.median.points.size Numeric. Size of median points in RLE plots. The default is set to 1.
#' @param rle.median.points.color Character. Color of median points in RLE plots. The default is set to `red`.
#' @param rle.geom.hline.color Character. Color of the horizontal line in RLE plots. The default is set to `cyan`.
#' @param rle.plot.ncol Numeric. Number of columns when putting more than one RLE plots in a grid. The default is set to 2.
#' @param rle.plot.nrow Numeric. Number of rows when putting more than one RLE plots in a grid. The default is set to 3.
#' @param rle.var.plot.ncol Numeric. Number of columns when putting more than one scatter plot of RLE medians or QR against
#' a variable in a grid. The default is set to 3.
#' @param rle.var.plot.nrow Numeric. Number of rows when putting more than one scatter plot of RLE medians or QR against
#' a variable in a grid. The default is set to 3.
#' @param rle.colors List. List of colors for variables in RLE plots. The default is set to `NULL`, the the function will
#' use the default colors.
#' @param fast.pca Logical. This specifies whether to use fast algorithm to do PCA ot not. The default is set to `TRUE`
#' @param compute.nb.pcs Numeric. Number of PCs to compute for PCA. The default is set to 10.
#' @param nb.pcs.toplot.pca Numeric. Number of PCs to plot against each other in PCA plots. The default is set to 3.
#' This must be ≤ `compute.nb.pcs`.
#' @param center Logical. Whether to center data before applying PCA. The default is set to `TRUE`.
#' @param scale Logical. Whether to scale data before SVD. The default is set to `FALSE`.
#' @param svd.bsparam A `BiocParallelParam` object for parallelization. The default is set to `bsparam()`. See `runSVD()`
#' in BiocSingular for more details.
#' @param pca.variables.colors List. List of colors for categorical variables in PCA plots. The default is set to `NULL`,
#' then the function will use the default colors.
#' @param color.palette Character. Name of the color palette to use for plots. The default is set to `nrc`.
#' @param pca.plot.ncol Numeric. Number of columns when putting more than one PCA plot in a grid. The default is set to 2.
#' @param pca.plot.nrow Numeric. Number of rows when putting more than one PCA plot in a grid. The default is set to 3.
#' @param pca.var.plot.ncol Numeric. Number of columns in a grid when more than one  scatter plots of PCs against continuous
#' variables are generated. The default is set to 3.
#' @param pca.var.plot.nrow Numeric. Number of rows in a grid when more than one  scatter plots of PCs against continuous
#' variables are generated. The default is set to 3.
#' @param pca.stroke.size Numeric. Stroke size of points in PCA plots. The default is set to 0.1.
#' @param pca.stroke.color Character. Color of stroke of points in PCA plots. The default is set to `gray`.
#' @param pca.points.alpha Numeric. Transparency of points in PCA plot. The default is set to 0.5.
#' @param pca.densities.alpha Numeric. Transparency of densities in PCA plots. The default is set to 0.5.
#' @param pca.legend.position Character. Position of legend in PCA plots: `top`, `bottom`, `left`, `right`, or `none`.
#' The default is set to `right`.
#' @param sil.dist.measure Character. A character string that specifies how to measure the distance for silhouette. The options
#' are `euclidean`, `maximum`, `manhattan`, `canberra`, `binary`, or `minkowski`. The default is set to `euclidean`.
#' @param sli.nb.pcs Numeric. Number of PCs to use for silhouette calculation. The default is set to 3. The is must be
#' ≤ `compute.nb.pcs`.
#' @param ari.clustering.method Character. A character string that specifies which clustering method to sue for ARI calculation.
#' The options are: `mclust` or `hclust`. The default is set to `hclust`.
#' @param ari.hclust.method Character. A character string that specifies the agglomeration method for `hclust`. The
#' options are `ward.D`, `ward.D2`, `single`, `complete`, `average`, `mcquitty`, `median`, `centroid`. The default is set
#' to `complete`.
#' @param ari.hclust.dist.measure Character. A character string that specifies the distance measure for `hclust`. The options
#' are `euclidean`, `maximum`, `manhattan`, `canberra`, `binary`, `minkowski`. The default is set to `euclidean`.
#' @param ari.nb.pcs Numeric. Number of PCs to use for ARI calculation. The default is set to 3. This must be ≤ `compute.nb.pcs`.
#' @param vca.nb.pcs Numeric. Number of PCs to use for vector correlation calculation. The default is set to 3
#' @param lra.nb.pcs Numeric. Number of PCs to use for linear regression calculation.. The default is set to 3.
#' @param vca.nb.pcs.to.plot Numeric. Number of PCs to include in vector correlation plots. The default is set to 10.
#' @param lra.nb.pcs.to.plot Numeric. Number of PCs to include in linear regression plots. The default is set to 10.
#' @param kbet.nb.pcs Numeric. Number of PCs to use for kBET analysis. The default is set to 10.
#' @param k0 number of nearest neighbors to test on (neighborhood size). The default is set to `NULL`. Refer to the
#' `kBET` function form the kBET package for more information.
#' @param knn an n x k matrix of nearest neighbors for each cell (optional). The default is set to `NULL`. Refer to the
#' `kBET` function form the kBET package for more information.
#' @param lisi.nb.pcs Numeric. Number of PCs to use for LISI analysis. The default is set to 3.
#' @param perplexity Numeric. The effective number of each cell's neighbors. The default is set to 10.
#' @param nn.eps Numeric. Error bound for nearest neighbor search with RANN:nn2(). The default is set to 0.0, implies exact
#' nearest neighbor search.
#' @param corr.method Character. A character string that specifies which correlation method should be used. The options are
#' `pearson`, `kendall`, `spearman`. The default is set to `spearman`.
#' @param a Numeric. Significance level for confidence intervals in correlation. The default is set to 0.05.
#' @param rho Numeric. Hypothesized correlation for testing. The default is set to 0.
#' @param correlation.plot.ncol Numeric. Number of columns of a gird plot when putting more than one correlation plot .
#' The default is set to 3.
#' @param correlation.plot.nrow  Numeric. Number of rows of a gird plot when putting more than one correlation plot .
#' The default is set to
#' @param anova.method Character. A character string that specifies which correlation method should be used. The options
#' are `aov` or `welch`. The default is set to `aov`.
#' @param anova.plot.ncol Numeric. Number of columns of a gird plot when putting more than one ANOVA plot. The default
#' is set to 3
#' @param anova.plot.nrow Numeric. Number of columns of a gird plot when putting more than one ANOVA plot .The default
#' is set to 3.
#' @param pcorr.method Character. A character string that specifies which correlation method should be used for partial
#' correlation. The options are `pearson`, `kendall`, `spearman`. The default is set to `spearman`.
#' @param pcorr.genes Vector. A logical vector that specify genes for pairwise correlation analysis. The default is set
#' to `NULL`, then all genes will be selected.
#' @param pcorr.select.genes Logical. If `TRUE`, the function will compute the correlation between individual genes and the
#' specified variable, then select a subset of genes based on the "corr.coff.cutoff" for downstream analysis. The default
#' is set to `TRUE`. This will speed up the computational time.
#' @param pcorr.reference.data  Character. A character string specifying the name of the data to be used for selecting genes
#' based on the correlation analysis. The default is se to `NULL`, which means all specified assays will be used.
#' @param pcorr.corr.cutoff Numeric. A numeric value used as a cutoff for selecting genes. The default is set to 0.7.
#' @param pcorr.filter.genes Logical. Whether to filter pairwise correlations before plotting. The default is set to `TRUE`.
#' @param pcorr.corr.dif.cutoff Numeric. A cutoff for difference in correlation coefficients between ordinary correlation
#' and partial correlation.  The default is set to 0.3. This means if a gene shows 0.3 difference between the two correlation
#' analyses, it will be selected.
#' @param pcorr.plot.ncol Numeric. Number of columns of a gird plot when putting more than one partial-pair wise plot.
#' The default is set to 3.
#' @param pcorr.plot.nrow Numeric. Number of rows of a gird plot when putting more than one partial-pair wise plot.
#' The default is set to 2.
#' @param deg.method Character. A character string that specifies which differential expression method should be used. The
#' default is set to `limma`.
#' @param deg.plot.ncol Numeric. Number of columns of a gird plot when putting more than one p-value histograms. The default
#' is set to 1.
#' @param deg.plot.nrow Numeric. Number of rows of a gird plot when putting more than one p-value histograms. The default
#' is set to 1.
#' @param gene.set.score.reference.data Character. A character string that specifies the gene set score whose data should
#' be used as reference, so all other scores from other datasets will be plotted against the reference data.
#' If set to `NULL`, all possible pair-wise plots of the scores will be generated.
#' @param gene.set.score.regress.out.variables Character. A character string or a vector specifying the name(s) of column(s)
#' in the sample annotation to regress out from the dataset(s) before performing gene set scoring analysis.
#' The default is set to `NULL`.
#' @param gene.set.score.list List. A list of gene sets used for enrichment analysis. All gene signature in the list will
#' be used to for scoring analysis separately.
#' @param gene.set.score.normalization Character. A character string that specifies which normalization should be applied
#' before applying the gene set scoring analysis.The default is set to `NULL`, so the dataset(s) will be used without any
#' normalization.
#' @param check.se.obj Logical. Whether to check the structure of the SummarizedExperiment object. The default is set to
#' `TRUE`
#' @param remove.na Character. Whether to remove `NA` values from assays. Options: `assays`, `none`. The default is set to
#' `none`.
#' @param override.check Logical. If set to `TRUE`, skip recalculating metrics if already present in the metadata of the
#' SummarizedExperiment object. The default is set to `FALSE`.
#' @param verbose Logical. Print progress messages while function is running. The default is set to `TRUE`.
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
        assay.names                          = 'all',
        bio.variables,
        uv.variables,
        assessment.level                     = 'L2',
        plots.to.exclude                     = NULL,
        apply.log                            = TRUE,
        pseudo.count                         = 1,
        general.points.size                  = 1.5,
        rle.iqr.width                        = 2,
        rle.median.points.size               = 1,
        rle.median.points.color              = 'red',
        rle.geom.hline.color                 = "cyan",
        rle.plot.ncol                        = 2,
        rle.plot.nrow                        = 3,
        rle.var.plot.ncol                    = 3,
        rle.var.plot.nrow                    = 3,
        rle.colors                           = NULL,
        fast.pca                             = TRUE,
        compute.nb.pcs                       = 10,
        nb.pcs.toplot.pca                    = 3,
        center                               = TRUE,
        scale                                = FALSE,
        svd.bsparam                          = bsparam(),
        pca.variables.colors                 = NULL,
        color.palette                        = 'nrc',
        pca.plot.nrow                        = 2,
        pca.plot.ncol                        = 3,
        pca.var.plot.ncol                    = 1,
        pca.var.plot.nrow                    = 3,
        pca.stroke.size                      = 0.05,
        pca.stroke.color                     = 'grey',
        pca.points.alpha                     = 0.5,
        pca.densities.alpha                  = 0.5,
        pca.legend.position                  = 'bottom',
        sil.dist.measure                     = 'euclidian',
        sli.nb.pcs                           = 3,
        ari.clustering.method                = "hclust",
        ari.hclust.method                    = "complete",
        ari.hclust.dist.measure              = "euclidian",
        ari.nb.pcs                           = 3,
        vca.nb.pcs                           = 3,
        lra.nb.pcs                           = 3,
        vca.nb.pcs.to.plot                   = 10,
        lra.nb.pcs.to.plot                   = 10,
        kbet.nb.pcs                          = 3,
        k0                                   = NULL,
        knn                                  = NULL,
        lisi.nb.pcs                          = 3,
        perplexity                           = 10,
        nn.eps                               = 0,
        corr.method                          = 'spearman',
        a                                    = 0.05,
        rho                                  = 0,
        correlation.plot.ncol                = 3,
        correlation.plot.nrow                = 3,
        anova.method                         = 'aov',
        anova.plot.ncol                      = 3,
        anova.plot.nrow                      = 3,
        pcorr.method                         = 'spearman',
        pcorr.genes                          = NULL,
        pcorr.select.genes                   = FALSE,
        pcorr.reference.data                 = NULL,
        pcorr.corr.cutoff                    = 0.6,
        pcorr.filter.genes                   = TRUE,
        pcorr.corr.dif.cutoff                = 0.1,
        pcorr.plot.ncol                      = 2,
        pcorr.plot.nrow                      = 2,
        deg.method                           = 'limma',
        deg.plot.ncol                        = 1,
        deg.plot.nrow                        = 1,
        gene.set.score.reference.data        = NULL,
        gene.set.score.regress.out.variables = NULL,
        gene.set.score.list,
        gene.set.score.normalization         = NULL,
        check.se.obj                         = TRUE,
        remove.na                            = 'none',
        override.check                       = FALSE,
        verbose                              = TRUE
        ){
    printColoredMessage(message = '------------The assessVariation function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the inputs of function ####
     if (!is.vector(assay.names)  | is.logical(assay.names) ){
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
    if (is.null(bio.variables) & is.null(uv.variables)){
        stop('Both "bio.variables" and "uv.variables" cannot be empty or NULL.')
    }
    if (is.logical(plots.to.exclude)){
        stop('The "to.exclude" must be a vector or NULL.')
    }
    if (!assessment.level %in% c('L1', 'L2')){
        stop('The "assessment.level" must be one of the "L1" or "L2".')
    }
    if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE".')
    }
    if (isTRUE(apply.log)){
         if (length(pseudo.count) > 1 | pseudo.count < 0 | is.null(pseudo.count))
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
            se.obj      = se.obj,
            assay.names = levels(assay.names),
            variables   = c(bio.variables, uv.variables),
            remove.na   = remove.na,
            verbose     = verbose)
    }
    # Getting all possible metrics for each variable #####
    printColoredMessage(
        message = '-- Finding all possible assessment metrics:',
        color = 'magenta',
        verbose = verbose
        )
    se.obj <- getAssessmentMetrics(
        se.obj = se.obj,
        uv.variables = uv.variables,
        bio.variables = bio.variables,
        plot.output = FALSE,
        save.se.obj = TRUE
        )
    metrics.table <- se.obj@metadata$AssessmentMetrics$metrics.table
    metrics.table$new.col <- paste0(
        metrics.table$Metrics,
        '_',
        metrics.table$Assessments
        )
    if (assessment.level == 'L1'){
        bio.metrics.table <- lapply(
            bio.variables,
            function(x){
                if (is.numeric(se.obj[[x]])){
                    bio.metrics.table <- metrics.table[metrics.table$Variables == x, ]
                    keep <- bio.metrics.table$new.col %in% c(
                        'LRA_averageRseq',
                        'Correlation_corrCoeff',
                        'PartialCorrelation_corrCoeff',
                        'LISI_meanScore'
                        )
                    bio.metrics.table <- bio.metrics.table[keep , ]
                } else {
                    bio.metrics.table <- metrics.table[metrics.table$Variables == x, ]
                    keep <- bio.metrics.table$new.col %in% c(
                        'VCA_averageCorr',
                        'ARI_ari',
                        'Silhouette_silhouetteCoeff',
                        'LISI_meanScore'
                        )
                    bio.metrics.table <- bio.metrics.table[keep , ]
                }
            })
        bio.metrics.table <- do.call(rbind, bio.metrics.table)
        uv.metrics.table <- lapply(
            uv.variables,
            function(x){
                if (is.numeric(se.obj[[x]])){
                    uv.metrics.table <- metrics.table[metrics.table$Variables == x, ]
                    keep <- uv.metrics.table$new.col %in% c(
                        'LRA_averageRseq',
                        'Correlation_corrCoeff',
                        'PartialCorrelation_corrCoeff',
                        'LISI_meanScore'
                    )
                    uv.metrics.table <- uv.metrics.table[keep , ]
                } else {
                    uv.metrics.table <- metrics.table[metrics.table$Variables == x, ]
                    keep <- uv.metrics.table$new.col %in% c(
                        'VCA_averageCorr',
                        'ARI_ari',
                        'Silhouette_silhouetteCoeff',
                        'DGE_pvalueNull',
                        'ANOVA_pvalueNull',
                        'ANOVA_fStat',
                        'LISI_meanScore',
                        'KBET_meanScore'
                    )
                    uv.metrics.table <- uv.metrics.table[keep , ]
                }
            })
        uv.metrics.table <- do.call(rbind, uv.metrics.table)
        metrics.table <- rbind(bio.metrics.table, uv.metrics.table, metrics.table[metrics.table$Variables == 'General' , ])
    }
    if (assessment.level == 'L2'){
        metrics.table <- metrics.table
    }
    # Excluding metrics and plots #####
     if (!is.null(plots.to.exclude)){
        printColoredMessage(
            message = paste0('- Excluding all the specified etrics and plots.'),
            color = 'blue',
            verbose = verbose
            )
        plots.to.exclude <- unique(plots.to.exclude)
        if (sum(plots.to.exclude %in% metrics.table$Code) != length(plots.to.exclude)){
            stop('All or some of the "plots.to.exclude" cannot be found in the assessment table.')
        }
        metrics.table <- metrics.table[!metrics.table$Code %in% plots.to.exclude, ]
    }
    printColoredMessage(
        message = paste0(
            '- In total, ',
            nrow(metrics.table),
            ' assessment plots will be generated.'
            ),
        color = 'blue',
        verbose = verbose
        )
    # RLE #####
    ## compute rle #####
     if ('RLE' %in% metrics.table$Metrics){
         if ('rlePlot' %in% metrics.table$PlotTypes){
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
     if ('rlePlot' %in% metrics.table$PlotTypes){
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
     if ('coloredRLEplot' %in% metrics.table$PlotTypes){
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
     if ('rleMedians' %in% metrics.table$Factors){
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
     if ('rleIqr' %in% metrics.table$Factors){
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
     if (sum(c('PCA', 'LRA', 'VCA', 'ARI', 'Silhouette')  %in% metrics.table$Metrics) > 0 ) {
        se.obj <- RUVprps::computePCA(
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
     if ('PCA' %in% metrics.table$Metrics & 'scatterPlot' %in% metrics.table$PlotTypes){
        pca.scatter.vars <- metrics.table$Metrics == 'PCA' & metrics.table$PlotTypes == 'scatterPlot'
        pca.scatter.vars <- metrics.table$Variables[pca.scatter.vars]
        for(i in pca.scatter.vars){
            if (is.numeric(se.obj[[i]])){
                plot.ncol = pca.var.plot.ncol
                plot.nrow = pca.var.plot.nrow
            } else {
                plot.ncol = pca.plot.ncol
                plot.nrow = pca.plot.nrow
            }
            se.obj <- RUVprps::plotPCA(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = nb.pcs.toplot.pca,
                plot.type = "scatter",
                variable.colors = pca.variables.colors,
                color.palette = color.palette,
                points.size = general.points.size,
                stroke.color = pca.stroke.color,
                stroke.size = pca.stroke.size,
                points.alpha = pca.points.alpha,
                densities.alpha = pca.densities.alpha,
                legend.position = pca.legend.position,
                plot.ncol = plot.ncol,
                plot.nrow = plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = TRUE)
        }
     }
    ## boxplot pca ####
     if ('PCA' %in% metrics.table$Metrics & 'boxPlot' %in% metrics.table$PlotTypes){
        pca.boxplot.vars <- metrics.table$Metrics == 'PCA' & metrics.table$PlotTypes == 'boxPlot'
        pca.boxplot.vars <- metrics.table$Variables[pca.boxplot.vars]
        for(i in pca.boxplot.vars){
            se.obj <- RUVprps::plotPCA(
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
                plot.ncol = pca.var.plot.ncol,
                plot.nrow = pca.var.plot.nrow,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = TRUE)
        }
     }
    # Vector correlation ####
    ## compute vector correlation ####
    if (vca.nb.pcs.to.plot > vca.nb.pcs){
        vca.nb.pcs.initial <- vca.nb.pcs
        vca.nb.pcs <- vca.nb.pcs.to.plot
    }
    if ('VCA' %in% metrics.table$Metrics){
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
                verbose = verbose
                )
        }
    }
    ## plot vector correlation ####
    if ('VCA' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'VCA'
        pc.vec.corr.vars <- metrics.table$Variables[index]
        for(i in pc.vec.corr.vars){
            se.obj <- plotPCVariableCorrelation(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = vca.nb.pcs.to.plot,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose
                )
        }
    }
    if (vca.nb.pcs.to.plot > vca.nb.pcs.initial){
        if ('VCA' %in% metrics.table$Metrics){
            pc.vec.corr.vars <- metrics.table$Metrics == 'VCA'
            pc.vec.corr.vars <- metrics.table$Variables[pc.vec.corr.vars]
            for(i in pc.vec.corr.vars){
                se.obj <- computePCVariableCorrelation(
                    se.obj = se.obj,
                    assay.names = levels(assay.names),
                    variable = i,
                    fast.pca = fast.pca,
                    nb.pcs = vca.nb.pcs.initial,
                    save.se.obj = TRUE,
                    verbose = verbose
                )
            }
        }
    }
    # Linear regression ####
    ## compute linear regression ####
    if (lra.nb.pcs.to.plot > lra.nb.pcs){
        lra.nb.pcs.initial <- lra.nb.pcs
        lra.nb.pcs <- lra.nb.pcs.to.plot
    }
    if ('LRA' %in% metrics.table$Metrics){
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
     if ('LRA' %in% metrics.table$Metrics){
        pc.reg.vars <- metrics.table$Metrics == 'LRA'
        pc.reg.vars <- metrics.table$Variables[pc.reg.vars]
        for(i in pc.reg.vars){
            se.obj <- plotPCVariableRegression(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                nb.pcs = lra.nb.pcs.to.plot,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
     }
    if (lra.nb.pcs.to.plot > lra.nb.pcs.initial){
        if ('LRA' %in% metrics.table$Metrics){
            pc.reg.vars <- metrics.table$Metrics == 'LRA'
            pc.reg.vars <- metrics.table$Variables[pc.reg.vars]
            for(i in pc.reg.vars){
                se.obj <- computePCVariableRegression(
                    se.obj = se.obj,
                    assay.names = levels(assay.names),
                    variable = i,
                    fast.pca = fast.pca,
                    nb.pcs = lra.nb.pcs.initial,
                    save.se.obj = TRUE,
                    verbose = verbose)
            }
        }
    }

    # Silhouette coefficient ####
    ## compute silhouette coefficients ####
     if ('Silhouette' %in% metrics.table$Metrics){
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
     if ('Silhouette' %in% metrics.table$Metrics & 'barPlot' %in% metrics.table$PlotTypes){
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
     if ('Silhouette' %in% metrics.table$Metrics & 'combinedPlot' %in% metrics.table$PlotTypes){
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
     if ('ARI' %in% metrics.table$Metrics){
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
     if ('ARI' %in% metrics.table$Metrics & 'barPlot' %in% metrics.table$PlotTypes){
        index <- metrics.table$Metrics == 'ARI' & metrics.table$PlotTypes == 'barPlot'
        ari.single.vars <- metrics.table$Variables[index]
         if (ari.clustering.method == 'mclust'){
            ari.method <- 'mclust'
        } else ari.method <- paste0('hclust.', ari.hclust.method, '.', ari.hclust.dist.measure)
        for(i in ari.single.vars){
            se.obj <- plotARI(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = i,
                ari.method = ari.method,
                plot.type = 'single.plot',
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
     }
    ## combined adjusted rand index  ####
     if ('ARI' %in% metrics.table$Metrics & 'combinedPlot' %in% metrics.table$PlotTypes){
        index <- metrics.table$Metrics == 'ARI' & metrics.table$PlotTypes == 'combinedPlot'
        ari.combined.vars <- metrics.table$Variables[index]
         if (ari.clustering.method == 'mclust'){
            ari.method <- 'mclust'
        } else ari.method <- paste0('hclust.', ari.hclust.method, '.', ari.hclust.dist.measure)
        for(i in ari.combined.vars){
            se.obj <- plotARI(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variables = strsplit(x = i, split = '&')[[1]],
                ari.method = ari.method,
                plot.type = 'combined.plot',
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose)
        }
     }

    # LISI ####
    ## compute local inverse Simpson's index  ####
    if ('LISI' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'LISI'
        lisi.vars <- unique(metrics.table$Variables[index])
        for (i in lisi.vars){
            se.obj <- computeLisi(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                perplexity = perplexity,
                nn.eps = nn.eps,
                fast.pca = fast.pca,
                nb.pcs = lisi.nb.pcs,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    ## plot local inverse Simpson's index  ####
    if ('LISI' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'LISI'
        lisi.vars <- unique(metrics.table$Variables[index])
        for (i in lisi.vars){
            se.obj <- plotLisi(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose
                )
        }
    }

    # kBET ####
    ## k-nearest neighbor batch effect test ####
    if ('KBET' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'KBET'
        kebt.vars <- unique(metrics.table$Variables[index])
        for (i in kebt.vars){
            se.obj <- computeKbet(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                k0 = k0,
                knn = knn,
                fast.pca = fast.pca,
                nb.pcs = kbet.nb.pcs,
                save.se.obj = TRUE,
                verbose = verbose
                )
        }
    }
    ## plot k-nearest neighbor batch effect test ####
    if ('KBET' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'KBET'
        kebt.vars <- unique(metrics.table$Variables[index])
        for (i in kebt.vars){
            se.obj <- plotKbet(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                fast.pca = fast.pca,
                plot.output = FALSE,
                save.se.obj = TRUE,
                verbose = verbose
                )
        }
    }


    # Gene variable correlation ####
    ## compute gene variable correlations ####
     if ('Correlation' %in% metrics.table$Metrics){
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
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                plot.top.genes = FALSE,
                nb.top.genes = NULL,
                apply.round = TRUE,
                check.se.obj = FALSE,
                remove.na = 'none',
                override.check = override.check,
                save.se.obj = TRUE,
                verbose = verbose)
        }
    }
    ## plot gene variable correlations ####
     if ('Correlation' %in% metrics.table$Metrics){
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
     if ('ANOVA' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'ANOVA'
        gene.var.anova.vars <- unique(metrics.table$Variables[index])
        for(i in gene.var.anova.vars){
            se.obj <- computeGenesVariableAnova(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                variable = i,
                method = anova.method,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                plot.top.genes = FALSE,
                nb.top.genes = NULL,
                apply.round = TRUE,
                check.se.obj = FALSE,
                remove.na = 'none',
                override.check = override.check,
                save.se.obj = TRUE)
        }
    }
    ## plot gene variable anova ####
     if ('ANOVA' %in% metrics.table$Metrics){
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
     if ('DGE' %in% metrics.table$Metrics){
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
    if ('DGE' %in% metrics.table$Metrics){
        deg.plot.ncol.initial <- deg.plot.ncol
        deg.plot.nrow.initial <- deg.plot.nrow
        index <- metrics.table$Metrics == 'DGE'
        dge.vars <- unique(metrics.table$Variables[index])
        for(i in dge.vars){
            possible.groups <- combn(x = unique(se.obj[[i]]), m = 2)
            if (is.null(deg.plot.ncol) & is.null(deg.plot.nrow)){
                if (ncol(possible.groups) > 10){
                    deg.plot.ncol = 1
                    deg.plot.nrow = 1
                }
                if (ncol(possible.groups) > 5 && ncol(possible.groups) < 10){
                    deg.plot.ncol = 2
                    deg.plot.nrow = 2
                }
                if (ncol(possible.groups) < 5){
                    deg.plot.ncol = 3
                    deg.plot.nrow = 3
                }
            }
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
            deg.plot.ncol = deg.plot.ncol.initial
            deg.plot.nrow = deg.plot.nrow.initial
        }
    }
    # Partial correlation ####
    ## compute partial correlation ####
     if ('PartialCorrelation' %in% metrics.table$Metrics){
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
     if ('PartialCorrelation' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'PartialCorrelation'
        pcorr.vars <- unique(metrics.table$Variables[index])
        for(i in pcorr.vars){
            index <- metrics.table$Metrics == 'PartialCorrelation' & metrics.table$Variables == i
            plot.types <- metrics.table$PlotTypes[index]
                for(j in plot.types){
                     if (j == 'scatterPlot'){
                        plot.type = 'scatter.plot'
                    } else  if (j == 'barPlot'){
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
     if ('GeneSetScore' %in% metrics.table$Metrics){
        index <- metrics.table$Metrics == 'GeneSetScore'
        gene.set.vars <- unique(metrics.table$Variables[index])
        gene.set.vars <- gene.set.vars[gene.set.vars %in% names(gene.set.score.list)]
        if (length(gene.set.vars) > 0){
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
                    gene.set.name = i,
                    check.se.obj = FALSE,
                    save.se.obj = TRUE,
                    verbose = verbose
                )
            }
        }
    }

    ## plot gene set scores ####
     if ('GeneSetScore' %in% metrics.table$Metrics){
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

