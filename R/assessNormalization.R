#' Assess the performance of RNA-seq normalization methods.
#'
#' @author Ramyar Molania
#'
#' @references
#' **Molania R., ..., Speed, T. P., A new normalization for Nanostring nCounter gene expression data, Nucleic Acids
#' Research, 2019.
#' **Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @description
#' This function summarizes a range of global and gene level metrics obtained using the `assessVariation` function for
#' individual biological and unwanted variables. The functions returns numerical assessments in a table to assess the
#' performance of difference normalization. Refer to details for more information.
#'
#' @details
#' It is essential to assess the performance of a normalization as a remover of unwanted variation and a preserver of
#' biological variation in the data. An ideal normalization should preserve all known sources of biological variation
#' and reduce or remove the impact of all known sources of unwanted variation. To rank methods according to their
#' performance, we aggregate all evaluation metrics into a overall score as follows.
#' Several assessment will be performed:
#' For each categorical variable:
#' - PCA plot of the categorical variable.
#' - Silhouette and ARI computed on the categorical variable.
#' - Differential analysis based ANOVA between the gene expression and the categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and the categorical variable.
#' For each continous variable:
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Correlation between gene expression and continuous variable.
#'
#' It will output the following plots:
#' - PCA plot of each categorical variable.
#' - Boxplot of the F-test distribution from ANOVA between the gene expression and each categorical variable.
#' - Vector correlation between the first cumulative PCs of the gene expression and each categorical variable.
#' - Combined Silhouette plot of the combined pair of all categorical variables.
#' - Linear regression between the first cumulative PC and continuous variable.
#' - Boxplot of the correlation between gene expression and continuous variable.
#' - It will also output the RLE plot distribution.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names character or character vector. One or more names of assays to select from the SummarizedExperiment
#' object.
#' The default is 'all', which selects all assays in the object.
#' @param bio.variables character or character vector. One or more column names indicating known biological variables
#' (categorical or continuous) in the SummarizedExperiment object.
#' @param uv.variables character or character vector. One or more column names indicating unwanted variables
#' (categorical or continuous) in the SummarizedExperiment object.
#' @param assessment.level TTT
#' @param assessments.to.exclude character or character vector. Names of assessment metrics to exclude, as returned by
#' the `getAssessmentMetrics()` function. Default is NULL.
#' @param fast.pca Logical. Whether to use fast PCA. The default is to `TRUE`.
#' @param select.top.ruv TTTT
#' @param sil.dist.measure character. Distance measure to use for silhouette analysis. Options: 'euclidean', 'maximum',
#' 'manhattan', 'canberra', 'binary', or 'minkowski'. Default is 'euclidean'.
#' @param ari.clustering.method character. Clustering method for ARI computation. Options: 'mclust' or 'hclust'. Default
#' is 'hclust'.
#' @param ari.hclust.method character. Method used in `hclust` for ARI. Options: 'ward.D', 'ward.D2', 'single', 'complete',
#' 'average', 'mcquitty', 'median', or 'centroid'. Default is 'complete'.
#' @param ari.hclust.dist.measure character. Distance measure for `hclust` in ARI analysis. Options:'euclidean', 'maximum',
#' 'manhattan', 'canberra', 'binary', or 'minkowski'. Default is 'euclidean'.
#' @param corr.method character. Correlation method for gene-variable correlation. Options: 'spearman' or 'pearson'.
#' Default is 'spearman'.
#' @param corr.cutoff Numeric. Threshold for selecting genes with an absolute correlation lower than this value with
#' continuous variables. The default is set to 0.2.
#' @param anova.method character. ANOVA method for gene-variable comparison. Options: 'aov' or 'welch'. Default is 'aov'.
#' @param fvalue.cutoff Numeric. Threshold for selecting genes with an F-statistic lower than this value with categorical
#' variables. Te default is set to 1.
#' @param pcorr.method character. Correlation method for gene-gene ordinary and partial correlations. Options: 'spearman'
#' or 'pearson'. The default is set to `spearman`.
#' @param pcorr.cutoff Numeric. Threshold for selecting gene-gene correlations where the difference between ordinary and
#' partial correlations is less than this value.
#' @param bio.weight Numeric. Weight for the biological preservation score. Default is 0.6. See details for more information.
#' @param uv.weight Numeric. Weight for the removal of unwanted variation score. Default is 0.4. See details for more
#' information.
#' @param plot.output Logical. If TRUE, displays the final assessment plot.
#' @param save.se.obj Logical. If TRUE, saves the results in the SummarizedExperiment object. Default is TRUE.
#' @param output.name character. File name for the results stored in the SummarizedExperiment. If NULL, the function creates
#' a default name using `paste0('RUVprps_assessNormalization', length(assays(se.obj)), '_assays.')`.
#' @param verbose Logical. If TRUE, displays messages during function execution.
#' @param sli.nb.pcs Numeric. Number of principal components to use for silhouette analysis. Default is 3. Must not exceed
#' `compute.nb.pcs`.
#' @param ari.nb.pcs Numeric. Number of principal components to use for ARI analysis. Default is 3. Must not exceed
#' `compute.nb.pcs`.
#' @param vca.nb.pcs Numeric. Number of principal components to use for variance component analysis. Default is 3. Must
#' not exceed `compute.nb.pcs`.
#' @param lra.nb.pcs Numeric. Number of principal components to use for latent representation analysis. Default is 3.
#' Must not exceed `compute.nb.pcs`.
#' @return A SummarizedExperiment object containing the assessment matrix, plot, and table, or a list containing all results.
#'
#' @importFrom ggh4x strip_nested elem_list_text elem_list_rect facet_nested
#' @importFrom dplyr summarise_at case_match row_number
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom SummarizedExperiment assays colData
#' @importFrom stats ks.test IQR fligner.test coef
#' @importFrom gridExtra grid.arrange grid.table
#' @importFrom ggforestplot geom_stripes
#' @importFrom graphics plot.new text
#' @importFrom stats kruskal.test
#' @importFrom qvalue qvalue
#' @import RColorBrewer
#' @export

assessNormalization <- function(
        se.obj,
        assay.names = 'all',
        bio.variables,
        uv.variables,
        assessment.level = 'L1',
        assessments.to.exclude = NULL,
        select.top.ruv = TRUE,
        fast.pca = TRUE,
        sil.dist.measure = 'euclidian',
        sli.nb.pcs = 3,
        ari.clustering.method = 'hclust',
        ari.hclust.method = 'complete',
        ari.hclust.dist.measure = 'euclidian',
        ari.nb.pcs = 3,
        corr.method = 'spearman',
        corr.cutoff = 0.2,
        anova.method = 'aov',
        fvalue.cutoff = 1,
        vca.nb.pcs = 3,
        lra.nb.pcs = 3,
        pcorr.method = 'spearman',
        pcorr.cutoff = 0.3,
        bio.weight = 0.6,
        uv.weight = 0.4,
        plot.output = TRUE,
        save.se.obj = TRUE,
        output.name = NULL,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The assessNormalization function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the inputs of function ####
    if (length(assay.names) == 1 && assay.names != 'all') {
        if (!assay.names %in% names(assays(se.obj)))
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    if (length(assay.names) > 1) {
        if (length(setdiff(assay.names, names(assays(se.obj)))) > 0)
            stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    if (is.null(bio.variables) | is.null(uv.variables)){
        stop('To performe "overall.performance" both "uv.variables" and "bio.variables" must be provided.')
    }
    if (!is.null(corr.cutoff)){
        if (!is.list(corr.cutoff)){
            stop('The "corr.cutoff" must be a list of correlation cutoff for each continuous variables.')
        }
    }

    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else assay.names <- factor(x = assay.names, levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Get all possible metrics for each variable #####
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

    printColoredMessage(
        message = paste0(
            '- The totall of ',
            sum(metrics.table$Assessments!='DA'),
            ' assessment metrics will be used for normalization performance assessment.'),
        color = 'blue',
        verbose = verbose
        )

    # Filter assessment #####
    if (!is.null(assessments.to.exclude)){
        printColoredMessage(
            message = '-- Filtering assessment metrics:',
            color = 'magenta',
            verbose = verbose
        )
        metrics.table <- metrics.table[!metrics.table$Code %in% assessments.to.exclude, ]
        printColoredMessage(
            message = paste0(
                '- After filteration, the totall of ',
                sum(metrics.table$Assessments!='DA'),
                ' assessment metrics will be used for normalization performance assessment.'),
            color = 'blue',
            verbose = verbose
        )

    } else metrics.table <- metrics.table

    # Summarizing all the selected assessments ####
    printColoredMessage(
        message = '-- Summarizing all the selected assessment metrics:',
        color = 'magenta',
        verbose = verbose
        )
    cols.names <- c('data', 'variable', 'test', 'measurements')
    ## RLE general scores #####
    ### rle medians ####
    if ('General_rleMed' %in% paste(metrics.table$Variables, metrics.table$Assessments, sep = '_')){
        printColoredMessage(
            message = '- Summarize the RLE medians of the general RLE plots:',
            color = 'orange',
            verbose = verbose
            )
        #### check the RLE medians ####
        printColoredMessage(
            message = '* check to see all the RLE medians are computed.',
            color = 'blue',
            verbose = verbose
            )
        check.out <- lapply(
            levels(assay.names),
            function(x){
                if (is.null(se.obj@metadata$Metrics[[x]]$global.level$RLE$gene.median.center$general$data$rle.med)){
                    stop(paste0('The RLE medians for the ',
                                x ,
                                ' data cannot be found in the SummarizedExperiment object. ',
                                'Please run either the "compuetRLE" or "assessVariation" functions, and  ',
                                'make sure to compute RLE medians.'))
                }
            })
        ### compute rle medians scores ####
        printColoredMessage(
            message = '- Computing the RLE medians scores for each data.',
            color = 'blue',
            verbose = verbose
            )
        general.rle.med.scores <- sapply(
            levels(assay.names),
            function(x)
                abs(se.obj@metadata$Metrics[[x]]$global.level$RLE$gene.median.center$general$data$rle.med)
            )
        max.rle.meds <- max(general.rle.med.scores) * ncol(se.obj)
        general.rle.med.scores <- sapply(
            levels(assay.names),
            function(x){
                med <- general.rle.med.scores[ , x]
                1 - sum(med)/max.rle.meds
            })
        names(general.rle.med.scores) <- levels(assay.names)

        ### put all together ####
        general.rle.med.scores <- as.data.frame(x = general.rle.med.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            dplyr::mutate(variable = 'RLE') %>%
            dplyr::mutate(test = 'RLE medians') %>%
            data.frame(.)
        row.names(general.rle.med.scores) <- c(1:nrow(general.rle.med.scores))
        colnames(general.rle.med.scores)[1] <- 'measurements'
        general.rle.med.scores <- general.rle.med.scores[ , cols.names]
    } else general.rle.med.scores <- NULL

    ### rle iqr ####
    if ('General_rleIqr' %in% paste(metrics.table$Variables, metrics.table$Assessments, sep = '_')){
        printColoredMessage(
            message = '- Summarize the RLE IQR of the general RLE plots:',
            color = 'orange',
            verbose = verbose
        )
        #### check the RLE iqr  ####
        printColoredMessage(
            message = '* check to see all the RLE IQRs are computed.',
            color = 'blue',
            verbose = verbose
        )
        check.out <- lapply(
            levels(assay.names),
            function(x){
                if (is.null(se.obj@metadata$Metrics[[x]]$global.level$RLE$gene.median.center$general$data$rle.iqr))
                    stop(paste0('The RLE IQRs for the ', x, ' data cannot be found in the SummarizedExperiment object.',
                                'Please run either the "compuetRLE" or "assessVariation" functions, and ',
                                'make sure to compute RLE IQR.'))
            })

        ### compute rle iqr scores ####
        printColoredMessage(
            message = '* compute the RLE IQR scores.',
            color = 'blue',
            verbose = verbose
            )
        general.rle.iqr.scores <- sapply(
            levels(assay.names),
            function(x){
                iqr <- se.obj@metadata$Metrics[[x]]$global.level$RLE$gene.median.center$general$data$rle.iqr
                iqr.qun_0.1_0.5_0.9 <- quantile(x = iqr, probs = c(0.1, 0.5, 0.9))
                1 - c(iqr.qun_0.1_0.5_0.9[[3]] - iqr.qun_0.1_0.5_0.9[[1]])/iqr.qun_0.1_0.5_0.9[[2]]
            })
        names(general.rle.iqr.scores) <- levels(assay.names)
        ### put all together ####
        general.rle.iqr.scores <- as.data.frame(x = general.rle.iqr.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            dplyr::mutate(variable = 'RLE') %>%
            dplyr::mutate(test = 'RLE IQRs') %>%
            data.frame(.)
        row.names(general.rle.iqr.scores) <- c(1:nrow(general.rle.iqr.scores))
        colnames(general.rle.iqr.scores)[1] <- 'measurements'
        general.rle.iqr.scores <- general.rle.iqr.scores[ , cols.names]
    } else general.rle.iqr.scores <- NULL

    ## RLE and variable association scores ####
    ### correlation between rle medians and variable #####
    if ('rleMedians_scatterPlot' %in% paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_')){
        printColoredMessage(
            message = '- Summarizing the association between the rle medians and the continuous variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check the RLE medians ####
        printColoredMessage(
            message = '- Checking to see all the RLE medians are computed.',
            color = 'blue',
            verbose = verbose
            )
        check.out <- lapply(
            levels(assay.names),
            function(x){
                if (is.null(se.obj@metadata$Metrics[[x]]$global.level$RLE$gene.median.center$general$data$rle.med))
                    stop(paste0('The RLE medians for the ',
                                x ,
                                ' data cannot be found in the SummarizedExperiment object.',
                                'Please run either the "compuetRLE" or "assessVariation" functions and ',
                                'make sure to compute RLE medians.'))
            })
        ### compute correlation between rle medians and variable ####
        printColoredMessage(
            message = '- Computing the association scores.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_') == 'rleMedians_scatterPlot'
        selected.vars <- metrics.table$Variables[selected.vars]
        # rle.med.var.corr.scores <- lapply(
        #     selected.vars,
        #     function(x){
        #         rle.med.var.corr.scores <- sapply(
        #             levels(assay.names),
        #             function(y){
        #                 abs(suppressWarnings(stats::cor.test(
        #                     x = se.obj@metadata$Metrics[[y]]$global.level$RLE$gene.median.center$general$data$rle.med,
        #                     y = colData(se.obj)[[x]], method = corr.method)[[4]])[[1]])
        #             })
        #     })
        # names(rle.med.var.corr.scores) <- selected.vars
        rle.med.var.corr.scores <- lapply(
            selected.vars,
            function(x){
                rle.med.var.corr.scores <- sapply(
                    levels(assay.names),
                    function(y){
                            a <- se.obj@metadata$Metrics[[y]]$global.level$RLE$gene.median.center$general$data$rle.med
                            b <- colData(se.obj)[[x]]
                            slope <- stats::coef(lm(b ~ a))[2]
                            scaled <- 1 / (1 + exp(-slope))
                            unname(scaled)
                    })
            })
        names(rle.med.var.corr.scores) <- selected.vars

        ### put all together ####
        rle.med.var.corr.scores <- as.data.frame(x = rle.med.var.corr.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Association with RLE medians') %>%
            data.frame(.)
        rle.med.var.corr.scores <- rle.med.var.corr.scores[ , cols.names]
        } else rle.med.var.corr.scores <- NULL

    ### rle iqr and variable correlation #####
    if ('rleIqr_scatterPlot' %in% paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_')){
        printColoredMessage(
            message = '- Summarizing the association between the rle IQR and the continuous variable(s):',
            color = 'orange',
            verbose = verbose
        )
        #### check the RLE medians ####
        printColoredMessage(
            message = '- Checking to see all the RLE IQRs are computed.',
            color = 'blue',
            verbose = verbose
        )
        check.out <- lapply(
            levels(assay.names),
            function(x){
                if (is.null(se.obj@metadata$Metrics[[x]]$global.level$RLE$gene.median.center$general$data$rle.iqr)){
                    stop(paste0('The RLE IQRs for the ',
                                x,
                                ' data cannot be found in the SummarizedExperiment object.',
                                'Please run the "compuetRLE" or "assessVariation" functions and ',
                                'make sure to compute the RLE IQRs.'))
                }
            })

        ### compute correlation between rle medians and variable ####
        printColoredMessage(
            message = '- Computing the association scores.',
            color = 'blue',
            verbose = verbose
        )
        selected.vars <- paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_') == 'rleIqr_scatterPlot'
        selected.vars <- metrics.table$Variables[selected.vars]
        # rle.iqr.var.corr.scores <- lapply(
        #     selected.vars,
        #     function(x){
        #         rle.med.var.corr.scores <- sapply(
        #             levels(assay.names),
        #             function(y){
        #                 abs(suppressWarnings(stats::cor.test(
        #                     x = se.obj@metadata$Metrics[[y]]$global.level$RLE$gene.median.center$general$data$rle.iqr,
        #                     y = colData(se.obj)[[x]], method = corr.method)[[4]][[1]]))
        #             })
        #     })
        # names(rle.iqr.var.corr.scores) <- selected.vars
        rle.iqr.var.corr.scores <- lapply(
            selected.vars,
            function(x){
                rle.med.var.corr.scores <- sapply(
                    levels(assay.names),
                    function(y){
                            a <- se.obj@metadata$Metrics[[y]]$global.level$RLE$gene.median.center$general$data$rle.iqr
                            b <- colData(se.obj)[[x]]
                            slope <- stats::coef(lm(b ~ a))[2]
                            scaled <- 1 / (1 + exp(-slope))
                            unname(scaled)
                    })
            })
        names(rle.iqr.var.corr.scores) <- selected.vars

        ### put all together ####
        rle.iqr.var.corr.scores <- as.data.frame(x = rle.iqr.var.corr.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Association with RLE IQRs') %>%
            data.frame(.)
        rle.iqr.var.corr.scores <- rle.iqr.var.corr.scores[ , cols.names]
    } else rle.iqr.var.corr.scores <- NULL

    ### rle medians and variable association #####
    if ('rleMedians_boxPlot' %in% paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_')){
        printColoredMessage(
            message = '- Summarizing the association between the rle medians and the categorical variables:',
            color = 'orange',
            verbose = verbose
            )
        #### check the RLE medians ####
        printColoredMessage(
            message = '- Checking to see all the RLE medians are computed.',
            color = 'blue',
            verbose = verbose
            )
        check.out <- lapply(
            levels(assay.names),
            function(x){
                if (is.null(se.obj@metadata$Metrics[[x]]$global.level$RLE$gene.median.center$general$data$rle.med)){
                    stop(paste0('The RLE mediasns for the ', x , ' data cannot be found in the SummarizedExperiment object.',
                                'Please run the "compuetRLE" or "assessVariation" functions and ',
                                'make sure to compute the RLE medians.'))
                }
            })

        ### compute association between rle medians and variable ####
        printColoredMessage(
            message = '- Computing the association scores.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_') == 'rleMedians_boxPlot'
        selected.vars <- metrics.table$Variables[selected.vars]
        rle.med.var.association.scores <- lapply(
            selected.vars,
            function(x){
                rle.med.var.corr.scores <- sapply(
                    levels(assay.names),
                    function(y){
                        suppressWarnings(
                            p.value <- kruskal.test(
                            x = se.obj@metadata$Metrics[[y]]$global.level$RLE$gene.median.center$general$data$rle.med,
                            g = colData(se.obj)[[x]])$p.value
                            )
                        if (p.value < 0.01) p.value = 0.01
                        p.value
                    })
            })
        names(rle.med.var.association.scores) <- selected.vars

        ### put all together ####
        rle.med.var.association.scores <- as.data.frame(x = rle.med.var.association.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Association with RLE medians') %>%
            data.frame(.)
        rle.med.var.association.scores <- rle.med.var.association.scores[ , cols.names]
    } else rle.med.var.association.scores <- NULL

    ### rle iqr and variable association #####
    if ('rleIqr_boxPlot' %in% paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_')){
        printColoredMessage(
            message = '- Summarizing the association between the rle IQRs and variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check the RLE medians ####
        printColoredMessage(
            message = '- Checking to see all the RLE IQRs are computed.',
            color = 'blue',
            verbose = verbose
            )
        check.out <- lapply(
            levels(assay.names),
            function(x){
                if (is.null(se.obj@metadata$Metrics[[x]]$global.level$RLE$gene.median.center$general$data$rle.iqr)){
                    stop(paste0('The RLE IQRs for the ',
                                x,
                                ' data cannot be found in the SummarizedExperiment object.',
                                'Please run the "compuetRLE" or "assessVariation" functions.'))
                }
            })
        ### compute association between rle IQRs and variable ####
        printColoredMessage(
            message = '* compute the association scores.',
            color = 'blue',
            verbose = verbose
        )
        selected.vars <- paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_') == 'rleIqr_boxPlot'
        selected.vars <- metrics.table$Variables[selected.vars]
        rle.iqr.var.association.scores <- lapply(
            selected.vars,
            function(x){
                rle.med.var.corr.scores <- sapply(
                    levels(assay.names),
                    function(y){
                        suppressWarnings(
                            p.value <- kruskal.test(
                                x = se.obj@metadata$Metrics[[y]]$global.level$RLE$gene.median.center$general$data$rle.iqr,
                                g = colData(se.obj)[[x]])$p.value)
                        if (p.value < 0.01) p.value = 0.01
                        p.value
                    })
            })
        names(rle.iqr.var.association.scores) <- selected.vars
        ### put all together ####
        rle.iqr.var.association.scores <- as.data.frame(x = rle.iqr.var.association.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Association with RLE IQRs') %>%
            data.frame(.)
        rle.iqr.var.association.scores <- rle.iqr.var.association.scores[ , cols.names]
    } else rle.iqr.var.association.scores <- NULL

    ## Vector correlation scores ####
    if ('VCA' %in% metrics.table$Metrics){
        printColoredMessage(
            message = '- Summarizing the vector correlations for each categorical variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check the vector correlations ####
        printColoredMessage(
            message = '- Checking to see all the vector correlations are computed.',
            color = 'blue',
            verbose = verbose
        )
        selected.vars <- metrics.table$Variables[metrics.table$Metrics == 'VCA']
        if (isTRUE(fast.pca)){
            svd.method <- 'fast.svd'
        } else svd.method <- 'ordinary.svd'
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$global.level$VCA[[svd.method]][[x]]$vector.correlations)){
                            stop(paste0('The vector correlation for the ',
                                        x, ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computePCVariableCorrelation" or "assessVariation" functions.'))
                        }
                    })
            })

        ## average vector correlations ####
        printColoredMessage(
            message = '- Computing the vector correlation scorss.',
            color = 'blue',
            verbose = verbose
            )
        pc.vec.corr.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        mean(
                            se.obj@metadata$Metrics[[y]]$global.level$VCA[[svd.method]][[x]]$vector.correlations[1:vca.nb.pcs]
                        )
                    })
            })
        names(pc.vec.corr.scores) <- selected.vars

        ### pcs vector correlation ####
        pc.vec.corr.scores <- as.data.frame(x = pc.vec.corr.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'VCA') %>%
            data.frame(.)
        pc.vec.corr.scores <- pc.vec.corr.scores[ , cols.names]
    } else pc.vec.corr.scores <- NULL

    ## Linear regression scores ####
    if ('LRA' %in% metrics.table$Metrics){
        printColoredMessage(
            message = '- Summarizing the R.squared form the linear regression analysis for each continous variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check the linear regression ####
        printColoredMessage(
            message = '- Checking to see all the R.squared of the linear regression analysis are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- metrics.table$Variables[metrics.table$Metrics == 'LRA']
        if (isTRUE(fast.pca)){
            svd.method <- 'fast.svd'
        } else svd.method <- 'ordinary.svd'
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$global.level$LRA[[svd.method]][[x]]$r.squared)){
                            stop(paste0('The vector correlation for the ',
                                        x,
                                        ' variable of the ' , y ,
                                        ' data cannot be found
                                        in the SummarizedExperiment object.',
                                        'Please run the "computePCVariableRegression" or "assessVariation" functions.'))
                        }
                    })
            })
        ## average r squared ####
        printColoredMessage(
            message = '- Computing the linear regression scores.',
            color = 'blue',
            verbose = verbose
            )
        pc.lin.reg.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y)
                        mean(se.obj@metadata$Metrics[[y]]$global.level$LRA[[svd.method]][[x]]$r.squared[1:lra.nb.pcs])
                    )
            })
        names(pc.lin.reg.scores) <- selected.vars

        ### Put all together ####
        pc.lin.reg.scores <- as.data.frame(x = pc.lin.reg.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'LRA') %>%
            data.frame(.)
        pc.lin.reg.scores <- pc.lin.reg.scores[ , cols.names]
    } else pc.lin.reg.scores <- NULL

    ## Silhouette scores ####
    if ('silhouetteCoeff_barPlot' %in% paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_')){
        printColoredMessage(
            message = '- Summarizing the silhouette coefficients analysis for each categorical variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check the linear regression ####
        printColoredMessage(
            message = '- Checking to see all the silhouette coefficients are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_') == 'silhouetteCoeff_barPlot'
        selected.vars <- metrics.table$Variables[selected.vars]
        silhouette.method <- paste0('sil.', sil.dist.measure)
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$global.level$Silhouette[[silhouette.method]][[x]]$silhouette.coeff)){
                            stop(paste0('The Silhouette for the ',
                                        x,
                                        ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeSilhouette" or "assessVariation" functions.'))
                        }
                    })
            })
        #### compute the silhouette scores ####
        printColoredMessage(
            message = '- Computing the silhouette scores.',
            color = 'blue',
            verbose = verbose
        )
        sil.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y)
                        c(se.obj@metadata$Metrics[[y]]$global.level$Silhouette$sil.euclidian[[x]]$silhouette.coeff + 1)/2
                    )
            })
        names(sil.scores) <- selected.vars

        ### Put all together ####
        sil.scores <- as.data.frame(x = sil.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Silhouette') %>%
            data.frame(.)
        sil.scores <- sil.scores[ , cols.names]
    } else sil.scores <- NULL

    ## ARI scores ####
    if ('ariCoeff_barPlot' %in% paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_')){
        printColoredMessage(
            message = '- Summarizing the adjusted rand index for each categorical variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check the linear regression ####
        printColoredMessage(
            message = '- Checking to see all the adjusted rand index are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Factors, metrics.table$PlotTypes, sep = '_') == 'ariCoeff_barPlot'
        selected.vars <- metrics.table$Variables[selected.vars]
        if (ari.clustering.method == 'mclust'){
            ari.method <- 'mclust'
        } else ari.method <- paste0('hclust.', ari.hclust.method, '.', ari.hclust.dist.measure)
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$global.level$ARI[[ari.method]][[x]]$ari)){
                            stop(paste0('The adjusted rand index for the ',
                                        x,
                                        ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeARI" or "assessVariation" functions.'))
                        }
                    })
            })
        #### obtain silhouette ####
        printColoredMessage(
            message = '- Computing the adjusted rand index scores.',
            color = 'blue',
            verbose = verbose
            )
        ari.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y)
                        c(se.obj@metadata$Metrics[[y]]$global.level$ARI$hclust.complete.euclidian[[x]]$ari + 1)/2
                )
            })
        names(ari.scores) <- selected.vars

        ### Put all together ####
        ari.scores <- as.data.frame(x = ari.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'ARI') %>%
            data.frame(.)
        ari.scores <- ari.scores[ , cols.names]
    } else ari.scores <- NULL

    ## Gene variable correlations scores ####
    ### number of genes with a correlation cutoff ####
    if ('Correlation_corrCoeff' %in% paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') ){
        printColoredMessage(
            message = '- Summarizing the correlation coefficients of gene-variable correlation analysis for each continuous:',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the gene-variable correlations are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') == 'Correlation_corrCoeff'
        selected.vars <- metrics.table$Variables[selected.vars]
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$Correlation[[corr.method]][[x]]$correlations.pvalues)){
                            stop(paste0('The correlation for the ',
                                        x,
                                        ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeGenesVariableCorrelation" or "assessVariation" functions.'))
                        }
                    })
            })
        ### count number of genes ####
        printColoredMessage(
            message = '- Computing the correlation coefficients scores.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars.bio <- selected.vars[selected.vars %in% bio.variables]
        if (length(selected.vars.bio) > 0){
            gene.var.corr.coef.scores.bio <- lapply(
                selected.vars.bio,
                function(x){
                    sapply(
                        levels(assay.names),
                        function(y){
                            corr.results <- abs(se.obj@metadata$Metrics[[y]]$gene.level$Correlation$spearman[[x]]$correlations.pvalues)
                            sum(corr.results[ ,'correlation'] > corr.cutoff[[x]])/nrow(se.obj)
                        })
                })
            names(gene.var.corr.coef.scores.bio) <- selected.vars.bio
        } else gene.var.corr.coef.scores.bio <- NULL


        selected.vars.uv <- selected.vars[selected.vars %in% uv.variables]
        if (length(selected.vars.uv) > 0){
            gene.var.corr.coef.scores.uv <- lapply(
                selected.vars.uv,
                function(x){
                    sapply(
                        levels(assay.names),
                        function(y){
                            corr.results <- abs(se.obj@metadata$Metrics[[y]]$gene.level$Correlation$spearman[[x]]$correlations.pvalues)
                            sum(corr.results[ ,'correlation'] < corr.cutoff[[x]])/nrow(se.obj)
                        })
                })
            names(gene.var.corr.coef.scores.uv) <- selected.vars.uv
        } else gene.var.corr.coef.scores.uv <- NULL
        gene.var.corr.coef.scores <- c(gene.var.corr.coef.scores.bio, gene.var.corr.coef.scores.uv)

        ### Put all together ####
        gene.var.corr.coef.scores <- as.data.frame(x = gene.var.corr.coef.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level correlation (correlation cutoff)') %>%
            data.frame(.)
        gene.var.corr.coef.scores <- gene.var.corr.coef.scores[ , cols.names]
        } else gene.var.corr.coef.scores <- NULL

    ### p-value distributions  ####
    if ('Correlation_pvalueDis' %in% paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') ){
        printColoredMessage(
            message = '- Summarizing the p-value distribution of gene-variable correlation analysis for each continuous variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the gene-variable correlations are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') == 'Correlation_pvalueDis'
        selected.vars <- metrics.table$Variables[selected.vars]
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$Correlation$spearman[[x]]$correlations.pvalues)){
                            stop(paste0('The correlation for the ',
                                        x, ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeGenesVariableCorrelation" or "assessVariation" functions.'))
                        }
                    })
            })

        #### assess p-value distribution ####
        printColoredMessage(
            message = '- Computing the p-value distribution scores.',
            color = 'blue',
            verbose = verbose
            )
        gene.var.corr.pvalue.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        corr.results <- se.obj@metadata$Metrics[[y]]$gene.level$Correlation[[corr.method]][[x]]$correlations.pvalues
                        selected.pvalues <- corr.results[, 'p-value'][corr.results[, 'p-value'] > 0.05 ]
                        suppressWarnings(1 - ks.test(x = selected.pvalues, y = "punif")$statistic[[1]])
                    })
            })
        names(gene.var.corr.pvalue.scores) <- selected.vars

        #### Put all together ####
        gene.var.corr.pvalue.scores <- as.data.frame(x = gene.var.corr.pvalue.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level correlation (p-value distribution)') %>%
            data.frame(.)
        gene.var.corr.pvalue.scores <- gene.var.corr.pvalue.scores[ , cols.names]
    } else gene.var.corr.pvalue.scores <- NULL

    ### null p-value ####
    if ('Correlation_pvalueNull' %in% paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') ){
        printColoredMessage(
            message = '- Summarizing the Null p-values of the gene-variable correlation analysis for continuous variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the gene-variable correlations are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') == 'Correlation_pvalueNull'
        selected.vars <- metrics.table$Variables[selected.vars]
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$Correlation[[corr.method]][[x]]$correlations.pvalues)){
                            stop(paste0('The correlation for the ',
                                        x,
                                        ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeGenesVariableCorrelation" or "assessVariation" functions.'))
                        }
                    })
            })
        #### count number of genes ####
        printColoredMessage(
            message = '* compute the null p-value scores.',
            color = 'blue',
            verbose = verbose
        )
        gene.var.corr.qvalue.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        corr.results <- se.obj@metadata$Metrics[[y]]$gene.level$Correlation$spearman[[x]]$correlations.pvalues
                        suppressWarnings(qvalue::qvalue(p = corr.results[, 'p-value'])$pi0)
                    })
            })
        names(gene.var.corr.qvalue.scores) <- selected.vars

        #### Put all together ####
        gene.var.corr.qvalue.scores <- as.data.frame(x = gene.var.corr.qvalue.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level correlation (null p-value)') %>%
            data.frame(.)
        gene.var.corr.qvalue.scores <- gene.var.corr.qvalue.scores[ , cols.names]
    } else gene.var.corr.qvalue.scores <- NULL

    ## Gene variable anova scores ####
    ### number of genes with a f-stat cutoff ####
    if ('ANOVA_fStat' %in% paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') ){
        printColoredMessage(
            message = '- Summarizing the F-statistics of gene-variable ANOVA analysis for each categorical variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the ANOVA are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') == 'ANOVA_fStat'
        selected.vars <- metrics.table$Variables[selected.vars]
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$ANOVA[[anova.method]][[x]]$fstatistics.pvalues)){
                            stop(paste0('The ANOVA for the ', x, ' variable of the ' , y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeGenesVariableAnova" or "assessVariation" functions.'))
                        }
                    })
            })
        ### count number of genes ####
        printColoredMessage(
            message = '- Computing the F-statistics scores.',
            color = 'blue',
            verbose = verbose
        )
        gene.var.anova.fstat.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        corr.results <- se.obj@metadata$Metrics[[y]]$gene.level$ANOVA[[anova.method]][[x]]$fstatistics.pvalues
                        sum(corr.results[, 'statistic'] < fvalue.cutoff)/nrow(se.obj)
                    })
            })
        names(gene.var.anova.fstat.scores) <- selected.vars

        ### put all together ####
        gene.var.anova.fstat.scores <- as.data.frame(x = gene.var.anova.fstat.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level ANOVA (F-values cutoff)') %>%
            data.frame(.)
        gene.var.anova.fstat.scores <- gene.var.anova.fstat.scores[ , cols.names]
    } else gene.var.anova.fstat.scores <- NULL

    ### p-value distributions  ####
    if ('ANOVA_pvalueDis' %in% paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') ){
        printColoredMessage(
            message = '- Summarizing the p-value distribution of gene-variable ANOVA analysis for each categorical variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the ANOVA are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') == 'ANOVA_pvalueDis'
        selected.vars <- metrics.table$Variables[selected.vars]
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$ANOVA[[anova.method]][[x]]$fstatistics.pvalues)){
                            stop(paste0('The ANOVA for the ',
                                        x,
                                        ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeGenesVariableAnova" or "assessVariation" functions.'))
                        }
                    })
            })
        #### count number of genes ####
        printColoredMessage(
            message = '- Computing the p-value distribution scores.',
            color = 'blue',
            verbose = verbose
            )
        gene.var.anova.pvalue.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        anova.results <- se.obj@metadata$Metrics[[y]]$gene.level$ANOVA$aov[[x]]$fstatistics.pvalues
                        selected.pvalues <- anova.results[, 'pvalue'][anova.results[, 'pvalue'] > 0.05 ]
                        suppressWarnings(1-ks.test(x = selected.pvalues, y = "punif")$statistic[[1]])
                    })
            })
        names(gene.var.anova.pvalue.scores) <- selected.vars

        #### put all together ####
        gene.var.anova.pvalue.scores <- as.data.frame(x = gene.var.anova.pvalue.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level ANOVA (p-value distribution)') %>%
            data.frame(.)
        gene.var.anova.pvalue.scores <- gene.var.anova.pvalue.scores[ , cols.names]
    } else gene.var.anova.pvalue.scores <- NULL

    ### null p-value ####
    if ('ANOVA_pvalueNull' %in% paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') ){
        printColoredMessage(
            message = '- Summarizing the null p-values of gene-variable ANOVA analysis for each categorical variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the ANOVA are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') == 'ANOVA_pvalueNull'
        selected.vars <- metrics.table$Variables[selected.vars]
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$ANOVA$aov[[x]]$fstatistics.pvalues)){
                            stop(paste0('The ANOVA for the ', x, ' variable of the ' , y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeGenesVariableAnova" or "assessVariation" functions.'))
                        }
                    })
            })
        #### count number of genes ####
        printColoredMessage(
            message = '- Computing the null p-value scores.',
            color = 'blue',
            verbose = verbose
        )
        gene.var.anova.qvalue.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        corr.results <- se.obj@metadata$Metrics[[y]]$gene.level$ANOVA$aov[[x]]$fstatistics.pvalues
                        pvalues <- corr.results[, 'pvalue']
                        pvalues <- pvalues[!is.na(pvalues) & !is.infinite(pvalues)]
                        result <- tryCatch({
                            qvalue.results <- qvalue::qvalue(p = pvalues, pi0 = NULL)
                            list(pi0 = qvalue.results$pi0)
                        }, warning = function(w) {
                            message("Warning occurred: ", conditionMessage(w))
                            NULL
                        }, error = function(e) {
                            message("Error occurred: ", conditionMessage(e))
                            NULL
                        })
                        # Check results
                        if (!is.null(result)) {
                            return(result$pi0)
                        } else {
                            message("Failed to compute q-values. q-values will be 0.")
                            return(0)
                        }
                    })
            })
        names(gene.var.anova.qvalue.scores) <- selected.vars

        #### put all together ####
        gene.var.anova.qvalue.scores <- as.data.frame(x = gene.var.anova.qvalue.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'Gene-level ANOVA (null p-value)') %>%
            data.frame(.)
        gene.var.anova.qvalue.scores <- gene.var.anova.qvalue.scores[ , cols.names]
    } else gene.var.anova.qvalue.scores <- NULL

    ## Partial correlation scores ####
    ### number of genes ####
    if ('PartialCorrelation' %in% paste(metrics.table$Metrics) ){
        printColoredMessage(
            message = '- Summarizing the correlation differences in gene-gene partial correlation analysis for each continuous variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the gene-gene partial correlation are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- unique(metrics.table$Variables[metrics.table$Metrics == 'PartialCorrelation'])
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$PPcorr[[pcorr.method]][[x]]$correlations)){
                            stop(paste0('The gene-gene partial correlation for the ',
                                        x, ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "computeGenesPartialCorrelation" or "assessVariation" functions.'))
                        }
                    })
            })
        #### count number of genes ####
        printColoredMessage(
            message = '- Computing the correlation differences scores.',
            color = 'blue',
            verbose = verbose
            )
        ppcorr.corr.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        corr.results <- se.obj@metadata$Metrics[[y]]$gene.level$PPcorr$spearman[[x]]$correlations
                        diff.corr <- abs(corr.results[ , 'pp.cor'] - corr.results[ , 'p.cor'])
                        sum(diff.corr < pcorr.cutoff)/nrow(corr.results)
                     })
            })
        names(ppcorr.corr.scores) <- selected.vars

        #### put all together  ####
        ppcorr.corr.scores <- as.data.frame(x = ppcorr.corr.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'PPcorr (correlation differences cutoff)') %>%
            data.frame(.)
        ppcorr.corr.scores <- ppcorr.corr.scores[ , cols.names]
    } else ppcorr.corr.scores <- NULL

    ## DGE scores ####
    ### p-value distributions  ####
    if ('DGE_pvalueDis' %in% paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') ){
        printColoredMessage(
            message = '- Summarize the p-value distribution of differentail gene expressin analysis for each categorical variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '* check to see all the differentail gene expression are computed.',
            color = 'blue',
            verbose = verbose
        )
        selected.vars <- paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') == 'DGE_pvalueDis'
        selected.vars <- metrics.table$Variables[selected.vars]
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$DGE$limma[[x]]$p.values)){
                            stop(paste0('The differentail gene expressin for the ',
                                        x, ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "compuetDGE" or "assessVariation" functions.'))
                        }
                    })
            })
        #### compute the p-values scores ####
        printColoredMessage(
            message = '- Computing the p-value distribution scores.',
            color = 'blue',
            verbose = verbose
            )
        dge.pvalue.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        mean(sapply(
                            names(se.obj@metadata$Metrics[[y]]$gene.level$DGE$limma[[x]]$p.values),
                            function(z){
                                corr.results <- se.obj@metadata$Metrics[[y]]$gene.level$DGE$limma[[x]]$p.values[[z]]
                                selected.pvalues <- corr.results[, 'pvalue'][corr.results[, 'pvalue'] > 0.05 ]
                                suppressWarnings(1 - ks.test(x = selected.pvalues, y = "punif")$statistic[[1]])
                            }))
                    })
            })
        names(dge.pvalue.scores) <- selected.vars

        #### put all together  ####
        dge.pvalue.scores <- as.data.frame(x = dge.pvalue.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'DGE (p-value distribution)') %>%
            data.frame(.)
        dge.pvalue.scores <- dge.pvalue.scores[ , cols.names]
    } else dge.pvalue.scores <- NULL

    ### q-value ####
    if ('DGE_pvalueNull' %in% paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') ){
        printColoredMessage(
            message = '- Summarizing the null p-values of differentail gene expression analysis for each categorical variable(s):',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the differentail gene expression are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- paste(metrics.table$Metrics, metrics.table$Assessments, sep = '_') == 'DGE_pvalueNull'
        selected.vars <- metrics.table$Variables[selected.vars]
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$gene.level$DGE$limma[[x]]$p.values))
                            stop(paste0('The DGE for the ',  x, ' variable of the ' , y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "compuetDGE" or "assessVariation" functions.'))
                    })
            })
        #### compute the null p-value ####
        printColoredMessage(
            message = '* compute the null p-value scores.',
            color = 'blue',
            verbose = verbose
        )
        dge.qvalue.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        mean(sapply(
                            names(se.obj@metadata$Metrics[[y]]$gene.level$DGE$limma[[x]]$p.values),
                            function(z){
                                pvalue.results <- se.obj@metadata$Metrics[[y]]$gene.level$DGE$limma[[x]]$p.values[[z]]
                                suppressWarnings(qvalue::qvalue(p = pvalue.results[, 'pvalue'])$pi0)
                            }))
                    })
            })
        names(dge.qvalue.scores) <- selected.vars

        #### put all together ####
        dge.qvalue.scores <- as.data.frame(x = dge.qvalue.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'DGE (null p-value)') %>%
            data.frame(.)
        dge.qvalue.scores <- dge.qvalue.scores[ , cols.names]
    } else dge.qvalue.scores <- NULL

    ## Lisi scores ####
    if ('LISI' %in% metrics.table$Metrics){
        printColoredMessage(
            message = '- Summarizing the LISI scores for different variables:',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the LISI are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- metrics.table$Metrics == 'LISI'
        selected.vars <- metrics.table$Variables[selected.vars]
        if (isTRUE(fast.pca)){
            svd.method <- 'fast.svd'
        } else svd.method <- 'ordinary.svd'
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$global.level$LISI[[svd.method]][[x]]$lisi))
                            stop(paste0('The LISI for the ',  x, ' variable of the ' , y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "compuetLisi" or "assessVariation" functions.'))
                    })
            })
        #### compute the null p-value ####
        printColoredMessage(
            message = '- Computing the average of LISI scores.',
            color = 'blue',
            verbose = verbose
            )
        lisi.average.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        vec <- se.obj@metadata$Metrics[[y]]$global.level$LISI[[svd.method]][[x]]$lisi
                        mean((vec - min(vec)) / (max(vec) - min(vec)))
                    })
            })
        names(lisi.average.scores) <- selected.vars

        #### put all together ####
        lisi.average.scores <- as.data.frame(x = lisi.average.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'LISI (average scores)') %>%
            data.frame(.)
        lisi.average.scores <- lisi.average.scores[ , cols.names]
    } else lisi.average.scores <- NULL
    ## KBET scores ####
    if ('KBET' %in% metrics.table$Metrics){
        printColoredMessage(
            message = '- Summarizing the KBET scores for different variables:',
            color = 'orange',
            verbose = verbose
            )
        #### check correlation data ####
        printColoredMessage(
            message = '- Checking to see all the KBET are computed.',
            color = 'blue',
            verbose = verbose
            )
        selected.vars <- metrics.table$Metrics == 'KBET'
        selected.vars <- metrics.table$Variables[selected.vars]
        if (isTRUE(fast.pca)){
            svd.method <- 'fast.svd'
        } else svd.method <- 'ordinary.svd'
        check.out <- lapply(
            selected.vars,
            function(x){
                lapply(
                    levels(assay.names),
                    function(y){
                        if (is.null(se.obj@metadata$Metrics[[y]]$global.level$KBET[[svd.method]][[x]]$kBET))
                            stop(paste0('The KBET for the ',
                                        x,
                                        ' variable of the ' ,
                                        y ,
                                        ' data cannot be found in the SummarizedExperiment object.',
                                        'Please run the "compuetLisi" or "assessVariation" functions.'))
                    })
            })
        #### compute the null p-value ####
        printColoredMessage(
            message = '- Computing the average of KBET scores.',
            color = 'blue',
            verbose = verbose
            )
        kbet.average.scores <- lapply(
            selected.vars,
            function(x){
                sapply(
                    levels(assay.names),
                    function(y){
                        1- mean(se.obj@metadata$Metrics[[y]]$global.level$KBET[[svd.method]][[x]]$kBET)
                    })
            })
        names(kbet.average.scores) <- selected.vars

        #### put all together ####
        kbet.average.scores <- as.data.frame(x = kbet.average.scores, check.names = FALSE) %>%
            dplyr::mutate(data = row.names(.)) %>%
            pivot_longer(-data, names_to = 'variable', values_to = 'measurements') %>%
            dplyr::mutate(test = 'kBET (average scores)') %>%
            data.frame(.)
        kbet.average.scores <- kbet.average.scores[ , cols.names]
    } else kbet.average.scores <- NULL

    # Preparing the final matrix ####
    printColoredMessage(
        message = '-- Preparing final matrix to plot:',
        color = 'magenta',
        verbose = verbose
        )
    ## putting all measurements together ####
    printColoredMessage(
        message = '-- Putting all the computed score together.',
        color = 'blue',
        verbose = verbose
        )
    all.measurements <- rbind(
        general.rle.med.scores,
        general.rle.iqr.scores,
        rle.med.var.corr.scores,
        rle.iqr.var.corr.scores,
        rle.med.var.association.scores,
        rle.iqr.var.association.scores,
        pc.vec.corr.scores,
        pc.lin.reg.scores,
        sil.scores,
        ari.scores,
        gene.var.corr.coef.scores,
        gene.var.corr.pvalue.scores,
        gene.var.corr.qvalue.scores,
        gene.var.anova.fstat.scores,
        gene.var.anova.pvalue.scores,
        gene.var.anova.qvalue.scores,
        dge.pvalue.scores,
        dge.qvalue.scores,
        ppcorr.corr.scores,
        lisi.average.scores,
        kbet.average.scores
        )
    all.measurements[ ,'data'] <- gsub(':', ".", all.measurements[ , 'data'])
    assay.names <- gsub(':', ".", assay.names)

    ## adding the labels of both biological and unwanted variables ####
    printColoredMessage(
        message = '- Adding the biological or unwanted labels to the variables.',
        color = 'blue',
        verbose = verbose
        )
    all.measurements$group <- NULL
    if (is.null(bio.variables)){
        all.measurements$group <- 'Removal of unwanted variation'
    } else if (is.null(uv.variables)){
        all.measurements$group <- 'Preservation of biological variation'
    } else if (!is.null(bio.variables) & !is.null(uv.variables)){
        all.measurements$group[all.measurements$variable %in% bio.variables] <- 'Preservation of biological variation'
        all.measurements$group[all.measurements$variable %in% uv.variables] <- 'Removal of unwanted variation'
        all.measurements$group[all.measurements$variable %in% 'RLE'] <- 'Removal of unwanted variation'
    }
    if (sum(is.na(all.measurements$group)) > 0){
        stop('Please check the variable names.')
    }

    ## modifying the direction of some scores ####
    printColoredMessage(
        message = '- Modifying the direction of several scores for unwanted and biological variables.',
        color = 'orange',
        verbose = verbose
        )
    ### unwanted variables ####
    if (!is.null(uv.variables)) {
        vars.class <- sapply(
            uv.variables,
            function(x) class(colData(se.obj)[[x]])
            )
        uv.categorical.vars <- names(vars.class[vars.class %in% c('character', 'factor')])
        uv.continuous.vars <- names(vars.class[vars.class %in% c('numeric', 'integer')])
        #### continuous variables ####
        if (!is.null(uv.continuous.vars)){
            for(a in uv.continuous.vars){
                # LRA
                index.a <- all.measurements$variable == a & all.measurements$test == 'LRA'
                if (sum(index.a) > 0){
                    all.measurements$measurements[index.a] <- 1 - all.measurements$measurements[index.a]
                }
                # correlation with RLE median
                index.b <- all.measurements$variable == a & all.measurements$test == 'Association with RLE medians'
                if (sum(index.b) > 0){
                    all.measurements$measurements[index.b] <- 1 - all.measurements$measurements[index.b]
                }
                # correlation with RLE IQR
                index.c <- all.measurements$variable == a & all.measurements$test == 'Association with RLE IQRs'
                if (sum(index.c) > 0){
                    all.measurements$measurements[index.c] <- 1 - all.measurements$measurements[index.c]
                }
                rm(a)
            }
        }
        ### categorical variables ####
        if (!is.null(uv.categorical.vars)){
            for(b in uv.categorical.vars){
                # Silhouette
                index.a <- all.measurements$variable == b & all.measurements$test == 'Silhouette'
                if (sum(index.a) >0){
                    all.measurements$measurements[index.a] <- 1 - all.measurements$measurements[index.a]
                }
                # ARI
                index.b <- all.measurements$variable == b & all.measurements$test == 'ARI'
                if (sum(index.b) > 0){
                    all.measurements$measurements[index.b] <- 1 - all.measurements$measurements[index.b]
                }
                # VCA
                index.c <- all.measurements$variable == b & all.measurements$test == 'VCA'
                if (sum(index.c)){
                    all.measurements$measurements[index.c] <- 1 - all.measurements$measurements[index.c]
                }
                rm(b)
            }
        }
    }
    ### biological variables ####
    if (!is.null(bio.variables)) {
        vars.class <- sapply(
            bio.variables,
            function(x) class(colData(se.obj)[[x]])
            )
        bio.categorical.vars <- names(vars.class[vars.class %in% c('character', 'factor')])
        bio.continuous.vars <- names(vars.class[vars.class %in% c('numeric', 'integer')])
        #### continuous variables ####
        if (!is.null(bio.continuous.vars)){
            for(a in bio.continuous.vars){
                # correlation p-values
                index.a <- all.measurements$variable == a & all.measurements$test == 'Gene-level correlation (null p-value)'
                if (sum(index.a)){
                    all.measurements$measurements[index.a] <- 1 - all.measurements$measurements[index.a]
                }
                # LISI
                index.b <- all.measurements$variable == a & all.measurements$test == 'LISI (average scores)'
                if (sum(index.b) > 0){
                    all.measurements$measurements[index.b] <- 1 - all.measurements$measurements[index.b]
                }
                # correlation cutoff F-values
                # index.b <- all.measurements$variable == a & all.measurements$test == 'Gene-level correlation (correlation cutoff)'
                # if (sum(index.b) > 0){
                #     all.measurements$measurements[index.b] <- 1 - all.measurements$measurements[index.b]
                # }
                rm(a)
            }
        }
        ### categorical variables ####
        if (!is.null(bio.categorical.vars)){
            for(b in bio.categorical.vars){
                # ANOVA null p-values
                index.a <- all.measurements$variable == b & all.measurements$test == 'Gene-level ANOVA (null p-value)'
                if (sum(index.a) > 0){
                    all.measurements$measurements[index.a] <- 1 - all.measurements$measurements[index.a]
                }
                # ANOVA F-values
                index.b <- all.measurements$variable == b & all.measurements$test == 'Gene-level ANOVA (F-values cutoff)'
                if (sum(index.b) > 0){
                    all.measurements$measurements[index.b] <- 1 - all.measurements$measurements[index.b]
                }
                # DGE null pvalues
                index.c <- all.measurements$variable == b & all.measurements$test == 'DGE (null p-value)'
                if (sum(index.c) > 0){
                    all.measurements$measurements[index.c] <- 1 - all.measurements$measurements[index.c]
                }
                # PPcorr
                index.d <- all.measurements$variable == b & all.measurements$test == 'PPcorr (correlation differences cutoff)'
                if (sum(index.d) > 0){
                    all.measurements$measurements[index.d] <- 1 - all.measurements$measurements[index.d]
                }
                # LISI
                index.e <- all.measurements$variable == b & all.measurements$test == 'LISI (average scores)'
                if (sum(index.e) > 0){
                    all.measurements$measurements[index.e] <- 1 - all.measurements$measurements[index.e]
                }
            }
        }
    }

    ## average of biological and unwanted variables scores ######
    printColoredMessage(
        message = '- Averaging the scores for each biological and unwanted variables.',
        color = 'blue',
        verbose = verbose
        )
    bio.uv.overall.scores <- all.measurements %>%
        group_by(data, group) %>%
        summarise_at(vars("measurements"), mean) %>%
        dplyr::mutate(test = 'Score') %>%
        dplyr::mutate(variable = 'Score') %>%
        data.frame(.)
    bio.uv.overall.scores <- bio.uv.overall.scores [ , colnames(all.measurements)]
    ## calculate final scores ####
    printColoredMessage(
        message = '- Calculating final scores for each data.',
        color = 'blue',
        verbose = verbose
        )
    final.overall.scores <- lapply(
        assay.names,
        function(x){
            if (is.null(bio.variables)){
                bio.uv.overall.scores$measurements
            } else if (is.null(uv.variables)){
                bio.uv.overall.scores$measurements
            } else if (!is.null(bio.variables) & !is.null(uv.variables)){
                index.a <- bio.uv.overall.scores$data == x & bio.uv.overall.scores$group == 'Removal of unwanted variation'
                index.b <- bio.uv.overall.scores$data == x & bio.uv.overall.scores$group == 'Preservation of biological variation'
                c(uv.weight*bio.uv.overall.scores$measurements[index.a] + bio.weight*bio.uv.overall.scores$measurements[index.b])
            }
        })
    names(final.overall.scores) <- assay.names
    final.overall.scores <- final.overall.scores %>%
        data.frame(.) %>%
        pivot_longer(everything(), names_to = 'data', values_to = 'measurements') %>%
        dplyr::mutate(variable = 'Score') %>%
        dplyr::mutate(group = 'Final performance') %>%
        dplyr::mutate(test = 'Score') %>%
        data.frame()
    final.overall.scores <- final.overall.scores[ , colnames(all.measurements)]

    ## generate the plot ####
    printColoredMessage(
        message = '- Generating the final assessmet plot.',
        color = 'blue',
        verbose = verbose
        )
    if (isTRUE(select.top.ruv)){
        none.ruv <- assay.names[!assay.names %in% grep('RUVIIIprps', assay.names, value = TRUE)]
        final.overall.scores.ruv <- final.overall.scores[final.overall.scores$data %in% grep('RUVIIIprps', assay.names, value = TRUE) , ]
        final.overall.scores.ruv <- final.overall.scores.ruv[order(-final.overall.scores.ruv$measurements) , ]
        final.overall.scores.ruv <- final.overall.scores.ruv[1 , , drop = FALSE]
        final.overall.scores.none.ruv <- final.overall.scores[final.overall.scores$data %in% none.ruv , ]
        final.overall.scores <- rbind(final.overall.scores.ruv, final.overall.scores.none.ruv)
        all.measurements <- all.measurements[all.measurements$data %in% final.overall.scores$data , ]
        bio.uv.overall.scores <- bio.uv.overall.scores[bio.uv.overall.scores$data %in% final.overall.scores$data , ]
        assay.names <- assay.names[assay.names %in% bio.uv.overall.scores$data ]
    }
    all.measurements <- rbind(
        all.measurements,
        bio.uv.overall.scores,
        final.overall.scores
        )
    index.u <- all.measurements$variable == 'Score' & all.measurements$group == 'Removal of unwanted variation'
    all.measurements$variable[index.u] <- 'UV'
    all.measurements$group[index.u] <- 'Final performance'

    index.b <- all.measurements$variable == 'Score' & all.measurements$group == 'Preservation of biological variation'
    all.measurements$variable[index.b] <- 'Bio'
    all.measurements$group[index.b] <- 'Final performance'

    all.measurements$variable[all.measurements$variable == 'Score'] <- 'Final'
    final.overall.scores <- final.overall.scores[order(final.overall.scores$measurements, decreasing = FALSE) , ]
    all.measurements$data <- factor(
        x = all.measurements$data,
        levels = final.overall.scores$data
        )
    all.measurements$variable <- factor(
        x = all.measurements$variable,
        levels = unique(all.measurements$variable)
        )
    all.measurements$group.test.variable <- paste0(
        all.measurements$group,
        all.measurements$test,
        all.measurements$variable
        )
    all.measurements <- all.measurements %>%
        group_by(group.test.variable) %>%
        mutate(rank = row_number(-measurements)) %>%
        data.frame(.)
    all.measurements$rank <- as.factor(all.measurements$rank)
    ## plot
    strip_background <- strip_nested(
        text_x = elem_list_text(colour = "white", face = "bold"),
        background_x =
            elem_list_rect(
                fill = c(
                    # level1 colors
                    case_match(
                        unique(all.measurements$group),
                        "Final performance" ~ "darkgreen",
                        "Removal of unwanted variation" ~ "darkred",
                        .default = "orange3"
                    ),
                    # level2 colors
                    case_match(
                        all.measurements$group,
                        "RLE" ~ "grey",
                        .default = "grey")
                ))
    )
    num.colors <- length(assay.names)
    original.palette <- RColorBrewer::brewer.pal(n = 8, name = 'BrBG')[8:1]
    interpolated.colors <- grDevices::colorRampPalette(original.palette)(num.colors)
    all.measurements$group[all.measurements$group == 'Final performance'] <- 'Final scores'

    assessment.plot <- ggplot(data = all.measurements, aes(x = test, y = data)) +
        geom_point(aes(size = measurements, color = rank)) +
        scale_color_manual(values = interpolated.colors, name = 'Rank') +
        scale_size_continuous(name = "Measurements") +
        theme_bw()  +
        xlab('') +
        ylab('') +
        facet_nested(
            .~ group + variable,
            strip = strip_background,
            space = 'free_x',
            scales = 'free') +
        theme(
              panel.background = element_blank(),
              axis.line = element_line(colour = 'black', linewidth = 1),
              strip.background = element_rect(color = "grey30", fill = "grey90"),
              strip.text = element_text(size = c(12)),
              panel.border = element_rect(color = "grey90"),
              axis.line.y = element_line(colour = 'white', linewidth = 1),
              axis.text.x = element_text(size = 10, angle = 35, hjust = 1),
              axis.text.y = element_text(size = 12),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank()) +
        geom_stripes(odd = "#22222222", even = "#66666666") +
        guides(color = guide_legend(override.aes = list(size = 5), ncol = 3), size = guide_legend(ncol = 3))
    if (isTRUE(plot.output)) print(assessment.plot)

    # Saving results ####
    ### add results to the SummarizedExperiment object ####
    if (is.null(output.name))
        output.name <- paste0(
            'RUVprps_assessNormalization',
            length(assays(se.obj)),
            '_assays'
            )
    if (isTRUE(save.se.obj)){
        printColoredMessage(
            message = '-- Saving the selected bio genes to the metadata of the SummarizedExperiment object.',
            color = 'magenta',
            verbose = verbose
            )
        ## check if metadata metric already exist
        if (!'NormAssessment' %in% names(se.obj@metadata)) {
            se.obj@metadata[['NormAssessment']] <- list()
        }
        ## check if metadata metric already exist for this assay
        if (!output.name %in% names(se.obj@metadata[['NormAssessment']])) {
            se.obj@metadata[['NormAssessment']][[output.name]] <- list()
        }
        ## check if metadata metric already exist for this assay and this metric
        if (!'AssessmentTable' %in% names(se.obj@metadata[['NormAssessment']][[output.name]])) {
            se.obj@metadata[['NormAssessment']][[output.name]][['AssessmentTable']] <- list()
        }
        se.obj@metadata[['NormAssessment']][[output.name]][['AssessmentTable']] <- all.measurements

        if (!'AssessmentPlot' %in% names(se.obj@metadata[['NormAssessment']][[output.name]])) {
            se.obj@metadata[['NormAssessment']][[output.name]][['AssessmentPlot']] <- list()
        }
        se.obj@metadata[['NormAssessment']][[output.name]][['AssessmentPlot']] <- assessment.plot

        if (!'AssessmentTable' %in% names(se.obj@metadata[['NormAssessment']][[output.name]])) {
            se.obj@metadata[['NormAssessment']][[output.name]][['AssessmentTable']] <- list()
        }
        se.obj@metadata[['NormAssessment']][[output.name]][['AssessmentTable']] <- metrics.table

        printColoredMessage(
            message = '- The assessment table and plot are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The assessNormalization function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ### export results as list ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = '- The assessment table and plot are outputed as list.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The assessNormalization function finished.',
            color = 'white',
            verbose = verbose
            )
        return(list(
            AssessmentTable = all.measurements,
            AssessmentPlot = assessment.plot,
            AssessmentTable = metrics.table)
            )
    }
    return(se.obj)
}



