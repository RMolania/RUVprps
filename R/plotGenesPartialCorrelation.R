#' Plot gene-gene partial correlation results.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function computes all possible gene-gene pairwise ordinary and partial correlation of the data sets in the
#' SummarizedExperiment object. All genes or a subset of genes can specified.
#'
#' @details
#' Partial correlation is used to estimate correlation between two variables while controlling for third
#' variables.
#'
#' @param se.obj A `SummarizedExperiment` object.
#' @param assay.names Character. A character string or a vector of character strings specifying the names of the data sets
#' in the `SummarizedExperiment` object for which the computed correlations are to be obtained. The default is set to
#' "all", indicating that all data sets in the object will be selected.
#' @param variable Character. A character string indicating the column name in the `SummarizedExperiment` object
#' used to compute partial correlation.
#' @param method Character. A character string specifying the correlation method used for the variable. Options are
#' `pearson` or `spearman`. The default is set to `spearman`.
#' @param plot.type Character. A character string specifying the type of plot to generate for the partial pairwise gene
#' correlation. Options are: "scatter plot", "bar plot", and "histogram". The default is set to "bar plot".
#' @param filter.genes Logical. Indicates whether to filter genes before generating the plots. The default is set to `TRUE`.
#' @param corr.dif.cutoff Numeric. A numeric cutoff value used to filter genes. Any gene pair with a difference between the
#' ordinary and partial correlation less than absolute value of `corr.dif.cutoff` will be excluded. The default threshold
#' is set to 0.4.
#' @param plot.ncol Numeric. The number of columns in the plot grid. When more than one assay is selected, the function
#' arranges the plots in a grid accordingly. The default is set to 2.
#' @param plot.nrow Numeric. The number of rows in the plot grid. When more than three assays are selected, the function
#' arranges the plots in a grid accordingly. The default is set to 3.
#' @param plot.output Logical. If `TRUE`, displays the histogram or boxplot of partial and ordinary correlation coefficients.
#' The default is set to `FALSE`.
#' @param save.se.obj Logical. If `TRUE`, saves the results in the metadata of the `SummarizedExperiment` object.
#' If `FALSE`, outputs the result as a list. Tghe default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, displays messages for the different steps of the function.
#'
#' @return A `SummarizedExperiment` object containing the correlation plots, or a list of plots for individual assays.
#'
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

plotGenesPartialCorrelation <- function(
        se.obj,
        assay.names = "all",
        variable = NULL,
        method = 'spearman',
        plot.type = 'barplot',
        filter.genes = TRUE,
        corr.dif.cutoff = 0.4,
        plot.ncol = 2,
        plot.nrow = 3,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotGenesPartialCorrelation function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking function inputs ####
    if (is.null(assay.names) | is.logical(variable)) {
        stop('The "assay.names" cannot be NULL or logical.')
    }
    if (is.logical(variable) | is.null(variable)){
        stop('The "variable" cannot be NULL or logical.')
    }
    if (length(variable) > 1){
        stop('The "variable" must contain only one variable.')
    }
    if (!class(colData(se.obj)[[variable]]) %in% c('numeric', 'integer')){
        stop('The "variable" must be a continous variable.')
    }
    if (!is.logical(filter.genes)){
        stop('The "filter.genes" must be logical.')
    }
    if (isTRUE(filter.genes)){
        if (corr.dif.cutoff < 0 | corr.dif.cutoff > 1){
            stop('The "filter.genes" value must be postive value between 0 and 1.')
        }
    }
    if (plot.ncol < 1 | plot.nrow < 1){
        stop('The "plot.ncol" and "plot.nrow"  must be postive value.')
    }
    if (!is.logical(plot.output)){
        stop('The "plot.output" must be logical.')
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical.')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical.')
    }
    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Obtaining the correlations data ####
    printColoredMessage(
        message = '-- Obtaining the computed correlation data from the SummarizedExperiment object:',
        color = 'magenta',
        verbose = verbose
        )
    all.corr.data <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = levels(assay.names),
        assessment = 'PPcorr',
        assessment.type = 'gene.level',
        method = method,
        variables = variable,
        file.name = 'correlations',
        sub.file.name = NULL,
        required.function = 'computeGenesPartialCorrelation',
        message.to.print = 'partial correlation'
        )

    # Filter the correlation
    if (isTRUE(filter.genes)){
        all.corr.data <- lapply(
            levels(assay.names),
            function(x) {
                corr.dif <- all.corr.data[[x]]$p.cor - all.corr.data[[x]]$pp.cor
                select.genes <- abs(corr.dif) > corr.dif.cutoff
                if (sum(select.genes) == 0){
                    printColoredMessage(
                        message = paste0(
                            '- All gene-gene correlations are filtered for the "',
                            x,
                            '" data.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    all.corr.data[[x]] <- NULL
                } else if (sum(select.genes) != 0){
                    printColoredMessage(
                        message = paste0(
                            '- ',
                            nrow(all.corr.data[[x]]) - sum(select.genes)  ,
                            ' gene-gene correlation are filtered for "',
                            x,
                            '" data.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    all.corr.data[[x]][select.genes , ]
                }
            })
        names(all.corr.data) <- levels(assay.names)
        all.corr.data <- Filter(Negate(is.null), all.corr.data)
        if (length(all.corr.data) == 0)
            stop('Any assay contain gene-gene correlations for plotting.')
        assay.names <- droplevels(assay.names[assay.names %in% names(all.corr.data)])
    }

    # Generating different plots of partial and ordinary correlations ####
    ## generating histograms for each data ####
    if (plot.type == 'histogram'){
        printColoredMessage(
            message = '-- Generating the histograms of the correlations :',
            color = 'magenta',
            verbose = verbose
        )
        all.corr.hist <- lapply(
            levels(assay.names),
            function(x){
                corr.data <- as.data.frame(all.corr.data[[x]][ , c('p.cor', 'pp.cor')])
                colnames(corr.data) <- c('Correlation', 'Partial correlation')
                corr.data <- pivot_longer(corr.data, everything(), names_to = 'corr.type', values_to = 'corr')
                corr <- corr.type <- NULL
                p.joy <- ggplot(data = corr.data, aes(x = corr, y = corr.type)) +
                    geom_density(alpha = 0.4, color = "gray60") +
                    coord_cartesian(xlim = c(-1, 1)) +
                    ggtitle(x) +
                    xlab('Correlation coefficients') +
                    ylab('') +
                    labs(caption = paste0(
                        'Analysis: ',
                        'histograms of the partial and ordinary pair wise gene correlations\n',
                        "Variable: ",
                        variable)) +
                    theme_minimal() +
                    theme(
                        panel.grid.major.y = element_line(color = "black", linewidth = 0.5),
                        plot.caption = element_text(hjust = 0, vjust = 0),
                        plot.title = element_text(size = 16),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14, angle = 65),
                    )
                if (isTRUE(plot.output) & length(levels(assay.names)) == 1) print(p.joy)
                return(p.joy)
            })
        names(all.corr.hist) <- levels(assay.names)
        ### overall histograms of all data sets ####
        if (length(levels(assay.names)) > 1){
            for(i in levels(assay.names) ){
                all.corr.hist[[i]]$labels$caption <- NULL
            }
            overall.hist.plots <- ggarrange(
                plotlist = all.corr.hist ,
                ncol = plot.ncol,
                nrow = plot.nrow
            )
            if (class(overall.hist.plots)[[1]] == 'list'){
                plot.list <- lapply(
                    seq(length(overall.hist.plots)),
                    function(x){
                        annotate_figure(
                            p = overall.hist.plots[[x]],
                            top = text_grob(
                                label = "Partial pairwise correlation",
                                color = "orange",
                                face = "bold",
                                size = 18),
                            bottom = text_grob(
                                label = paste0(
                                    'Analysis: ',
                                    'histograms of the partial and ordinary pair wise gene correlations\n',
                                    "Variable: ",
                                    variable),
                                color = "black",
                                hjust = 1,
                                x = 1,
                                size = 10))
                    })
                overall.hist.plots <- ggarrange(
                    plotlist = plot.list,
                    ncol = 1,
                    nrow = 1)
            } else {
                overall.hist.plots <- annotate_figure(
                    p = overall.hist.plots,
                    top = text_grob(
                        label = "Partial pairwise correlation",
                        color = "black",
                        face = "bold",
                        size = 18
                    ),
                    bottom = text_grob(
                        label = paste0(
                            'Analysis: ',
                            'histograms of the partial and ordinary pair wise gene correlations\n',
                            "Variable: ",
                            variable),
                        color = "black",
                        hjust = 1,
                        x = 1,
                        size = 10))
            }
            printColoredMessage(
                message = '- All the scatter plots of all assays are combined into one plot.',
                color = 'blue',
                verbose = verbose
            )
            if (isTRUE(plot.output)) print(overall.hist.plots)
        }
    }

    ## generating scatter plots for all data sets  ####
    if (plot.type == 'scatter.plot'){
        printColoredMessage(
            message = '-- Generating scatter plots of the correlations :',
            color = 'magenta',
            verbose = verbose
        )
        all.corr.scatter <- lapply(
            levels(assay.names),
            function(x){
                corr.data <- as.data.frame(all.corr.data[[x]][ , c('p.cor', 'pp.cor')])
                Correlation <- Partial.correlation <- NULL
                colnames(corr.data) <- c('Correlation', 'Partial.correlation')
                p.hex <- ggplot(data = corr.data, aes(x = Correlation, y = Partial.correlation)) +
                    geom_hex(alpha = 0.6) +
                    xlim(-1,1) +
                    ylim(-1,1) +
                    ggtitle(x) +
                    geom_abline(intercept = 0, slope = 1, color = "red") +
                    xlab('Correlations') +
                    ylab('Partial correlations') +
                    labs(caption = paste0(
                        'Analysis: ',
                        'scatter plots of the partial and ordinary pair wise gene correlations\n',
                        "Variable: ",
                        variable)
                        ) +
                    theme_minimal() +
                    theme(
                        plot.caption = element_text(hjust = 0, vjust = 0),
                        plot.title = element_text(size = 16),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 12))
                if (isTRUE(plot.output) & length(levels(assay.names)) == 1) print(p.hex)
                return(p.hex)
            })
        names(all.corr.scatter) <- levels(assay.names)

        ### overall scatter plots of all data sets ####
        if (length(levels(assay.names)) > 1){
            for (i in levels(assay.names) ){
                all.corr.scatter[[i]]$labels$caption <- NULL
            }
            overall.scatter.plots <- ggarrange(
                plotlist = all.corr.scatter ,
                common.legend = TRUE,
                ncol = plot.ncol,
                nrow = plot.nrow
            )
            if (class(overall.scatter.plots)[[1]] == 'list'){
                plot.list <- lapply(
                    seq(length(overall.scatter.plots)),
                    function(x){
                        annotate_figure(
                            p = overall.scatter.plots[[x]],
                            top = text_grob(
                                label = "Partial pairwise correlation",
                                color = "orange",
                                face = "bold",
                                size = 18),
                            bottom = text_grob(
                                label = paste0(
                                    'Analysis: ',
                                    'histograms of the partial and ordinary pair wise gene correlations\n',
                                    "Variable: ",
                                    variable),
                                color = "black",
                                hjust = 1,
                                x = 1,
                                size = 10))
                    })
                overall.scatter.plots <- ggarrange(
                    plotlist = plot.list,
                    ncol = 1,
                    nrow = 1)
            } else {
                overall.scatter.plots <- annotate_figure(
                    p = overall.scatter.plots,
                    top = text_grob(
                        label = "Partial pairwise correlation",
                        color = "black",
                        face = "bold",
                        size = 18
                    ),
                    bottom = text_grob(
                        label = paste0(
                            'Analysis: ',
                            'scatter plots of the partial and ordinary pair wise gene correlations\n',
                            "Variable: ",
                            variable),
                        color = "black",
                        hjust = 1,
                        x = 1,
                        size = 10))
            }
            printColoredMessage(
                message = '- All the scatter plots of all assays are combined into one plot.',
                color = 'blue',
                verbose = verbose
            )
            if (isTRUE(plot.output)) print(overall.scatter.plots)
        }
    }
    ## generating barplot for each data set ####
    if (plot.type == 'barplot'){
        printColoredMessage(
            message = '-- Generating the barplot of the correlations :',
            color = 'magenta',
            verbose = verbose
        )
        all.corr.barplots <- lapply(
            levels(assay.names),
            function(x) {
                corr.dif <- all.corr.data[[x]]$p.cor - all.corr.data[[x]]$pp.cor
                select.genes <- abs(corr.dif) > corr.dif.cutoff
                p.barplot <- ggplot() +
                    geom_col(aes(x = 1, y = sum(select.genes))) +
                    ggtitle(x) +
                    xlab(x) +
                    ylab('Number of genes') +
                    labs(caption = paste0(
                        'Analysis: ',
                        'barplot of the number of genes that show difference between partial and ordinary correlations\n',
                        "Variable: ",
                        variable)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        plot.title = element_text(size = 16),
                        axis.text.x = element_text(size = 0),
                        axis.text.y = element_text(size = 12))
                if (isTRUE(plot.output) & length(levels(assay.names)) == 1) print(p.barplot)
                return(p.barplot)
            })
        names(all.corr.barplots) <- levels(assay.names)

        ### generate barplot for all data sets ####
        datasets <- corr.dif.group <- dd <- NULL
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = paste0('-- Putting all the barplots togather:'),
                color = 'magenta',
                verbose = verbose)
            all.corr.data <- lapply(
                levels(assay.names),
                function(x) {
                    temp.data <- all.corr.data[[x]][ , c('p.cor', 'pp.cor')]
                    temp.data$datasets <- rep(x, nrow(temp.data))
                    return(temp.data)
                })
            all.corr.data <- do.call(rbind, all.corr.data)
            all.corr.data$corr.dif <- all.corr.data$p.cor - all.corr.data$pp.cor
            all.corr.data$corr.dif.group <- ifelse(test = all.corr.data$corr.dif > corr.dif.cutoff, TRUE, FALSE)
            all.corr.data <- group_by(.data = all.corr.data, datasets)
            all.corr.data <- summarise(.data = all.corr.data, dd = sum(corr.dif.group))
            all.corr.data <- as.data.frame(all.corr.data)
            all.corr.data$datasets <- factor(
                x =  all.corr.data$datasets,
                levels = levels(assay.names))
            overall.barplots <- ggplot(all.corr.data, aes(x = datasets, y = dd)) +
                geom_col() +
                xlab('Datasets') +
                ylab('Number of genes') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    plot.title = element_text(size = 15),
                    axis.text.x = element_text(
                        size = 10,
                        angle = 25,
                        hjust = 1),
                    axis.text.y = element_text(size = 12))
            overall.barplots <- annotate_figure(
                p = overall.barplots,
                top = text_grob(
                    label = "Partial pairwise correlation",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ',
                        'number of paired genes that show at least ',
                        corr.dif.cutoff,
                        ' difference between partial and ordinary correlation\n',
                        "Variable: ", variable),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10)
                )
            printColoredMessage(
                message = '- The individual barplot are combined into one.',
                color = 'blue',
                verbose = verbose
                )
            if (isTRUE(plot.output))
                suppressMessages(print(overall.barplots))
        }
    }

    # Saving the plots ####
    ## adding results to the SummarizedExperiment object ####
    printColoredMessage(
        message = '-- Saving all the plots :',
        color = 'magenta',
        verbose = verbose)
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving all the plots to the "metadata" of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose)
        # for all assays

        if (plot.type == 'barplot'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'gene.level',
                assessment = 'PPcorr',
                method = method,
                variables = variable,
                file.name = 'barplot',
                results.data = all.corr.barplots
            )
        }
        if (plot.type == 'histogram'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'gene.level',
                assessment = 'PPcorr',
                method = method,
                variables = variable,
                file.name = 'histogram',
                results.data = all.corr.hist
            )
        }
        if (plot.type == 'scatter.plot'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'gene.level',
                assessment = 'PPcorr',
                method = method,
                variables = variable,
                file.name = 'scatter.plot',
                results.data = all.corr.scatter
            )
        }
        printColoredMessage(
            message = '- The plots of individual assay is saved to the metadata@metric in SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        # Overall plot
        if (length(levels(assay.names)) > 1){
            if (plot.type == 'barplot'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'gene.level',
                    assessment = 'PPcorr',
                    method = method,
                    variables = variable,
                    file.name = 'barplot',
                    plot.data = overall.barplots
                )
            }
            if (plot.type == 'histogram'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'gene.level',
                    assessment = 'PPcorr',
                    method = method,
                    variables = variable,
                    file.name = 'histogram',
                    plot.data = overall.hist.plots
                )
            }
            if (plot.type == 'scatter.plot'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'gene.level',
                    assessment = 'PPcorr',
                    method = method,
                    variables = variable,
                    file.name = 'scatter.plot',
                    plot.data = overall.scatter.plots
                )
            }
            printColoredMessage(
                message = paste0(
                    '- The combined plots of all assays are saved to the',
                    ' "se.obj@metadata$plot$RLE" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(
            message = '------------The plotGenesPartialCorrelation function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj = se.obj)
    }

    ## outputting the plots as alist ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = paste0('- The plots of individual assay is outputed as list.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The plotGenesPartialCorrelation function finished.',
            color = 'white',
            verbose = verbose
            )
        if (plot.type == 'barplot'){
            if (length(levels(assay.names)) > 1){
                return(list(
                    all.corr.barplots = all.corr.barplots,
                    overall.barplots = overall.barplots)
                    )
            } else return(all.corr.barplots = all.corr.barplots)
        }
        if (plot.type == 'histogram'){
            if (length(levels(assay.names)) > 1){
                return(list(
                    all.corr.hist = all.corr.hist,
                    overall.hist.plots = overall.hist.plots)
                    )
            } else return(all.corr.hist = all.corr.hist)
        }
        if (plot.type == 'scatter.plot'){
            if (length(levels(assay.names)) > 1){
                return(list(
                    all.corr.scatter = all.corr.scatter,
                    overall.scatter.plots = overall.scatter.plots)
                    )
            } else return(all.corr.scatter = all.corr.scatter)
        }
    }
}


