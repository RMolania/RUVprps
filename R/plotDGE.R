#' Plot histograms of p-values from differential gene expression analysis
#'
#' @author Ramyar Molania
#'
#' @description
#' Generates histograms of p-values obtained from differential gene expression
#' (DGE) analyses to assess the presence and impact of unwanted sources of
#' variation, such as batch effects, across multiple assays.
#'
#' @details
#' Differential gene expression (DGE) analysis is performed using the Wilcoxon
#' signed-rank test on log-transformed expression values (e.g., raw counts,
#' normalized data, or batch-corrected data). To evaluate the influence of
#' unwanted variation, DGE analysis is conducted across batches or other
#' categorical variables. In the absence of systematic bias or batch effects,
#' the histogram of unadjusted p-values is expected to follow a uniform
#' distribution. Deviations from uniformity may indicate the presence of
#' technical artifacts or confounding factors.
#'
#' @param se.obj A \code{SummarizedExperiment} object containing expression
#'   data and sample-level metadata.
#' @param assay.names Character vector specifying the name(s) of assay(s) in
#'   the \code{SummarizedExperiment} object to be used for DGE analysis.
#'   The default is \code{"all"}, indicating that all available assays will be
#'   evaluated.
#' @param variable Character string specifying the column name in
#'   \code{colData(se.obj)} containing the categorical variable of interest
#'   (e.g., batch, condition, or sample type).
#' @param method Character string specifying the differential expression method
#'   to be used.
#' @param plot.ncol Numeric. Number of columns in the plot grid when multiple
#'   datasets are provided (default: 1).
#' @param plot.nrow Numeric. Number of rows in the plot grid when multiple
#'   datasets are provided (default: 1).
#' @param plot.output Logical. If \code{TRUE}, p-value histograms are displayed
#'   during execution. Default is \code{FALSE}.
#' @param save.se.obj Logical. If \code{TRUE}, the results are stored in the
#'   metadata of the input \code{SummarizedExperiment} object and returned.
#'   If \code{FALSE}, the function returns the results as a separate object.
#'   Default is \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, informative messages describing
#'   function progress are printed.
#'
#' @return A \code{SummarizedExperiment} object with appended metadata or a
#'   list containing the computed assessment results.
#'
#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#'
#' @export

plotDGE <- function(
        se.obj,
        assay.names = 'all',
        variable,
        method = 'limma',
        plot.ncol = 1,
        plot.nrow = 2,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotDGE function starts:',
                        color = 'white',
                        verbose = verbose)
    # Check the inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be empty.')
    } else if (is.null(variable)) {
        stop('The "variable" cannot be empty.')
    }

    # Check the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Obtain computed p-values for each DE contrasts ####
    printColoredMessage(
        message = paste0('-- Obtaining computed p-values of all contrasts of the DGE analysis for the "', variable, '" variable.') ,
        color = 'magenta',
        verbose = verbose
    )
    all.de.tests <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = levels(assay.names),
        assessment = 'DGE',
        assessment.type = 'gene.level',
        method = method,
        variables = variable,
        file.name = 'p.values',
        sub.file.name = NULL,
        required.function = 'computeDGE',
        message.to.print = 'DGE'
    )
    # Generate p-values histograms for each contrast ####
    ## specified ylim for histograms ####
    breaks <- seq(from = 0, to = 1, by = .1)
    ylim.pvalue <- sapply(
        levels(assay.names),
        function(x) {
            sapply(
                names(all.de.tests[[x]]),
                function(y){
                    binned <- cut(
                        x = all.de.tests[[x]][[y]][['pvalue']],
                        breaks = breaks,
                        include.lowest = TRUE)
                    frequency <- table(binned)[1]
                })
        })
    ylim.pvalue <- ceiling(x = max(ylim.pvalue))

    ## generate p-values histograms for each assay  ####
    all.pval.histograms <- lapply(
        levels(assay.names),
        function(x){
            printColoredMessage(
                message = paste0('- Generate p-values histograms for the ', x, ' data.'),
                color = 'blue',
                verbose = verbose
            )
            if (length(unique(colData(se.obj)[[variable]])) == 2 ){
                pval.data <- do.call(rbind, all.de.tests[[x]])
                pval.data$contrasts <- rep(names(all.de.tests[[x]]), each = nrow(se.obj))
                pval.plot <- ggplot(pval.data, aes(x = pvalue)) +
                    geom_histogram(binwidth = 0.1) +
                    scale_y_continuous(
                        labels = function(x) format(x / 1000, scientific = F),
                        limits = c(0, ylim.pvalue)) +
                    ggtitle(x) +
                    xlab('p-values') +
                    ylab(expression('Frequency'~10^3)) +
                    theme(
                        panel.background = element_blank(),
                        plot.title = element_text(size = 12),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8)
                        )
            } else if (length(unique(colData(se.obj)[[variable]])) > 2 ) {
                pval.data <- do.call(rbind, all.de.tests[[x]])
                pval.data$contrasts <- rep(names(all.de.tests[[x]]), each = nrow(se.obj))
                pval.plot <- ggplot(pval.data, aes(x = pvalue)) +
                    geom_histogram(binwidth = 0.1) +
                    scale_y_continuous(labels = function(x) format(x /1000, scientific = F), limits = c(0, ylim.pvalue)) +
                    ggtitle(x) +
                    xlab('p-values') +
                    ylab(expression('Frequency'~10^3)) +
                    facet_wrap(~contrasts) +
                    theme(
                        panel.background = element_blank(),
                        plot.title = element_text(size = 12),
                        strip.text = element_text(size = 8),
                        strip.text.x.top = element_text(size = 8),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 10),
                        axis.title.y = element_text(size = 10),
                        axis.text.x = element_text(size = 8),
                        axis.text.y = element_text(size = 8))
            }
            if (isTRUE(plot.output) & length(assay.names) == 1) print(pval.plot)
            return(pval.plot)
        })
    names(all.pval.histograms) <- levels(assay.names)

    ## put all the p-values histograms into one  ####
    printColoredMessage(
        message = '- Putting all the p-values histograms together:' ,
        color = 'magenta',
        verbose = verbose
        )
    overall.pval.histograms <- ggpubr::ggarrange(
        plotlist = all.pval.histograms,
        ncol = plot.ncol,
        nrow = plot.nrow
        )
    if (class(overall.pval.histograms)[[1]] == 'list'){
        plot.list <- lapply(
            seq(length(overall.pval.histograms)),
            function(x){
                annotate_figure(
                    p = overall.pval.histograms[[x]],
                    top = text_grob(
                        label = "Differential gene expression analysis.",
                        color = "black",
                        face = "bold",
                        size = 18),
                    bottom = text_grob(
                        label = paste0(
                            'Analysis: ',
                            'differential gene expression analysis using Wilcoxon for all possible contrasts\n',
                            "Variable: ",
                            variable),
                        color = "black",
                        hjust = 1,
                        x = 1,
                        size = 10))
            })
        overall.pval.histograms <- ggarrange(
            plotlist = plot.list,
            ncol = 1,
            nrow = 1
            )
    } else {
        overall.pval.histograms <- annotate_figure(
            p = overall.pval.histograms,
            top = text_grob(
                label = "Differential gene expression analysis.",
                color = "black",
                face = "bold",
                size = 18),
            bottom = text_grob(
                label = paste0(
                    'Analysis: ',
                    'differential gene expression analysis using Wilcoxon for all possible contrasts\n',
                    "Variable: ",
                    variable),
                color = "black",
                hjust = 1,
                x = 1,
                size = 10))
    }
    if (isTRUE(plot.output)) suppressMessages(print(overall.pval.histograms))

    # Save the plots ####
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'gene.level',
            assessment = 'DGE',
            method = method,
            variables = variable,
            file.name = 'plot',
            results.data = all.pval.histograms
            )
        printColoredMessage(
            message = 'The Wilcoxon results for indiviaul assay are saved to metadata@metric',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(message = '------------The genesDEA function finished.',
                            color = 'white',
                            verbose = verbose)
        if (length(assay.names) > 1){
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'gene.level',
                assessment = 'DGE',
                method = method,
                variables = variable,
                file.name = 'histogram',
                plot.data = overall.pval.histograms
            )
            printColoredMessage(
                message = 'The p-value histograms of all assays are saved to metadata@plot',
                color = 'blue',
                verbose = verbose
            )
        }
        printColoredMessage(
            message = '------------The plotDGE function finished.',
            color = 'white',
            verbose = verbose
        )
        return(se.obj = se.obj)
    }
    ## return the results as a list ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = 'All the p-value histograms are outputed as list.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The plotDGE function finished.',
            color = 'white',
            verbose = verbose
            )
        if (length(assay.names) == 1){
            return(all.pval.histograms = all.pval.histograms)
        } else {
            return(list(
                all.pval.histograms = all.pval.histograms,
                overall.pval.histograms = overall.pval.histograms)
            )
        }
    }
}




