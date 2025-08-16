#' Computes the average Silhouette coefficients.
#'
#' @author Ramyar Molania
#'
#' @references
#' *Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @description
#' This functions generates barplots of average Silhouette coefficients for individual assays. If two variables are
#' provided, the function creates combined scatter plots of the average silhouette coefficients of each variables for the
#' individual assays.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or vector of character strings specifying the name(s) of the data set(s)
#' in the SummarizedExperiment object to generate barplots or scatter plots of the computed Silhouette coefficients.
#' By default, all assays in the SummarizedExperiment object will be selected.
#' @param variables Character. A character or character vectors of the column names in the SummarizedExperiment object
#' that the silhouette coefficient has been calculated by the `computeSilhouette` function.
#' @param silhouette.method Character. Indicates which computed Silhouette method should be used for plotting. The method
#' must match one of those defined in the `computeSilhouette()` function. The default is set to `sil.euclidian`, which is
#' also the default of the `computeSilhouette()` function. Refer to the `computeSilhouette()` function for more details.
#' @param plot.type Character. Indicates how to plot the silhouette coefficients. Options are "single.plot" and `combined.plot`.
#' If a single variable is provided, `plot.type` must be set to `single.plot`, which generates a barplot of the Silhouette
#' coefficients. If two variables are provided and `plot.type` is set to `combined.plot`, the function generates a scatter
#' plot of the Silhouette coefficients for each variable against each other. The default is set to `single.plot`.
#' @param plot.output Logical. If `TRUE,` the individual barplots or scatter plots will be printed during function execution.
#' The default is set to `TRUE`.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the SummarizedExperiment object
#' or to return the results as a list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, process messages will be displayed.
#'
#' @return A SummarizedExperiment object or a list containing all the plots of the computed average Silhouette coefficients
#' for the categorical variables.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

plotSilhouette <- function(
        se.obj,
        assay.names = 'all',
        variables,
        silhouette.method = 'sil.euclidian',
        plot.type = 'single.plot',
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotSilhouette function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the function inputs ####
    if (is.null(assay.names)) {
        stop('The "assay.names" cannot be NULL')
    }
    if (is.list(assay.names) | !is.character(assay.names)){
        stop('The "assay.names" must be a vector of the assay names(s) or "assay.names == all".')
    }
    if (is.null(variables)) {
        stop('The "variables" cannot be NULL')
    }
    if (!is.character(variables)){
        stop('The "variables" must a character or a vector of characters.')
    }
    if (!plot.type %in% c('single.plot', 'combined.plot')) {
        stop('The "plot.type" must be one of the "single.plot" or "combined.plot".')
    }
    if (plot.type == 'combined.plot') {
        if (length(variables) == 1)
            stop('To plot combined Silhouette, two variables must be provided.')
    }

    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Plotting the Silhouette coefficients ####
    ## Single plot ####
    if (plot.type == 'single.plot') {
        ### Obtaining silhouette ####
        printColoredMessage(
            message = paste0('-- Obtaining the computed Silhouette coefficient from the SummarizedExperiment object:'),
            color = 'magenta',
            verbose = verbose
        )
        all.silhouette <- getMetricFromSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment = 'Silhouette',
            assessment.type = 'global.level',
            method = silhouette.method,
            variables = variables,
            file.name = 'silhouette.coeff',
            sub.file.name = NULL,
            required.function = 'computeSilhouette',
            message.to.print = 'silhouette coefficient'
            )

        ## Plotting the Silhouette coefficient for individual assay ####
        printColoredMessage(
            message = paste0('-- Generating barplot of the Silhouette coefficient for the individual data set(s):'),
            color = 'magenta',
            verbose = verbose
            )
        all.single.silhouette.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0(
                        '- Creating barplot of the silhouette coefficient for the "'
                        , x,
                        '" data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                p.silhouette <- ggplot() +
                    geom_col(aes(y = all.silhouette[[x]], x = 1)) +
                    ylab('Silhouette ') +
                    xlab(x) +
                    ggtitle(variables) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        plot.title = element_text(size = 16),
                        axis.text.x = element_text(size = 0),
                        axis.text.y = element_text(size = 12))
                if (plot.output & length(assay.names) == 1) print(p.silhouette)
                return(p.silhouette)
            })
        names(all.single.silhouette.plots) <- levels(assay.names)
        everything <- datasets <- silhou.coff <- NULL

        ## Putting all the plots of individual assays ####
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = '- Putting all the silhouette coefficient togather:',
                color = 'magenta',
                verbose = verbose
                )
            overall.single.silhouette.plot <- as.data.frame(all.silhouette) %>%
                tidyr::pivot_longer(
                    everything(),
                    names_to = 'datasets',
                    values_to = 'silhou.coff'
                    )
            overall.single.silhouette.plot$datasets <- factor(
                x = overall.single.silhouette.plot$datasets,
                levels = levels(assay.names))
            overall.single.silhouette.plot <- ggplot(overall.single.silhouette.plot,
                       aes(x = datasets, y = silhou.coff)) +
                geom_col() +
                ggtitle(variables) +
                xlab('Datasets') +
                ylab('Silhouette coefficient ') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    plot.title = element_text(size = 15),
                    axis.text.x = element_text(
                        size = 12,
                        angle = 25,
                        hjust = 1),
                    axis.text.y = element_text(size = 12))
            overall.single.silhouette.plot <- annotate_figure(
                p = overall.single.silhouette.plot,
                top = text_grob(
                    label = "Average silhouette coefficient",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ',
                        'average silhouette coefficient using the first PCs and the ',
                        variables,
                        ' variable.'),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10)
            )
            printColoredMessage(
                message = '- The individual assay silhouette coefficient barplot are combined into one.',
                color = 'blue',
                verbose = verbose)
            if (isTRUE(plot.output))
                suppressMessages(print(overall.single.silhouette.plot))
        }

    }
    ## Combined plot ####
    if (plot.type == 'combined.plot') {
        printColoredMessage(
            message = paste0('-- Obtain computed silhouette for the from the SummarizedExperiment object:'),
            color = 'magenta',
            verbose = verbose
        )
        all.silhouette <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('- Obtain silhouette for the "', x, '" data.'),
                    color = 'blue',
                    verbose = verbose
                )
                if (!silhouette.method %in% names(se.obj@metadata[['Metrics']][[x]][['global.level']][['Silhouette']])) {
                    stop(paste0('The ', silhouette.method ,'has not been computed yet for the ', variables,' variable and the ', x,' data.'))
                }
                for (i in variables) {
                    if (!i %in% names(se.obj@metadata[['Metrics']][[x]][['global.level']][['Silhouette']][[silhouette.method]])) {
                        stop(paste0('The ', silhouette.method ,' has not been computed yet for the ', i, ' variable and the ', x,' data.'))
                    }
                }
                silhouette <- c()
                for (i in 1:length(variables))
                    silhouette[i] <-
                    se.obj@metadata[['Metrics']][[x]][['global.level']][['Silhouette']][[silhouette.method]][[variables[i]]]$silhouette.coeff
                return(silhouette)
            })
        names(all.silhouette) <- levels(assay.names)

        ### Individual plots ####
        datasets <- NULL
        all.combined.silhouette.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0('- Plot combined Silhouettes for the "', x, '" data.'),
                    color = 'blue',
                    verbose = verbose
                )
                all.silhouettes <- as.data.frame(t(as.data.frame(all.silhouette[[x]])))
                row.names(all.silhouettes) <- x
                colnames(all.silhouettes) <- variables
                all.silhouettes$datasets <- row.names(all.silhouettes)
                p.combined <- ggplot(all.silhouettes, aes_string(x = sym(variables[1]), y = sym(variables[2]) ) ) +
                    geom_point() +
                    ggrepel::geom_text_repel(
                        aes(label = datasets),
                        hjust = 0,
                        vjust = 0) +
                    xlab(paste0('Silhouette (', variables[1], ')')) +
                    ylab(paste0('Silhouette (', variables[2], ')')) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 18),
                        axis.title.y = element_text(size = 18),
                        plot.title = element_text(size = 15),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 12)
                    )
                if (isTRUE(plot.output) & length(assay.names) == 1) print(p.combined)
                return(p.combined)
            })
        names(all.combined.silhouette.plots) <- levels(assay.names)

        ## Overall plots ####
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = paste0('- Plot Silhouette for all the assays(s)'),
                color = 'blue',
                verbose = verbose
            )
            all.silhouette <- as.data.frame(t(as.data.frame(all.silhouette)))
            colnames(all.silhouette) <- variables
            all.silhouette$datasets <- row.names(all.silhouette)
            overall.combined.silhouette.plot <- ggplot(all.silhouette, aes_string(x = sym(variables[1]), y = sym(variables[2]))) +
                geom_point() +
                ggrepel::geom_text_repel(aes(label = datasets),
                                hjust = 0,
                                vjust = 0) +
                xlab(paste0('Silhouette (', variables[1], ')')) +
                ylab(paste0('Silhouette (', variables[2], ')')) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    plot.title = element_text(size = 15),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12))
            overall.combined.silhouette.plot <- annotate_figure(
                p = overall.combined.silhouette.plot,
                top = text_grob(
                    label = "Average silhouette coefficient",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ', 'average silhouette coefficient using the first PCs and the ', variables, ' variable.'),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10)
            )
        }
    }

    # Saving the results ####
    printColoredMessage(
        message = '-- Saving the silhouette  barplots:',
        color = 'magenta',
        verbose = verbose
        )
    ## Adding results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '- Saving all the silhouette coefficient barplots to the "metadata" in the SummarizedExperiment object:',
            color = 'blue',
            verbose = verbose
            )
        if (plot.type == 'single.plot'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'global.level',
                assessment = 'Silhouette',
                method = silhouette.method,
                variables = variables,
                file.name = 'single.plot',
                results.data = all.single.silhouette.plots
            )
        } else if (plot.type == 'combined.plot'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = assay.names,
                assessment.type = 'global.level',
                assessment = 'Silhouette',
                method = silhouette.method,
                variables = paste0(variables, collapse = '&'),
                file.name = 'combined.plot',
                results.data = all.combined.silhouette.plots
            )
        }
        printColoredMessage(
            message = paste0(
                '- The Silhouette barplot of the individual assay(s) is saved to the ',
                ' "se.obj@metadata$metric$AssayName$ARI" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )

        if (length(assay.names) > 1) {
            if (plot.type == 'single.plot'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'global.level',
                    assessment = 'Silhouette',
                    method = silhouette.method,
                    variables = variables,
                    file.name = 'single.plot',
                    plot.data = overall.single.silhouette.plot
                )
            } else if (plot.type == 'combined.plot'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'global.level',
                    assessment = 'Silhouette',
                    method = silhouette.method,
                    variables = paste0(variables, collapse = '&'),
                    file.name = 'combined.plot',
                    plot.data = overall.combined.silhouette.plot
                )
            }
            printColoredMessage(
                message = paste0(
                    '- The combined silhouette coefficient barplot all the assays is saved to the ',
                    ' "se.obj@metadata$plot$Silhouette" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        }
        printColoredMessage(
            message = '------------The plotSilhouette function finished.',
            color = 'white',
            verbose = verbose)
        return(se.obj = se.obj)

        ## Returning only the Silhouette coefficients results ####
    }
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = paste0('-All the Silhouette plots re saved as list.'),
            color = 'blue',
            verbose = verbose)
        printColoredMessage(
            message = '------------The plotSilhouette function finished.',
            color = 'white',
            verbose = verbose)
        if (plot.type == 'single.plot') {
            if (length(assay.names) == 1){
                return(all.silhouette.plots = list(all.single.silhouette.plots = all.single.silhouette.plots))
            } else {
                return(all.silhouette.plots = list(
                        all.single.silhouette.plots = all.single.silhouette.plots,
                        overall.single.silhouette.plot = overall.single.silhouette.plot ))
            }
        } else if (plot.type == 'combined.plot') {
            if (length(assay.names) == 1){
                return(all.silhouette.plots = list(all.combined.silhouette.plots = all.combined.silhouette.plots))
            } else{
                return(all.silhouette.plots = list(
                        all.combined.silhouette.plots = all.combined.silhouette.plots,
                        overall.combined.silhouette.plot = overall.combined.silhouette.plot))
            }
        }
    }
}
