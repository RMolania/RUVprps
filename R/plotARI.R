#' Generate barplot or scatter plot the Adjusted Rand Index (ARI).
#'
#' @author Ramyar Molania
#'
#' @references
#' Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @description
#' This function generates barplots of the Adjusted Rand Index (ARI) for individual datasets in the `SummarizedExperiment`
#' object. If two variables are provided, the function creates scatter plots comparing the ARIs of each variable across
#' the assays. The `computeARI()` function must applied before using the `plotARI()` function.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character or character vector. Specifies the name(s) of the assay(s) in the `SummarizedExperiment`
#' object to be used for generating barplots or scatter plots of the computed Adjusted Rand Index (ARI). By default,
#' all assays in the `SummarizedExperiment` object will be selected.
#' @param variables Character or character vector of length one or two. Indicates one or two column names in the
#' `SummarizedExperiment` object that contain categorical variables, such as sample subtypes or batch labels. If two
#' variables are provided, the function plots the ARIs against each other for all specified assays.
#' @param ari.method Character. Specifies the ARI computation method to use for plotting. This must match one of the
#' methods supported by the `computeARI` function. The default is `"hclust.complete.euclidean"`, which is also the default
#'  in `computeARI`. See the `computeARI` function documentation for more details.
#' @param plot.type Character. Specifies the plot type for the Adjusted Rand Index. Options are `"single.plot"` and
#' `"combined.plot"`. If one variable is provided in `variables`, `plot.type` must be set to `"single.plot"` to generate
#' a barplot. If two variables are provided and `plot.type` is set to `"combined.plot"`, the function generates a scatter
#' plot comparing the ARIs. The default is `"single.plot"`.
#' @param plot.output Logical. If `TRUE`, the individual barplots or scatter plots will be printed during function execution.
#' Default is `TRUE`.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the `SummarizedExperiment` object
#' or return them as a list. The default is `TRUE`.
#' @param verbose Logical. If `TRUE`, messages for different steps of the function will be displayed.
#'
#' @return A `SummarizedExperiment` object or a list containing all the plots of the computed ARIs for the categorical
#' variable(s).
#'
#' @importFrom SummarizedExperiment assays assay
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

plotARI <- function(
        se.obj,
        assay.names = 'all',
        variables,
        ari.method = 'hclust.complete.euclidian',
        plot.type = 'single.plot',
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotARI function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the inputs ####
    if (is.null(assay.names) | is.logical(assay.names)) {
        stop('The "assay.names" cannot be NUll or logical.')
    }
    if (is.null(variables) | is.logical(variables)) {
        stop('The "variables" cannot be NULL or logical.')
    }
    if (!plot.type %in% c('single.plot', 'combined.plot')) {
        stop('The "plot.type" must be one of the "single.plot" or "combined.plot".')
    }
    if (plot.type == 'combined.plot') {
        if (length(variables) == 1)
            stop('To plot combined ARI, two variables must be provided.')
    }
    if (plot.type == 'single.plot') {
        if (length(variables) > 1)
            stop('To plot "single.plot" ARI, only one variable must be provided.')
    }

    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }

    # Plotting the ARI values ####
    ## Single plot ####
    if (plot.type == 'single.plot') {
        ### obtain ari ####
        printColoredMessage(
            message = paste0('-- Obtaining computed ARI from the SummarizedExperiment object:'),
            color = 'magenta',
            verbose = verbose
            )
        all.ari <- getMetricFromSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = assay.names,
            assessment = 'ARI',
            assessment.type = 'global.level',
            method = ari.method,
            variables = variables,
            file.name = 'ari',
            sub.file.name = NULL,
            required.function = 'compuetARI',
            message.to.print = 'ARI'
            )
        ## plot for individual assay ####
        printColoredMessage(
            message = paste0('--Generating barplot of the ARI values for the individual data set(s):'),
            color = 'magenta',
            verbose = verbose
            )
        all.single.ari.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0(
                        '- Creating barplot of the ARI for the "',
                        x,
                        '" data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                ari.plot <- ggplot() +
                    geom_col(aes(y = all.ari[[x]], x = 1)) +
                    ggtitle(variables) +
                    xlab(x) +
                    ylab('Adjusted rand index (ARI) ') +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        plot.title = element_text(size = 16),
                        axis.text.x = element_text(size = 0),
                        axis.text.y = element_text(size = 12))
                if(isTRUE(plot.output) & length(assay.names) == 1) print(ari.plot)
                return(ari.plot)
            })
        names(all.single.ari.plots) <- levels(assay.names)

        ## putting all the plots of individual assays together ####
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = '-- Putting all the ARIs together:',
                color = 'magenta',
                verbose = verbose
                )
            overall.single.ari.plot <- as.data.frame(all.ari) %>%
                tidyr::pivot_longer(
                    everything(),
                    names_to = 'datasets',
                    values_to = 'ari'
                    )
            overall.single.ari.plot$datasets <- factor(
                x =  overall.single.ari.plot$datasets,
                levels = assay.names
                )
            overall.single.ari.plot <- ggplot(overall.single.ari.plot, aes(x = datasets, y = ari)) +
                geom_col() +
                ylab('Adjusted rand index (ARI)') +
                xlab('Datasets') +
                ggtitle(variables) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 16),
                    axis.title.y = element_text(size = 16),
                    plot.title = element_text(size = 18),
                    axis.text.x = element_text(
                        size = 14,
                        angle = 25,
                        hjust = 1),
                    axis.text.y = element_text(size = 12))
            overall.single.ari.plot <- annotate_figure(
                p = overall.single.ari.plot,
                top = text_grob(
                    label = "Adjusted rand index (ARI)",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ',
                        'Adjusted rand index between the first PCs and the ',
                        variables,
                        ' variable.'),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10)
            )
            printColoredMessage(
                message = '- The individual ARI barplots from each dataset have been combined into a single plot.',
                color = 'blue',
                verbose = verbose
                )
            if (isTRUE(plot.output))
                suppressMessages(print(overall.single.ari.plot))
        }
    }

    ## Combined plot ####
    if (plot.type == 'combined.plot') {
        printColoredMessage(
            message = '-- Obtain computed ARI for the from the SummarizedExperiment object:',
            color = 'magenta',
            verbose = verbose
        )
        all.ari <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0(
                        '- Obtaining ARI value for the ',
                        x,
                        ' data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                if (!ari.method %in% names(se.obj@metadata[['Metrics']][[x]][['global.level']][['ARI']])) {
                    stop(paste0('The ', ari.method ,'has not been computed yet for the ', variables, ' variable and the ', x, ' assay.'))
                }
                for (i in variables) {
                    if (!i %in% names(se.obj@metadata[['Metrics']][[x]][['global.level']][['ARI']][[ari.method]])) {
                        stop(paste0('The ', ari.method ,' has not been computed yet for the ', i, ' variable and the ', x, ' assay.'))
                    }
                }
                ari <- c()
                for (i in 1:length(variables))
                    ari[i] <- se.obj@metadata[['Metrics']][[x]][['global.level']][['ARI']][[ari.method]][[variables[i]]]$ari
                return(ari)
            })
        names(all.ari) <- levels(assay.names)
        ### Individual plots ####
        printColoredMessage(
            message = paste0('-- Plotting combined ARI'),
            color = 'magenta',
            verbose = verbose
            )
        datasets <- NULL
        all.combined.ari.plots <- lapply(
            levels(assay.names),
            function(x) {
                printColoredMessage(
                    message = paste0(
                        '- Plotting ARI for the ',
                        x,
                        ' data.'),
                    color = 'blue',
                    verbose = verbose
                    )
                all.aris <- as.data.frame(t(as.data.frame(all.ari[[x]])))
                row.names(all.aris) <- x
                colnames(all.aris) <- variables
                all.aris$datasets <- row.names(all.aris)
                p.combined <- ggplot(all.aris, aes_string(x = sym(variables[1]), y = sym(variables[2]))) +
                    geom_point() +
                    ggrepel::geom_text_repel(aes(label = datasets),
                                             hjust = 0,
                                             vjust = 0) +
                    xlab(paste0('Adjusted rand index (', variables[1], ')')) +
                    ylab(paste0('Adjusted rand index (', variables[2], ')')) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 18),
                        axis.title.y = element_text(size = 18),
                        plot.title = element_text(size = 15),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 12)
                    )
                if(isTRUE(plot.output) & length(assay.names) == 1) print(p.combined)
                return(p.combined)
            })
        names(all.combined.ari.plots) <- levels(assay.names)

        ## overall plots ####
        datasets <- NULL
        if (length(assay.names) > 1) {
            printColoredMessage(
                message = '-- Generating the combined ARI plot for all the data stes:',
                color = 'magenta',
                verbose = verbose
                )
            all.ari <- as.data.frame(t(as.data.frame(all.ari)))
            colnames(all.ari) <- variables
            all.ari$datasets <- row.names(all.ari)
            overall.combined.ari.plot <-
                ggplot(all.ari, aes_string(x = sym(variables[1]), y = sym(variables[2]))) +
                geom_point() +
                ggrepel::geom_text_repel(aes(label = datasets),
                                         hjust = 0,
                                         vjust = 0) +
                xlab(paste0('Adjusted rand index (', variables[1], ')')) +
                ylab(paste0('Adjusted rand index (', variables[2], ')')) +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 18),
                    axis.title.y = element_text(size = 18),
                    plot.title = element_text(size = 15),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12)
                )
            overall.combined.ari.plot <- annotate_figure(
                p = overall.combined.ari.plot,
                top = text_grob(
                    label = "Adjusted rand index",
                    color = "orange",
                    face = "bold",
                    size = 18),
                bottom = text_grob(
                    label = paste0(
                        'Analysis: ',
                        'Adjusted rand index between the first PCs and the ',
                        variables,
                        ' variable.'),
                    color = "black",
                    hjust = 1,
                    x = 1,
                    size = 10)
                )
        }
    }

    # Saving the results ####
    printColoredMessage(
        message = '-- Saving the ARI barplots:',
        color = 'magenta',
        verbose = verbose
        )
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '-- Save all the ARI barplots to the "metadata" in the SummarizedExperiment object:',
            color = 'blue',
            verbose = verbose
            )
        if(plot.type == 'single.plot'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = levels(assay.names),
                assessment.type = 'global.level',
                assessment = 'ARI',
                method = ari.method,
                variables = variables,
                file.name = 'single.plot',
                results.data = all.single.ari.plots
                )
        } else if (plot.type == 'combined.plot'){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                slot = 'Metrics',
                assay.names = levels(assay.names),
                assessment.type = 'global.level',
                assessment = 'ARI',
                method = ari.method,
                variables = paste0(variables, collapse = '&'),
                file.name = 'combined.plot',
                results.data = all.combined.ari.plots
                )
        }
        printColoredMessage(
            message = paste0(
                '- The ARI barplot of the individual assay(s) is saved to the ',
                ' "se.obj@metadata$metric$AssayName$ARI" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
            )
        if (length(assay.names) > 1) {
            if(plot.type == 'single.plot'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'global.level',
                    assessment = 'ARI',
                    method = ari.method,
                    variables = variables,
                    file.name = 'single.plot',
                    plot.data = overall.single.ari.plot
                    )
            } else if (plot.type == 'combined.plot'){
                se.obj <- addOverallPlotToSeObj(
                    se.obj = se.obj,
                    slot = 'Plots',
                    assessment.type = 'global.level',
                    assessment = 'ARI',
                    method = ari.method,
                    variables = paste0(variables, collapse = '&'),
                    file.name = 'combined.plot',
                    plot.data = overall.combined.ari.plot
                    )
            }
            printColoredMessage(
                message = paste0(
                    '- The combined ARI barplot all the assays is saved to the ',
                    ' "se.obj@metadata$plot$ARI" in the SummarizedExperiment object.'),
                color = 'blue',
                verbose = verbose)
        }
        printColoredMessage(message = '------------The plotARI function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj = se.obj)


    }
    ## return only the adjusted rand index results ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = paste0('- All the ARI plots re saved as list.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The plotARI function finished.',
            color = 'white',
            verbose = verbose
            )
        if (plot.type == 'single.plot') {
            if(length(assay.names) == 1){
                if(length(assay.names) > 1){
                    return( all.ari.plots = list(all.single.ari.plots = all.single.ari.plots))
                }else{
                    return( all.ari.plots = list(
                        all.single.ari.plots = all.single.ari.plots,
                        overall.single.ari.plot = overall.single.ari.plot))}
            }
        } else if (plot.type == 'combined.plot') {
            if(length(assay.names) == 1){
                return(all.ari.plots = list(all.combined.ari.plots = all.combined.ari.plots))
            }else{
                return(all.ari.plots = list(
                    all.combined.ari.plots = all.combined.ari.plots,
                    overall.combined.ari.plot = overall.combined.ari.plot))
            }
        }
    }
}
