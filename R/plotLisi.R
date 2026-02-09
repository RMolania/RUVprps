#' Generate barplots the local inverse Simpson's index (lisi).
#'
#' @author Ramyar Molania
#'
#' @description
#' This function generates barplots of the local inverse Simpson's index (lisi) for individual dataset(s) in the
#' `SummarizedExperiment` object. Before applying the function, the `computeLisi` function should have been applied on the
#' specified dataset(s).
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character or character vector. Specifies the name(s) of the dataset(s) in the `SummarizedExperiment`
#' object to be used for generating barplots or scatter plots of the computed LISI By default, all dataset(s) in the
#' `SummarizedExperiment` object will be selected.
#' @param variable Character. A character specifying a column name in the sample annotation of the SummarizedExperiment
#' object that corresponds to the variable for which LISI has been computed.
#' @param fast.pca Character. A character specifying which PCA approach has been used to compute LISI The default is set
#' to `fast.pca`. Refer to the `computeLisi` function for more information.
#' @param plot.output Logical. If `TRUE`, the individual barplots or scatter plots will be printed during function execution.
#' The default is set to `TRUE`.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the `SummarizedExperiment` object
#' or return them as a list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, messages for different steps of the function will be displayed.
#'
#' @return A `SummarizedExperiment` object or a list containing all the plots of the computed LISI scores for the
#' categorical variable(s).
#'
#' @importFrom SummarizedExperiment assays assay
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#'
#' @export

plotLisi <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The plotLisi function starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking the inputs ####
    if (is.null(assay.names) | is.logical(assay.names)) {
        stop('The "assay.names" cannot be NUll or logical.')
    }
    if (is.null(variable) | is.logical(variable)) {
        stop('The "variable" cannot be NULL or logical.')
    }
    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    if (isTRUE(fast.pca)){
        method = 'fast.svd'
    } else method = 'svd'
    # Plotting the ARI values ####
    ### obtain LISI ####
    printColoredMessage(
        message = paste0('-- Obtaining computed LISI from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
        )
    all.lisi <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = assay.names,
        assessment = 'LISI',
        assessment.type = 'global.level',
        method = method,
        variables = variable,
        file.name = 'lisi',
        sub.file.name = NULL,
        required.function = 'compuetLisi',
        message.to.print = 'LISI'
    )
    ## plot for individual assay ####
    printColoredMessage(
        message = paste0('-- Generating barplot of the LISI values for the individual data set(s):'),
        color = 'magenta',
        verbose = verbose
    )
    all.single.lisi.plots <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0(
                    '- Creating boxrplot of the LISI for the "',
                    x,
                    '" data.'),
                color = 'blue',
                verbose = verbose
            )
            ari.plot <- ggplot() +
                geom_boxplot(aes(y = all.lisi[[x]], x = 1)) +
                ggtitle(variable) +
                xlab(x) +
                ylab("Local inverse Simpson's index") +
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
    names(all.single.lisi.plots) <- levels(assay.names)

    ## putting all the plots of individual assays together ####
    if (length(assay.names) > 1) {
        printColoredMessage(
            message = '-- Putting all the LISI together:',
            color = 'magenta',
            verbose = verbose
        )
        overall.single.lisi.plot <- as.data.frame(all.lisi) %>%
            tidyr::pivot_longer(
                everything(),
                names_to = 'datasets',
                values_to = 'lisi'
            )
        overall.single.lisi.plot$datasets <- factor(
            x =  overall.single.lisi.plot$datasets,
            levels = assay.names
        )
        overall.single.lisi.plot <- ggplot(overall.single.lisi.plot, aes(x = datasets, y = lisi)) +
            geom_boxplot() +
            ylab("Local inverse Simpson's index") +
            xlab('Datasets') +
            ggtitle(variable) +
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
        overall.single.lisi.plot <- annotate_figure(
            p = overall.single.lisi.plot,
            top = text_grob(
                label = "Local inverse Simpson's index (LISI)",
                color = "orange",
                face = "bold",
                size = 18),
            bottom = text_grob(
                label = paste0(
                    'Analysis: ',
                    "Local inverse Simpson's index between the first PCs and the ",
                    variable,
                    ' variable.'),
                color = "black",
                hjust = 1,
                x = 1,
                size = 10)
        )
        printColoredMessage(
            message = '- The individual LISI boxplots from each dataset have been combined into a single plot.',
            color = 'blue',
            verbose = verbose
        )
        if (isTRUE(plot.output))
            suppressMessages(print(overall.single.lisi.plot))
    }
    # Saving the results ####
    printColoredMessage(
        message = '-- Saving the LISI barplots:',
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
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'global.level',
            assessment = 'LISI',
            method = method,
            variables = variable,
            file.name = 'boxplot',
            results.data = all.single.lisi.plots
            )
        printColoredMessage(
            message = paste0(
                '- The LISI barplot of the individual assay(s) is saved to the ',
                ' "se.obj@metadata$metric$AssayName$ARI" in the SummarizedExperiment object.'),
            color = 'blue',
            verbose = verbose
        )
        if (length(assay.names) > 1) {
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'global.level',
                assessment = 'LISI',
                method = method,
                variables = variable,
                file.name = 'boxplot',
                plot.data = overall.single.lisi.plot
                )
        }
        printColoredMessage(
            message = '------------The plotLisi function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj = se.obj)
    }
    ## return only the adjusted rand index results ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = paste0('- All the LISI plots re saved as list.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The plotLisi function finished.',
            color = 'white',
            verbose = verbose
            )
        if(length(assay.names) == 1){
            return(all.lisi.plots = list(all.single.lisi.plots = all.single.lisi.plots))
        }
        if (length(assay.names) > 1){
            return(all.lisi.plots = list(
                all.single.lisi.plots = all.single.lisi.plots,
                overall.single.lisi.plot = overall.single.lisi.plot
                ))
        }
    }
}
