#' Generate boxplots for the k-nearest neighbor batch effect test (kBET).
#'
#' @author Ramyar Molania
#'
#' @description
#' This function generates barplots of the k-nearest neighbor batch effect test (kBET) for individual dataset(s) in the
#' `SummarizedExperiment` object. Before applying the function, the `computeKbet` function should have been applied on the
#' specified dataset(s).
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character or character vector. Specifies the name(s) of the dataset(s) in the `SummarizedExperiment`
#' object to be used for generating barplots or scatter plots of the computed kBET. By default, all dataset(s) in the
#' `SummarizedExperiment` object will be selected.
#' @param variable Character. A character specifying a column name in the sample annotation of the SummarizedExperiment
#' object that corresponds to the variable for which kBET has been computed.
#' @param fast.pca Character. A character specifying which PCA approach has been used to compute kBET. The default is set
#' to `fast.pca`. Refer to the `computeKbet` function for more information.
#' @param plot.output Logical. If `TRUE`, the individual barplots or scatter plots will be printed during function execution.
#' The default is set to `TRUE`.
#' @param save.se.obj Logical. Indicates whether to save the plots in the metadata of the `SummarizedExperiment` object
#' or return them as a list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, messages for different steps of the function will be displayed.
#'
#' @return A `SummarizedExperiment` object or a list containing all the plots of the computed kBET scores for the
#' categorical variable(s).
#'
#' @importFrom SummarizedExperiment assays assay
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

plotKbet <- function(
        se.obj,
        assay.names = 'all',
        variable,
        fast.pca = TRUE,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotKbet function starts:',
                        color = 'white',
                        verbose = verbose)
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
        message = paste0('-- Obtaining computed kBET from the SummarizedExperiment object:'),
        color = 'magenta',
        verbose = verbose
        )
    all.kbet <- getMetricFromSeObj(
        se.obj = se.obj,
        slot = 'Metrics',
        assay.names = assay.names,
        assessment = 'KBET',
        assessment.type = 'global.level',
        method = method,
        variables = variable,
        file.name = 'kBET',
        sub.file.name = NULL,
        required.function = 'compuetKbet',
        message.to.print = 'kBET'
        )
    ## plot for individual assay ####
    printColoredMessage(
        message = paste0('-- Generating boxplot of the kBET values for the individual data set(s):'),
        color = 'magenta',
        verbose = verbose
        )
    all.single.kbet.plots <- lapply(
        levels(assay.names),
        function(x) {
            printColoredMessage(
                message = paste0(
                    '- Creating boxrplot of the kBET for the "',
                    x,
                    '" data.'),
                color = 'blue',
                verbose = verbose
                )
            kept.plot <- ggplot() +
                geom_boxplot(aes(y = all.kbet[[x]], x = 1)) +
                ggtitle(variable) +
                xlab(x) +
                ylab(" k-nearest neighbour batch effect (kBET)") +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    plot.title = element_text(size = 16),
                    axis.text.x = element_text(size = 0),
                    axis.text.y = element_text(size = 12))
            if(isTRUE(plot.output) & length(assay.names) == 1) print(kept.plot)
            return(kept.plot)
        })
    names(all.single.kbet.plots) <- levels(assay.names)

    ## putting all the plots of individual assays together ####
    if (length(assay.names) > 1) {
        printColoredMessage(
            message = '-- Putting all the kEBT together:',
            color = 'magenta',
            verbose = verbose
            )
        overall.single.kbet.plot <- as.data.frame(all.kbet) %>%
            tidyr::pivot_longer(
                everything(),
                names_to = 'datasets',
                values_to = 'kebt'
            )
        overall.single.kbet.plot$datasets <- factor(
            x =  overall.single.kbet.plot$datasets,
            levels = assay.names
            )
        overall.single.kbet.plot <- ggplot(overall.single.kbet.plot, aes(x = datasets, y = kebt)) +
            geom_boxplot() +
            ylab("k-nearest neighbour batch effect test (kBET)") +
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
                axis.text.y = element_text(size = 12)
                )
        overall.single.kbet.plot <- annotate_figure(
            p = overall.single.kbet.plot,
            top = text_grob(
                label = "k-nearest neighbour batch effect test (kBET)",
                color = "black",
                face = "bold",
                size = 18),
            bottom = text_grob(
                label = paste0(
                    'Analysis: ',
                    "k-nearest neighbour batch effect test (kBET) between the first PCs and the ",
                    variable,
                    ' variable.'),
                color = "black",
                hjust = 1,
                x = 1,
                size = 10)
            )
        printColoredMessage(
            message = '- The individual boxplots of the kBET scores from each dataset have been combined into a single plot.',
            color = 'blue',
            verbose = verbose
        )
        if (isTRUE(plot.output))
            suppressMessages(print(overall.single.kbet.plot))
    }
    # Saving the results ####
    printColoredMessage(
        message = '-- Saving the kBET barplots:',
        color = 'magenta',
        verbose = verbose
        )
    ## add results to the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = '-- Save all the kBET boxplots to the "metadata" in the SummarizedExperiment object:',
            color = 'blue',
            verbose = verbose
            )
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'global.level',
            assessment = 'KBET',
            method = method,
            variables = variable,
            file.name = 'boxplot',
            results.data = all.single.kbet.plots
            )
        if (length(assay.names) > 1) {
            se.obj <- addOverallPlotToSeObj(
                se.obj = se.obj,
                slot = 'Plots',
                assessment.type = 'global.level',
                assessment = 'KBET',
                method = method,
                variables = variable,
                file.name = 'boxplot',
                plot.data = overall.single.kbet.plot
                )
        }
        printColoredMessage(
            message = '------------The plotKbet function finished.',
            color = 'white',
            verbose = verbose
        )
        return(se.obj = se.obj)
    }
    ## return only the adjusted rand index results ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(
            message = paste0('- All the plots of the kBET scores are outputed as a list.'),
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The plotKbet function finished.',
            color = 'white',
            verbose = verbose
            )
        if(length(assay.names) == 1){
            return(all.kbet.plots = list(all.single.kbet.plots = all.single.kbet.plots))
        }
        if (length(assay.names) > 1){
            return(all.kbet.plots = list(
                all.single.kbet.plots = all.single.kbet.plots,
                overall.single.kbet.plot = overall.single.kbet.plot
            ))
        }
    }
}
