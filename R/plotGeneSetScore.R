#' Plots the score of the gene set enrichment analysis.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function plots the score of the gene set scoring analysis performed by the `computeGeneSetScore()` function.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character or character vector that specifies the name(s) of the data (assays) in the
#' SummarizedExperiment object. The specified `gene.set.name` must be calculated by the `computeGeneSetScore()` function
#' before. The default is set to `all`, indicating all data sets in the SummarizedExperiment object will be selected.
#' @param gene.set.name Character. A character string indicating the name of the gene set stored in the metadata of the
#' SummarizedExperiment object. If `NULL`, the function will select the gene set name generated using the default parameters
#' of the `computeGeneSetScore()` function. See `computeGeneSetScore()` for more details.
#' @param reference.score Character. A character string specifying the name of the data set whose gene set score will be
#' used as the reference. All scores from other data sets will be plotted against this reference. The default is set to
#' `NULL`, indicating that scores from all data sets will be plotted against each other.
#' @param plot.output Logical. Whether to plot the gene set score or not. The default is set to `TRUE`.
#' @param save.se.obj Logical. Indicates whether to save the plots of scores in the metadata of the SummarizedExperiment
#' object or to output the result as a list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, shows the messages of different steps of the function.
#'
#' @importFrom ggpubr ggarrange stat_cor stat_compare_means
#' @importFrom SummarizedExperiment assays
#' @importFrom tidyr pivot_longer
#' @importFrom GGally ggpairs
#' @import ggplot2
#' @export

plotGeneSetScore <- function(
        se.obj,
        assay.names = 'all',
        reference.score = NULL,
        gene.set.name,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The plotGeneSetScore function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the function inputs ####
    if (!is.character(assay.names)){
        stop('The "assay.names" must be a character or a vector of characters of names of the data in the SummarizedExperiment object.')
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
    # Plotting the gene set scores ####
    ## Plotting the gene set scores against a reference score ####
    if (!is.null(reference.score)){
        all.assay.score <- lapply(
            levels(assay.names),
            function(x){
                se.obj@metadata$Metrics[[x]]$global.level$GeneSetScore$singscore[[gene.set.name]]$score
            })
        colnames(all.assay.score) <- levels(assay.names)
        all.assay.score <- as.data.frame(all.assay.score)
        if (!reference.score %in% levels(assay.names)){
            all.assay.score$ref <- se.obj[[reference.score]]
        } else if (reference.score %in% levels(assay.names)){
            index <- colnames(all.assay.score) %in% reference.score
            colnames(all.assay.score)[index] <- 'ref'
        }
        scores <- NULL
        all.assay.score <- pivot_longer(data = all.assay.score, -ref, names_to = 'datasets', values_to = 'scores')
        all.assay.score$datasets <- factor(x = all.assay.score$datasets, levels = assay.names)
        p.all.scores.plot <- ggplot(data = all.assay.score, aes(x = scores, y = ref)) +
            geom_point(color = 'grey') +
            ggtitle('Comapre the gene set scores to the reference score.') +
            xlab('Scores') +
            ylab('Reference scores')  +
            facet_wrap(~datasets) +
            geom_smooth(formula = y ~ x, method = 'lm', colour = "darkgreen") +
            ggpubr::stat_cor(
                aes(label = r.label),
                color = "navy") +
            theme(panel.background = element_blank(),
                  axis.line = element_line(colour = 'black', linewidth = 1),
                  axis.title.x = element_text(size = 12),
                  axis.title.y = element_text(size = 12),
                  axis.text.x = element_text(size = 9),
                  axis.text.y = element_text(size = 9)
            )
        if (isTRUE(plot.output)) print(p.all.scores.plot)
    }
    ## Plotting the gene set scores against each other s####
    if (is.null(reference.score)){
        all.assay.score <- lapply(
            levels(assay.names),
            function(x){
                se.obj@metadata$Metrics[[x]]$global.level$GeneSetScore$singscore[[gene.set.name]]$score
            })
        all.assay.score <- do.call(cbind, all.assay.score)
        colnames(all.assay.score) <- levels(assay.names)
        all.assay.score <- as.data.frame(all.assay.score)
        p.all.scores.plot <- ggpairs(data = all.assay.score) +
            theme(axis.line = element_line(colour = 'black', linewidth = 1),
                  axis.text.x = element_text(size = 6),
                  axis.text.y = element_text(size = 6))
        if (isTRUE(plot.output)) print(p.all.scores.plot)
    }
    # Saving the results ####
    ## Adding the plots to the metadata of the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- Saving all the score plots into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        se.obj <- addOverallPlotToSeObj(
            se.obj = se.obj,
            slot = 'Plots',
            assessment.type = 'global.level',
            assessment = 'GeneSetSocore',
            method = 'singscore',
            variables = gene.set.name,
            file.name = 'general',
            plot.data = p.all.scores.plot
        )
        printColoredMessage(
            message = '------------The plotGeneSetScore function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    printColoredMessage(message = '------------The plotGeneSetScore function finished.',
                        color = 'white',
                        verbose = verbose)
    ## Saving the plots as list ####
    if (isFALSE(save.se.obj)){
        return(list(p.all.scores.plot = p.all.scores.plot))
    }
}
