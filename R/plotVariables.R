#' Plot variables in a SummarizedExperiment object.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function plots any variables in a SummarizedExperiment object
#'
#' @param se.obj TTT
#' @param x.variables TTT
#' @param y.variables TTTT
#' @param cat.test TTTT
#' @param cont.test TTTT
#' @param generate.heatmap TTTT
#' @param legend.title TTTT
#' @param legend.labels TTTT
#' @param nb.ncol TTTT
#' @param x.lab TTTT
#' @param y.lab TTTTT
#' @param plot.out.put TTTTT
#' @param color.palette TTTTT
#' @param check.se.object TTTTT
#' @param remove.na TTTTT
#' @param verbose TTTTT
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom dplyr any_of mutate
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#'
#' @export

plotVariables <- function(
        se.obj,
        x.variables,
        y.variables,
        cat.test = 'anova',
        cont.test = 'spearman',
        generate.heatmap = FALSE,
        legend.title = NULL,
        legend.labels = NULL,
        nb.ncol = 4,
        x.lab = NULL,
        y.lab = NULL,
        plot.out.put = TRUE,
        color.palette = 'nrc',
        check.se.object = FALSE,
        remove.na = 'sample.annotation',
        verbose = TRUE
        ){
    # Generating box plots ####
    if (!is.null(x.variables) & !is.null(y.variables)){
        p.x.y <- colData(se.obj)[, c(x.variables, y.variables)] %>%
            data.frame(check.names = FALSE) %>%
            pivot_longer(-any_of(y.variables), names_to = 'variables', values_to = 'values') %>%
            mutate(
                variables = factor(variables, levels = x.variables),
                y_value = .data[[y.variables]] ) %>%
            ggplot(aes(x = values, y = y_value)) +
            geom_boxplot() +  # or any appropriate geom
            facet_wrap(~ variables, scales = "free_x") +
            geom_boxplot() +
            facet_wrap(~ variables, ncol = nb.ncol) +
            stat_compare_means(
                method = 'anova',
                label = "p.format",
                color = 'navy',
                label.y.npc = "top",
                label.x.npc = "center") +
            xlab(x.lab) +
            ylab(y.lab) +
            theme(
                panel.background = element_blank(),
                strip.text = element_text(size = 12),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 14)
            )
        if (isTRUE(plot.out.put)) print(p.x.y)
    } else p.x.y <- NULL
    # Generating heatmap ####
    if (isTRUE(generate.heatmap) & !is.null(x.variables)){
        if (is.null(legend.labels)){
            legend.labels <- unique(unlist(lapply(droplevels(colData(se.obj)[ , c(x.variables)]), as.character)))
            legend.labels <- sort(legend.labels)
        }
        variable.matrix <- as.matrix(droplevels(colData(se.obj)[ , c(x.variables)]))
        selected.colores <- selectColors(
            nb.color = 1:length(unique(as.vector(legend.labels))) ,
            group = color.palette
            )
        h.plot <- ComplexHeatmap::Heatmap(
            matrix = variable.matrix,
            show_row_names = FALSE,
            column_names_rot = 35,
            column_names_gp = grid::gpar(fontsize = 7),
            col = selected.colores,
            heatmap_legend_param = list(
                at = 1:length(unique(as.vector(variable.matrix))),
                labels = legend.labels,
                title = legend.title,
                color_bar = "discrete",
                legend_height = unit(4, "cm"))
        )
        if (isTRUE(plot.out.put)) print(h.plot)
    } else h.plot <- NULL
    return(list(p.x.y = p.x.y, h.plot = h.plot))
}








