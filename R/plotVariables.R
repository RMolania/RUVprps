#' Plot variables in a SummarizedExperiment object.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function plots variables stored in a `SummarizedExperiment` object, including categorical and continuous
#' variables. It can generate scatter plots, boxplots, or heatmaps depending on the input and options selected.
#'
#' @param se.obj A `SummarizedExperiment` object containing the data and sample annotations to be plotted.
#' @param x.variables Character. The variable(s) from the sample annotations to be used as x-axis in the plot(s).
#' @param y.variables Character. The variable(s) from the sample annotations to be used as y-axis in the plot(s).
#' @param cat.test Character. The statistical test to apply when comparing categorical variables (e.g., `"chisq"`, `"fisher"`).
#' @param cont.test Character. The statistical test to apply when comparing continuous variables (e.g., `"t.test"`, `"wilcox"`).
#' @param generate.heatmap Logical. If `TRUE`, a heatmap of the specified variables will be generated in addition to plots.
#' The default is `FALSE`.
#' @param legend.title Character. A string specifying the title of the legend.
#' @param legend.labels Character vector. Labels to be used for categories in the legend.
#' @param nb.ncol Numeric. Number of columns to arrange plots in when multiple plots are generated. The default is 1.
#' @param x.lab Character. A string specifying the x-axis label.
#' @param y.lab Character. A string specifying the y-axis label.
#' @param plot.out.put Logical. If `TRUE`, the plots are displayed during function execution. If `FALSE`, the plots
#' are returned as a list without being displayed. The default is `TRUE`.
#' @param color.palette Character or function. Specifies the color palette to use for plots. Defaults to `ggplot2`
#' color scales if not provided.
#' @param check.se.object Logical. Indicates whether to validate the structure of the `SummarizedExperiment` object
#' before plotting. The default is `TRUE`.
#' @param remove.na Character. Specifies how to handle missing values. Options are: `"all"` (remove all rows with NAs),
#' `"pairwise"` (remove only for specific variable comparisons), or `"none"`. The default is `"all"`.
#' @param verbose Logical. If `TRUE`, messages describing progress and steps of the function will be displayed. The
#' default is `TRUE`.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom dplyr any_of mutate
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#'
#' @return A list of plots (scatter plots, boxplots, or heatmaps) depending on the input options.
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








