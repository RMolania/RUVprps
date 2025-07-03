#' Plots variables in a SummarizedExperiment object.
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
#' @param check.se.object TTTTT
#' @param remove.na TTTTT
#' @param verbose TTTTT

plotVariables <- function(
        se.obj,
        x.variables,
        y.variables,
        cat.test = 'anova',
        cont.test = 'spearman',
        generate.heatmap = FALSE,
        check.se.object = TRUE,
        remove.na = 'sample.annotation',
        verbose = TRUE
        ){
    p1 <- colData(se.obj)[ , c(x.variables, y.variables)] %>%
        data.frame(check.names = FALSE) %>%
        pivot_longer(-any_of(y.variables), names_to = 'cms', values_to = 'estimate') %>%
        mutate(cms = factor(cms , levels = x.variables)) %>%
        ggplot(., aes(x = estimate, y = Library.size)) +
        geom_boxplot() +
        facet_wrap(~ cms, ncol = 4) +
        stat_compare_means(method = 'anova', label.x  = 2, label.y = 26.3, color = 'navy') +
        xlab('') +
        ylab(expression(Log[2]~ 'library size')) +
        theme(
            panel.background = element_blank(),
            strip.text = element_text(size = 12),
            axis.line = element_line(colour = 'black', linewidth = 1),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
            axis.text.y = element_text(size = 10)
        )
    print(p1)
    if (isTRUE(generate.heatmap)){
        # Selecting colores ####
        selected.colores <-  c(
            c("#E7B800", "#2E9FDF", 'red4'),
            RColorBrewer::brewer.pal(8, "Dark2")[-5],
            RColorBrewer::brewer.pal(10, "Paired")
            )
        cms.labels <- as.matrix(droplevels(colData(se.obj)[ , c(x.variables)]))
        h.plot <- ComplexHeatmap::Heatmap(
            matrix = cms.labels,
            show_row_names = FALSE,
            column_names_rot = 35,
            column_names_gp = grid::gpar(fontsize = 7),
            col = selected.colores[1:6],  # use named vector
            heatmap_legend_param = list(
                at = 1:6,
                labels = c("CMS1", "CMS2", "CMS3", "CMS4", "Normal", "Not classified"),
                title = "CMS",
                color_bar = "discrete",
                legend_height = unit(4, "cm"))
        )
        print(h.plot)
    }
}








