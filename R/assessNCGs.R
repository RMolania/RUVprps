#' Assess performance of NCG.
#'
#' @author Ramyar Molania
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character or a vector of characters specifying the name(s) of the data (assays) in the
#' SummarizedExperiment object to be selected. These data stes will be log2-transformed with a pseudo count. The default
#' is se to `all`, which indicates that all data sets in the SummarizedExperiment object will be selected.
#' @param variables.to.assess.ncg TTT
#' @param ncg TTTT
#' @param apply.log TTTT
#' @param pseudo.count Numeric. A numeric value as pseudo count value to be added to all measurements in the selected
#' data sets before applying the log2 transformation, to avoid `-Inf` values for zero measurements. The default is set to 1.
#' @param nb.pcs TTT
#' @param svd.bsparam TTT
#' @param center TTTT
#' @param scale TTT
#' @param plot.output TTT
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. The default it is
#'  set to `TRUE`.
#' @param remove.na Character. A Character that indicates whether to remove NA or missing values from the data sets or
#' not. The options are: `assays` or `none`. The default is set to `assays`.  Refer to the `checkSeObj()` function for more
#' details.
#' @param verbose Logical. Indicates whether to display output messages during function execution. The default is set to
#' `TRUE`.
#'
#' @return The function returns a log2 transformed of all specified data sets as a list object.

assessNCGs <- function(
        se.obj,
        assay.name,
        variables.to.assess.ncg,
        ncg,
        apply.log = TRUE,
        pseudo.count = 1,
        nb.pcs,
        svd.bsparam = bsparam(),
        center = TRUE,
        scale = FALSE,
        plot.output = TRUE,
        check.se.obj = TRUE,
        remove.na = 'none',
        verbose = TRUE
        ){
    printColoredMessage(
        message = '-- Assessing the performance of selected NCG set:',
        color = 'magenta',
        verbose = verbose
        )
    printColoredMessage(
        message = '- Performing PCA using only the selected genes as NCGs.',
        color = 'blue',
        verbose = verbose
        )
    ### Applying log2 + pseudo count transformation ####
    if (isTRUE(apply.log)){
        expr.data <- applyLog(
            se.obj = se.obj,
            assay.names = assay.name,
            pseudo.count = pseudo.count,
            check.se.obj = check.se.obj,
            remove.na = remove.na,
            verbose = verbose
            )[[assay.name]]
    }
    if (isFALSE(apply.log)) {
        printColoredMessage(
            message = paste0(
                '- The ',
                assay.name,
                ' data will be used without any log transformation.'),
            color = 'blue',
            verbose = verbose
        )
        expr.data <- assay(x = se.obj, i = assay.name)
    }
    pca.data <- BiocSingular::runSVD(
        x = t(expr.data[ncg, ]),
        k = nb.pcs,
        BSPARAM = svd.bsparam,
        center = center,
        scale = scale
    )
    d <- pca.data$d
    n <- ncol(expr.data)  # number of samples (after transpose)
    pc.var <- (d^2) / (n - 1)
    centered.data <- scale(t(expr.data[ncg , ]), center = center, scale = scale)
    total.var <- sum(colVars(centered.data))
    percentage <- round(x = c(pc.var / total.var) * 100, digits = 2)
    ## regression and vector correlations ####
    printColoredMessage(
        message = paste0(
            '- Exploring the association of the first ',
            nb.pcs,
            '  PCs with the ',
            paste0(variables.to.assess.ncg, collapse = ' & '),
            ' variables.'),
        color = 'blue',
        verbose = verbose
    )
    all.corr <- lapply(
        variables.to.assess.ncg,
        function(x) {
            if (class(se.obj[[x]]) %in% c('numeric', 'integer')) {
                rSquared <- sapply(
                    1:nb.pcs,
                    function(y) summary(lm(se.obj[[x]] ~ pca.data$u[, 1:y]))$r.squared)
            } else if (class(se.obj[[x]]) %in% c('factor', 'character')) {
                catvar.dummies <- dummy_cols(se.obj[[x]])
                catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                cca.pcs <- sapply(
                    1:nb.pcs,
                    function(y) {
                        cca <- cancor(x = pca.data$u[, 1:y, drop = FALSE], y = catvar.dummies)
                        1 - prod(1 - cca$cor ^ 2)
                    })
            }
        })
    names(all.corr) <- variables.to.assess.ncg
    pca.ncg <- as.data.frame(do.call(cbind, all.corr))
    pca.ncg['pcs'] <- c(1:nb.pcs)
    pca.ncg <- tidyr::pivot_longer(
        data = pca.ncg, -pcs,
        names_to = 'Groups',
        values_to = 'ls'
    )
    assess.ncg.plot <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
        geom_line(aes(color = Groups), size = 1) +
        geom_point(aes(color = Groups), size = 2) +
        xlab('PCs') +
        ylab (expression("Correlations")) +
        ggtitle('Assessment of the NCGs') +
        scale_x_continuous(breaks = (1:nb.pcs),labels = c('PC1', paste0('PC1:', 2:nb.pcs))) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0, 1)) +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', linewidth = 1),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(
                size = 10,
                angle = 25,
                hjust = 1),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 10),
            plot.title = element_text(size = 16)
        )
    p.pca.percentage <- data.frame(var = percentage, no = 1:nb.pcs) %>%
        ggplot(., aes(x = no, y = percentage)) +
        geom_point(size = 3) +
        ylab('Variation(%)') +
        ylim(c(0,100)) +
        geom_line() +
        theme(
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black', linewidth = 1),
            axis.line.x  = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 14),
            plot.title = element_text(size = 16),
            plot.margin = unit(c(0, 0, 3, 0), "pt")
        )
    assess.ncg.plot <- assess.ncg.plot / p.pca.percentage + plot_layout(heights = c(3, 1))
    if (isTRUE(plot.output)) print(assess.ncg.plot)
    return(assess.ncg.plot)
}
