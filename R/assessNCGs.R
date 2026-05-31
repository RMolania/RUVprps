#' Assess the performance of NCG.
#'
#' @author Ramyar Molania
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character specifying the name of the data (assay) in the SummarizedExperiment object
#' to be selected. This data be will be used to assess the performance of NCGs. The data should the same one used for the
#' identification of the NCG.
#' @param variables.to.assess.ncg Character or character vector. Specifies the column name(s) in the sample annotation
#' of the SummarizedExperiment object corresponding to variables for which NCG performance should be assessed.
##' @param ncg Logical. A logical vector specifying a set of NCG whose performance will be assessed.
#' @param apply.log Logical. Indicates whether to apply log2 transformation to the selected data set. The default is set
#' to `TRUE`.
#' @param pseudo.count Numeric. A numeric value to be added to all measurements in the selected
#' data sets before applying the log2 transformation to avoid `-Inf` values for zero measurements. The default is set to 1.
#' @param nb.pcs Numeric. Specifies the number of principal components to retain when performing PCA on the data. The PCs
#' will be used for the assessment. The default is set to 3.
#' @param svd.bsparam Character. A BiocParallelParam object specifying how palatalization should be performed for performing
#' PCA using SVD. The default is set to `bsparam()`. We refer to the `runSVD()` function from the **BiocSingular** R package
#' for further details.
#' @param center Logical. Indicates whether to center the data or not before calculating PCs. If center is `TRUE`, then
#' centering is done by subtracting the column means of the assay from their corresponding columns. The default is set
#' to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data or not before calculating PCs. If scale is set to `TRUE`, then
#' scaling is done by dividing the (centered) columns of the assays by their standard deviations if center is `TRUE`, and
#' the root mean square otherwise. The default is set to `FALSE`.
#' @param plot.output Logical. Indicates whether to generate and display performance plots of NCG. Default is `TRUE`.
#' @param check.se.obj Logical. Indicates whether to validate the SummarizedExperiment object before analysis. The default
#' is set to `TRUE`.
#' @param remove.na Character. Indicates whether to remove NA or missing values from the data sets. Options are: `assays`
#' or `none`. Default is `assays`. Refer to the `checkSeObj()` function for more details.
#' @param verbose Logical. Indicates whether to display messages during function execution. The default is set `TRUE`.
#'
#' @return A line-dot plot show the association between the first `nb.pcs` PCs and the variables specified by `variables.to.assess.ncg`.
#'
#' @export

assessNCGs <- function(
        se.obj,
        assay.name,
        variables.to.assess.ncg,
        ncg,
        apply.log = TRUE,
        pseudo.count = 1,
        nb.pcs = 3,
        svd.bsparam = bsparam(),
        center = TRUE,
        scale = FALSE,
        plot.output = TRUE,
        check.se.obj = TRUE,
        remove.na = 'none',
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The assessNCGs function starts:',
        color = 'white',
        verbose = verbose
        )
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
    printColoredMessage(
        message = '------------The assessNCGs function finished.',
        color = 'white',
        verbose = verbose
        )
    return(assess.ncg.plot)
}
