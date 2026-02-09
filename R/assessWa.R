#' Assesses the Wa matrix of RUV-III.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function assesses the association between the W matrix from RUV-III normalized data and known variables. It calculates
#' the correlation between the columns of the W matrix and both the specified biological and unwanted variables.

#' @details
#' The function performs linear regression between each individual continuous variable specified by users and the columns
#' of the W matrix in a cumulative manner, and then computes R^2 .For categorical variables, the function applies vector
#' correlation between the variables and the columns of the W matrix. Ideally, the columns of the W matrix should show a
#' strong correlation with unwanted variation and no or a weak correlation with biological variation. To compare the
#' performance of different $k$ values, the correlations obtained for biological variables are subtracted from 1, and then
#' the correlations for biological and unwanted variables are averaged separately. The final score for each set of $k$
#' values is sum of the half of the average scores. A higher score indicates a better value for $k$. It is important to
#' note that the assessment of W, should be used in combination with all other normalization assessments to select suitable
#' RUV-III normalized data.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param variables Character. A character or character vector of variables within the SummarizedExperiment object. These
#' can be categorical, continuous, or a combination of both. If oit is specified, a line-dot plots between correlations
#' and the columns of W , in cumulative way, will be generated. If set to `NULL`, both `bio.variables` and `bio.variables`
#' must be provided.
#' @param bio.variables Character. A character or character vector of biological variables within the SummarizedExperiment
#' object. These can be categorical, continuous, or a combination of both. The defeat is set to `NULL`.
#' @param uv.variables Character. A character or character vector of unwanted variables within the SummarizedExperiment
#' object. These can be categorical, continuous, or a combination of both. The defeat is set to `NULL`.
#' @param compare.wa Logical. If `TRUE` and both 'bio.variables' and 'uv.variables' are provided, the function generates
#' performance scores for each W matrix of RUV-III normalized data. The default is set to `FALSE`.
#' @param nb.pcs Numeric. A numeric value indicates the number of PCs to be used.
#' @param plot.output Logical. If `TRUE`, the function will show the output plots while is running. The default is set
#' to `TRUE`.
#' @param save.se.obj Logical. If `TRUE`, the plots will be saved to the metadata of the SummarizedExperiment object. The
#' default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, the function will display progress messages. The default is set to `TRUE`.
#'
#' @importFrom dplyr bind_rows summarise mutate group_by
#' @importFrom fastDummies dummy_cols
#' @importFrom tidyr pivot_longer
#' @importFrom gtools mixedorder
#' @export

assessWa <- function(
        se.obj ,
        variables,
        bio.variables = NULL,
        uv.variables = NULL,
        compare.wa = FALSE,
        nb.pcs = 10,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    # Checking the function inputs ####
    if (!is.null(variables)){
        if (!is.character(variables)){
            stop('The "variables" must be a charachter or a vector of charachters.')
        }
        if (!is.null(bio.variables) | !is.null(uv.variables)){
            stop('One of the "variables" or "bio.variables" and "uv.variables" must be specified.')
        }
        if (sum(variables %in% colnames(colData(se.obj))) != length(variables)){
            stop('All or some of the "variables" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (!is.null(bio.variables)){
        if (!is.null(variables)){
            stop('One of the "variables" or "bio.variables" and "uv.variables" must be specified.')
        }
        if (is.null(uv.variables)){
            stop('To generate a pefroamnce socre, both "bio.variables" and "uv.variables" must provided.')
        }
        if (!is.character(bio.variables)){
            stop('The "bio.variables" must be a charachter or a vector of charachters.')
        }
        if (!is.character(uv.variables)){
            stop('The "uv.variables" must be a charachter or a vector of charachters.')
        }
        if (sum(bio.variables %in% colnames(colData(se.obj))) != length(bio.variables)){
            stop('All or some of the "bio.variables" cannot be found in the SummarizedExperiment object.')
        }
        if (sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables)){
            stop('All or some of the "uv.variables" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (!is.null(uv.variables)){
        if (!is.null(variables)){
            stop('One of the "variables" or "bio.variables" and "uv.variables" must be specified.')
        }
        if (is.null(bio.variables)){
            stop('To generate a pefroamnce socre, both "bio.variables" and "uv.variables" must provided.')
        }
        if (!is.character(bio.variables)){
            stop('The "bio.variables" must be a charachter or a vector of charachters.')
        }
        if (!is.character(uv.variables)){
            stop('The "uv.variables" must be a charachter or a vector of charachters.')
        }
        if (sum(bio.variables %in% colnames(colData(se.obj))) != length(bio.variables)){
            stop('All or some of the "bio.variables" cannot be found in the SummarizedExperiment object.')
        }
        if (sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables)){
            stop('All or some of the "uv.variables" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (isTRUE(compare.wa)){
        if (is.null(uv.variables) | is.null(bio.variables)){
            stop('To compare the different W matrix, both the "bio.variables" and "uv.variables" must be provided.')
        }
    }
    # Comparing W  ####
    if (isTRUE(compare.wa)){
        printColoredMessage(
            message = '-- Comparing the performance of different W values:',
            color = 'magenta',
            verbose = verbose
            )
        ## Retrieving all the W matrix ####
        printColoredMessage(
            message = '- Retrieving all the W matrix of each RUV-III normalized data:',
            color = 'blue',
            verbose = verbose
            )
        data.names <- names(se.obj@metadata$RUVIII$Wa)
        names.order <- mixedorder(data.names)
        data.names <- data.names[names.order]
        all.wa <- lapply(
            data.names,
            function(x) se.obj@metadata$RUVIII$Wa[[x]]
            )
        names(all.wa) <- data.names

        # Applying PCA on WA data ####
        all.wa.pca <- lapply(
            data.names,
            function(x){
                pca.res <- irlba::svdr(
                    x = t(all.wa[[x]]),
                    k = nb.pcs,
                    center = TRUE
                    )
            })
        names(all.wa.pca) <- data.names
        # Finding the class of variables ####
        printColoredMessage(
            message = '- Finding the class of variables:',
            color = 'blue',
            verbose = verbose
            )
        all.vars <- c(bio.variables, uv.variables)
        class.all.vars <- sapply(
            all.vars,
            function(x) class(colData(x = se.obj)[[x]])
            )
        cat.vars <- names(class.all.vars[class.all.vars %in% c('factor', 'character')])
        cont.vars <- names(class.all.vars[class.all.vars %in% c('numeric', 'integer')])

        # Performing linear regression analysis ####
        if (length(cont.vars) > 0){
            cont.vars.r.squareds <- lapply(
                names(all.wa),
                function(x) {
                    pca.res <- all.wa.pca[[x]]$u
                    r.squareds <- sapply(
                        cont.vars,
                        function(y) {
                            lm.reg <- summary(lm(colData(x = se.obj)[[y]] ~ pca.res))
                            lm.reg$r.squared
                        })
                    names(r.squareds) <- cont.vars
                    r.squareds
                })
            names(cont.vars.r.squareds) <- names(all.wa)
            cont.vars.r.squareds <- as.data.frame(cont.vars.r.squareds)
        } else cont.vars.r.squareds <- NULL

        # Performing vector correlation analysis ####
        if (length(cat.vars) > 0){
            cat.vars.vec.corr <- lapply(
                names(all.wa),
                function(x) {
                    pca.res <- all.wa.pca[[x]]$u
                    vec.corr <- sapply(
                        cat.vars,
                        function(y) {
                            catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[y]])
                            catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                            cca <- cancor(x = pca.res, y = catvar.dummies)
                            1 - prod(1 - cca$cor ^ 2)
                        })
                    names(vec.corr) <- cat.vars
                    vec.corr
                })
            names(cat.vars.vec.corr) <- names(all.wa)
            cat.vars.vec.corr <- as.data.frame(cat.vars.vec.corr)
        } else cat.vars.vec.corr <- NULL

        ### Putting all together
        all <- bind_rows(cont.vars.r.squareds, cat.vars.vec.corr) %>%
            round(digits = 3)
        all.a <- mutate(.data = all, var = row.names(all)) %>%
            pivot_longer(cols = -var, values_to = 'corr', names_to = 'data') %>%
            mutate(data =  gsub('\\.', ':', data)) %>%
            mutate(data = factor(data, levels = names(all.wa)))

        p.w.1 <- ggplot(data = all.a, aes(x = data, y = corr, group = var)) +
            geom_line(aes(color = var), linewidth = 1) +
            geom_point(aes(color = var), size = 3) +
            scale_color_manual(
                values = selectColors(nb.color = 1:length(all.vars), group = "pan.selection.a"),
                name = 'Variables') +
            xlab('W (estimated unwanted factors)') +
            ylab('Correlations') +
            ylim(c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(size = 13, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 13),
                legend.title = element_text(size = 18),
                legend.text  = element_text(size = 14),
                legend.key.size = unit(1.5, "lines")
            )
        all <- mutate(.data = all, var = row.names(all)) %>%
            pivot_longer(cols = -var, values_to = 'corr', names_to = 'data') %>%
            mutate(groups = 'unwanted') %>%
            mutate(groups = ifelse(var %in% bio.variables, 'wanted', groups)) %>%
            mutate(corr = ifelse(groups == 'wanted', 1 - corr, corr)) %>%
            group_by(data, groups) %>%
            summarise(corr = mean(corr)) %>%
            summarise(assess = corr[groups == 'wanted']/2 + corr[groups == 'unwanted']/2) %>%
            mutate(data =  gsub('\\.', ':', data)) %>%
            mutate(data = factor(data, levels = names(all.wa)))
        p.w.2 <- ggplot(all, aes(x = data, y = assess)) +
            geom_bar(stat = 'identity', fill = 'grey') +
            xlab('W (estimated unwanted factors)') +
            ylab('Summarized correlations') +
            ylim(c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(size = 13, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 13)
            )
        p.w <- ggarrange(p.w.1, p.w.2 , ncol = 2)
        if (isTRUE(plot.output)) print(p.w)
    }
    #  Assessing W ####
    if (isFALSE(compare.wa)){
        printColoredMessage(
            message = '-- Assessing the performance of different W values:',
            color = 'magenta',
            verbose = verbose
            )
        ## Retrieving all the W matrix ####
        printColoredMessage(
            message = '- Retrieving all the W matrix of each RUV-III normalized data:',
            color = 'blue',
            verbose = verbose
            )
        data.names <- gsub('\\.', '_', names(se.obj@metadata$RUVIII$Wa))
        data.names <- mixedorder(data.names)
        all.wa <- lapply(
            data.names,
            function(x) se.obj@metadata$RUVIII$Wa[[x]]
            )
        names(all.wa) <- gsub('\\.', '_', names(se.obj@metadata$RUVIII$Wa))

        # Applying PCA on WA data ####
        all.wa.pca <- lapply(
            data.names,
            function(x){
                pca.res <- irlba::svdr(
                    x = t(all.wa[[x]]),
                    k = nb.pcs,
                    center = TRUE
                )
            })
        names(all.wa.pca) <- data.names

        # Finding the class of variables ####
        class.all.vars <- sapply(
            variables,
            function(x) class(colData(x = se.obj)[[x]])
            )
        cat.vars <- names(class.all.vars[class.all.vars %in% c('factor', 'character')])
        cont.vars <- names(class.all.vars[class.all.vars %in% c('numeric', 'integer')])

        # Peforming linear regression ####
        if (length(cont.vars) > 0){
            cont.vars.r.squareds <- lapply(
                names(all.wa),
                function(x) {
                    pca.res <- all.wa.pca[[x]]$u
                    r.squareds <- sapply(
                        cont.vars,
                        function(y) {
                            lm.reg <- summary(lm(colData(x = se.obj)[[y]] ~ pca.res))
                            lm.reg$r.squared
                        })
                    names(r.squareds) <- cont.vars
                    r.squareds
                })
            names(cont.vars.r.squareds) <- names(all.wa)
            cont.vars.r.squareds <- as.data.frame(cont.vars.r.squareds)
        } else cont.vars.r.squareds <- NULL

        # Performing vector correlation analysis ####
        if (length(cat.vars) > 0){
            cat.vars.vec.corr <- lapply(
                names(all.wa),
                function(x) {
                    pca.res <- all.wa.pca[[x]]$u
                    vec.corr <- sapply(
                        cat.vars,
                        function(y) {
                            catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[y]])
                            catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                            cca <- cancor(x = pca.res, y = catvar.dummies)
                            1 - prod(1 - cca$cor ^ 2)
                        })
                    names(vec.corr) <- cat.vars
                    vec.corr
                })
            names(cat.vars.vec.corr) <- names(all.wa)
            cat.vars.vec.corr <- as.data.frame(cat.vars.vec.corr)
        } else cat.vars.vec.corr <- NULL

        # Putting all together ####
        all.corrs <- bind_rows(cont.vars.r.squareds, cat.vars.vec.corr) %>%
            round(digits = 3)
        all.corrs$variable <- row.names(all.corrs)
        all.corrs <- pivot_longer(
            data = all.corrs,
            cols = -variable,
            values_to = 'corr',
            names_to = 'data') %>%
            mutate(data = factor(data, levels = gsub('\\.', '_', names(se.obj@metadata$RUVIII$Wa))))
        all.corrs <- data.frame(all.corrs)
        p.w <- ggplot(all.corrs, aes(x = data, y = corr, group = variable)) +
            geom_line(aes(color = variable), linewidth = 1) +
            geom_point(aes(color = variable), size = 3) +
            ylab('Correlations') +
            xlab('W') +
            theme(
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12)
            )
        if (isTRUE(plot.output)) print(p.w)
    }
    # Saving the results ####
    printColoredMessage(
        message = '-- Saving the results:',
        color = 'magenta',
        verbose = verbose
    )
    if (isTRUE(save.se.obj)){
        if (length(se.obj@metadata) == 0) {
            se.obj@metadata[['RUVIII']] <- list()
        }
        # check if RUVIII already exists in the metadata
        if (!'RUVIII' %in% names(se.obj@metadata)) {
            se.obj@metadata[['RUVIII']] <- list()
        }
        ## check if W already exists in the RUVIII
        if (!'CompareWa' %in% names(se.obj@metadata[['RUVIII']])) {
            se.obj@metadata[['RUVIII']][['CompareWa']] <- list()
        }
        se.obj@metadata[['RUVIII']][['CompareWa']] <- p.w
        return(se.obj)
    }
    if (isFALSE(save.se.obj)){
        return(p.w = p.w)
    }
}

