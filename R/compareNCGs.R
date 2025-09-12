#' Assesses and compares the performance of NCGs gene sets.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function globally assesses and compares the performance of various sets of negative control genes (NCGs) based on
#' their ability to capture unwanted variation while remaining unassociated with biological variation. This is essential
#' to select a suitable set of NCGs for RUV-III normalization.
#'
#' @details
#' For each NCG set, the function performs PCA on the specified data using only the genes in the set. It then evaluates
#' the association of the top `nb.pcs` principal components with variables as defined by `variables`. Linear regression
#' is used for continuous variables, and vector correlation is used for categorical variables. For each NCG
#' set, $R^2$ values (from regression) and vector correlation coefficients are computed to quantify these associations.
#' The $R^2$ and vector correlation values calculated for individual variables are plotted against the cumulative PCs.
#' Ideally, a suitable set of NCGs should show high and low correlation values, respectively, for the unwanted and biological
#' variables.\
#' If both `bio.variables` and `ub.variables` are specified, then the function generates a performance score between 0 and
#' 1 is for NCG each set, reflecting how effectively the NCGs capture unwanted variation while minimizing association
#' with biological factors. To compute such scores, the correlations obtained for biological variables are subtracted
#' from 1, and then the correlations for biological and unwanted variables are averaged, separately. The final score for
#' each NCG set is calculated as the sum of half the average scores. A higher score indicates a better set of NCGs.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character that specifies the name of the data in the SummarizedExperiment object. This
#' data should one that will be used as the input data for RUV-III-PRPS normalization.
#' @param variables Character. One or more variable names in the SummarizedExperiment object. These will be used to assess
#' the performance of the NCGs. If this is set to `NULL`, both `bio.variables` and `uv.variables` must be specified.
#' @param bio.variables Character. One or more variable names representing biological variables (e.g., cancer subtypes,
#' tumor purity) within the SummarizedExperiment object. This can be a vector of categorical, continuous, or mixed variable
#' types. The default is set to `NULL`. This must be specified if the performance scores is need to be calculated.
#' @param uv.variables Character. One or more variable names representing unwanted variables (e.g., batch effects,
#' library size) within the SummarizedExperiment object. This can also be a vector of categorical, continuous, or mixed
#' variable types. The default is set to `NULL`. This must be specified if the performance scores is need to be calculated.
#' @param ncg.type Character. A character that specifies the type of NCGs to be selected and assessed their performnce.
#' Options are `'supervised'`, `'un.supervised'`, `pre.selected` and `all.genes`.
#' @param  ncg.group.name Character. A character that specifies the name of the group of NCGs stored under the 'supervised' or
#' 'un.supervised' slots in the SummarizedExperiment object. Refer to supervised and unsupervised functions for finding
#' NCGS for more details.
#' @param ncg.set.names Character. Specifies the exact name of an NCG set under the given ` ncg.group.name` in either the
#' 'supervised' or 'un.supervised' slot or in the `pre.selected` slot. The default is se to `'all'`, which selects all
#' available NCG sets.
#' @param nb.pcs Numeric. A numeric value that indicates the number of principal components (PCs) to compute and use for
#' the assessment analysis. The default is set to `10`.
#' @param center Logical. Indicates whether to center the data prior to applying PCA suing SVD. If `TRUE`, centering is
#' done by subtracting the column means from the respective columns. The default is set to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data prior to applying PCA using SVD. If `TRUE`, the centered
#' columns are divided by their standard deviations (if `center = TRUE`), or by the root mean square
#' (if `center = FALSE`). The default is set to `FALSE`.
#' @param apply.log Logical. Indicates whether to apply a log transformation to the data before computing the SVD. The
#' default is set to `TRUE`. Data should be in log scale before applying SVD.
#' @param pseudo.count Numeric. A numeric value as a pseudo count to be added to all measurements before applying the
#' log transformation. The default is set to `1`.
#' @param plot.output Logical. Indicates whether to generate the assessment plot during function execution. The default
#' is set to `TRUE`.
#' @param output.name Character. Specifies the name under which the assessment plot will be saved in the
#' SummarizedExperiment object. If `NULL`, a name will be generated using:
#' `paste0('comparison_', paste0(names(ncgs), collapse = '&'))`.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not before applying
#' nationalizations. The default is set  to `TRUE`. We refer to the `checkSeObj()` function for more details.
#' @param remove.na Character. A character that indicates whether to remove NA or missing values from the data sets or
#' not. The options are `'assays'` or `'none'`. The default is set to `'assays'`. Refer to the `checkSeObj()` function
#' for more details.
#' @param save.se.obj Logical. Indicates whether to save the assessment plot metadata of the SummarizedExperiment object.
#' The default is set to `TRUE`. The plot will be saved to:
#' @param verbose Logical. Indicates whether to display output messages during function execution. The default is set
#' to `TRUE`.
#'
#' @return Either a SummarizedExperiment object containing the NCG assessment plots, or out theses results as list if
#' `save.se.obj` is set to `FALSE`.
#'
#' @importFrom dplyr bind_rows summarise mutate group_by
#' @importFrom fastDummies dummy_cols
#' @importFrom gtools mixedorder
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @export

compareNCGs <- function(
        se.obj ,
        assay.name,
        variables,
        bio.variables = NULL,
        uv.variables = NULL,
        ncg.type,
        ncg.group.name,
        ncg.set.names = 'all',
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        apply.log = TRUE,
        pseudo.count = 1,
        plot.output = TRUE,
        output.name = NULL,
        check.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    # Checking the function inputs ####
    if (is.logical(assay.name) | is.null(assay.name)){
        stop ('The "assay.name" cannot be NULL or logical.')
    }
    if (length(assay.name) > 1 | assay.name == 'all') {
        stop('The "assay.name" must be a name of a data (assay) in the SummarizedExperiment object.')
    }
    if (!is.null(bio.variables) & !is.null(uv.variables)){
        variables <- NULL
    }
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
    if (is.null(ncg.type) | is.null( ncg.group.name)){
        stop('The "ncg.type" and " ncg.group.name" must be provided.')
    }
    if (!ncg.type %in% c('supervised', 'un.supervised', 'pre.selected', 'all.genes')){
        stop('The "ncg.type" must be one of the "supervised", "un.supervised", "pre.selected" or "all.genes".')
    }
    if (!is.character(ncg.type) | length(ncg.type) > 1){
        stop('The "ncg.type" must be a single character.')
    }
    if (!is.character( ncg.group.name) | length( ncg.group.name) > 1){
        stop('The "ncg.type" must be a single character.')
    }
    if (!is.character(ncg.set.names)){
        stop('The "ncg.set.names" must be a character.')
    }
    if (!is.numeric(nb.pcs) | nb.pcs < 0){
        stop('The "nb.pcs" must be a psotive numeric value.')
    }
    if (!is.logical(center) | !is.logical(scale)){
        stop('Both "center" and "scale" must be logical.')
    }
    if (!is.logical(apply.log)){
        stop('The "apply.log" and "scale" must be logical.')
    }
    if (!is.numeric(pseudo.count) | pseudo.count < 0){
        stop('The "pseudo.count" must be a psotive numeric value.')
    }
    if (!is.logical(plot.output)){
        stop('The "plot.output" and "scale" must be logical.')
    }
    if (!is.null(output.name)){
        if (!is.character(output.name)){
            stop('The "output.name" must be character.')
        }
        if (length(output.name) > 1){
            stop('The "output.name" must be a single character.')
        }
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical.')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical.')
    }

    # Checking SummarizedExperiment object ####
    if (isTRUE(check.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(variables, bio.variables, uv.variables),
            remove.na = remove.na,
            verbose = verbose
        )
    }

    # Obtaining the NCGs sets ####
    printColoredMessage(
        message = '-- Obtaining the specified NCG sets:',
        color = 'magenta',
        verbose = verbose
        )
    ## Supervised NCGs selection ####
    if (ncg.type == 'supervised'){
        printColoredMessage(
            message = '- Obtaining the specified NCG sets found by the supervised approaches:',
            color = 'blue',
            verbose = verbose
            )
        if (sum(ncg.set.names == 'all') == 1){
            ncg.set.names <- names(se.obj@metadata$NCG$supervised[[ncg.group.name]])
            ncg.set.names <- ncg.set.names[!ncg.set.names %in% 'assessment.plot']
            ncgs <- lapply(
                ncg.set.names,
                function(x){
                    se.obj@metadata$NCG$supervised[[ncg.group.name]][[x]]$ncg.set
                })
            names(ncgs) <- ncg.set.names
        } else {
            if (sum(ncg.set.names %in% names(se.obj@metadata$NCG$supervised[[ncg.group.name]])) != length(ncg.set.names)){
                stop('All or some of the "ncg.set.names" cannot be found in the metadata of the SummarizedExperiment obejct.')
            } else {
                ncgs <- lapply(
                    ncg.set.names,
                    function(x){
                        se.obj@metadata$NCG$supervised[[ncg.group.name]][[x]]$ncg.set
                    })
                names(ncgs) <- ncg.set.names
            }
        }
        printColoredMessage(
            message = paste0(
                '- The total number of ',
                length(ncg.set.names),
                ' NCG sets are obtained:'),
            color = 'blue',
            verbose = verbose
            )
        print(sapply(ncgs, sum))

    }
    ## Un supervised NCGs selection ####
    if (ncg.type == 'un.supervised'){
        printColoredMessage(
            message = '- Obtaining the specified NCG sets found by the un-supervised approaches:',
            color = 'blue',
            verbose = verbose
            )
        if (ncg.set.names == 'all'){
            ncg.set.names <- names(se.obj@metadata$NCG$un.supervised[[ncg.group.name]])
            ncg.set.names <- ncg.set.names[!ncg.set.names %in% 'assessment.plot']
            ncgs <- lapply(
                ncg.set.names,
                function(x){
                    se.obj@metadata$NCG$un.supervised[[ncg.group.name]][[x]]$ncg.set
                })
            names(ncgs) <- ncg.set.names

        } else {
            if (sum(ncg.set.names %in% names(se.obj@metadata$NCG$un.supervised[[ ncg.group.name]]$ncg.set)) != length(ncg.set.names) ){
                stop('All or some of the "ncg.set.names" cannot be found in the SummarizedExperiment.')
            } else{
                ncgs <- lapply(
                    ncg.set.names,
                    function(x){
                        se.obj@metadata$NCG$un.supervised[[ncg.group.name]][[x]]$ncg.set
                    })
                names(ncgs) <- ncg.set.names
            }
        }
        printColoredMessage(
            message = paste0(
                '- The total number of ',
                length(ncg.set.names),
                ' NCG sets are obtained:'),
            color = 'blue',
            verbose = verbose
        )
        print(sapply(ncgs, sum))
    }
    ## Pre selected sets ####
    if (ncg.type == 'pre.selected'){
        printColoredMessage(
            message = '- Obtaining specified NCG sets provided as pre selection:',
            color = 'blue',
            verbose = verbose
            )
        if (ncg.set.names == 'all'){
            ncgs <- se.obj@metadata[['NCG']][['pre.selected']]
        } else {
            if (sum(ncg.set.names %in% names(se.obj@metadata[['NCG']][['pre.selected']])) != length(ncg.set.names) ){
                stop('All or some of the "ncg.set.names" cannot be found in the metadata of the SummarizedExperiment object.')
            } else{
                ncgs <- lapply(
                    ncg.set.names,
                    function(x){
                        se.obj@metadata$NCG$un.supervised[[ ncg.group.name]]$ncg.set[[x]]
                    })
                names(ncgs) <- ncg.set.names
            }
        }
    }
    ## All genes ####
    if (ncg.type == 'all.genes'){
        printColoredMessage(
            message = '- Using all genes as NCGs:',
            color = 'blue',
            verbose = verbose
        )
        ncgs <- list(all.genes = row.names(se.obj))
    }
    # Comparing and assessing the performance of NCGs ####
    printColoredMessage(
        message = '-- Assessing the performance of the specified NCG sets:',
        color = 'magenta',
        verbose = verbose
        )
    printColoredMessage(
        message = '- Performing PCA using only selected NCGs:',
        color = 'blue',
        verbose = verbose
        )
    ## Performing PCA on NCGs ####
    all.pcs.on.ncgs <- lapply(
        names(ncgs),
        function(x){
            pcs <- computePCA(
                se.obj = se.obj[ncgs[[x]] , ],
                assay.names = assay.name,
                fast.pca = TRUE,
                nb.pcs = nb.pcs,
                center = center,
                scale = scale,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                svd.bsparam = bsparam(),
                check.se.obj = FALSE,
                remove.na = 'none',
                save.se.obj = FALSE
                )
            pcs[[assay.name]]$svd$u
        })
    names(all.pcs.on.ncgs) <- names(ncgs)
    ## Comparing the performance of different NCG sets ####
    if (!is.null(bio.variables) & !is.null(uv.variables)){
        variables <- NULL
    }
    if (is.null(variables)) {
        printColoredMessage(
            message = '- Comparing the performance of different NCG sets:',
            color = 'blue',
            verbose = verbose
            )
        #### Finding the class of the variables ####
        printColoredMessage(
            message = '- Finding the class of the variable:',
            color = 'blue',
            verbose = verbose
            )
        all.vars <- c(bio.variables, uv.variables)
        class.all.vars <- sapply(
            all.vars,
            function(x) class(colData(x = se.obj)[[x]])
            )
        cont.vars <- names(class.all.vars[class.all.vars %in% c('numeric', 'integer')])
        cat.vars <- names(class.all.vars[class.all.vars %in% c('factor', 'character')])
        ## Linear regression analysis ####
        if (length(cont.vars) > 1){
            cont.vars.r.squareds <- lapply(
                names(ncgs),
                function(x) {
                    pcs <- all.pcs.on.ncgs[[x]]
                    r.squareds <- sapply(
                        cont.vars,
                        function(y) {
                            lm.reg <- summary(lm(colData(x = se.obj)[[y]] ~ pcs))
                            lm.reg$r.squared
                        })
                    names(r.squareds) <- cont.vars
                    r.squareds
                })
            names(cont.vars.r.squareds) <- names(ncgs)
            cont.vars.r.squareds <- as.data.frame(cont.vars.r.squareds)
        } else cont.vars.r.squareds <- NULL
        ## Vector correlation analysis ####
        if (length(cat.vars) > 1){
            cat.vars.vec.corr <- lapply(
                names(ncgs),
                function(x) {
                    pcs <- all.pcs.on.ncgs[[x]]
                    vec.corr <- sapply(
                        cat.vars,
                        function(y) {
                            catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[y]])
                            catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                            cca <- cancor(x = pcs, y = catvar.dummies)
                            1 - prod(1 - cca$cor ^ 2)
                        })
                    names(vec.corr) <- cat.vars
                    vec.corr
                })
            names(cat.vars.vec.corr) <- names(ncgs)
            cat.vars.vec.corr <- as.data.frame(cat.vars.vec.corr)
        } else cat.vars.vec.corr <- NULL
        ## Putting all together ####
        all <- bind_rows(cont.vars.r.squareds, cat.vars.vec.corr) %>%
            round(digits = 3)
        all <- mutate(.data = all, var = row.names(all)) %>%
            pivot_longer(cols = -var, values_to = 'corr', names_to = 'geneset') %>%
            mutate(groups = 'unwanted') %>%
            mutate(groups = ifelse(var %in% bio.variables, 'wanted', groups)) %>%
            mutate(corr = ifelse(groups == 'wanted', 1 - corr, corr)) %>%
            group_by(geneset, groups) %>%
            summarise(corr = mean(corr)) %>%
            summarise(assess = corr[groups == 'wanted']/2 + corr[groups == 'unwanted']/2) %>%
            arrange(desc(assess)) %>%
            mutate(data = factor(geneset, levels = geneset))
        wrapped.labels <- gsub("\\.", "\n", levels(all$data))
        p.comparing.ncgs <- ggplot(all, aes(x = data, y = assess)) +
            geom_bar(stat = 'identity', fill = 'grey') +
            ylim(c(0, 1)) +
            ggtitle('') +
            scale_x_discrete(labels = wrapped.labels) +
            xlab('') +
            ylab('Summarized correlations') +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 18),
                axis.title.y = element_text(size = 18),
                plot.title = element_text(size = 15),
                axis.text.x = element_text(
                    size = 12,
                    angle = 35,
                    hjust = 0.9,
                    vjust = 1
                ),
                axis.text.y = element_text(size = 12)
            )
        if (isTRUE(plot.output)) print(p.comparing.ncgs)
    } else p.comparing.ncgs <- NULL
    ## Assessing the performance of NCG sets ####
    if (!is.null(variables)){
        printColoredMessage(
            message = '- Assessing the performance of NCG sets:',
            color = 'blue',
            verbose = verbose
            )
        ### Performing linear regression and vector correlation analysis ####
        all.ncg.assess.plots <- lapply(
            names(all.pcs.on.ncgs),
            function(n){
                all.corr <- lapply(
                    variables,
                    function(x){
                        if (class(se.obj[[x]]) %in% c('numeric', 'integer')){
                            rSquared <- sapply(
                                1:nb.pcs,
                                function(y) summary(lm(se.obj[[x]] ~ all.pcs.on.ncgs[[n]][, 1:y]))$r.squared
                            )
                        } else if (class(se.obj[[x]]) %in% c('factor', 'character')){
                            catvar.dummies <- dummy_cols(se.obj[[x]])
                            catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                            cca.pcs <- sapply(
                                1:nb.pcs,
                                function(y){ cca <- cancor(
                                    x = all.pcs.on.ncgs[[n]][, 1:y, drop = FALSE],
                                    y = catvar.dummies)
                                1 - prod(1 - cca$cor^2)
                                })
                        }
                    })
                names(all.corr) <- variables
                pca.ncg <- as.data.frame(do.call(cbind, all.corr))
                pca.ncg['pcs'] <- c(1:nb.pcs)
                pca.ncg <- tidyr::pivot_longer(
                    data = pca.ncg,
                    -pcs,
                    names_to = 'Groups',
                    values_to = 'ls'
                    )
                p.assess.ncg <- ggplot(pca.ncg, aes(x = pcs, y = ls, group = Groups)) +
                    geom_line(aes(color = Groups), linewidth = 1) +
                    geom_point(aes(color = Groups), size = 2) +
                    xlab('PCs') +
                    ylab (expression("Correlations")) +
                    ggtitle(n) +
                    scale_x_continuous(breaks = (1:nb.pcs),labels = c('PC1', paste0('PC1:', 2:nb.pcs)) ) +
                    scale_y_continuous(breaks = scales::pretty_breaks(n = nb.pcs), limits = c(0,1)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14),
                        axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
                        axis.text.y = element_text(size = 12),
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 14),
                        plot.title = element_text(size = 16)
                    )
            })
        p.assess.ncg <- ggarrange(
            plotlist = all.ncg.assess.plots,
            common.legend = TRUE,
            ncol = 2
            )
        if (isTRUE(plot.output)) print(p.assess.ncg)
    } else p.assess.ncg <- NULL
    ## Saving the results ####
    printColoredMessage(
        message = '-- Saving the results:',
        color = 'magenta',
        verbose = verbose
        )
    if (is.null(output.name)){
        output.name <- paste0(
            'comparsion_',
            paste0(names(ncgs), collapse = '&')
            )
    }
    if (isTRUE(save.se.obj)){
        if ( ncg.type %in% c('supervised', 'un.supervised')){
            if (length(se.obj@metadata) == 0) {
                se.obj@metadata[['NCG']] <- list()
            }
            if (!ncg.type %in% names(se.obj@metadata[['NCG']])){
                se.obj@metadata[['NCG']][[ncg.type]] <- list()
            }
            if (! ncg.group.name %in% names(se.obj@metadata[['NCG']][[ncg.type]])){
                se.obj@metadata[['NCG']][[ncg.type]][[ ncg.group.name]] <- list()
            }
            if (!'assessment.plot' %in% names(se.obj@metadata[['NCG']][[ncg.type]][[ ncg.group.name]])){
                se.obj@metadata[['NCG']][[ncg.type]][[ ncg.group.name]][['assessment.plot']] <- list()
            }
            if (!output.name %in% names(se.obj@metadata[['NCG']][[ncg.type]][[ ncg.group.name]][['assessment.plot']] )){
                se.obj@metadata[['NCG']][[ncg.type]][[ ncg.group.name]][['assessment.plot']][[output.name]] <- list()
            }
            if (!is.null(p.assess.ncg)){
                se.obj@metadata[['NCG']][[ncg.type]][[ ncg.group.name]][['assessment.plot']][[output.name]]<- p.assess.ncg
            }
            if (!is.null(p.comparing.ncgs)){
                se.obj@metadata[['NCG']][[ncg.type]][[ncg.group.name]][['assessment.plot']][[output.name]]<- p.comparing.ncgs
            }
            return(se.obj)
        }
        if ( ncg.group.name == 'pre.selected'){
            if (length(se.obj@metadata$NCG) == 0 ) {
                se.obj@metadata[['NCG']] <- list()
            }
            if (!'pre.selected' %in% names(se.obj@metadata[['NCG']])){
                se.obj@metadata[['NCG']][['pre.selected']] <- list()
            }
            if (!'assessment.plot' %in% names(se.obj@metadata[['NCG']][['pre.selected']])){
                se.obj@metadata[['NCG']][['pre.selected']][['assessment.plot']] <- list()
            }
            if (!output.name %in% names(se.obj@metadata[['NCG']][['pre.selected']][['assessment.plot']])){
                se.obj@metadata[['NCG']][['pre.selected']][['gene.set']][['assessment.plot']][[output.name]] <- list()
            }
            if (!is.null(p.assess.ncg)){
                se.obj@metadata[['NCG']][[ncg.type]][[ ncg.group.name]][['assessment.plot']][[output.name]]<- p.assess.ncg
            }
            if (!is.null(p.comparing.ncgs)){
                se.obj@metadata[['NCG']][[ncg.type]][[ ncg.group.name]][['assessment.plot']][[output.name]]<- p.comparing.ncgs
            }
            se.obj@metadata[['NCG']][['pre.selected']][['gene.set']][['assessment.plot']][[output.name]] <- p.assess.ncg
        }
        if ( ncg.group.name == 'all.genes'){
            if (length(se.obj@metadata$NCG) == 0 ) {
                se.obj@metadata[['NCG']] <- list()
            }
            if (!'all.genes' %in% names(se.obj@metadata[['NCG']])){
                se.obj@metadata[['NCG']][['all.genes']] <- list()
            }
            if (!'assessment.plot' %in% names(se.obj@metadata[['NCG']][['pre.selected']])){
                se.obj@metadata[['NCG']][['all.genes']][['assessment.plot']] <- list()
            }
            if (!output.name %in% names(se.obj@metadata[['NCG']][['all.genes']][['assessment.plot']])){
                se.obj@metadata[['NCG']][['all.genes']][['assessment.plot']][[output.name]] <- list()
            }
            if (!is.null(p.assess.ncg)){
                se.obj@metadata[['NCG']][['all.genes']][['assessment.plot']][[output.name]]<- p.assess.ncg
            }
            if (!is.null(p.comparing.ncgs)){
                se.obj@metadata[['NCG']][['all.genes']][['assessment.plot']][[output.name]]<- p.comparing.ncgs
            }
        }
    }
    if (isFALSE(save.se.obj)){
        return(list(p.assess.ncg = p.assess.ncg))
    }
}


