#' Computes sample level enrichment score of a gene set.

#' @author Ramyar Molania

#' @description
#' The function uses the `singscore` function from the R/Bioconductor **singscore** package to calculate sample-wise
#' scores of a gene set across all data(assay) in a SummarizedExperiment object.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character or a vector of characters specifying the name(s) of the data set(s) in the
#' SummarizedExperiment object for which to calculate sample-wise scores. The default is set to `all`, indicating that all
#' assays in the SummarizedExperiment object will be selected.
#' @param upset.genes Character vector. A vector of gene names/IDs that are up-regulated in the gene set.
#' @param downset.genes Character vector. A vector of gene names/IDs that are down-regulated in the gene set. The default
#' is set to `NULL`.
#' @param normalization Character. A character string indicating the type of normalization to be used before computing
#' the scores. The option are: `CPM`, `TMM`, `VST`, `full`, `upper` and `median`. The default is set to `NULL`.
#' Refer to the `applyOtherNormalizations()` function for more details.
#' @param regress.out.variables Character. A characer or vector of characters indicating the name(s) of the column(s) in
#' the SummarizedExperiment object to be regressed out from the data before computing the scores. The default is set
#' to `NULL`.
#' @param apply.log Logical. Indicates whether to apply a log transformation to the data. The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value as pseudo-count to be added to all measurements of the data set(s) before
#' applying log transformation to avoid `-Inf` for measurements equal to 0. The default is set to 1.
#' @param assess.score Logical. If `TRUE`, the association between the computed scores and specified variables will be
#' assessed. The default is set to `FALSE`. See the details for more information.
#' @param variables.to.assess Character. A character string or vector of strings indicating the name(s) of the column(s)
#' in the SummarizedExperiment object to be used for assessing the computed scores. These can be continuous or categorical
#' variables. The default is set to `NULL`.
#' @param corr.method Character. A character string indicating which correlation method should be used for the correlation
#' analysis of the computed scores with the specified continuous variable(s). Options include `pearson`, `kendall` amd
#' `spearman`. The default is set to `spearman`.
#' @param gene.set.name Character. A character string indicating the name to be used to save the score in the metadata of
#' the SummarizedExperiment object. If `NULL`, the function will select a name as follows:
#' gene.set.name <- paste0('singscore|', length(c(upset.genes, downset.genes)), 'genes')
#' @param plot.output Logical. If `TRUE`, the assessment plot will be printed while running the function. The default is
#' se to `TRUE`.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. The default is se to `TRUE`.
#' Refer to the `checkSeObj()` function for more details.
#' @param save.se.obj Logical. Indicates whether to save the score results in the metadata of the SummarizedExperiment
#' object or to output the result as a list. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, messages describing the steps of the function will be shown.

#' @importFrom SummarizedExperiment assays assay
#' @importFrom singscore rankGenes simpleScore
#' @importFrom fastDummies dummy_cols
#' @importFrom stats cor.test cancor
#' @import ggplot2
#' @export

computeGeneSetScore <- function(
        se.obj,
        assay.names = 'all',
        upset.genes,
        downset.genes = NULL,
        normalization = NULL,
        regress.out.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        assess.score = FALSE,
        variables.to.assess = NULL,
        corr.method = 'spearman',
        gene.set.name = NULL,
        plot.output = TRUE,
        check.se.obj = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The computeGeneSetScore function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the function inputs ####
    if (is.logical(assay.names) | is.null(assay.names)){
        stop('The "assay.names" cannot be NULL or logical.')
    }
    if (!is.null(upset.genes)){
        if (!is.character(upset.genes)){
            stop('The "upset.genes" must be character.')
        }
    }
    if (!is.null(downset.genes)){
        if (!is.character(downset.genes)){
            stop('The "downset.genes" must be character.')
        }
    }
    if (!is.null(normalization)){
        if (!is.character(normalization)){
            stop('The "normalization" must be a character.')
        }
        if (!normalization %in% c('CPM', 'TMM', 'upper', 'median', 'full', 'VST')){
            stop('The "normalization" must be one of the "CPM", "TMM", "upper", "median", "full", or "VST".')
        }
    }
    if (!is.null(regress.out.variables)){
        if (!is.character(regress.out.variables)){
            stop('The "regress.out.variables" must be character.')
        }
        if (sum(regress.out.variables %in% colnames(colData(se.obj))) !=length(regress.out.variables) ){
            stop('All or some of the "regress.out.variables" variable cannot be found in the SummarizedExperiment object.')
        }
    }
    if (!is.logical(apply.log)){
        stop('The "apply.log" must be logical.')
    }
    if (isTRUE(apply.log)){
        if(!is.numeric(pseudo.count) | pseudo.count < 0){
            stop('The "pseudo.count" must be a positive numeric value.')
        }
    }
    if (!is.logical(assess.score)){
        stop('The "assess.score" must be logical.')
    }
    if (!is.null(variables.to.assess)){
        if (!is.character(variables.to.assess)){
            stop('The "variables.to.assess" must be character.')
        }
        if (sum(variables.to.assess %in% colnames(colData(se.obj))) !=length(variables.to.assess)){
            stop('All or some of the "variables.to.assess" variable cannot be found in the SummarizedExperiment object.')
        }
    }
    if (!corr.method %in% c('spearman', 'pearson', 'kendall')){
        stop ('The "corr.method" must be one of the "spearman", "pearson" or "kendall".')
    }
    if (!is.null(gene.set.name)){
        if (!is.character(gene.set.name) | length(gene.set.name) > 1){
            stop('The "gene.set.name" must be character.')
        }
    }
    if (!is.logical(plot.output)){
        stop('The "plot.output" must be logical.')
    }
    if (!is.logical(check.se.obj)){
        stop('The "check.se.obj" must be logical.')
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
    # Data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Data transformation and normalization:',
        color = 'magenta',
        verbose = verbose
        )
    all.assays <- lapply(
        levels(assay.names),
        function(x){
            ## Applying log ####
            if (is.null(normalization) & is.null(regress.out.variables)){
                do.log <- TRUE
            } else if (is.null(normalization) & !is.null(regress.out.variables)){
                do.log <- TRUE
            }
            if (isTRUE(do.log)){
                if (isTRUE(apply.log) & !is.null(pseudo.count)){
                    printColoredMessage(
                        message = paste0(
                            '- applying log2 + ',
                            pseudo.count,
                            ' (pseudo.count) on the ',
                            x,
                            ' data.'),
                        color = 'blue',
                        verbose = verbose)
                    expr.data <- log2(assay(x = se.obj, i = x) + pseudo.count)
                } else if (isTRUE(apply.log) & is.null(pseudo.count)){
                    printColoredMessage(
                        message = paste0(
                            '- Applying log2 on the "',
                            x,
                            '" data.'),
                        color = 'blue',
                        verbose = verbose)
                    expr.data <- log2(assay(x = se.obj, i = x))
                } else if (isFALSE(apply.log)) {
                    printColoredMessage(
                        message = paste0(
                            '- The "',
                            x,
                            '" data will be used without any log transformation.'),
                        color = 'blue',
                        verbose = verbose)
                    expr.data <- assay(x = se.obj, i = x)
                }
            }
            ## Normalization ####
            if (!is.null(normalization)) {
                printColoredMessage(
                    message = '-- Data normalization:',
                    color = 'magenta',
                    verbose = verbose
                    )
                expr.data <- applyOtherNormalizations(
                    se.obj = se.obj,
                    assay.name = x,
                    method = normalization,
                    pseudo.count = pseudo.count,
                    apply.log = apply.log,
                    check.se.obj = FALSE,
                    save.se.obj = FALSE,
                    remove.na = 'none',
                    verbose = verbose
                )
            }
            # Regressing  out variables ####
            if (!is.null(regress.out.variables)){
                printColoredMessage(
                    message = '-- Regressing out unwanted or biological variables:',
                    color = 'magenta',
                    verbose = verbose
                )
                printColoredMessage(
                    message = paste0(
                        'The ', paste0(regress.out.variables, collapse = ' & '),
                        ' will be regressed out from the data,', ' please make sure your data is log transformed.'),
                    color = 'blue',
                    verbose = verbose
                    )
                printColoredMessage(
                    message = paste0(
                        'We do not recommend regressing out ',
                        paste0(regress.out.variables, collapse = ' & '),
                        ' if they are largely associated with the ',
                        paste0(regress.out.variables, collapse = ' & '),
                        ' variables.'),
                    color = 'red',
                    verbose = verbose
                    )
                expr.data <- t(expr.data)
                uv.variables.all <- paste('se.obj', regress.out.variables, sep = '$')
                expr.data <- lm(as.formula(paste(
                    'expr.data',
                    paste0(uv.variables.all, collapse = '+') ,
                    sep = '~')))
                expr.data <- t(expr.data$residuals)
                colnames(expr.data) <- colnames(se.obj)
                row.names(expr.data) <- row.names(se.obj)
            }
            return(expr.data)
            }
        )
    names(all.assays) <- levels(assay.names)

    # Computing gene set scoring analysis ####
    all.scores <- lapply(
        levels(assay.names),
        function(x){
            # Ranking the data ####
            rank.data <- singscore::rankGenes(expreMatrix = all.assays[[x]])
            # Applying the singscore with only upset gene set ####
            if (is.null(downset.genes)){
                gene.set.score <- singscore::simpleScore(
                    rankData = rank.data,
                    upSet = upset.genes,
                    downSet = NULL,
                    subSamples = NULL,
                    centerScore = TRUE,
                    dispersionFun = mad,
                    knownDirection = TRUE
                )
            }
            # Applying the singscore with only downset gene set ####
            if (is.null(upset.genes)){
                gene.set.score <- singscore::simpleScore(
                    rankData = rank.data,
                    upSet = NULL,
                    downSet = downset.genes,
                    subSamples = NULL,
                    centerScore = TRUE,
                    dispersionFun = mad,
                    knownDirection = TRUE
                )
            }
            if (is.null(downset.genes) & !is.null(upset.genes)){
                gene.set.score <- singscore::simpleScore(
                    rankData = rank.data,
                    upSet = upset.genes,
                    downSet = downset.genes,
                    subSamples = NULL,
                    centerScore = TRUE,
                    dispersionFun = mad,
                    knownDirection = TRUE
                )
            }
            gene.set.score <- gene.set.score$TotalScore
            return(gene.set.score)
        })
    rm(all.assays)
    names(all.scores) <- levels(assay.names)

    #  Assessing the association between the scores and the specified variables ####
    if (isTRUE(assess.score)){
        printColoredMessage(
            message = '-- Assessing the association between the scores and the specified variables.',
            color = 'magenta',
            verbose = verbose
        )
        ## Finding the class of variables ####
        class.variables <- sapply(
            variables.to.assess,
            function(x) class(se.obj[[x]]))
        continuous.variables <- variables.to.assess[class.variables %in% c('numeric', 'integer')]
        categorical.variables <- variables.to.assess[class.variables %in% c('factor', 'charachter')]

        # Assessing the  association ####
        all.assessment.plots <- lapply(
            levels(assay.names),
            function(x){
                ## Applying correlation for continuous variable ####
                if (length(continuous.variables) > 1){
                    corr.continuous <- sapply(
                        continuous.variables,
                        function(y){
                            cor.test(x = all.scores[[x]], y = se.obj[[y]], method = corr.method)[[4]][[1]]
                        })
                }
                ## Applying vector correlation for categorical variable ####
                if (length(categorical.variables) > 1){
                    corr.categorical <- sapply(
                        categorical.variables,
                        function(y){
                            catvar.dummies <- fastDummies::dummy_cols(se.obj@colData[[y]])
                            catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                            cca <- cancor(x = all.scores[[x]], y = catvar.dummies)
                            1 - prod(1 - cca$cor ^ 2)
                        })
                }
                ## Plotting the results #####
                corr <- variable.name <- NULL
                all.corr <- data.frame(
                    corr = c(corr.continuous, corr.categorical),
                    variable.name = names(c(corr.continuous, corr.categorical))
                )
                p.assess.geneset <- ggplot(data = all.corr, aes(x = variable.name, y = corr)) +
                    geom_col() +
                    ggtitle('Assessment of the gene set scoring') +
                    xlab('Variables') +
                    ylab('Correlations') +
                    ylim(c(-1,1)) +
                    theme(
                        panel.background = element_blank(),
                        axis.line = element_line(colour = 'black', linewidth = 1),
                        plot.title = element_text(size = 18),
                        axis.title.x = element_text(size = 16),
                        axis.title.y = element_text(size = 16),
                        axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
                        axis.text.y = element_text(size = 12))
                if (isTRUE(plot.output)) print(p.assess.geneset)
                return(p.assess.geneset)
            })
        names(all.assessment.plots) <- levels(assay.names)
    }
    # Saving the results ####
    if (isTRUE(save.se.obj)){
        printColoredMessage(
            message = '- Saving all the gene set enrichment score into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        if (is.null(gene.set.name)){
            gene.set.name <- paste0(
                'singscore|',
                length(c(upset.genes, downset.genes)),
                'genes')
        }
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            slot = 'Metrics',
            assessment = 'GeneSetScore',
            assessment.type = 'global.level',
            method = 'singscore',
            file.name = 'score',
            variables = gene.set.name,
            results.data = all.scores
        )
        if (isTRUE(assess.score)){
            se.obj <- addMetricToSeObj(
                se.obj = se.obj,
                assay.names = levels(assay.names),
                slot = 'Metrics',
                assessment = 'GeneSetScore',
                assessment.type = 'global.level',
                method = 'singscore',
                file.name = 'assessment.plot',
                variables = gene.set.name,
                results.data = all.assessment.plots
            )
        }
        printColoredMessage(message = '------------The computeGeneSetScore function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## save the data as list ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(message = '------------The computeGeneSetScore function finished.',
                            color = 'white',
                            verbose = verbose)
        if (isTRUE(assess.score)){
            return(list(all.scores = all.scores, all.assessment.plots = all.assessment.plots))
        } else return(list(all.scores = all.scores))
    }
}
