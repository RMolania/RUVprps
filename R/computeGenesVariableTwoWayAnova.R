#' Computes ANOVA between individual gene expression and a categorical variable
#'
#' @author Ramyar Molania
#'
#' @description
#' This function calculates the ANOVA between individual gene expression of the assay(s) in a SummarizedExperiment object
#' and a categorical variable as factor.
#'
#' @details
#' ANOVA enables us to assess the effects of a given qualitative variable (which we call a factor) on gene expression
#' measurements across any set of groups (labeled by the levels of the factor) under study. We use ANOVA F-statistics
#' to summarize the effects of a qualitative source of unwanted variation (for example, batches) on the expression levels
#' of individual genes, where genes having large F-statistics are deemed to be affected by the unwanted variation.
#' We also use ANOVA tests (the aov() function in R) to assign P values to the association between tumor purity and
#' molecular subtypes.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character. A character string or vector of character strings specifying the name(s) of the assay(s)
#' in the SummarizedExperiment object to compute the ANOVA. By default is set to 'all', which means all assays of the
#' SummarizedExperiment object will be selected.
#' @param bio.variables Character. A character string indicating a column name in the sample annotation of the SummarizedExperiment
#' object that contains a categorical variable, such as experimental batches, etc.
#' @param uv.variables Character. A character string indicating a column name in the sample annotation of the SummarizedExperiment
#' object that contains a categorical variable, such as experimental batches, etc.
#' @param nb.bio.clusters TTT
#' @param bio.clustering.method TTTT
#' @param nb.uv.clusters TTT
#' @param uv.clustering.method TTT
#' @param samples.to.use TTT
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing ANOVA. The
#' default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value representing a pseudo count to be added to all measurements before applying
#' the log transformation. The default is set to 1.
#' @param nb.cores TTT
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object. The default is set to `TRUE`.
#' This means the function will apply the `checkSeObj()` function.
#' @param remove.na Character. A character string specifying whether to eliminate missing values from `assays`, `sample.annotation`,
#' `both`, or `none`. When 'assays' is chosen, genes with missing values will be omitted. If 'sample.annotation' is selected,
#' samples with NA or missing values for each 'variable' will be excluded. The default is 'both'.
#' @param override.check Logical. When set to TRUE, the function checks whether ANOVA has already been computed for the
#' current parameters on the SummarizedExperiment object. If it has, the metric will not be recalculated. The default is
#' set to `FALSE`.
#' @param save.se.obj Logical. Indicates whether to save the results, ANOVA F-statistics, and p-values in the metadata
#' of the SummarizedExperiment object or to output these results as a list or vector. The default is set to `TRUE`.
#' @param verbose Logical. If `TRUE`, displays the messages of different steps of the function.
#' @return Either a SummarizedExperiment object containing the log2 F-statistics and p-values of ANOVA for the continuous
#' variable or a list of these results.
#'
#' @importFrom car Anova
#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer %>%
#' @importFrom dplyr mutate
#' @import ggplot2
#' @export

computeGenesVariableTwoWayAnova <- function(
        se.obj,
        assay.names,
        bio.variables,
        uv.variables,
        nb.bio.clusters = 3,
        bio.clustering.method = 'kmeans',
        nb.uv.clusters = 3,
        uv.clustering.method = 'kmeans',
        samples.to.use = 'all',
        apply.log = TRUE,
        pseudo.count = 1,
        nb.cores = 1,
        check.se.obj = TRUE,
        remove.na = 'none',
        override.check = FALSE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    # Checking samples to use ####
    if (is.logical(samples.to.use)){
        se.obj.initial <- se.obj
        se.obj <- se.obj[ , samples.to.use]
        }

    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
        } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if (!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
        }

    # Assessing the SummarizedExperiment object ####
    if (isTRUE(check.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            variables = c(bio.variables, uv.variables),
            remove.na = remove.na,
            verbose = verbose
            )
        }
    # Data log transformation ####
    if (isTRUE(apply.log)){
        printColoredMessage(
            message = '-- Applying log transformation on all the specified assay(s):',
            color = 'magenta',
            verbose = verbose
            )
        all.assays <- applyLog(
            se.obj = se.obj,
            assay.names = levels(assay.names),
            pseudo.count = pseudo.count,
            check.se.obj = check.se.obj,
            verbose = verbose
            )
        }
    if (isFALSE(apply.log)){
        printColoredMessage(
            message = '-- The specified assay(s) will be used for LMM, without applying log transformation.',
            color = 'blue',
            verbose = verbose
            )
        all.assays <- lapply(
            levels(assay.names),
            function(x) assay(x = se.obj, i = x)
            )
        names(all.assays) <- levels(assay.names)
        }

    # Creating all possible homogeneous sample groups ####
    ### biological groups ####
    printColoredMessage(
        message = '-- Creating all possible groups with respect to the specified sources of biological variation:',
        color = 'magenta',
        verbose = verbose
        )
    all.bio.groups <- createHomogeneousBioGroups(
        se.obj = se.obj,
        bio.variables = bio.variables,
        nb.clusters = nb.bio.clusters,
        clustering.method = bio.clustering.method,
        check.se.obj = FALSE,
        save.se.obj = FALSE,
        remove.na = 'none',
        verbose = verbose
        )

    ### unwanted groups ####
    printColoredMessage(
        message = '-- Creating all possible groups with respect to the specified sources of unwanted variation:',
        color = 'magenta',
        verbose = verbose
        )
    all.uv.groups <- createHomogeneousUVGroups(
        se.obj = se.obj,
        uv.variables = uv.variables,
        nb.clusters = nb.uv.clusters,
        clustering.method = uv.clustering.method,
        check.se.obj = FALSE,
        save.se.obj = FALSE,
        remove.na = 'none',
        verbose = verbose
        )

    ## Applying two way ANOVA ####
    printColoredMessage(
        message = '-- Performing two way ANOVA:',
        color = 'magenta',
        verbose = verbose
        )

    printColoredMessage(
        message = paste0(
            '- This is between all individual gene-level expression',
            'and considering both biological and unwanted variables created above as factors.'),
        color = 'blue',
        verbose = verbose
        )
    ### two way ANOVA ####
    all.twa.results <- mclapply(
        levels(assay.names),
        function(d){
            expr.data <- all.assays[[d]]
            all.twa <- lapply(
                1:nrow(expr.data),
                function(x){
                    sub.data <- data.frame(batch = all.uv.groups, bio = all.bio.groups, gene = expr.data[x, ])
                    lm.fit <- lm(gene ~ batch + bio , sub.data)
                    result.twa <- car::Anova(lm.fit, type = "II")
                    total.sum.sq <- sum(result.twa$`Sum Sq`)
                    data.frame(
                        uv.pct = result.twa$`Sum Sq`[1] / total.sum.sq,
                        bio.pct = result.twa$`Sum Sq`[2] / total.sum.sq,
                        residual.pct = result.twa$`Sum Sq`[3] / total.sum.sq,
                        uv.fvalue = result.twa$`F value`[1],
                        bio.fvalue = result.twa$`F value`[2]
                        )
                })
            all.twa <- do.call(rbind, all.twa)
            all.twa <- round(x = all.twa, digits = 3)
            row.names(all.twa) <- row.names(se.obj)
            return(all.twa)
        }, mc.cores = nb.cores
        )
    names(all.twa.results) <-  levels(assay.names)
    variables <- paste(
        paste0('bio:', paste0(bio.variables, collapse = '||')),
        paste0('uv:', paste0(uv.variables, collapse = '||')),
        sep = '_'
        )
    if (isTRUE(save.se.obj)){
        se.obj <- addMetricToSeObj(
            se.obj = se.obj,
            slot = 'Metrics',
            assay.names = levels(assay.names),
            assessment.type = 'gene.level',
            assessment = 'TwoWayAnova',
            method = 'aov',
            variables = variables,
            file.name = 'gene.variance.fvalue',
            results.data = all.twa.results
        )
        return(se.obj)
    }
    if(isFALSE(save.se.obj)){
        return(all.twa.results)
    }
}
