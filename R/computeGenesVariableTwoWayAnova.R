#' Compute two-way ANOVA across all genes in RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function calculates the two-way ANOVA between individual gene expression of the dataset(s) in a `SummarizedExperiment`
#' object and either categorical, continuous or a combination of these variable as a factor.
#'
#' @details
#' ANOVA enables assessment of the effect of a given qualitative variable (factor) on gene expression measurements
#' across groups labeled by the factor's levels. F-statistics, effect size and p-values summarize the effects of a
#' qualitative  variation (e.g., batch) on gene expression, with larger F-statistics indicating stronger association. P-values
#' are assigned using R's `aov()` function to quantify associations with variables such as tumor purity or molecular subtypes.
#'
#' @param se.obj A `SummarizedExperiment` object.
#' @param assay.names Character. A character string or vector specifying the name(s) of the assay(s) to compute ANOVA.
#' Defaults to `"all"` to select all assays.
#' @param bio.variables Character. Column name(s) in `colData(se.obj)` representing biological variables for ANOVA.
#' @param uv.variables Character. Column name(s) in `colData(se.obj)` representing unwanted variation variables.
#' @param nb.bio.clusters Numeric. Number of clusters to use for grouping biological variables if clustering is applied.
#' @param bio.clustering.method Character. Method used for clustering biological variables (e.g., `"kmeans"`, `"hierarchical"`).
#' @param nb.uv.clusters Numeric. Number of clusters to use for unwanted variation variables if clustering is applied.
#' @param uv.clustering.method Character. Method used for clustering unwanted variation variables (e.g., `"kmeans"`, `"hierarchical"`).
#' @param samples.to.use Character or NULL. Vector of sample IDs to include in the analysis. If `NULL`, all samples are used.
#' @param apply.log Logical. Indicates whether to log-transform the assay data before ANOVA. Default is `TRUE`.
#' @param pseudo.count Numeric. Pseudo-count added to assay values before log transformation. Default is 1.
#' @param nb.cores Numeric. Number of CPU cores to use for parallel computation. Default is 1.
#' @param check.se.obj Logical. If `TRUE`, validates the `SummarizedExperiment` object before analysis. Default is `TRUE`.
#' @param remove.na Character. Specifies how to handle missing values: `"assays"`, `"sample.annotation"`, `"both"`, or `"none"`.
#' Default is `"both"`.
#' @param override.check Logical. If `TRUE`, recalculates ANOVA even if results already exist in the object metadata. The
#' default is set to `FALSE`.
#' @param save.se.obj Logical. If `TRUE`, saves F-statistics and p-values in the metadata of the `SummarizedExperiment` object;
#' otherwise returns results as a list. Default is `TRUE`.
#' @param verbose Logical. If `TRUE`, displays messages for progress and steps. Default is `TRUE`.
#'
#' @return Either a `SummarizedExperiment` object containing log2-transformed F-statistics and p-values for ANOVA
#' or a list of these results.
#'
#' @importFrom car Anova
#' @importFrom SummarizedExperiment assays assay
#' @importFrom tidyr pivot_longer %>%
#' @importFrom dplyr mutate
#' @import ggplot2
#'
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
