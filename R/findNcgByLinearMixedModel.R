#' Finds NCGs using two-way ANOVA.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function utilizes two-way ANOVA to identify a set of suitable genes as negative control genes (NCG) for RUV-III
#' normalization. Both biological and unwanted variation sources is necessary and should be specified.
#'
#' @details
#' The function begins by creating all possible sample groups based on biological and unwanted variation separately.
#' Subsequently, these groups are used as factors in two-way ANOVA to identify genes highly influenced by biological and
#' unwanted variation. Finally, the function selects genes with the possible highest F-statistics for unwanted variation and
#' lowest F-statistics for biological variation. Various approaches are employed for the final gene selection; please refer
#' to the details for more information.
#' The function uses 5 ways to summarize two gene-level F-statistics obtained for the biological and unwanted variation.
#' The function uses either the values or the ranks of F-statistics for NCGs selection. The function ranks the
#' negative of F-statistics values for unwanted variation. The lower the ranks, the greater the impact of unwanted
#' variation on genes. The function ranks the F-statistics for biological variation. The higher the ranks, the greater
#' the impact of biological variation on genes. The options are `prod`, `sum`, `average`, `auto`, `non.overlap` and
#' `quantile`.
#'
#' If `prod`, `sum` and `average` is set:
#'
#' * The product, sum or average of ranks of F-statistics is calculated. Then, the function selects `nb.ncg` numbers of
#' genes as negative control genes that have the lowest ranks.
#'
#' If `non.overlap` is selected:
#' \enumerate{
#'    \item The function selects the top `top.rank.bio.genes` genes that have the highest ranks of F-statistics
#'    for biological variation.
#'    \item The function selects the top `top.rank.uv.genes` genes that have the lowest ranks of F-statistics for
#'    unwanted variation.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1. This will be a set of genes as
#'    negative control genes.
#' }
#'
#' If `auto` is selected:
#' \enumerate{
#'    \item The function selects the top `top.rank.bio.genes` genes that have the highest ranks of F-statistics for
#'    biological variation.
#'    \item The function selects the top `top.rank.uv.genes` genes that have the lowest ranks of F-statistics for
#'    unwanted variation.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1.
#'    \item If the number of selected genes is larger or smaller than the specified `nb.ncg`, the function applies an
#'    auto search to find approximate `nb.ncg` of genes as negative control genes as follow. The auto search will either
#'    decrease or increase the values of either `top.rank.bio.genes` or `top.rank.uv.genes` or both till to find
#'    approximate `nb.ncg` of genes as negative control genes.
#' }
#'
#' If `quantile` is selected:
#' \enumerate{
#'    \item The function selects the `bio.percentile` percentile of F-statistics for biological variation. Then, selects
#'    all the genes that have F-statistics larger the calculated percentile.
#'    \item The function selects the `uv.percentile` percentile of F-statistics for unwanted variation. Then, selects
#'    all the genes that have F-statistics larger the calculated percentile.
#'    \item The function excludes all genes obtained in 2 from the ones obtained 1.
#' }
#'
#' Assess the performance of NCGS:
#' * The function can assess the initial performance of selected NCGs. This analysis involves principal component analysis
#' on only the selected NCG and then explore the R^2 or vector correlation between the `nb.pcs` first principal components
#' and with the specified variables. Ideal NCGS, should show high and low R^2 or vector correlation for unwanted and
#' biological variation respectively.
#'
#' @references
#' * Gandolfo L. C. & Speed, T. P., RLE plots: visualizing unwanted variation in high dimensional data. PLoS ONE, 2018.
#' * Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character string indicating the name of the data (assay) in the `SummarizedExperiment`
#' object. This data should be the one that will be used as input data for the RUV-III normalization.
#' @param bio.variables Character. A character string or vector of strings indicating the column name(s) of the biological
#' variable(s) in the SummarizedExperiment object. These variable can be categorical or continuous or a combination. This
#' argument cannot be `NULL`.
#' @param uv.variables Character. A character string or vector of strings indicating the column name(s) of the unwanted
#' variable(s) in the SummarizedExperiment object. These variable can be categorical or continuous or a combination.This
#'  argument cannot be `NULL`.
#' @param ncg.selection.method Character. A character that indicates how to summarize different statistics and select a
#' set of genes as negative control genes. The options are: `prod`, `average`, `sum`, `non.overlap`, `auto`, and `quantile`.
#' The default is set to `non.overlap`. For more information, refer to the details of the function.
#' @param nb.ncg Numeric. A numeric value that specifies the number of genes to be chosen as negative control genes (NCG)
#' when the `ncg.selection.method` parameter is set to `auto`. This value, `nb.ncg`, corresponds to a fraction of the total
#' genes in the SummarizedExperiment object. The default is set to 0.1.
#' @param top.rank.bio.genes Numeric. A numeric value that indicates the percentage of top-ranked genes that are highly
#' affected by biological variation. This is required to be specified when the `ncg.selection.method` is either `auto`
#' or `non.overlap`. The default is set to 0.2.
#' @param top.rank.uv.genes Numeric. A numeric value that indicates the percentage of top-ranked genes that are highly
#' affected by unwanted variables. This is required to be specified when the `ncg.selection.method` is either `auto` or
#' `non.overlap`. The default is set to 0.2.
#' @param bio.percentile Numeric. A numeric value that specifies the percentile cut-off to select genes that are highly
#' affected by biological variation. This is required to be specified when the `ncg.selection.method` is set to `quantile`.
#' The default is set to 0.8.
#' @param uv.percentile Numeric. A numeric value that specifies the percentile cut-off to select genes that are highly
#' affected by unwanted variation. This is required to be specified when the `ncg.selection.method` is set to `quantile`.
#' The default is set to 0.8.
#' @param grid.group Character. A character that indicates whether the grid search should be performed on biological,
#' unwanted, or both factors when the `ncg.selection.method` is set to `auto`. The options are `bio`, `uv`, or `both`.
#' The default is set to `uv`. For more details, refer to the function documentation.
#' @param grid.direction Character. A character that indicates whether the grid search should be performed in decreasing
#' or increasing order when the `ncg.selection.method` is set to `auto`. The options are: `increase` and `decrease`. The
#' default is set to `decrease`.
#' @param grid.nb Numeric. A numeric value that indicates the number of genes for grid search when the `ncg.selection.method`
#' is set to `auto`. In the `auto` approach, the grid search starts with the initial `top.rank.uv.genes` and
#' `top.rank.bio.genes` values and adds or drops the `grid.nb` in each loop to find `nb.ncg` of genes as negative control
#' genes. The default is set to 20.
#' @param create.ncg.rank.plot Logical. Indicates whether to generate a heatmap that shows the rank of the all genes
#' with respect to their biological and unwanted variation effects. The default is set to `FALSE`.
#' @param plot.ncg.rank Logical. Indicates whether to plot a heatmap that shows the rank of the all genes
#' with respect to their biological and unwanted variation effects, while function is running. The default is set to `FALSE`.
#' @param bio.clustering.method Character. A character Indicating which clustering methods should be used to group continuous
#' sources of biological variation if any is provided. The options are: `kmeans`, `cut`, and `quantile`. The default is
#' set to `kmeans` clustering. Refer to the `createHomogeneousBioGroups()` function for more details.
#' @param nb.bio.clusters Numeric. A numeric value indicating the number of clusters for each continuous source of biological
#' variation. The default is set to 2.
#' @param uv.clustering.method Character. Indicates which clustering methods should be used to group continuous sources
#' of unwanted variation if any is provided. The options are: `kmeans`, `cut`, and `quantile`. The default is
#' set to `kmeans` clustering. Refer to the `createHomogeneousUvGroups()` function for more details.
#' @param nb.uv.clusters Numeric. A numeric that indicates the number of clusters for each continuous source of unwanted
#' variation. The default is set to 2.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing any statistical
#' analysis. The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value to be added as a pseudo count to all measurements before log transformation
#' .The default is se to 1.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected genes as negative control or not.
#' This analysis involves principal component analysis on the selected genes, followed by exploration of the R^2 or vector
#' correlation between the first `nb.pcs` principal components and the biological and unwanted variables. The default is
#' set to `TRUE`.
#' @param variables.to.assess.ncg Character. A character string or vector of strings indicating the column names in sample
#' annotation of of the  SummarizedExperiment object that contain variables whose association with the selected genes as
#' NCG needs to be evaluated. The default is set to `NULL`. This means all the variables specified in the `bio.variables`
#' and `uv.variables` will be assessed.
#' @param nb.pcs Numeric. A numeric value that indicates the number of the first principal components of selected negative
#' control genes to be used to assess their performance. The default is set to 10.
#' @param center Logical. Indicates whether to scale the data before applying SVD. If `TRUE`, centering is done by subtracting
#' the column means of the assay from their corresponding columns. The default is set to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data before applying SVD. If `TRUE`, scaling is done by dividing
#' the(centered) columns of the assays by their standard deviations if centering is `TRUE`, and by the root mean square
#' otherwise. The default is set to `FALSE`.
#' @param plot.ncg.assessment Logical. Indicates whether to plot the output of the NCG assessment while function is running
#' . The default is set to `TRUE`.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object before any analysis. If `TRUE`,
#'  the function `checkSeObj()` will be used. The default is set to `TRUE`.
#' @param remove.na Character set. Indicates whether to remove NA or missing values from the SummarizedExperiment object
#' The options are: `assays`, the `sample.annotation`, `both`, or `none`. If `assays` is selected, genes containing NA or
#' missing values will be excluded. If `sample.annotation` is selected, the samples containing NA or missing values for
#' any `bio.variables` or `uv.variables` will be excluded. The default is set to `none`.
#' @param ncg.group.name Character. A character to be used as name of the group of NCG. The default is set to `NULL`, then
#' the function create a names as following: `paste0('ncg|unsupervised')`. We refer to the details of the function for
#' more details.
#' @param ncg.set.name Character. A character to be used as name of the NCG set based on current variables and parameters
#' The default is set to `NULL`, then the function create a names as following:
#' `paste0(sum(ncg.selected),'|',paste0(bio.variables, collapse = '&'),'|',paste0(uv.variables, collapse = '&'),'|AnoCorrAs:',
#' ncg.selection.method,'|',assay.name)`.We refer to the details of the function for more details.
#' @param save.imf Logical. Indicates whether to save the intermediate file. If `TRUE`, the function saves the results
#' of the statistical analyses in the metadata of the SummarizedExperiment object. If users want to change the parameters
#' including `nb.ncg`, `ncg.selection.method`, `top.rank.bio.genes`, and `top.rank.uv.genes`, the analyses will not be
#' re-calculated. The default is set to `FALSE`.
#' @param imf.name Character string. A name to save the intermediate file. If `NULL`, the function generates a name.
#' @param use.imf Logical. Indicates whether to use the intermediate file. The default is set to `FALSE`.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the SummarizedExperiment
#' object or output the result. The default is `TRUE`.
#' @param verbose Logical. If `TRUE`, shows messages of different steps of the function.
#'
#' @importFrom variancePartition fitExtractVarPartModel
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocParallel SnowParam
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @export






# se.obj = read.se.obj
# assay.name = 'RawCount'
# bio.variables
# uv.variables
# form =  ~ (1|Time.interval) + Library.size + (1|CMS) + Tumour.purity
# use.rank = FALSE
# ncg.selection.method = 'quantile'
# adjust.data = FALSE
# adjustment.method = 'lmm'
# adjustment.variables = 'uv'
# samples.to.use = 'all'
# nb.ncg = 0.1
# top.rank.bio.genes = 0.7
# top.rank.uv.genes = 0.7
# bio.percentile = 0.2
# uv.percentile = 0.2
# grid.group = 'uv'
# grid.direction = 'increase'
# grid.nb = 20
# assess.ncg = TRUE
# apply.log = TRUE
# pseudo.count = 1
# variables.to.assess.ncg = NULL
# nb.pcs = 5
# center = TRUE
# scale = FALSE
# plot.ncg.assessment = TRUE
# nb.cores = 10
# check.se.obj = TRUE
# remove.na = 'none'
# ncg.group.name = NULL
# ncg.set.name = NULL
# save.imf = FALSE
# imf.name = NULL
# use.imf = FALSE
# save.se.obj = TRUE
# verbose = TRUE

findNcgByLinearMixedModel <- function(
        se.obj,
        assay.name,
        bio.variables,
        uv.variables,
        form,
        use.rank = FALSE,
        ncg.selection.method = 'quantile',
        adjust.data = FALSE,
        adjustment.method = 'lmm',
        adjustment.variables = 'uv',
        samples.to.use = 'all',
        nb.ncg = 0.1,
        top.rank.bio.genes = 0.7,
        top.rank.uv.genes = 0.7,
        bio.percentile = 0.2,
        uv.percentile = 0.2,
        grid.group = 'uv',
        grid.direction = 'increase',
        grid.nb = 20,
        assess.ncg = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        plot.ncg.assessment = TRUE,
        nb.cores = 1,
        check.se.obj = TRUE,
        remove.na = 'none',
        ncg.group.name = NULL,
        ncg.set.name = NULL,
        save.imf = FALSE,
        imf.name = NULL,
        use.imf = FALSE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    if (is.logical(samples.to.use)){
        se.obj.all <- se.obj
        se.obj <- se.obj[ , samples.to.use]
    }
    if (isFALSE(adjust.data)){
        all.lmm <- computeGenesVarianceComposition(
            se.obj = se.obj,
            assay.names = assay.name,
            form = form,
            adjust.data = FALSE,
            adjustment.form = NULL,
            adjustment.method = adjustment.method,
            samples.to.use = 'all',
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            nb.cores = nb.cores,
            plot.output = FALSE,
            output.name = NULL,
            check.se.obj = FALSE,
            remove.na = 'none',
            save.se.obj = FALSE,
            verbose = verbose
            )
    }
    if (isTRUE(adjust.data)){
        if (adjustment.variables == 'uv'){
            all.lmm <- computeGenesVarianceComposition(
                se.obj = se.obj,
                assay.names = assay.name,
                form = form,
                adjust.data = FALSE,
                adjustment.form = NULL,
                adjustment.method = adjustment.method,
                samples.to.use = 'all',
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                nb.cores = nb.cores,
                plot.output = FALSE,
                output.name = NULL,
                check.se.obj = FALSE,
                remove.na = 'none',
                save.se.obj = FALSE,
                verbose = verbose
                )
            adjustment.form <- changeLmmFormula(form = form, out.put = 'sub.formula', sub.set = uv.variables)
            form <- changeLmmFormula(form = form, out.put = 'sub.formula', sub.set = bio.variables)
            all.lmm.adjusted <- computeGenesVarianceComposition(
                se.obj = se.obj,
                assay.names = assay.name,
                form = form,
                adjust.data = adjust.data,
                adjustment.form = adjustment.form,
                adjustment.method = adjustment.method,
                samples.to.use = 'all',
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                nb.cores = nb.cores,
                plot.output = FALSE,
                output.name = NULL,
                check.se.obj = FALSE,
                remove.na = 'none',
                save.se.obj = FALSE,
                verbose = verbose
                )
            for(i in bio.variables){
                all.lmm[[assay.name]][[i]] <- all.lmm.adjusted[[assay.name]][[i]]
            }
        }
        if (adjustment.variables == 'both'){
            ## UV
            adjustment.form <- changeLmmFormula(form = form, out.put = 'sub.formula', sub.set = uv.variables)
            form.bio <- changeLmmFormula(form = form, out.put = 'sub.formula', sub.set = bio.variables)
            all.lmm.uv.adjusted <- computeGenesVarianceComposition(
                se.obj = se.obj,
                assay.names = assay.name,
                form = form.bio,
                adjust.data = adjust.data,
                adjustment.form = adjustment.form,
                adjustment.method = adjustment.method,
                samples.to.use = 'all',
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                nb.cores = nb.cores,
                plot.output = FALSE,
                output.name = NULL,
                check.se.obj = FALSE,
                remove.na = 'none',
                save.se.obj = FALSE,
                verbose = verbose
                )
            ## bio
            adjustment.form <- changeLmmFormula(form = form, out.put = 'sub.formula', sub.set = bio.variables)
            form.uv <- changeLmmFormula(form = form, out.put = 'sub.formula', sub.set = uv.variables)
            all.lmm.bio.adjusted <- computeGenesVarianceComposition(
                se.obj = se.obj,
                assay.names = assay.name,
                form = form.uv,
                adjust.data = adjust.data,
                adjustment.form = adjustment.form,
                adjustment.method = adjustment.method,
                samples.to.use = 'all',
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                nb.cores = nb.cores,
                plot.output = FALSE,
                output.name = NULL,
                check.se.obj = FALSE,
                remove.na = 'none',
                save.se.obj = FALSE,
                verbose = verbose
                )
            all.lmm <- cbind(all.lmm.bio.adjusted[[assay.name]], all.lmm.uv.adjusted[[assay.name]])
        }
    }
    all.lmm <- all.lmm[[assay.name]]
    all.lmm$uv <- rowSums(all.lmm[, uv.variables, drop = FALSE])
    all.lmm$bio <- rowSums(all.lmm[, bio.variables, drop = FALSE])
    set.seed(4455)
    all.lmm$uv.rank <- rank(x = -all.lmm$uv, ties.method = 'random')
    set.seed(4455)
    all.lmm$bio.rank <- rank(x = all.lmm$bio, ties.method = 'random')

    # Selecting a set of genes as NCG ####
    printColoredMessage(
        message = '-- Selecting a set of genes as NCG:',
        color = 'magenta',
        verbose = verbose
        )
    ## Applying ratio approach ####
    if (isFALSE(use.rank)){
        var.partition <- data.frame(
            gene = rownames(all.lmm),
            tech.var = all.lmm$uv,
            bio.var = all.lmm$bio,
            ratio = all.lmm$uv / (all.lmm$bio + 1e-6)  # Avoid division by zero
        )
        high.tech.genes <- var.partition[order(-var.partition$ratio), ]
        # high.tech.genes <- high.tech.genes[high.tech.genes$tech.var > top.rank.uv.genes , ]
        # nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
        # ncg.selected <- row.names(se.obj) %in% row.names(high.tech.genes)[1:nb.ncg]
        if (ncg.selection.method == 'quantile'){
            high.tech.genes <- high.tech.genes[high.tech.genes$tech.var > quantile(x = high.tech.genes$tech.var, probs = uv.percentile) , ]
            high.tech.genes <- high.tech.genes[high.tech.genes$bio.var < quantile(x = high.tech.genes$bio.var, probs = bio.percentile) , ]
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% high.tech.genes$gene[1:nb.ncg]
        } else {
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% high.tech.genes$gene[1:nb.ncg]
        }
    }
    ## Applying product, average and sum of ranks ####
    if (ncg.selection.method %in% c('prod', 'average', 'sum') & isTRUE(use.rank)){
        ### Product of ranks ####
        if (isTRUE(ncg.selection.method == 'prod')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the product of ranks.',
                color = 'blue',
                verbose = verbose)
            all.lmm$all.rank <- all.lmm$bio.rank * all.lmm$uv.rank
            if (sum(is.infinite(all.lmm$all.rank)) > 0)
                stop('The product of ranks results in infinity values.')
        }
        ## Average of ranks ####
        if (isTRUE(ncg.selection.method == 'average')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the average of ranks.',
                color = 'blue',
                verbose = verbose)
            all.lmm$all.rank <- rowMeans(all.lmm[ , c('bio.rank', 'uv.rank')])
        }
        ## Sum of ranks ####
        if (isTRUE(ncg.selection.method == 'sum')){
            printColoredMessage(
                message = '- A set of genes will be selected as NCGs based on the sum of ranks.',
                color = 'blue',
                verbose = verbose
            )
            all.lmm$all.rank <- rowSums(all.lmm[ , c('bio.rank', 'uv.rank')])
        }
        ## Selecting top genes as NCGS ####
        nb.ncg <- round(x = nb.ncg * nrow(se.obj), digits = 0)
        printColoredMessage(
            message = paste0(
                '- Selecting ',
                nb.ncg ,
                ' genes as NCGS.'),
            color = 'blue',
            verbose = verbose
        )
        all.lmm <- all.lmm[order(all.lmm$all.rank, decreasing = FALSE), ]
        ncg.selected <- row.names(all.lmm)[1:nb.ncg]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
    }

    ## Applying quantile approach ####
    if (ncg.selection.method == 'quantile' & isTRUE(use.rank)){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "quantile" approach.',
            color = 'blue',
            verbose = verbose
        )
        ### Selecting biological percentile ####
        bio.quan <- quantile(x = all.lmm$bio.rank , probs = bio.percentile)
        top.bio.genes <- row.names(all.lmm)[all.lmm$bio.rank > bio.quan]

        ## Selecting UV percentile ####
        uv.quan <- quantile(x = all.lmm$uv.rank , probs = uv.percentile)
        top.uv.genes <- row.names(all.lmm)[all.lmm$uv.rank < uv.quan]
        printColoredMessage(
            message = paste0(
                '- Selecting ',
                length(top.uv.genes),
                ' genes with the uv ranked F-statistics lower than ',
                uv.quan,
                ' (' ,
                uv.percentile* 100,
                '% percentile), and exclude any genes presents in ',
                length(top.bio.genes),
                ' genes with the biological ranked F-statistics higher than ',
                bio.quan,
                ' (' ,
                bio.percentile* 100,
                '% percentile).'),
            color = 'blue',
            verbose = verbose)
        ## Selecting top genes as NCGS ####
        top.uv.genes <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        if (isTRUE(length(top.uv.genes) == 0)) stop('No NCGs can be found based on the current parameters.')
        ncg.selected <- row.names(se.obj) %in% top.uv.genes
    }
    ## Applying non overlap approach ####
    if (ncg.selection.method == 'non.overlap' & isTRUE(use.rank)){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "non.overlap" approach.',
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '- Selecting top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation, and then exclude top ',
                top.rank.bio.genes *100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
        )
        ### Selecting genes highly affected by biological variation ####
        top.rank.bio.genes.nb <- round(c(1-top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(all.lmm)[all.lmm$bio.rank > top.rank.bio.genes.nb]
        ## Selecting genes highly affected by unwanted variation ####
        top.rank.uv.genes <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        top.uv.genes <- row.names(all.lmm)[all.lmm$uv.rank <  top.rank.uv.genes]
        ## Selecting top genes as NCGS ####
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        ncg.selected <- row.names(se.obj) %in% ncg.selected
        if (isTRUE(sum(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
    }
    ## Applying auto approach ####
    if (ncg.selection.method == 'auto' & isTRUE(use.rank)){
        printColoredMessage(
            message = '- A set of genes will be selected as NCGs based on the "auto" approach.',
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '- Selecting top ',
                top.rank.uv.genes * 100,
                '% of highly affected genes by the unwanted variation, and then exclude all genes in top ',
                top.rank.bio.genes * 100,
                '% of highly affected genes by the bioloigcal variation.'),
            color = 'blue',
            verbose = verbose
        )
        ### Selecting genes affected by biological variation ####
        nb.ncg <- round(nb.ncg * nrow(se.obj), digits = 0)
        top.rank.bio.genes.nb <- round(c(1 - top.rank.bio.genes) * nrow(se.obj), digits = 0)
        top.bio.genes <- row.names(all.lmm)[all.lmm$bio.rank > top.rank.bio.genes.nb]
        ## Selecting genes affected by unwanted variation ####
        top.rank.uv.genes.nb <- round(top.rank.uv.genes * nrow(se.obj), digits = 0)
        top.uv.genes <- row.names(all.lmm)[all.lmm$uv.rank < top.rank.uv.genes.nb]
        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
        # if (isTRUE(length(ncg.selected) == 0)) stop('NCGs cannot be found based on the current parameters.')
        printColoredMessage(
            message = paste0(
                '- ',
                length(ncg.selected),
                ' genes are found.'),
            color = 'blue',
            verbose = verbose
        )
        ncg.ranges <- round(x = 0.01 *nb.ncg, digits = 0)
        if (length(ncg.selected) > c(nb.ncg + ncg.ranges) | length(ncg.selected) < c(nb.ncg - ncg.ranges) ){
            if (isTRUE(nb.ncg > length(ncg.selected))){
                con <- parse(text = paste0("nb.ncg", ">", "length(ncg.selected)"))
                printColoredMessage(
                    message = paste0(
                        '- The number of selected genes ',
                        length(ncg.selected),
                        ' is less than the number (',
                        nb.ncg ,
                        ') of specified genes ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose)

            }
            if (isTRUE(nb.ncg < length(ncg.selected))){
                con <- parse(text = paste0("length(ncg.selected)", ">", "nb.ncg"))
                printColoredMessage(
                    message = paste0(
                        '- The number of selected genes ',
                        length(ncg.selected),
                        ' is larger than ',
                        nb.ncg ,
                        ', which is the specified number of NCG ',
                        'by "nb.ncg". A grid search will be performed.'),
                    color = 'blue',
                    verbose = verbose)
            }
            #### Applying grid search ####
            ##### grid group: both bio and uv variable ####
            if (grid.group == 'both'){
                printColoredMessage(
                    message = '- The grid search will be applied on both biological and unwanted factors. ',
                    color = 'blue',
                    verbose = verbose
                )
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                    )
                    lo <- min(
                        nrow(se.obj) - top.rank.uv.genes.nb,
                        top.rank.bio.genes.nb
                    )
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while (eval(con) & top.rank.uv.genes.nb < nrow(se.obj) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if (top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        top.uv.genes <- row.names(all.lmm)[all.lmm$uv.rank <  top.rank.uv.genes.nb]
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if (top.rank.bio.genes.nb < 1) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(all.lmm)[all.lmm$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of both "top.rank.uv.genes" and "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                    )
                    lo <- min(top.rank.uv.genes.nb, c(nrow(se.obj) - top.rank.bio.genes.nb))
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while (eval(con) & top.rank.uv.genes.nb > 1 & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if (top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        top.uv.genes <- row.names(all.lmm)[all.lmm$uv.rank <  top.rank.uv.genes.nb]
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if (top.rank.bio.genes.nb > nrow(se.obj)) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(all.lmm)[all.lmm$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                # genes selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                ## bio
                top.rank.bio.genes <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes/nrow(se.obj) * 100, digits = 2)
                if (top.rank.bio.genes >= 100) top.rank.bio.genes = 100
                ## uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if (top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Updating the selection. Selecting top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then excluding any genes in top ',
                        top.rank.bio.genes,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
            }
            ##### grid group: bio ####
            if (grid.group == 'bio'){
                printColoredMessage(
                    message = '- The grid search will be applied on biological factor. ',
                    color = 'blue',
                    verbose = verbose
                )
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the number of "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                    )
                    lo <- top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb - grid.nb
                        if (top.rank.bio.genes.nb < 1 ) top.rank.bio.genes.nb = 1
                        top.bio.genes <- row.names(all.lmm)[all.lmm$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the number of "top.rank.bio.genes". ',
                        color = 'blue',
                        verbose = verbose
                    )
                    lo <- nrow(se.obj) - top.rank.bio.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.bio.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # bio genes
                        top.rank.bio.genes.nb <- top.rank.bio.genes.nb + grid.nb
                        if (top.rank.bio.genes.nb > nrow(se.obj) ) top.rank.bio.genes.nb = nrow(se.obj)
                        top.bio.genes <- row.names(all.lmm)[ all.lmm$bio.rank > top.rank.bio.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
                }
                # gene selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # bio
                top.rank.bio.genes.nb <- nrow(se.obj) - top.rank.bio.genes.nb
                top.rank.bio.genes <- round(top.rank.bio.genes.nb/nrow(se.obj) * 100, digits = 2)
                if (top.rank.bio.genes >= 100) top.rank.bio.genes = 100

                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Updating the selection. Select top ',
                        top.rank.uv.genes * 100,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
            }
            ##### grid group: uv ####
            if (grid.group == 'uv'){
                printColoredMessage(
                    message = '- The grid search will be applied on unwanted factor. ',
                    color = 'blue',
                    verbose = verbose
                )
                ###### increasing order ####
                if (grid.direction == 'increase'){
                    printColoredMessage(
                        message = '- The grid search will increase the value of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose
                    )
                    lo <- nrow(se.obj) - top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb < nrow(se.obj)){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb + grid.nb
                        if (top.rank.uv.genes.nb > nrow(se.obj)) top.rank.uv.genes.nb = nrow(se.obj)
                        top.uv.genes <- row.names(all.lmm)[all.lmm$uv.rank <  top.rank.uv.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('NCGs cannot be found based on the current parameters.')
                }
                ##### decreasing order ####
                if (grid.direction == 'decrease'){
                    printColoredMessage(
                        message = '- The grid search will decrease the value of "top.rank.uv.genes". ',
                        color = 'blue',
                        verbose = verbose
                    )
                    lo <- top.rank.uv.genes.nb
                    pro.bar <- progress_estimated(round(lo/grid.nb, digits = 0) + 2)
                    while(eval(con) & top.rank.uv.genes.nb > 1){
                        pro.bar$pause(0.1)$tick()$print()
                        # uv genes
                        top.rank.uv.genes.nb <- top.rank.uv.genes.nb - grid.nb
                        if (top.rank.uv.genes.nb < 1) top.rank.uv.genes.nb = 1
                        top.uv.genes <- row.names(all.lmm)[all.lmm$uv.rank <  top.rank.uv.genes.nb]
                        ncg.selected <- top.uv.genes[!top.uv.genes %in% top.bio.genes]
                    }
                    if (length(ncg.selected) == 0)
                        stop('No NCGs can be found based on the current parameters.')
                }
                # gene selection
                ncg.selected <- row.names(se.obj) %in% ncg.selected
                ##### update numbers ####
                # uv
                top.rank.uv.genes <- round(top.rank.uv.genes.nb/nrow(se.obj) * 100, digits = 2)
                if (top.rank.uv.genes >= 100) top.rank.uv.genes = 100
                message(' ')
                printColoredMessage(
                    message = paste0(
                        '- Updating the selection. Selecting top ',
                        top.rank.uv.genes,
                        '% of highly affected genes by the unwanted variation, and then exclude any genes in top ',
                        top.rank.bio.genes * 100,
                        '% of highly affected genes by the bioloigcal variation.'),
                    color = 'blue',
                    verbose = verbose)
            }
        } else {
            printColoredMessage(
                message = paste0(
                    length(ncg.selected),
                    ' genes are selected as NCGs.'),
                color = 'blue',
                verbose = verbose)
        }
    }
    printColoredMessage(
        message = paste0(
            '- ',
            sum(ncg.selected),
            ' genes are selected as negative control genes.'),
        color = 'blue',
        verbose = verbose
    )

    # Plotting the F-statistics ####
    if (isTRUE(create.ncg.rank.plot)){
        if (isTRUE(use.fvalues)){
            all.lmm$ncg <- factor(x = ncg.selected, levels = c('TRUE', 'FALSE'))
            p.fvalues <- ggplot(all.lmm, aes(x = log2(bio.fvalue), y = log2(uv.fvalue), color =  ncg)) +
                geom_point(alpha = .1) +
                scale_color_manual(values = c('darkgreen', 'grey10')) +
                xlab('Biology (log2 of F-values)') +
                ylab('Unwanted variation (log2 of F-values)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            p.fvalues.rank <- ggplot(all.lmm, aes(x = bio.rank, y = uv.rank.plot , color =  ncg)) +
                geom_point(alpha = .1) +
                scale_color_manual(values = c('darkgreen', 'grey10')) +
                xlab('Biology (rank of F-values)') +
                ylab('Unwanted variation (rank of F-values)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            all.plots <- ggarrange(p.fvalues + p.fvalues.rank)
            print(all.plots)
            rm(p.fvalues, p.fvalues.rank)
        }
        if (isFALSE(use.fvalues)){
            all.lmm$ncg <- factor(x = ncg.selected, levels = c('TRUE', 'FALSE'))
            p.fvalues <- ggplot(all.lmm, aes(x = bio.pct, y = uv.pct, color =  ncg)) +
                geom_point(alpha = .1) +
                scale_color_manual(values = c('darkgreen', 'grey10')) +
                xlab('Biology (percentage of variation)') +
                ylab('Unwanted variation (percentage of variation)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            p.fvalues.rank <- ggplot(all.lmm, aes(x = bio.rank, y = uv.rank.plot, color =  ncg)) +
                geom_point(alpha = .1) +
                scale_color_manual(values = c('darkgreen', 'grey10')) +
                xlab('Biology (rank of percentage of variation') +
                ylab('Unwanted variation (rank of percentage of variation)') +
                theme(
                    panel.background = element_blank(),
                    axis.line = element_line(colour = 'black', linewidth = 1),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    legend.title = element_text(size = 14),
                    plot.title = element_text(size = 16))
            all.plots <- ggarrange(p.fvalues + p.fvalues.rank)
            print(all.plots)
            rm(p.fvalues, p.fvalues.rank)
        }
    } else all.plots = NULL
    # Plot variance explained (optional, can be toggled if needed)
    # Assessment of selected set of NCG  ####
    ## pca ####
    if (isTRUE(assess.ncg)) {
        printColoredMessage(
            message = '-- Assessing the performance of selected NCG set:',
            color = 'magenta',
            verbose = verbose
        )
        if (is.null(variables.to.assess.ncg)) {
            all.variables <- c(bio.variables, uv.variables)
        } else all.variables <- variables.to.assess.ncg
        printColoredMessage(
            message = '- Performing PCA using only the selected genes as NCGs.',
            color = 'blue',
            verbose = verbose
        )
        ### Applying log2 + pseudo count transformation ####
        if (isTRUE(apply.log) & !is.null(pseudo.count)){
            printColoredMessage(
                message = paste0(
                    '- Applying log2 + ',
                    pseudo.count,
                    ' (pseudo.count) on the ',
                    assay.name,
                    ' data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        } else if (isTRUE(apply.log) & is.null(pseudo.count)){
            printColoredMessage(
                message = paste0(
                    '- Applying log2 on the ',
                    assay.name,
                    ' data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log)) {
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
            x = t(expr.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM =  bsparam(),
            center = center,
            scale = scale)$u

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
            all.variables,
            function(x) {
                if (class(se.obj[[x]]) %in% c('numeric', 'integer')) {
                    rSquared <- sapply(
                        1:nb.pcs,
                        function(y) summary(lm(se.obj[[x]] ~ pca.data[, 1:y]))$r.squared)
                } else if (class(se.obj[[x]]) %in% c('factor', 'character')) {
                    catvar.dummies <- dummy_cols(se.obj[[x]])
                    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca.pcs <- sapply(
                        1:nb.pcs,
                        function(y) {
                            cca <- cancor(x = pca.data[, 1:y, drop = FALSE], y = catvar.dummies)
                            1 - prod(1 - cca$cor ^ 2)
                        })
                }
            })
        names(all.corr) <- all.variables
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
            scale_y_continuous(breaks = scales::pretty_breaks(n = nb.pcs), limits = c(0, 1)) +
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
        if (isTRUE(plot.ncg.assessment)) print(assess.ncg.plot)
    }
    if (is.null(ncg.set.name)){
        ncg.set.name <- paste0(
            sum(ncg.selected),
            '|',
            paste0(bio.variables, collapse = '&'),
            '|',
            paste0(uv.variables, collapse = '&'),
            '|LMM',
            '|',
            assay.name)
    }
    if (is.null(ncg.group.name)){
        ncg.group.name <- 'LinearMixedModel'
    }
    ### Adding the results to the SummarizedExperiment object ####
    if (is.logical(samples.to.use)){
        se.obj <- se.obj.all
    }
    if (isTRUE(save.se.obj)){
        ## Check if metadata NCG already exists
        if (length(se.obj@metadata$NCG) == 0 ) {
            se.obj@metadata[['NCG']] <- list()
        }
        if (!'supervised' %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][['supervised']] <- list()
        }
        if (!ncg.group.name %in% names(se.obj@metadata[['NCG']][['supervised']])){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]] <- list()
        }
        if (!ncg.set.name %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]] )){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]] <- list()
        }
        if (!'ncg.set' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]])){
            se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- list()
        }
        se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- ncg.selected
        if (isTRUE(assess.ncg)){
            if (!'ranl.plot' %in% names(se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]])){
                se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- list()
            }
            se.obj@metadata[['NCG']][['supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- assess.ncg.plot
        }
        printColoredMessage(
            message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = '------------The findNcgByTwoWayAnova function finished.',
            color = 'white',
            verbose = verbose
        )
    }
    return(se.obj)
}
