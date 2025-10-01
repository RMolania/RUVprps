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
#' * Molania R., ..., Speed, T. P., Removing unwanted variation from large-scale RNA sequencing data with PRPS,
#' Nature Biotechnology, 2023
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character string indicating the name of the data (assay) in the `SummarizedExperiment`
#' object. This data should be the one that will be used as input data for the RUV-III normalization.
#' @param uv.variables Character. A character string or vector of strings indicating the column name(s) of the unwanted
#' variable(s) in the SummarizedExperiment object. These variables can be categorical, continuous, or a combination. This
#' argument cannot be `NULL`.
#' @param form Formula. A formula describing the relationship between biological and unwanted variation to be modeled
#' (e.g., `~ bio.variables + uv.variables`).
#' @param ncg.idenfitication.approach Character. A character that indicates how to summarize different statistics and select a
#' set of genes as negative control genes. The options are: `prod`, `average`, `sum`, `non.overlap`, `auto`, and `quantile`.
#' The default is set to `non.overlap`. For more information, refer to the details of the function.
#' @param nb.bio.pcs Numeric. A numeric value specifying the number of principal components of biological variation to be
#' included in the model. The default is set to 10.
#' @param nb.uv.pcs Numeric. A numeric value specifying the number of principal components of unwanted variation to be
#' included in the model. The default is set to 10.
#' @param ncg.selection.method Character. A character that indicates how to summarize different statistics and select a
#' set of genes as negative control genes. The options are: `prod`, `average`, `sum`, `non.overlap`, `auto`, and `quantile`.
#' The default is set to `non.overlap`. For more information, refer to the details of the function.
#' @param use.rank Logical. Indicates whether to use rank-based statistics for selecting negative control genes. The
#' default is set to `FALSE`.
#' @param regress.out.rle.med Logical. Indicates whether to regress out the median of the RLE (Relative Log Expression)
#' per sample after normalization. The default is set to `FALSE`.
#' @param samples.to.use Logical or numeric. A logical vector or numeric index specifying which samples to include in the
#' analysis. If `NULL`, all samples will be used.
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
#' @param filter.ncgs Logical. Indicates whether to apply filtering to the selected negative control genes (e.g., removing
#' genes with missing values or zero variance). The default is set to `TRUE`.
#' @param common.hk Logical. Indicates whether to restrict negative control gene selection to common housekeeping genes.
#' The default is set to `FALSE`.
#' @param nb.stable.genes Numeric. A numeric value specifying the number of stable genes (e.g., housekeeping) to include
#' as negative controls. The default is set to 50.
#' @param hk.group Character. A character string specifying the group of genes to be treated as housekeeping or stable
#' controls. The default is set to `NULL`.
#' @param create.ncg.rank.plot Logical. Indicates whether to generate a heatmap that shows the rank of all genes
#' with respect to their biological and unwanted variation effects. The default is set to `FALSE`.
#' @param assess.ncg Logical. Indicates whether to assess the performance of selected genes as negative controls or not.
#' This analysis involves principal component analysis on the selected genes, followed by exploration of the R^2 or vector
#' correlation between the first `nb.pcs` principal components and the biological and unwanted variables. The default is
#' set to `TRUE`.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data before performing any statistical
#' analysis. The default is set to `TRUE`.
#' @param pseudo.count Numeric. A numeric value to be added as a pseudo count to all measurements before log transformation.
#' The default is set to 1.
#' @param variables.to.assess.ncg Character. A character string or vector of strings indicating the column names in the sample
#' annotation of the SummarizedExperiment object that contain variables whose association with the selected genes as
#' NCG needs to be evaluated. The default is set to `NULL`. This means all the variables specified in the `bio.variables`
#' and `uv.variables` will be assessed.
#' @param nb.pcs Numeric. A numeric value that indicates the number of the first principal components of selected negative
#' control genes to be used to assess their performance. The default is set to 10.
#' @param center Logical. Indicates whether to center the data before applying SVD. If `TRUE`, centering is done by subtracting
#' the column means of the assay from their corresponding columns. The default is set to `TRUE`.
#' @param scale Logical. Indicates whether to scale the data before applying SVD. If `TRUE`, scaling is done by dividing
#' the (centered) columns of the assays by their standard deviations if centering is `TRUE`, and by the root mean square
#' otherwise. The default is set to `FALSE`.
#' @param svd.bsparam List. A list of parameters controlling the SVD or basis decomposition, such as the rank to compute or
#' the algorithm to use. The default is set to `NULL`.
#' @param plot.ncg.assessment Logical. Indicates whether to plot the output of the NCG assessment while the function is running.
#' The default is set to `TRUE`.
#' @param use.fvalues Logical. Indicates whether to use F-statistics (instead of correlation) for assessing the relationship
#' between negative control genes and biological/unwanted variables. The default is set to `FALSE`.
#' @param nb.cores Numeric. A numeric value specifying the number of CPU cores to use for parallel computation. The default
#' is set to 1.
#' @param check.se.obj Logical. Indicates whether to assess the SummarizedExperiment object before any analysis. If `TRUE`,
#' the function `checkSeObj()` will be used. The default is set to `TRUE`.
#' @param remove.na Character. Indicates whether to remove NA or missing values from the SummarizedExperiment object.
#' The options are: `assays`, `sample.annotation`, `both`, or `none`. If `assays` is selected, genes containing NA or
#' missing values will be excluded. If `sample.annotation` is selected, the samples containing NA or missing values for
#' any `bio.variables` or `uv.variables` will be excluded. The default is set to `none`.
#' @param ncg.group.name Character. A character to be used as the name of the group of NCGs. The default is set to `NULL`, in
#' which case the function creates a name as follows: `paste0('ncg|unsupervised')`. We refer to the details of the function
#' for more information.
#' @param ncg.set.name Character. A character to be used as the name of the NCG set based on current variables and parameters.
#' The default is set to `NULL`, in which case the function creates a name as follows:
#' `paste0(sum(ncg.selected), '|', paste0(bio.variables, collapse = '&'), '|', paste0(uv.variables, collapse = '&'),
#' '|AnoCorrAs:', ncg.selection.method, '|', assay.name)`. We refer to the details of the function for more information.
#' @param save.imf Logical. Indicates whether to save the intermediate file. If `TRUE`, the function saves the results
#' of the statistical analyses in the metadata of the SummarizedExperiment object. If users want to change the parameters
#' including `nb.ncg`, `ncg.selection.method`, `top.rank.bio.genes`, and `top.rank.uv.genes`, the analyses will not be
#' re-calculated. The default is set to `FALSE`.
#' @param imf.name Character. A character string specifying the name to save the intermediate file. If `NULL`, the function
#' generates a name automatically.
#' @param use.imf Logical. Indicates whether to use the intermediate file. The default is set to `FALSE`.
#' @param save.se.obj Logical. Indicates whether to save the result of the function in the metadata of the SummarizedExperiment
#' object or output the result. The default is `TRUE`.
#' @param verbose Logical. If `TRUE`, shows messages for different steps of the function.
#'
#' @importFrom variancePartition fitExtractVarPartModel
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocParallel SnowParam
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' @export

findNcgUnSupervisedByLinearMixedModel <- function(
        se.obj,
        assay.name,
        uv.variables,
        form,
        ncg.idenfitication.approach = 'LMM.BioUvAdjustment',
        nb.bio.pcs = 3,
        nb.uv.pcs = 3,
        ncg.selection.method = 'quantile',
        use.rank = FALSE,
        regress.out.rle.med = FALSE,
        samples.to.use = 'all',
        nb.ncg = 0.05,
        top.rank.bio.genes = 0.7,
        top.rank.uv.genes = 0.7,
        bio.percentile = 0.2,
        uv.percentile = 0.7,
        grid.group = 'uv',
        grid.direction = 'increase',
        grid.nb = 20,
        filter.ncgs = FALSE,
        common.hk = 'cancer',
        nb.stable.genes = 2000,
        hk.group = 'micorarray.hk.genes',
        create.ncg.rank.plot = FALSE,
        assess.ncg = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        variables.to.assess.ncg = NULL,
        nb.pcs = 5,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        plot.ncg.assessment = TRUE,
        use.fvalues = TRUE,
        nb.cores = NULL,
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
    printColoredMessage(
        message = '------------The findNcgUnSupervisedByLinearMixedModel function starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking the  functions inputs ####
    if (length(assay.name) > 1 | !is.character(assay.name)){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    }
    if (!is.character(uv.variables)){
        stop('The "uv.variables" must be a character vector of variable(s) in the SummarizedExperiment object.')
    }
    if (sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables)){
        stop('Some or all the "uv.variables" cannot be found in the SummarizedExperiment object.')
    }
    if (!is.logical(use.rank)) {
        stop('The "use.rank" must be "TRUE" or "FALSE.')
    }
    if (nb.ncg >= 1 | nb.ncg <= 0){
        stop('The "nb.ncg" should be a positve value  0 < nb.ncg < 1.')
    }
    if (!ncg.idenfitication.approach %in% c('LMM', 'LMM.UvAdjustment','LMM.BioUvAdjustment')){
        stop('The "ncg.idenfitication.approach" must be one of "LMM", "LMM.UvAdjustment" or "LMM.BioUvAdjustment".')
    }
    if (!ncg.selection.method %in% c('auto', 'non.overlap','quantile')){
        stop('The "ncg.selection.method" must be one of "auto", "quantile" or "non.overlap".')
    }
    if (top.rank.bio.genes > 1 | top.rank.bio.genes <= 0){
        stop('The "top.rank.bio.genes" msut be a positve value  0 < top.rank.bio.genes < 1.')
    }
    if (top.rank.uv.genes > 1 | top.rank.uv.genes <= 0){
        stop('The "top.rank.uv.genes" must be a positve value  0 < top.rank.uv.genes < 1.')
    }
    if (grid.nb < 1 | grid.nb > nrow(se.obj)){
        stop(paste0('The "grid.nb" must be a positve value  0 < grid.nb < ', nrow(se.obj), '.'))
    }
    if (!grid.group %in% c('bio', 'uv', 'both')){
        stop('The "grid.group" must be one of "bio", "uv" or "non.overlap".')
    }
    if (!grid.direction %in% c('increase', 'decrease', 'auto')){
        stop('The "grid.direction" must be one of "increase", "decrease" or "auto".')
    }
    if (isFALSE(is.logical(assess.ncg))){
        stop('The "assess.ncg" must be "TRUE" or "FALSE.')
    }
    if (length(nb.bio.pcs) > 1 | nb.bio.pcs < 0){
        stop('The "nb.bio.pcs" must be a postive integer value.')
    }
    if (length(nb.uv.pcs) > 1 | nb.uv.pcs < 0){
        stop('The "nb.uv.pcs" must be a postive integer value.')
    }
    if (length(nb.pcs) > 1 | nb.pcs < 0){
        stop('The "nb.pcs" must be a postive integer value.')
    }
    if (isFALSE(is.logical(scale))) {
        stop('The "scale" must be "TRUE" or "FALSE.')
    }
    if (isFALSE(is.logical(center))) {
        stop('The "center" must be "TRUE" or "FALSE.')
    }
    if (isFALSE(is.logical(apply.log))) {
        stop('The "apply.log" must be "TRUE" or "FALSE.')
    }
    if (length(pseudo.count) > 1){
        stop('The "pseudo.count" must be 0 or a postive integer value.')
    }
    if (pseudo.count < 0){
        stop('The "pseudo.count" must be 0 or a postive integer value.')
    }
    if (isFALSE(is.logical(check.se.obj))) {
        stop('The "check.se.obj" must be "TRUE" or "FALSE.')
    }
    if (is.null(check.se.obj)) {
        if (isTRUE(sum(uv.variables %in% colnames(colData(se.obj))) != length(uv.variables))) {
            stop('All or some of "uv.variables" cannot be found in the SummarizedExperiment object.')
        } else if (!is.null(variables.to.assess.ncg)) {
            if (isTRUE(sum(variables.to.assess.ncg %in% colnames(colData(se.obj))) != length(variables.to.assess.ncg))) {
                stop('All or some of "variables.to.assess.ncg" cannot be found in the SummarizedExperiment object.')
            }
        }
    }
    if (isTRUE(ncg.selection.method == 'quantile') & isTRUE(use.rank)){
        if (is.null(bio.percentile) | is.null(uv.percentile))
            stop('The "bio.percentile" or "uv.percentile" cannot be NULL.')
        if (bio.percentile > 1 | bio.percentile < 0)
            stop('The "bio.percentile" must be a postive value between 0 and 1.')
        if (uv.percentile > 1 | uv.percentile < 0)
            stop('The "uv.percentile" must be a postive value between 0 and 1.')
    }
    # Selecting samples to use for down-stream analysis ####
    if (is.logical(samples.to.use)){
        if (length(samples.to.use) != ncol(se.obj)){
            stop('The length of the "samples.to.use" must be the same as the number of columns in the SummarizedExperiment object.')
        }
        se.obj.all <- se.obj
        se.obj <- se.obj[ , samples.to.use]
    }

    # Specifying cores
    if (is.null(nb.cores)){
        if (.Platform$OS.type == "windows") {
            nb.cores <- as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS", unset = 1))
        } else {
            # macOS or Unix
            nb.cores <- as.numeric(system("sysctl -n hw.ncpu", intern = TRUE)) - 1
            if (is.na(nb.cores) || length(nb.cores) == 0) {
                nb.cores <- 1
            }
        }
    }
    if (ncg.idenfitication.approach == 'LMM'){
        if (isTRUE(apply.log)){
            expr.data <- applyLog(
                se.obj = se.obj,
                assay.names = assay.name,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = TRUE
            )[[assay.name]]
        } else expr.data <- assay(x = se.obj, i = assay.name)

        sample.annotation <- as.data.frame(colData(se.obj))
        gene.var.part <- fitExtractVarPartModel(
            exprObj = expr.data,
            formula = form,
            data = sample.annotation,
            BPPARAM = MulticoreParam(workers = nb.cores)
            )
        new.form <- changeLmmFormula(
            form = form,
            out.put = "character",
            sub.set = NULL
            )
        uv.var <- rowSums(gene.var.part[, new.form, drop = FALSE])
        bio.var <- gene.var.part$Residuals
        all.lmm <- data.frame(
            gene = rownames(gene.var.part),
            uv = uv.var,
            bio = bio.var,
            ratio = uv.var / (bio.var + 1e-6)  # Avoid divide-by-zero
            )
        set.seed(4455)
        all.lmm$uv.rank <- rank(x = -all.lmm$uv, ties.method = 'random')
        set.seed(4455)
        all.lmm$bio.rank <- rank(x = all.lmm$bio, ties.method = 'random')
    }
    if (ncg.idenfitication.approach == 'LMM.UvAdjustment'){
        if (isTRUE(apply.log)){
            expr.data <- applyLog(
                se.obj = se.obj,
                assay.names = assay.name,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = TRUE
            )[[assay.name]]
        } else expr.data <- assay(x = se.obj, i = assay.name)
        sample.annotation <- as.data.frame(colData(se.obj))
        new.form <- changeLmmFormula(
            form = form,
            out.put = "character",
            sub.set = uv.variables
            )
        adjusted.data <- lmFit(
            object = expr.data,
            design = model.matrix(
                as.formula(paste0('~ ',  paste0(new.form, collapse = '+'))),
                sample.annotation)
            )
        adjusted.data <- residuals(adjusted.data, expr.data)
        if (isTRUE(regress.out.rle.med)){
            rle.med <- colMedians(adjusted.data - rowMedians(adjusted.data))
            lm.fit.data <- lmFit(
                object = adjusted.data,
                design = model.matrix(~rle.med)
            )
            adjusted.data <- residuals(lm.fit.data, adjusted.data)
        }
        sv.dec <- BiocSingular::runSVD(
            x = t(adjusted.data),
            k = nb.bio.pcs,
            BSPARAM = svd.bsparam,
            center = center,
            scale = scale
            )
        rownames(sv.dec$u) <- colnames(se.obj)
        rownames(sv.dec$v) <- row.names(se.obj)
        percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
        percentage <- sapply(
            seq_along(percentage),
            function(i) round(percentage [i], 1)
            )
        bio.variables <- c()
        for(i in 1:nb.bio.pcs){
            col.nam <- paste0('Bio.pc', i)
            sample.annotation[[col.nam]] <- sv.dec$u[ , i]
            bio.variables <- c(bio.variables, col.nam)
            }
        new.form <- paste(
            paste0(form, collapse = ''),
            '+',
            paste0(bio.variables, collapse = ' + ')
            )
        gene.var.part <- fitExtractVarPartModel(
            exprObj = expr.data,
            formula = new.form,
            data = sample.annotation,
            BPPARAM = MulticoreParam(workers = nb.cores)
            )
        uv.var <- changeLmmFormula(
            form = form,
            out.put = "character",
            sub.set = uv.variables
            )
        uv.var <- rowSums(gene.var.part[, uv.var, drop = FALSE])
        bio.var <- rowSums(gene.var.part[, bio.variables, drop = FALSE])
        all.lmm <- data.frame(
            gene = rownames(gene.var.part),
            uv = uv.var,
            bio = bio.var,
            ratio = uv.var / (bio.var + 1e-6)
        )
        set.seed(4455)
        all.lmm$uv.rank <- rank(x = -all.lmm$uv, ties.method = 'random')
        set.seed(4455)
        all.lmm$bio.rank <- rank(x = all.lmm$bio, ties.method = 'random')
    }
    if (ncg.idenfitication.approach == 'LMM.BioUvAdjustment'){
        if (isTRUE(apply.log)){
            expr.data <- applyLog(
                se.obj = se.obj,
                assay.names = assay.name,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = TRUE
            )[[assay.name]]
        } else expr.data <- assay(x = se.obj, i = assay.name)
        sample.annotation <- as.data.frame(colData(se.obj))
        new.form <- changeLmmFormula(
            form = form,
            out.put = "character",
            sub.set = uv.variables
        )
        adjusted.data <- lmFit(
            object = expr.data,
            design = model.matrix(
                as.formula(paste0('~ ',  paste0(new.form, collapse = '+'))),
                sample.annotation)
        )
        adjusted.data <- residuals(adjusted.data, expr.data)
        if (isTRUE(regress.out.rle.med)){
            rle.med <- colMedians(adjusted.data - rowMedians(adjusted.data))
            lm.fit.data <- lmFit(
                object = adjusted.data,
                design = model.matrix(~rle.med)
            )
            adjusted.data <- residuals(lm.fit.data, adjusted.data)
        }
        sv.dec <- BiocSingular::runSVD(
            x = t(adjusted.data),
            k = nb.bio.pcs,
            BSPARAM = svd.bsparam,
            center = center,
            scale = scale
        )
        rownames(sv.dec$u) <- colnames(se.obj)
        rownames(sv.dec$v) <- row.names(se.obj)
        percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
        percentage <- sapply(
            seq_along(percentage),
            function(i) round(percentage [i], 1)
            )
        bio.variables <- c()
        for(i in 1:nb.bio.pcs){
            col.nam <- paste0('Bio.pc', i)
            sample.annotation[[col.nam]] <- sv.dec$u[ , i]
            bio.variables <- c(bio.variables, col.nam)
        }
        new.form <- paste(
            paste0(form, collapse = ''),
            '+',
            paste0(bio.variables, collapse = ' + ')
            )
        gene.var.part <- fitExtractVarPartModel(
            exprObj = expr.data,
            formula = new.form,
            data = sample.annotation,
            BPPARAM = MulticoreParam(workers = nb.cores)
            )
        uv.var <- changeLmmFormula(
            form = form,
            out.put = "character",
            sub.set = uv.variables
            )
        uv.var <- rowSums(gene.var.part[, uv.var, drop = FALSE])
        bio.var <- rowSums(gene.var.part[, bio.variables, drop = FALSE])
        all.var.result <- data.frame(
            gene = rownames(gene.var.part),
            uv.var = uv.var,
            bio.var = bio.var,
            ratio = uv.var / (bio.var + 1e-6)  # Avoid divide-by-zero
            )
        all.var.result <- all.var.result[order(-all.var.result$ratio), ]
        all.var.result <- all.var.result[all.var.result$uv.var > .7 , ]
        nb.ncg.a <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
        selected.ncg <- row.names(se.obj) %in% row.names(all.var.result)[1:nb.ncg.a]

        sv.dec <- BiocSingular::runSVD(
            x = t(expr.data[selected.ncg, ]),
            k = nb.bio.pcs,
            BSPARAM = svd.bsparam,
            center = center,
            scale = scale
            )
        rownames(sv.dec$u) <- colnames(se.obj)
        rownames(sv.dec$v) <- row.names(se.obj)[selected.ncg]
        percentage <- sv.dec$d ^ 2 / sum(sv.dec$d ^ 2) * 100
        percentage <- sapply(
            seq_along(percentage),
            function(i) round(percentage [i], 1))
        uv.variables.a <- c()
        for(i in 1:nb.bio.pcs){
            col.nam <- paste0('Uv.pc', i)
            sample.annotation[[col.nam]] <- sv.dec$u[ , i]
            uv.variables.a <- c(uv.variables.a, col.nam)
        }
        new.form <- paste(
            '~',
            paste0(uv.variables.a, collapse = ' + '),
            ' + ',
            paste0(bio.variables, collapse = ' + ')
            )
        gene.var.part <- fitExtractVarPartModel(
            exprObj = expr.data,
            formula = new.form,
            data = sample.annotation,
            BPPARAM = MulticoreParam(workers = nb.cores)
            )
        uv.var <- rowSums(gene.var.part[, uv.variables.a, drop = FALSE])
        bio.var <- rowSums(gene.var.part[, bio.variables, drop = FALSE])
        all.lmm <- data.frame(
            gene = rownames(gene.var.part),
            uv = uv.var,
            bio = bio.var,
            ratio = uv.var / (bio.var + 1e-6)
            )
        set.seed(4455)
        all.lmm$uv.rank <- rank(x = -all.lmm$uv, ties.method = 'random')
        set.seed(4455)
        all.lmm$bio.rank <- rank(x = all.lmm$bio, ties.method = 'random')
    }
    # Selecting a set of genes as NCG ####
    printColoredMessage(
        message = '-- Selecting a set of genes as NCG:',
        color = 'magenta',
        verbose = verbose
        )
    ## Applying ratio approach ####
    if (isFALSE(use.rank)){
        all.lmm <- all.lmm[order(-all.lmm$ratio), ]
        if (!is.null(uv.percentile) & !is.null(bio.percentile)){
            all.lmm <- all.lmm[all.lmm$uv > quantile(x = all.lmm$uv, probs = uv.percentile) , ]
            if (quantile(x = all.lmm$bio, probs = bio.percentile) == 0){
                all.lmm <- all.lmm
            } else {
                all.lmm <- all.lmm[all.lmm$bio < quantile(x = all.lmm$bio, probs = bio.percentile) , ]
            }
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% row.names(all.lmm)[1:nb.ncg]
        } else if (!is.null(uv.percentile) & is.null(bio.percentile)){
            all.lmm <- all.lmm[all.lmm$uv > quantile(x = all.lmm$uv, probs = uv.percentile) , ]
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% all.lmm$gene[1:nb.ncg]
        } else if (is.null(uv.percentile) & !is.null(bio.percentile)){
            all.lmm <- all.lmm[all.lmm$bio < quantile(x = all.lmm$bio, probs = bio.percentile) , ]
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% all.lmm$gene[1:nb.ncg]
        } else if (is.null(uv.percentile) & is.null(bio.percentile)){
            nb.ncg <- round(x = nrow(se.obj) * nb.ncg, digits = 0)
            ncg.selected <- row.names(se.obj) %in% all.lmm$gene[1:nb.ncg]
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

    # Filtering selected negative control genes ######
    if (isTRUE(filter.ncgs)){
        printColoredMessage(
            message = '- Filtering the selected NCGs based on publicly available stable genes:',
            color = 'blue',
            verbose = verbose
        )
        if (common.hk == 'cancer'){
            common.hks <- singscore::getStableGenes(n_stable = nb.stable.genes)
            common.hks <- intersect(common.hks, row.names(se.obj)[ncg.selected])
            ncg.selected <- row.names(se.obj) %in% common.hks
        }
        if (common.hk == 'non.cancer'){
            common.hks <- row.names(se.obj)[rowData(se.obj)[[hk.group]]]
            common.hks <- intersect(common.hks, row.names(se.obj)[ncg.selected])
            ncg.selected <- row.names(se.obj) %in% common.hks
        }
    }

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
            all.variables <- c(uv.variables, bio.variables)
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
        }
        if (isTRUE(apply.log) & is.null(pseudo.count)){
            printColoredMessage(
                message = paste0(
                    '- Applying log2 on the ',
                    assay.name,
                    ' data.'),
                color = 'blue',
                verbose = verbose
                )
            expr.data <- log2(assay(x = se.obj, i = assay.name))
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
            x = t(expr.data[ncg.selected, ]),
            k = nb.pcs,
            BSPARAM =  bsparam(),
            center = center,
            scale = scale
            )
        pc.var <- (pca.data$d^2) / (ncol(expr.data) - 1)
        centered.data <- scale(
            t(expr.data[ncg.selected , ]),
            center = center,
            scale = scale
            )
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
            all.variables,
            function(x) {
                if (class(sample.annotation[[x]]) %in% c('numeric', 'integer')) {
                    rSquared <- sapply(
                        1:nb.pcs,
                        function(y) summary(lm(sample.annotation[[x]] ~ pca.data$u[, 1:y]))$r.squared)
                } else if (class(sample.annotation[[x]]) %in% c('factor', 'character')) {
                    catvar.dummies <- dummy_cols(sample.annotation[[x]])
                    catvar.dummies <- catvar.dummies[, c(2:ncol(catvar.dummies))]
                    cca.pcs <- sapply(
                        1:nb.pcs,
                        function(y) {
                            cca <- cancor(x = pca.data$u[, 1:y, drop = FALSE], y = catvar.dummies)
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
            ggtitle('') +
            scale_x_continuous(breaks = (1:nb.pcs),labels = c('PC1', paste0('PC1:', 2:nb.pcs))) +
            scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) +
            theme(
                panel.background = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 14),
                axis.title.y = element_text(size = 14),
                axis.text.x = element_text(size = 10, angle = 25, hjust = 1),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 16),
                strip.text.x = element_text(size = 14)
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
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 16),
                plot.title = element_text(size = 16),
                plot.margin = unit(c(0, 0, 3, 0), "pt")
                )
        assess.ncg.plot <- assess.ncg.plot / p.pca.percentage + plot_layout(heights = c(3, 1))
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
        if (!'un.supervised' %in% names(se.obj@metadata[['NCG']])){
            se.obj@metadata[['NCG']][['un.supervised']] <- list()
        }
        if (!ncg.group.name %in% names(se.obj@metadata[['NCG']][['un.supervised']])){
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]] <- list()
        }
        if (!ncg.set.name %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]] )){
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]] <- list()
        }
        if (!'ncg.set' %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]])){
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- list()
        }
        se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['ncg.set']] <- ncg.selected
        if (isTRUE(assess.ncg)){
            if (!'ranl.plot' %in% names(se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]])){
                se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- list()
            }
            se.obj@metadata[['NCG']][['un.supervised']][[ncg.group.name]][[ncg.set.name]][['assessment.plot']] <- assess.ncg.plot
        }
        printColoredMessage(
            message = '- The NCGs are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = '------------The findNcgUnSupervisedByLinearMixedModel function finished.',
            color = 'white',
            verbose = verbose
        )
    }
    return(se.obj)
}
