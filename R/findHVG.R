#' Find highly variable genes on RNA-seq data.
#'
#' @author Ramyar Molania
#'
#' @param se.obj description
#' @param assay.name description
#' @param approach description
#' @param nb.hvg description
#' @param uv.variables description
#' @param bio.variables description
#' @param group.name description
#' @param hvg.selection  TTT
#' @param form description
#' @param adjust.data description
#' @param nb.bio.pcs description
#' @param bio.percentile TTT
#' @param span description
#' @param clustering.method description
#' @param nb.clusters description
#' @param samples.to.use description
#' @param regress.out.variables description
#' @param regress.out.rle.med description
#' @param normalization description
#' @param apply.log description
#' @param pseudo.count description
#' @param center TTT
#' @param scale TTT
#' @param svd.bsparam TTT
#' @param nb.cores TTT
#' @param hvg.group.name description
#' @param hvg.set.name description
#' @param check.se.obj TTT
#' @param remove.na description
#' @param save.se.obj description
#' @param verbose description
#'
#' @importFrom SummarizedExperiment assays colData
#' @importFrom matrixStats rowVars
#' @importFrom parallel mclapply
#' @importFrom stats loess
#' @importFrom limma lmFit
#' @export

findHVG <- function(
        se.obj,
        assay.name,
        approach = 'vst',
        nb.hvg = .1,
        uv.variables = NULL,
        bio.variables = NULL,
        group.name = NULL,
        hvg.selection = 'intersect',
        form = NULL,
        adjust.data = FALSE,
        nb.bio.pcs = 3,
        bio.percentile = .7,
        span = .3,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        samples.to.use = 'all',
        normalization = NULL,
        regress.out.variables = NULL,
        regress.out.rle.med = FALSE,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        apply.log = TRUE,
        pseudo.count = 1,
        nb.cores = 1,
        hvg.group.name = NULL,
        hvg.set.name = NULL,
        check.se.obj = TRUE,
        remove.na = 'none',
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The findHVG function starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking inputs ####
    if (!is.vector(assay.name) | length(assay.name) > 1 | is.logical(assay.name) | assay.name == 'all'){
        stop('The "assay.name" must be a single assay name in the SummarizedExperiment object.')
    }
    if (!approach %in% c('mad', 'cv', 'var', 'mad.cv' ,'mean.var', 'auto', 'vst', 'lmm')){
        stop('The "approach" muat be one of "mad", "cv", "var", "mad.cv","mean.var", "auto", "vst", "lmm".')
    }
    if (nb.hvg >= 1 | nb.hvg <= 0){
        stop('The "nb.hvg" must be a positve value 0 < nb.hvg < 1.')
    }
    if (!is.null(group.name)){
        if (!hvg.selection %in% c('intersect', 'union', 'max', 'average', 'median')){
            stop('The "hvg.selection" must be one of the "intersect", "union", "max", "average" or "median".')
        }
    }
    if (approach == 'vst'){
        if (span >= 1 | span <= 0){
            stop('The "span" must be a positve value 0 < span < 1.')
        }
    }
    if (approach == 'lmm'){
        if (is.null(bio.variables) & is.null(uv.variables)){
            stop('Both "bio.variables" and "uv.variables" cannot be NULL when approach is set to "lmm".')
        }
        if (is.null(form)){
            stop('The "form" cannot be NULL.')
        }
    }
    if (!is.logical(adjust.data)) {
        stop('The "adjust.data" must be "TRUE" or "FALSE".')
    }
    if (!is.logical(regress.out.rle.med)) {
        stop('The "regress.out.rle.med" must be "TRUE" or "FALSE".')
    }
    if (!is.logical(apply.log)) {
        stop('The "apply.log" must be "TRUE" or "FALSE".')
    }
    if (pseudo.count < 0){
        stop('The "pseudo.count" must be a postive numeric value.')
    }
    if (!is.logical(check.se.obj)) {
        stop('The "check.se.obj" must be "TRUE" or "FALSE".')
    }
    if (!is.logical(save.se.obj)) {
        stop('The "save.se.obj" must be "TRUE" or "FALSE".')
    }
    if (!is.logical(verbose)) {
        stop('The "verbose" must be "TRUE" or "FALSE".')
    }

    # Checking samples to use ####
    ### Adding the results to the SummarizedExperiment object ####
    if (is.logical(samples.to.use)){
        se.obj.initial <- se.obj
        se.obj <- se.obj[ , samples.to.use]
    }

    # Assessing the SummarizedExperiment object ####
    if (isTRUE(check.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(bio.variables, uv.variables, group.name),
            remove.na = remove.na,
            verbose = verbose
            )
    }

    # Assessing the group variable ####
    if (!is.null(group.name)){
        if (is.numeric(se.obj[[group.name]])){
            initial.group <- se.obj[[group.name]]
            se.obj[[group.name]] <- groupContinuousVariable(
                se.obj = se.obj,
                variable = group,
                nb.clusters = nb.clusters,
                clustering.method = clustering.method,
                perfix = '_',
                plot.output = TRUE,
                verbose = verbose
            )
        }
    }

    # Applying data transformation and normalization ####
    printColoredMessage(
        message = '-- Applying data transformation and normalization:',
        color = 'magenta',
        verbose = verbose
        )
    expr.data <- preProcessData(
        se.obj = se.obj,
        assay.name = assay.name,
        normalization = normalization,
        regress.out.variables = regress.out.variables,
        regress.out.rle.med = FALSE,
        apply.log = apply.log,
        pseudo.count = pseudo.count,
        check.se.obj = FALSE,
        remove.na = 'none',
        verbose = verbose
        )
    ## Applying library size normalization ####
    ## Applying MAD approach
    if (approach == 'mad'){
        if (is.null(group.name)){
            if (isTRUE(regress.out.rle.med)){
                rle.med <- colMedians(expr.data - rowMedians(expr.data))
                lm.fit.data <- lmFit(
                    object = expr.data,
                    design = model.matrix(~rle.med)
                    )
                expr.data <- residuals(lm.fit.data, expr.data)
            }
            hvg <- matrixStats::rowMads(x = expr.data)
        }
        if (!is.null(group.name)){
            hvg <- sapply(
                unique(se.obj[[group.name]]),
                function(x){
                    temp.data <- expr.data[ , se.obj[[group.name]] == x, drop = FALSE]
                    if (isTRUE(regress.out.rle.med)){
                        rle.med <- colMedians(temp.data - rowMedians(temp.data))
                        lm.fit.data <- limma::lmFit(
                            object = temp.data,
                            design = model.matrix(~rle.med)
                        )
                        temp.data <- residuals(lm.fit.data, temp.data)
                    }
                    matrixStats::rowMads(temp.data)
                })
            row.names(hvg) <- row.names(se.obj)
            if (hvg.selection %in% c('max', 'average', 'median')){
                if (hvg.selection == 'max'){
                    hvg <- matrixStats::rowMaxs(hvg)
                }
                if (hvg.selection == 'average'){
                    hvg <- matrixStats::rowMeans2(hvg)
                }
                if (hvg.selection == 'median'){
                    hvg <- matrixStats::rowMedians(hvg)
                }
                hvg <- sort(hvg, decreasing = TRUE)
                hvg <- names(hvg)[1:(nrow(se.obj) * nb.hvg)]
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'union'){
                all.hvg <- lapply(
                    1:ncol(hvg),
                    function(x){
                        temp.hvg <- hvg[order(hvg[, x], decreasing = TRUE), ]
                        row.names(temp.hvg)[1:c((nrow(se.obj) * nb.hvg))]
                    })
                hvg <- unique(unlist(all.hvg))
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'intersect'){
                all.hvg <- lapply(
                    1:ncol(hvg),
                    function(x){
                        temp.hvg <- hvg[order(hvg[, x], decreasing = TRUE), ]
                        row.names(temp.hvg)[1:c((nrow(se.obj) * nb.hvg))]
                    })
                hvg <- Reduce(f = intersect, x = all.hvg)
                hvg <- row.names(se.obj) %in% hvg
            }
        }
    }
    ## Applying CV approach ####
    if (approach == 'cv'){
        if (is.null(group.name)){
            if (isTRUE(regress.out.rle.med)){
                rle.med <- colMedians(expr.data - rowMedians(expr.data))
                lm.fit.data <- lmFit(
                    object = expr.data,
                    design = model.matrix(~rle.med)
                )
                expr.data <- residuals(lm.fit.data, expr.data)
            }
            hvg <- Rfast::rowcvs(x = expr.data)
            hvg <- row.names(se.obj) %in% hvg
        }
        if (!is.null(group.name)){
            hvg <- sapply(
                unique(se.obj[[group.name]]),
                function(x){
                    temp.data <- expr.data[ , se.obj[[group.name]] == x, drop = FALSE]
                    if (isTRUE(regress.out.rle.med)){
                        rle.med <- colMedians(temp.data - rowMedians(temp.data))
                        lm.fit.data <- limma::lmFit(
                            object = temp.data,
                            design = model.matrix(~rle.med)
                        )
                        temp.data <- residuals(lm.fit.data, temp.data)
                    }
                    Rfast::rowcvs(temp.data)
                })
            row.names(hvg) <- row.names(se.obj)
            if (hvg.selection %in% c('max', 'average', 'median')){
                if (hvg.selection == 'max'){
                    hvg <- matrixStats::rowMaxs(hvg)
                }
                if (hvg.selection == 'average'){
                    hvg <- matrixStats::rowMeans2(hvg)
                }
                if (hvg.selection == 'median'){
                    hvg <- matrixStats::rowMedians(hvg)
                }
                hvg <- sort(hvg, decreasing = TRUE)
                hvg <- names(hvg)[1:(nrow(se.obj) * nb.hvg)]
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'union'){
                all.hvg <- lapply(
                    1:ncol(hvg),
                    function(x){
                        temp.hvg <- hvg[order(hvg[, x], decreasing = TRUE), ]
                        row.names(temp.hvg)[1:c((nrow(se.obj) * nb.hvg))]
                    })
                hvg <- unique(unlist(all.hvg))
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'intersect'){
                all.hvg <- lapply(
                    1:ncol(hvg),
                    function(x){
                        temp.hvg <- hvg[order(hvg[, x], decreasing = TRUE), ]
                        row.names(temp.hvg)[1:c((nrow(se.obj) * nb.hvg))]
                    })
                hvg <- Reduce(f = intersect, x = all.hvg)
                hvg <- row.names(se.obj) %in% hvg
            }
        }
    }
    ## Applying Var approach ####
    if (approach == 'var'){
        if (is.null(group.name)){
            if (isTRUE(regress.out.rle.med)){
                rle.med <- colMedians(expr.data - rowMedians(expr.data))
                lm.fit.data <- lmFit(
                    object = expr.data,
                    design = model.matrix(~rle.med)
                )
                expr.data <- residuals(lm.fit.data, expr.data)
            }
            hvg <- matrixStats::rowVars(x = expr.data)
            hvg <- row.names(se.obj) %in% hvg
        }
        if (!is.null(group.name)){
            hvg <- sapply(
                unique(se.obj[[group.name]]),
                function(x){
                    temp.data <- expr.data[ , se.obj[[group.name]] == x, drop = FALSE]
                    if (isTRUE(regress.out.rle.med)){
                        rle.med <- colMedians(temp.data - rowMedians(temp.data))
                        lm.fit.data <- limma::lmFit(
                            object = temp.data,
                            design = model.matrix(~rle.med)
                        )
                        temp.data <- residuals(lm.fit.data, temp.data)
                    }
                    matrixStats::rowVars(temp.data)
                })
            row.names(hvg) <- row.names(se.obj)
            if (hvg.selection %in% c('max', 'average', 'median')){
                if (hvg.selection == 'max'){
                    hvg <- matrixStats::rowMaxs(hvg)
                }
                if (hvg.selection == 'average'){
                    hvg <- matrixStats::rowMeans2(hvg)
                }
                if (hvg.selection == 'median'){
                    hvg <- matrixStats::rowMedians(hvg)
                }
                hvg <- sort(hvg, decreasing = TRUE)
                hvg <- names(hvg)[1:(nrow(se.obj) * nb.hvg)]
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'union'){
                all.hvg <- lapply(
                    1:ncol(hvg),
                    function(x){
                        temp.hvg <- hvg[order(hvg[, x], decreasing = TRUE), ]
                        row.names(temp.hvg)[1:c((nrow(se.obj) * nb.hvg))]
                    })
                hvg <- unique(unlist(all.hvg))
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'intersect'){
                all.hvg <- lapply(
                    1:ncol(hvg),
                    function(x){
                        temp.hvg <- hvg[order(hvg[, x], decreasing = TRUE), ]
                        row.names(temp.hvg)[1:c((nrow(se.obj) * nb.hvg))]
                    })
                hvg <- Reduce(f = intersect, x = all.hvg)
                hvg <- row.names(se.obj) %in% hvg
            }
        }
    }
    ## Applying MAD and CV approach ####
    if (approach == 'mad.cv'){
        if (is.null(group.name)){
            if (isTRUE(regress.out.rle.med)){
                rle.med <- colMedians(expr.data - rowMedians(expr.data))
                lm.fit.data <- lmFit(
                    object = expr.data,
                    design = model.matrix(~rle.med)
                )
                expr.data <- residuals(lm.fit.data, expr.data)
            }
            hvg <- data.frame(
                mad = matrixStats::rowMads(x = expr.data),
                cv = Rfast::rowcvs(x = expr.data)
            )
            hvg <- matrixStats::rowMaxs(x = as.matrix(hvg))
            hvg <- sort(hvg, decreasing = TRUE)
            hvg <- names(hvg)[1:(nrow(se.obj) * nb.hvg)]
            hvg <- row.names(se.obj) %in% hvg
        }
        if (!is.null(group.name)){
            hvg <- sapply(
                unique(se.obj[[group.name]]),
                function(x){
                    temp.data <- expr.data[ , se.obj[[group.name]] == x, drop = FALSE]
                    if (isTRUE(regress.out.rle.med)){
                        rle.med <- colMedians(temp.data - rowMedians(temp.data))
                        lm.fit.data <- limma::lmFit(
                            object = temp.data,
                            design = model.matrix(~rle.med)
                        )
                        temp.data <- residuals(lm.fit.data, temp.data)
                    }
                    hvg <- data.frame(
                        mad = matrixStats::rowMads(x = temp.data),
                        cv = Rfast::rowcvs(x = temp.data)
                        )
                    matrixStats::rowMaxs(x = as.matrix(hvg))
                })
            if (hvg.selection %in% c('max', 'average', 'median')){
                if (hvg.selection == 'max'){
                    hvg <- matrixStats::rowMaxs(hvg)
                }
                if (hvg.selection == 'average'){
                    hvg <- matrixStats::rowMeans2(hvg)
                }
                if (hvg.selection == 'median'){
                    hvg <- matrixStats::rowMedians(hvg)
                }
                hvg <- sort(hvg, decreasing = TRUE)
                hvg <- names(hvg)[1:(nrow(se.obj) * nb.hvg)]
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'union'){
                all.hvg <- lapply(
                    1:ncol(hvg),
                    function(x){
                        temp.hvg <- hvg[order(hvg[, x], decreasing = TRUE), ]
                        row.names(temp.hvg)[1:c((nrow(se.obj) * nb.hvg))]
                    })
                hvg <- unique(unlist(all.hvg))
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'intersect'){
                all.hvg <- lapply(
                    1:ncol(hvg),
                    function(x){
                        temp.hvg <- hvg[order(hvg[, x], decreasing = TRUE), ]
                        row.names(temp.hvg)[1:c((nrow(se.obj) * nb.hvg))]
                    })
                hvg <- Reduce(f = intersect, x = all.hvg)
                hvg <- row.names(se.obj) %in% hvg
            }
        }
    }
    ## Applying mean-variance approach ####
    if (approach == 'mean.var'){
        if (is.null(group.name)){
            if (isTRUE(regress.out.rle.med)){
                rle.med <- colMedians(expr.data - rowMedians(expr.data))
                lm.fit.data <- lmFit(
                    object = expr.data,
                    design = model.matrix(~rle.med)
                )
                expr.data <- residuals(lm.fit.data, expr.data)
            }
            hvg <- scran::modelGeneVar(expr.data)
            hvg <- scran::getTopHVGs(stats = hvg, n = (nrow(se.obj) * nb.hvg))
            hvg <- row.names(se.obj) %in% hvg
        }
        if (!is.null(group.name)){
            if (isTRUE(regress.out.rle.med)){
                rle.med <- colMedians(expr.data - rowMedians(expr.data))
                lm.fit.data <- lmFit(
                    object = expr.data,
                    design = model.matrix(~rle.med)
                )
                expr.data <- residuals(lm.fit.data, expr.data)
            }
            hvg <- scran::modelGeneVar(
                x = expr.data,
                block = se.obj[[variable]]
                )
            hvg <- scran::getTopHVGs(stats = hvg, n = (nrow(se.obj) * nb.hvg ))
            hvg <- row.names(se.obj) %in% hvg
        }
    }
    ## Applying VST approach ####
    if (approach == 'vst'){
        expr.data <- assay(x = se.obj, i = assay.name)
        if (is.null(group.name)){
            nfeatures <- nrow(expr.data)
            hvf.info <- data.frame(
                mean = numeric(nfeatures),
                variance = numeric(nfeatures),
                variance.expected = numeric(nfeatures),
                variance.standardized = numeric(nfeatures),
                variable = logical(nfeatures),
                rank = rep(NA_integer_, nfeatures)
                )
            hvf.info$mean <- rowMeans(expr.data)
            hvf.info$variance <- rowVars(expr.data)
            not.const <- hvf.info$variance > 0
            fit <- stats::loess(
                log10(variance) ~ log10(mean),
                data = hvf.info[not.const, , drop = TRUE],
                span = span
                )
            hvf.info$variance.expected[not.const] <- 10^fit$fitted
            feature.mean <- hvf.info$mean
            feature.mean[feature.mean == 0] <- 0.1
            feature.sd <- sqrt(hvf.info$variance.expected)
            cap <- sqrt(ncol(expr.data)) * feature.sd + feature.mean
            expr.data <- pmin(expr.data, matrix(rep(cap, times = ncol(expr.data)), nrow = nfeatures))
            data.standard <- (expr.data - feature.mean) / feature.sd
            hvf.info$variance.standardized <- rowVars(data.standard)
            vf <- order(hvf.info$variance.standardized, decreasing = TRUE)[1:c(nrow(se.obj) * nb.hvg)]
            hvf.info$variable[vf] <- TRUE
            hvf.info$rank[vf] <- seq_along(vf)
            if (!is.null(rownames(data))) {
                rownames(hvf.info) <- rownames(data)
            }
            hvg <- row.names(expr.data)[hvf.info$variable]
            hvg <- row.names(se.obj) %in% hvg
        }
        if (!is.null(group.name)){
            all.hvg <- lapply(
                unique(se.obj[[group.name]]),
                function(x){
                    temp.expr <- expr.data[ , se.obj[[group.name]] == x, drop = FALSE]
                    nfeatures <- nrow(temp.expr)
                    hvf.info <- data.frame(
                        mean = numeric(nfeatures),
                        variance = numeric(nfeatures),
                        variance.expected = numeric(nfeatures),
                        variance.standardized = numeric(nfeatures),
                        variable = logical(nfeatures),
                        rank = rep(NA_integer_, nfeatures)
                        )
                    hvf.info$mean <- rowMeans(temp.expr)
                    hvf.info$variance <- rowVars(temp.expr)
                    not.const <- hvf.info$variance > 0
                    fit <- stats::loess(
                        log10(variance) ~ log10(mean),
                        data = hvf.info[not.const, , drop = TRUE],
                        span = span
                        )
                    hvf.info$variance.expected[not.const] <- 10^fit$fitted
                    feature.mean <- hvf.info$mean
                    feature.mean[feature.mean == 0] <- 0.1
                    feature.sd <- sqrt(hvf.info$variance.expected)
                    cap <- sqrt(ncol(temp.expr)) * feature.sd + feature.mean
                    temp.expr <- pmin(
                        temp.expr,
                        matrix(rep(cap, times = ncol(temp.expr)), nrow = nfeatures)
                        )
                    data.standard <- (temp.expr - feature.mean) / feature.sd
                    hvf.info$variance.standardized <- rowVars(data.standard)
                    vf <- order(hvf.info$variance.standardized, decreasing = TRUE)[1:c(nrow(se.obj) * nb.hvg)]
                    hvf.info$variable[vf] <- TRUE
                    hvf.info$rank[vf] <- seq_along(vf)
                    if (!is.null(rownames(data))) {
                        rownames(hvf.info) <- rownames(data)
                        }
                    hvg <- row.names(temp.expr)[hvf.info$variable]
                    hvg
                })
            if (hvg.selection == 'union'){
                hvg <- unique(unlist(all.hvg))
                hvg <- row.names(se.obj) %in% hvg
            }
            if (hvg.selection == 'intersect'){
                hvg <- Reduce(f = intersect, all.hvg)
                hvg <- row.names(se.obj) %in% hvg
            }
        }
    }
    ## Applying LMM approach ####
    if (approach == 'lmm'){
        ### Both biological and unwanted variable is known ####
        if (!is.null(bio.variables) & !is.null(uv.variables)){
            expr.data <- applyLog(
                se.obj = se.obj,
                assay.names = assay.name,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = TRUE
                )[[assay.name]]
            sample.annotation <- as.data.frame(colData(se.obj))
            gene.var.part <- fitExtractVarPartModel(
                exprObj = expr.data,
                formula = form,
                data = sample.annotation,
                BPPARAM = MulticoreParam(workers = 14)
                )
            new.form <- changeLmmFormula(
                form = form,
                out.put = "character",
                sub.set = NULL
                )
            uv.var <- rowSums(gene.var.part[ , uv.variables, drop = FALSE])
            bio.var <- rowSums(gene.var.part[ , bio.variables, drop = FALSE])
            all.lmm <- data.frame(
                genes = rownames(gene.var.part),
                uv = uv.var,
                bio = bio.var,
                ratio = bio.var / (uv.var + 1e-6)  # Avoid divide-by-zero
                )
            all.lmm <- all.lmm[all.lmm$bio > quantile(x = all.lmm$bio, probs = bio.percentile) , ]
            all.lmm <- all.lmm[order(all.lmm$ratio, decreasing = TRUE) , ]
            hvg <- all.lmm$genes[1:c(nrow(se.obj) * nb.hvg)]
            hvg <- row.names(se.obj) %in% hvg
        }
        ## Only unwanted variable is known ####
        if (is.null(bio.variables) & !is.null(uv.variables)){
            sample.annotation <- as.data.frame(colData(se.obj))
            if (isTRUE(adjust.data)){
                adjusted.data <- lmFit(
                    object = expr.data,
                    design = model.matrix(
                        as.formula(paste0('~ ',  paste0(uv.variables, collapse = '+'))),
                        sample.annotation)
                    )
                adjusted.data <- residuals(adjusted.data, expr.data)
                sv.dec <- BiocSingular::runSVD(
                    x = t(adjusted.data),
                    k = nb.bio.pcs,
                    BSPARAM = svd.bsparam,
                    center = center,
                    scale = scale
                    )
                bio.variables <- c()
                for(i in 1:nb.bio.pcs){
                    col.nam <- paste0('Bio.pc', i)
                    sample.annotation[[col.nam]] <- sv.dec$u[ , i]
                    bio.variables <- c(bio.variables, col.nam)
                    }
                new.form <- paste(
                    '~',
                    paste0(bio.variables, collapse = ' + ')
                    )
                gene.var.part <- fitExtractVarPartModel(
                    exprObj = adjusted.data,
                    formula = new.form,
                    data = sample.annotation,
                    BPPARAM = MulticoreParam(workers = nb.cores)
                    )
                bio.var <- rowSums(gene.var.part[, bio.variables, drop = FALSE])
                all.lmm <- data.frame(
                    gene = rownames(gene.var.part),
                    bio = bio.var,
                    residuals = gene.var.part$Residuals
                    )
                all.lmm <- all.lmm[order(all.lmm$bio, decreasing = TRUE) , ]
                hvg <- row.names(all.lmm)[1:c(nrow(se.obj) * nb.hvg)]
                hvg <- row.names(se.obj) %in% hvg
            }
            if (isFALSE(adjust.data)){
                adjusted.data <- lmFit(
                    object = expr.data,
                    design = model.matrix(
                        as.formula(paste0('~ ',  paste0(uv.variables, collapse = '+'))),
                        sample.annotation)
                    )
                adjusted.data <- residuals(adjusted.data, expr.data)
                new.form <- changeLmmFormula(
                    form = form,
                    out.put = "character",
                    sub.set = uv.variables
                    )
                sv.dec <- BiocSingular::runSVD(
                    x = t(adjusted.data),
                    k = nb.bio.pcs,
                    BSPARAM = svd.bsparam,
                    center = center,
                    scale = scale
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
                    ratio = bio.var / (uv.var + 1e-6)
                    )
                all.lmm <- all.lmm[all.lmm$bio > quantile(x = all.lmm$bio, probs = bio.percentile) , ]
                all.lmm <- all.lmm[order(all.lmm$ratio, decreasing = TRUE) , ]
                hvg <- row.names(all.lmm)[1:c(nrow(se.obj) * nb.hvg)]
                hvg <- row.names(se.obj) %in% hvg
            }
        }
        ## Only biological variable is known ####
        if (!is.null(bio.variables) & is.null(uv.variables)){
            sample.annotation <- as.data.frame(colData(se.obj))
            gene.var.part <- fitExtractVarPartModel(
                exprObj = expr.data,
                formula = new.form,
                data = sample.annotation,
                BPPARAM = MulticoreParam(workers = nb.cores)
                )
            bio.var <- rowSums(gene.var.part[, bio.variables, drop = FALSE])
            all.lmm <- data.frame(
                gene = rownames(gene.var.part),
                bio = bio.var,
                residuals = gene.var.part$Residuals
                )
            all.lmm <- all.lmm[order(all.lmm$bio, decreasing = TRUE) , ]
            hvg <- row.names(all.lmm)[1:c(nrow(se.obj) * nb.hvg)]
            hvg <- row.names(se.obj) %in% hvg
        }
    }
    # Saving results ####
    ### Adding the results to the SummarizedExperiment object ####
    if (is.logical(samples.to.use)){
        se.obj <- se.obj.initial
    }
    #### Selecting output name ####
    if (is.null(hvg.group.name)){
        hvg.group.name <- paste0('hvg_', approach)
    }

    if (is.null(hvg.set.name)){
        if (approach %in% c('mad', 'var', 'mad.cv', 'mean.var', 'vst')){
            if (is.null(group.name)){
                hvg.set.name <- paste0(
                    sum(hvg),
                    '|acrossAllSamples|',
                    assay.name)
            }
            if (!is.null(group.name)){
                hvg.set.name <- paste0(
                    sum(hvg),
                    '|per',
                    group.name,
                    '|',
                    assay.name)
            }
        }
        if (approach == 'lmm'){
            hvg.set.name <- paste0(
                sum(hvg),
                '|Bio:',
                paste0(bio.variables, collapse = '&'),
                '|Uv:',
                paste0(uv.variables, collapse = '&'),
                '|',
                assay.name)
        }
    }

    #### Saving all the results ####
    if (isTRUE(save.se.obj)){
        if (length(se.obj@metadata$HVG) == 0 ) {
            se.obj@metadata[['HVG']] <- list()
        }
        if (!hvg.group.name %in% names(se.obj@metadata[['HVG']])){
            se.obj@metadata[['HVG']][[hvg.group.name]] <- list()
        }
        if (!hvg.set.name %in% names(se.obj@metadata[['HVG']][[hvg.group.name]] )){
            se.obj@metadata[['HVG']][[hvg.group.name]][[hvg.set.name]] <- list()
        }
        if (!'hvg.set' %in% names(se.obj@metadata[['HVG']][[hvg.group.name]][[hvg.set.name]])){
            se.obj@metadata[['HVG']][[hvg.group.name]][[hvg.set.name]][['hvg.set']] <- list()
        }
        se.obj@metadata[['HVG']][[hvg.group.name]][[hvg.set.name]][['hvg.set']] <- hvg

        printColoredMessage(
            message = '- The hvg are saved to metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        printColoredMessage(
            message = '------------The findHVG function finished.',
            color = 'white',
            verbose = verbose
            )
        return(se.obj)
    }
    ### Export results as logical vector ####
    if (isFALSE(save.se.obj)){
        printColoredMessage(
            message = '-- The set of HVG is outpputed as a logical vector.',
            color = 'magenta',
            verbose = verbose
        )
        printColoredMessage(
            message = '------------The findHVG function finished.',
            color = 'white',
            verbose = verbose
        )
        return(hvg)
    }
}

