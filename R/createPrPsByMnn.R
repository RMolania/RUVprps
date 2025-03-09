#' Creates PRPS sets using k and mutual nearest neighbors in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function uses the k and mutual nearest neighbors approaches to create PRPS in the RNA-seq data. This function
#' can be used in situation that the biological variation are entirely unknown. The function applies the 'findKnn' function
#' to find similar samples per batch and then average them to create pseudo-samples. Then, function uses the 'findMnn' to
#' match up pseudo samples across batches to create pseudo-replicates.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character string indicating the name of the assay in the SummarizedExperiment object.
#' This assay will be used to first find k nearest neighbors and then mutual nearest neighbors data. This data must be
#' the one that will be used for the RUV-III normalization.
#' @param main.uv.variable Character. Indicates the name of a column in the sample annotation of the SummarizedExperiment
#' object. The 'main.uv.variable' can be either categorical or continuous. If 'main.uv.variable' is a continuous variable,
#' this will be divided into 'nb.clusters' groups using the 'clustering.method'.
#' @param clustering.method Character. A character string indicating the choice of clustering method for grouping the
#' 'main.uv.variable' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default
#' is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'main.uv.variable' is a
#' continuous variable. The default is 3.
#' @param filter.prps.sets Logical. If 'TRUE', the number of PRPS sets across each pair of batches will be filtered if
#' they are higher than the 'max.prps.sets' value. The default is 'TRUE'. A high number of PRPS sets will just increase
#' the computational time for the RUV-III normalization.
#' @param max.prps.sets Numeric. A numeric value specifying the maximum number for PRPS sets across each pair of batches.
#' The default is set to 10.
#' @param select.extreme.groups Character. Specifies whether to select extreme groups for analysis.
#' @param other.uv.variables Character. A character string or character vector representing the label of unwanted
#' variable(s), such as library size, etc., within the sample annotation (colData) of the SummarizedExperiment object.
#' This can comprise a vector containing either categorical, continuous, or a combination of both variables. These variables
#' will be considered when generating PRPS sets for the "main.uv.variable". This can help avoid potential contamination
#' from the "other.uv.variables" in the PRPS data of the "main.uv.variable". The default is set to 'NULL'.
#' @param other.uv.clustering.method Character. A character string indicating which clustering method should be used to
#' group each continuous unwanted variable, if specified in the "other.uv.variables". The options include 'kmeans', 'cut',
#' and 'quantile'. The default is 'kmeans'. Refer to the createHomogeneousUVGroups() function for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous unwanted
#' variable. The default is set to 3.
#' @param min.sample.for.ps Numeric. Specifies the minimum number of samples required for the PRPS sets.
#' @param min.mnn Numeric. The maximum number of nearest neighbors to compute. The default is set to 3.
#' @param min.batch.to.cover Numeric. The maximum number of nearest neighbors to compute. The default is set to 3.
#' @param hvg Character vector. A vector of the names of the highly variable genes. These genes will be used to find the
#' anchors across the batches. The default is NULL.
#' @param normalization Character. Indicates which normalization methods should be applied before finding the knn. The
#' default is 'cpm'. If set to NULL, no normalization will be applied.
#' @param apply.cosine.norm Logical. Specifies whether to apply cosine normalization.
#' @param regress.out.variables Character. Indicates the column names that contain biological variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by unwanted variation. The default is NULL, indicating that regression will not be applied.
#' @param check.prps.connectedness Logical. Indicates whether to assess the connectedness between the PRPS sets.
#' The default is set to TRUE. See the details for more information.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' Please note that any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value to be added as a pseudo count to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param mnn.bpparam BiocParallelParam. A BiocParallelParam object specifying how parallelization should be performed.
#' The default is SerialParam(). Refer to the 'findMutualNN' function from the BiocNeighbors R package for more details.
#' @param mnn.nbparam BiocParallelParam. A BiocParallelParam object specifying how parallelization should be performed to
#' find MNN. The default is KmknnParam(). Refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Character. Specifies whether to remove NA or missing values from the assays or not. The options are
#' 'assays' and 'none'. The default is "assays", so all the NA or missing values from the assay(s) will be removed before
#' computing RLE. See the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object
#' or to output the result as a list. By default, it is set to TRUE.
#' @param plot.output Logical. Specifies whether to generate and output plots.
#' @param output.name Character. A character string specifying the name of the output file. If set to 'NULL', the function
#' will select a name based on "paste0(uv.variable, '|', 'anchor', '|', assay.name)".
#' @param prps.group Character. A character string specifying the name of the output file. If set to 'NULL', the function
#' will select a name based on "paste0('prps_mnn_', uv.variable)".
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return The SummarizedExperiment object that contains all the PRPS data, knn, mnn, and plot results in the metadata, or
#' a list of the results.

#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocNeighbors findMutualNN
#' @importFrom batchelor cosineNorm
#' @importFrom RANN nn2
#' @export

createPrPsByMnn <- function(
        se.obj,
        assay.name,
        main.uv.variable,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        select.extereme.groups = FALSE,
        other.uv.variables = NULL,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        min.sample.for.ps = 3,
        min.mnn = 3,
        min.batches.to.cover = 'all',
        hvg = NULL,
        normalization = 'CPM',
        apply.cosine.norm = FALSE,
        regress.out.variables = NULL,
        check.prps.connectedness = TRUE,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        plot.output = TRUE,
        output.name = NULL,
        prps.group = NULL,
        verbose = TRUE
        ){
    printColoredMessage(
        message = '------------The createPrPsByMnn function starts:',
        color = 'white',
        verbose = verbose
        )
    # Checking the function inputs ####
    if (is.null(assay.name) | is.logical(assay.name)) {
        stop('The "assay.name" cannot be empty or logical.')
    }
    if (length(assay.name) > 1 | assay.name == 'all') {
        stop('The "assay.name" must be an assay name in the SummarizedExperiment object.')
    }
    if (isFALSE(assess.se.obj)){
        if (!assay.name %in% names(assays(se.obj))){
            stop('The "assay.name" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (length(main.uv.variable) > 1) {
        stop('The "main.uv.variable" must a categorical or continuous variable in the SummarizedExperiment object.')
    }
    if (is.null(main.uv.variable) | is.logical(main.uv.variable)) {
        stop('The "main.uv.variable" cannot be empty or logical(TRUE or FALSE).')
    }
    if (isFALSE(assess.se.obj)){
        if (!main.uv.variable %in% colnames(colData(se.obj))){
            stop('The "main.uv.variable" cannot be found in the SummarizedExperiment object.')
        }
    }
    if (is.numeric(colData(se.obj))){
        if (var(colData(se.obj)[[main.uv.variable]]) == 0){
            stop('The variance of the "main.uv.variable" is 0. No need to create PRPS for this variable.')
        }
    }
    if (!is.null(other.uv.variables)){
        if (main.uv.variable %in% other.uv.variables){
            stop('The "main.uv.variable" must not be in the "other.uv.variables".')
        }
        if (isFALSE(assess.se.obj)){
            if (sum(other.uv.variables %in% colnames(colData(se.obj))) != length(other.uv.variables)){
                stop('All or some of the "other.uv.variables" cannot be found in the SummarizedExperiment object.')
            }
        }
    }
    if (!is.logical(filter.prps.sets)){
        stop('The "filter.prps.sets" must be logical (TRUE or FALSE)')
    }
    if (isTRUE(filter.prps.sets)){
        if (!is.numeric(max.prps.sets) | max.prps.sets < 0){
            stop('The "max.prps.sets" must be postive numeric value.')
        }
    }
    if (min.sample.for.ps <= 1) {
        stop('The minimum value for the "min.sample.for.ps" is 2.')
    }
    if (!is.logical(apply.log)){
        stop('The "apply.log" must be logical (TRUE or FALSE).')
    }
    if (!is.logical(check.prps.connectedness)){
        stop('The "check.prps.connectedness" must be logical.')
    }
    if (!is.logical(assess.se.obj)){
        stop('The "assess.se.obj" must be logical (TRUE or FALSE).')
    }
    if (isTRUE(apply.log)){
        if (pseudo.count < 0){
            stop('The value for "pseudo.count" can not be negative.')
        }
    }
    if (!is.null(regress.out.variables)){
        if (isFALSE(assess.se.obj)){
            if (sum(regress.out.variables %in% colnames(colData(se.obj))) != length(regress.out.variables)){
                stop('All or some of the "regress.out.variables" cannot be found in the SummarizedExperiment object.')
            }
        }
    }
    if (!is.logical(plot.output)){
        stop('The "plot.output" must be logical (TRUE or FALSE).')
    }
    if (!is.logical(save.se.obj)){
        stop('The "save.se.obj" must be logical (TRUE or FALSE).')
    }
    if (is.logical(output.name)){
        stop('The "output.name" must be a character or NULL.')
    }
    if (is.logical(prps.group)){
        stop('The "prps.group" must be a character or NULL.')
    }
    if (!is.logical(verbose)){
        stop('The "verbose" must be logical (TRUE or FALSE).')
    }

    # Assessing the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(main.uv.variable, other.uv.variables, regress.out.variables),
            remove.na = remove.na,
            verbose = verbose
        )
    }

    # Applying data normalization and transformation and regression ####
    printColoredMessage(
        message = '-- Applying data normalization, transformation and regression:',
        color = 'magenta',
        verbose = verbose
        )

    ## normalization ####
    if (!is.null(normalization) & is.null(regress.out.variables)) {
        printColoredMessage(
            message = paste0(
                '- Applying the ',
                normalization,
                ' normalization on the data before finding MNN.'),
            color = 'blue',
            verbose = verbose
        )
        norm.data <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = TRUE
        )
    }
    ## normalization and regression ####
    if (!is.null(normalization) & !is.null(regress.out.variables)) {
        printColoredMessage(
            message = paste0(
                '- Applying the ',
                normalization,
                ' normalization and then regressing out the ',
                paste0(regress.out.variables, collapse = '&'),
                ' variable(s) from the data before finding MNN.'),
            color = 'blue',
            verbose = verbose
        )
        ### normalization ####
        norm.data <- applyOtherNormalizations(
            se.obj = se.obj,
            assay.name = assay.name,
            method = normalization,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = FALSE
        )
        ## regression ####
        sample.info <- as.data.frame(colData(se.obj)[regress.out.variables])
        colnames(sample.info) <- paste0('variable', 1:ncol(sample.info))
        regress.out.variables <- colnames(sample.info)
        norm.data <- t(norm.data)
        lm.formua <- paste('sample.info', regress.out.variables, sep = '$')
        norm.data <- lm(as.formula(paste(
            'norm.data',
            paste0(lm.formua, collapse = '+') ,
            sep = '~'
        )))
        norm.data <- t(norm.data$residuals)
        colnames(norm.data) <- colnames(norm.data)
        row.names(norm.data) <- row.names(norm.data)
    }
    ## regression ####
    if (is.null(normalization) & !is.null(regress.out.variables)){
        if(isTRUE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- Applying log2 transformation and then regressing out the ',
                    paste0(regress.out.variables, collapse = '&'),
                    ' variable(s) from the data before finding MNN. '),
                color = 'blue',
                verbose = verbose
            )
            if(!is.null(pseudo.count)){
                norm.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
            } else {
                norm.data <- log2(assay(x = se.obj, i = assay.name))
            }

        } else if (isFALSE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- regressing out the ',
                    paste0(regress.out.variables, collapse = '&'),
                    ' variable(s) from the data before finding MNN.'),
                color = 'blue',
                verbose = verbose
            )
            norm.data <- assay(x = se.obj, i = assay)
        }
        sample.info <- as.data.frame(colData(se.obj))
        norm.data <- t(norm.data)
        lm.formua <- paste('sample.info', regress.out.variables, sep = '$')
        norm.data <- lm(as.formula(paste(
            'norm.data',
            paste0(lm.formua, collapse = '+') ,
            sep = '~'
        )))
        norm.data <- t(norm.data$residuals)
        colnames(norm.data) <- colnames(norm.data)
        row.names(norm.data) <- row.names(norm.data)
        norm.data

    }
    ## log transformation ####
    if (is.null(normalization) & is.null(regress.out.variables)) {
        if (isTRUE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- Applying the log2 transformation on the data before finding MNN.'),
                color = 'blue',
                verbose = verbose
            )
            if (!is.null(pseudo.count)){
                norm.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
            } else {
                norm.data <- log2(assay(x = se.obj, i = assay.name))
            }

        } else if (isFALSE(apply.log)){
            printColoredMessage(
                message = paste0(
                    '- no library size normalization and transformation is applied on the data before finding MNN.'),
                color = 'blue',
                verbose = verbose
            )
            norm.data <- assay(x = se.obj, i = assay.name)
        }
    }

    # Assessing and grouping the main unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the main unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    initial.variable <- se.obj[[main.uv.variable]]
    if (is.numeric(initial.variable)){
        se.obj[[main.uv.variable]] <- groupContinuousVariable(
            se.obj = se.obj,
            variable = main.uv.variable,
            nb.clusters = nb.clusters,
            clustering.method = clustering.method,
            perfix = '_',
            verbose = verbose
        )
        if (isTRUE(select.extereme.groups)){
            max.group <- se.obj[[main.uv.variable]][initial.variable == max(initial.variable)]
            min.group <- se.obj[[main.uv.variable]][initial.variable == min(initial.variable)]
            selected.samples <- se.obj[[main.uv.variable]] %in% c(max.group, min.group)
            se.obj <- se.obj[ , selected.samples]
            se.obj[[main.uv.variable]] <- droplevels(se.obj[[main.uv.variable]])
            initial.variable <- initial.variable[selected.samples]
        }

    }
    if (!is.numeric(initial.variable)){
        if (length(unique(initial.variable)) == 1){
            stop('To create MNN, the "main.uv.variable" must have at least two groups/levels.')
        } else if (length(unique(initial.variable)) > 1){
            printColoredMessage(
                message = paste0(
                    '- The "',
                    main.uv.variable,
                    '" is a categorical variable with ',
                    length(unique(se.obj[[main.uv.variable]])),
                    ' levels.'),
                color = 'blue',
                verbose = verbose
            )
            se.obj[[main.uv.variable]] <- factor(x = se.obj[[main.uv.variable]])
        }
    }

    # Create PRPS data ####
    if (!is.null(other.uv.variables)){
        ## other uv variable is TRUE ####
        # assessing and grouping the other unwanted variable ####
        printColoredMessage(
            message = '- Assessing and grouping the other specified unwanted variable(s):',
            color = 'magenta',
            verbose = verbose
            )
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = other.uv.variables,
            nb.clusters = nb.other.uv.clusters,
            clustering.method = other.uv.clustering.method,
            assess.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
        all.uv.groups <- data.frame(
            main.uv = se.obj[[main.uv.variable]],
            other.uv = homo.uv.groups
            )
        covered.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                subgroups.size <- findRepeatingPatterns(
                    vec = all.uv.groups$main.uv[all.uv.groups$other.uv == x],
                    n.repeat = max(min.sample.for.ps, min.mnn)
                )
            })
        names(covered.batches) <- unique(all.uv.groups$other.uv)
        covered.batches.table <- as.data.frame(
            table(all.uv.groups$main.uv, all.uv.groups$other.uv)
            )
        covered.batches.table$selected <- covered.batches.table$Freq >= max(min.sample.for.ps, min.mnn)
        covered.batches.table <- ggplot(covered.batches.table, aes(x = Var2, y = Var1, color = selected)) +
            geom_point() +
            geom_text(aes(label = Freq , hjust = 0.5, vjust = 0.5), color = 'black') +
            xlab('Homogeneous groups (other unwanted variables)') +
            ylab('Main unwanted variable') +
            theme_bw() +
            theme(
                legend.key = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.y = element_text(size = 14),
                axis.text.x = element_text(size = 14, angle = 90, vjust = 1, hjust = 1),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 18),
                strip.text.y = element_text(size = 0)
            )
        if (isTRUE(plot.output)) print(covered.batches.table)
        selected.covered.batches <- lapply(
            1:length(covered.batches),
            function(x) length(covered.batches[[x]])
            )
        if (sum(selected.covered.batches == 1) == length(selected.covered.batches)){
            stop(paste0(
                ' Non of the sample groups with respect to the other unwanted variables that have at least ',
                max(min.sample.for.ps, min.mnn),
                ' samples across at least two sub-groups of the ',
                main.uv.variable,
                ' variable.'))
        }
        if (sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) == 0){
            printColoredMessage(
                message = paste0(
                    '- Non of the sample groups with respect to the other unwanted variables have at least ',
                    max(min.sample.for.ps, min.mnn),
                    ' samples across all the sub-groups of the "',
                    main.uv.variable,
                    '" variable.'),
                color = 'blue',
                verbose = verbose
            )
            if (isFALSE(check.prps.connectedness)){
                printColoredMessage(
                    message = '-- We recommend applying the "check.prps.connectedness"',
                    color = 'red',
                    verbose = verbose
                )
            } else if (isTRUE(check.prps.connectedness)){
                checkPRPSconnectedness(
                    data.input = table(all.uv.groups$main.uv, all.uv.groups$other.uv),
                    min.samples = c(min.mnn, min.sample.for.ps),
                    batch.name = main.uv.variable,
                    verbose = verbose
                    )
            }
        }
        if (sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) > 0  ){
            printColoredMessage(
                message = paste0(
                    '-- There are ',
                    sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) ,
                    ' groups with respect to the other unwanted variables that have at least ',
                    max(min.sample.for.ps, min.mnn),
                    ' samples across all sub-groups of the ',
                    main.uv.variable,
                    ' variable.'),
                color = 'blue',
                verbose = verbose
            )
        }

        # find mutual nearest neighbor ####
        printColoredMessage(
            message = '-- Creating PRPS data across all pairs of the subgroups:',
            color = 'magenta',
            verbose = verbose
        )
        all.possible.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                possible.batch <- findRepeatingPatterns(
                    vec = all.uv.groups[all.uv.groups$other.uv == x, ]$main.uv,
                    n.repeat = max(min.sample.for.ps, min.mnn)
                )
                if (length(possible.batch) > 1){
                    combn(x = possible.batch , m = 2)
                } else NA

            })
        names(all.possible.batches) <- unique(all.uv.groups$other.uv)
        all.possible.batches <- all.possible.batches[!is.na(all.possible.batches)]
        prps.data <- assay(x = se.obj, i = assay.name)
        sample.annotation <- as.data.frame(colData(x = se.obj))

        all.prps.data <- lapply(
            1:length(all.possible.batches),
            function(x){
                sub.sample.annotation <- sample.annotation[all.uv.groups$other.uv == names(all.possible.batches)[x] , ]
                pairs.batch <- all.possible.batches[[x]]
                sub.prps.data <- lapply(
                    1:ncol(pairs.batch),
                    function(y){
                        printColoredMessage(
                            message = paste0(
                                '- Creating PRPS data between the "',
                                pairs.batch[1 , y],
                                '" and "' ,
                                pairs.batch[2 , y],
                                '" subgroups within the "',
                                names(all.possible.batches)[x],
                                ' " subgroup.'),
                            color = 'orange',
                            verbose = verbose
                        )
                        ## sample annotation ####
                        sample.annot.a <- sub.sample.annotation[sub.sample.annotation[[main.uv.variable]] == pairs.batch[1 , y] , , drop = FALSE]
                        sample.annot.b <- sub.sample.annotation[sub.sample.annotation[[main.uv.variable]] == pairs.batch[2 , y] , , drop = FALSE]

                        ## highly variable genes ####
                        if (is.null(hvg)){
                            printColoredMessage(
                                message = '* highly variable are not specified, then using all genes.',
                                color = 'blue',
                                verbose = verbose
                            )
                            data.a <- norm.data[ , row.names(sample.annot.a)]
                            data.b <- norm.data[ , row.names(sample.annot.b)]
                        }
                        if (!is.null(hvg)){
                            printColoredMessage(
                                message = '* using the specified highly variable genes.',
                                color = 'blue',
                                verbose = verbose
                            )
                            data.a <- norm.data[hvg , row.names(sample.annot.a)]
                            data.b <- norm.data[hvg , row.names(sample.annot.b)]
                        }

                        ## cosine normalization ####
                        if (isTRUE(apply.cosine.norm)){
                            printColoredMessage(
                                message = '- Applying cosine normalization on the data:',
                                color = 'blue',
                                verbose = verbose
                            )
                            data.a <- cosineNorm(x = data.a, mode = 'matrix')
                            data.b <- cosineNorm(x = data.b, mode = 'matrix')
                        }

                        ## finding mnn for data b ####
                        printColoredMessage(
                            message = paste0(
                                '* Finding ',
                                min.sample.for.ps,
                                ' nearest neighbours for each sample of the "',
                                pairs.batch[2 , y],
                                '" using the "' ,
                                pairs.batch[1 , y],
                                '" subgroup:'),
                            color = 'blue',
                            verbose = verbose
                        )
                        knn.data.a.b <- RANN::nn2(
                            data = t(data.a),
                            query = t(data.b),
                            k = min.sample.for.ps
                        )
                        ## obtaining knn index ####
                        printColoredMessage(
                            message = '** obtaining the knn indexs.',
                            color = 'blue',
                            verbose = verbose
                        )
                        knn.data.a.b.index <- knn.data.a.b$nn.idx
                        colnames(knn.data.a.b.index) <- paste0('knn', seq(min.sample.for.ps))
                        row.names(knn.data.a.b.index) <- c(1:ncol(data.b))

                        ## obtaining the distance ####
                        printColoredMessage(
                            message = '** obtaining the distances.',
                            color = 'blue',
                            verbose = verbose
                        )
                        knn.data.a.b.distance <- as.data.frame(knn.data.a.b$nn.dists)
                        colnames(knn.data.a.b.distance) <- paste0('knn', seq(min.sample.for.ps))
                        row.names(knn.data.a.b.distance) <- c(1:ncol(data.b))
                        knn.data.a.b.distance$aver.dist <- rowMeans(knn.data.a.b.distance)

                        ## finding knn for data a ####
                        printColoredMessage(
                            message = paste0(
                                '* Finding ',
                                min.sample.for.ps,
                                ' nearest neighbours for each sample of the "',
                                pairs.batch[1 , y],
                                '" using "' ,
                                pairs.batch[2 , y],
                                '" subgroups:'),
                            color = 'blue',
                            verbose = verbose
                        )
                        knn.data.b.a <- RANN::nn2(
                            data = t(data.b),
                            query = t(data.a),
                            k = min.sample.for.ps
                        )
                        printColoredMessage(
                            message = '** Obtaining the knn indexs.',
                            color = 'blue',
                            verbose = verbose
                        )
                        knn.data.b.a.index <- knn.data.b.a$nn.idx
                        colnames(knn.data.b.a.index) <- paste0('knn', seq(min.sample.for.ps))
                        row.names(knn.data.b.a.index) <- c(1:ncol(data.a))
                        ### distance
                        printColoredMessage(
                            message = '** Obtaining the distances.',
                            color = 'blue',
                            verbose = verbose
                        )
                        knn.data.b.a.distance <- as.data.frame(knn.data.b.a$nn.dists)
                        colnames(knn.data.b.a.distance) <- paste0('knn', seq(min.sample.for.ps))
                        row.names(knn.data.b.a.distance) <- c(1:ncol(data.a))
                        knn.data.b.a.distance$aver.dist <- rowMeans(knn.data.b.a.distance)

                        ## final mnn ####
                        printColoredMessage(
                            message = paste0(
                                '* Finding MNN between the "',
                                pairs.batch[1 , y],
                                '" using "' ,
                                pairs.batch[2 , y],
                                '" subgroups:'),
                            color = 'blue',
                            verbose = verbose
                        )
                        all.mnn <- BiocNeighbors::findMutualNN(
                            data1 = t(data.a),
                            data2 = t(data.b),
                            k1 = min.mnn,
                            BNPARAM = mnn.nbparam,
                            BPPARAM = mnn.bpparam
                        )
                        printColoredMessage(
                            message = paste0(
                                '** ',
                                length(all.mnn$first),
                                ' MNN are found.'),
                            color = 'blue',
                            verbose = verbose
                        )
                        ## create prps index ####
                        printColoredMessage(
                            message = paste0('* creating PRPS indexs and score.'),
                            color = 'blue',
                            verbose = verbose
                        )
                        all.prps.index <- lapply(
                            1:length(all.mnn$first),
                            function(z){
                                sample.to.ave.b <- knn.data.b.a.index[all.mnn$first[z] , ]
                                sample.to.ave.a <- knn.data.a.b.index[all.mnn$second[z] , ]
                                prps.index <- data.frame(
                                    set.a = sample.to.ave.a,
                                    set.b = sample.to.ave.b
                                )
                                samples.ids <- data.frame(
                                    set.a = row.names(sample.annot.a)[sample.to.ave.a],
                                    set.b = row.names(sample.annot.b)[sample.to.ave.b]
                                )
                                avre.dist.a <- knn.data.b.a.distance[all.mnn$first[z] , 'aver.dist' ]
                                avre.dist.b <- knn.data.a.b.distance[all.mnn$second[z] , 'aver.dist' ]
                                list(
                                    prps.index = prps.index,
                                    aver.dist = c(avre.dist.a + avre.dist.b)/2,
                                    samples.ids = samples.ids
                                )
                            })
                        names(all.prps.index) <- paste0('prps.set', 1:length(all.mnn$first))
                        if (isTRUE(filter.prps.sets)){
                            printColoredMessage(
                                message = '* Filtering the number of PRPS sets.',
                                color = 'blue',
                                verbose = verbose
                            )
                            if (length(all.prps.index) >=  max.prps.sets){
                                printColoredMessage(
                                    message = '* The number of the PRPS set is larger than the specified "max.prps.sets", filtering PRPS sets.',
                                    color = 'blue',
                                    verbose = verbose
                                )
                                aver.dists <- sapply(
                                    1:length(all.prps.index),
                                    function(p) all.prps.index[[p]]$aver.dist
                                )
                                names(aver.dists) <- 1:length(all.prps.index)
                                aver.dists <- aver.dists[order(aver.dists, decreasing = FALSE)]
                                select.prps <- names(aver.dists)[1:max.prps.sets]
                                all.prps.index <- all.prps.index[as.numeric(select.prps)]
                                printColoredMessage(
                                    message = paste0(
                                        '** ',
                                        length(all.prps.index),
                                        ' PRPS set are kept.'),
                                    color = 'blue',
                                    verbose = verbose
                                )
                            }
                        }
                        ## prps data
                        # finding PRPS sets ####
                        printColoredMessage(
                            message = '* Creating the PRPS data:',
                            color = 'blue',
                            verbose = verbose
                        )
                        data.a <- prps.data[ , colnames(data.a)]
                        data.b <- prps.data[ , colnames(data.b)]
                        if (isTRUE(apply.log)){
                            printColoredMessage(
                                message = '** Applying log transformation on the data before creating PRPS. ',
                                color = 'blue',
                                verbose = verbose
                            )
                            data.a <- log2(data.a + pseudo.count)
                            data.b <- log2(data.b + pseudo.count)
                        } else {
                            printColoredMessage(
                                message = '** the data will be used without any transformation. ',
                                color = 'blue',
                                verbose = verbose
                            )
                        }
                        all.prps.data <- lapply(
                            1:length(all.prps.index),
                            function(a){
                                prps.set <- all.prps.index[[a]]
                                ps.a <- rowMeans(data.a[ , prps.set$prps.index$set.a, drop = FALSE])
                                ps.b <- rowMeans(data.b[ , prps.set$prps.index$set.b, drop = FALSE])
                                prps <- data.frame(ps.a, ps.b)
                                prps
                            })
                        all.prps.data <- do.call(cbind, all.prps.data)
                        colnames(all.prps.data) <- paste0(
                            pairs.batch[1 , y],
                            '||',
                            pairs.batch[2 , y],
                            '||',
                            rep(1:c(ncol(all.prps.data)/2), each = 2),
                            '||',
                            names(all.possible.batches[x])
                            )
                        all.prps.sample.annot <- lapply(
                            1:length(all.prps.index),
                            function(f){
                                prps.set <- all.prps.index[[f]]$samples.ids
                            })
                        names(all.prps.sample.annot) <- rep(
                            paste0(pairs.batch[1 , y], '_' , pairs.batch[2 , y]),
                            length(all.prps.sample.annot)
                        )
                        list(
                            all.prps.data = all.prps.data,
                            all.prps.sample.annot = all.prps.sample.annot
                        )
                    })
                sub.prps.epxr.data <- lapply(
                    1:length(sub.prps.data),
                    function(t) sub.prps.data[[t]]$all.prps.data
                )
                sub.prps.epxr.data <- do.call(
                    cbind,
                    sub.prps.epxr.data
                )
                sub.prps.info.data <- lapply(
                    1:length(sub.prps.data),
                    function(t) sub.prps.data[[t]]$all.prps.sample.annot
                )
                return(list(
                    sub.prps.epxr.data = sub.prps.epxr.data,
                    sub.prps.info.data = sub.prps.info.data,
                    groups = names(all.possible.batches[x]))
                )
            })

        all.prps.expr.data <- lapply(
            1:length(all.prps.data),
            function(x){
                all.prps.data[[x]]$sub.prps.epxr.data
            })
        all.prps.expr.data <- do.call(cbind, all.prps.expr.data)
        se.obj[[main.uv.variable]] <- initial.variable
        sample.annotation <- as.data.frame(colData(x = se.obj))

        ### plot PRPS map ####
        all.prps.plot <- lapply(
            1:length(all.prps.data),
            function(x){
                prps.index <- all.prps.data[[x]]$sub.prps.info.data
                expr.data <- lapply(
                    1:length(prps.index),
                    function(y){
                        m <- lapply(
                            1:length(prps.index[[y]]),
                            function(z){
                                data.frame(
                                    group1 = sample.annotation[prps.index[[y]][[z]]$set.a , ][[ main.uv.variable]],
                                    group2 =  sample.annotation[prps.index[[y]][[z]]$set.b , ][[ main.uv.variable]],
                                    set.name = rep(
                                        unique(names(prps.index[[y]])),
                                        length(prps.index[[y]][[z]]$set.b)
                                    ),
                                    group.name = rep(
                                        all.prps.data[[x]]$groups,
                                        length(prps.index[[y]][[z]]$set.b)
                                    )
                                )
                            })
                        m <- do.call(rbind, m )
                    })
                expr.data <- do.call(rbind, expr.data)
                expr.data$prps.set <- rep(1:c(nrow(expr.data)/3), each = min.sample.for.ps)
                expr.data
            })
        all.prps.plot <- do.call(rbind, all.prps.plot) %>%
            pivot_longer(-c(set.name, group.name, prps.set), values_to = 'var', names_to = 'rep')
        all.prps.plot$new.g <- paste(
            all.prps.plot$group.name,
            all.prps.plot$set.name,
            all.prps.plot$prps.set,
            sep = '||'
            )
        prps.map <- ggplot(all.prps.plot, aes(x = new.g, y = var)) +
            geom_point() +
            facet_grid(~new.g, scales = 'free', space = 'free') +
            scale_x_discrete(expand = c(0, 0.5)) +
            xlab('') +
            ylab(main.uv.variable) +
            # xlim(c(
            #     min(se.obj[[main.uv.variable]]),
            #     max(se.obj[[main.uv.variable]])
            # )) +
            # geom_hline(yintercept = c(
            #     min(se.obj[[main.uv.variable]]),
            #     max(se.obj[[main.uv.variable]])), color = 'gray70') +
            theme_bw() +
            theme(
                legend.key = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.y = element_text(size = 14),
                axis.text.x = element_text(size = 14, angle = 90, vjust = 1, hjust = 1),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 18),
                strip.text.y = element_text(size = 0),
                strip.text.x = element_text(size = 0)
            )
        if (isTRUE(plot.output)) print(prps.map)

        ### sample annotation ####
        all.prps.sample.annot <- lapply(
            1:length(all.prps.data),
            function(x){
                sample.annot <- all.prps.data[[x]]$sub.prps.info.data
                sample.annot <- sample.annot[[1]]
                sample.annot <- do.call(rbind, sample.annot)
                sample.annot$group <- all.prps.data[[x]]$groups
                sample.annot
            })
        all.prps.sample.annot <- do.call(rbind, all.prps.sample.annot)

    }

    ## other uv variable is FALSE ####
    if (is.null(other.uv.variables)){
        ## check the sample size of each group in the variable ####
        subgroups.size <- findRepeatingPatterns(
            vec = se.obj[[main.uv.variable]],
            n.repeat = max(min.sample.for.ps, min.mnn)
            )
        if (min.batches.to.cover == 'all') {
            if (length(subgroups.size) != length(unique(se.obj[[main.uv.variable]])) ){
                stop(paste0(
                    'Some sub-groups of the variable "',
                    main.uv.variable,
                    '" have less than ',
                    max(min.sample.for.ps, min.mnn),
                    ' samples. Then, MNN cannot be created across all batches.')
                    )
            }
        }
        if (is.numeric(min.batches.to.cover)) {
            if (length(subgroups.size) >= min.batches.to.cover){
                printColoredMessage(
                    message = paste0(
                        '- At least ',
                        min.batches.to.cover,
                        ' sub-groups of the variable ',
                        main.uv.variable,
                        ' have at least ',
                        max(min.sample.for.ps, min.mnn),
                        ' samples.'),
                    color = 'blue',
                    verbose = verbose
                )
            } else {
                stop(paste0(
                    'Some sub-groups of the variable "',
                    main.uv.variable,
                    '" have less than ',
                    max(min.sample.for.ps, min.mnn),
                    ' samples. Then, MNN cannot be created across all batches.')
                )
            }
        }
        # find mutual nearest neighbor ####
        printColoredMessage(
            message = '-- Creating PRPS data across all pairs of the subgroups:',
            color = 'magenta',
            verbose = verbose
        )
        pairs.batch <- combn(x = subgroups.size, m = 2)
        prps.data <- assay(x = se.obj, i = assay.name)
        sample.annotation <- as.data.frame(colData(x = se.obj))
        all.prps.data <- lapply(
            1:ncol(pairs.batch),
            function(y){
                printColoredMessage(
                    message = paste0(
                        '- Creating PRPS data between the "',
                        pairs.batch[1 , y],
                        '" and "' ,
                        pairs.batch[2 , y],
                        '" subgroups."'),
                    color = 'orange',
                    verbose = verbose
                )
                ## sample annotation ####
                sample.annot.a <- sample.annotation[sample.annotation[[main.uv.variable]] == pairs.batch[1 , y] , , drop = FALSE]
                sample.annot.b <- sample.annotation[sample.annotation[[main.uv.variable]] == pairs.batch[2 , y] , , drop = FALSE]

                ## highly variable genes ####
                if (is.null(hvg)){
                    printColoredMessage(
                        message = '* highly variable are not specified, then using all genes.',
                        color = 'blue',
                        verbose = verbose
                    )
                    data.a <- norm.data[ , row.names(sample.annot.a)]
                    data.b <- norm.data[ , row.names(sample.annot.b)]
                }
                if (!is.null(hvg)){
                    printColoredMessage(
                        message = '* using the specified highly variable genes.',
                        color = 'blue',
                        verbose = verbose
                    )
                    data.a <- norm.data[hvg , row.names(sample.annot.a)]
                    data.b <- norm.data[hvg , row.names(sample.annot.b)]
                }

                ## cosine normalization ####
                if (isTRUE(apply.cosine.norm)){
                    printColoredMessage(
                        message = '- Applying cosine normalization on the data:',
                        color = 'blue',
                        verbose = verbose
                    )
                    data.a <- cosineNorm(x = data.a, mode = 'matrix')
                    data.b <- cosineNorm(x = data.b, mode = 'matrix')
                }

                ## finding mnn for data b ####
                printColoredMessage(
                    message = paste0(
                        '* Finding ',
                        min.sample.for.ps,
                        ' nearest neighbours for each sample of the "',
                        pairs.batch[2 , y],
                        '" using the "' ,
                        pairs.batch[1 , y],
                        '" subgroup:'),
                    color = 'blue',
                    verbose = verbose
                )
                knn.data.a.b <- RANN::nn2(
                    data = t(data.a),
                    query = t(data.b),
                    k = min.sample.for.ps
                )
                ## obtaining knn index ####
                printColoredMessage(
                    message = '** obtaining the knn indexs.',
                    color = 'blue',
                    verbose = verbose
                )
                knn.data.a.b.index <- knn.data.a.b$nn.idx
                colnames(knn.data.a.b.index) <- paste0('knn', seq(min.sample.for.ps))
                row.names(knn.data.a.b.index) <- c(1:ncol(data.b))

                ## obtaining the distance ####
                printColoredMessage(
                    message = '** obtaining the distances.',
                    color = 'blue',
                    verbose = verbose
                )
                knn.data.a.b.distance <- as.data.frame(knn.data.a.b$nn.dists)
                colnames(knn.data.a.b.distance) <- paste0('knn', seq(min.sample.for.ps))
                row.names(knn.data.a.b.distance) <- c(1:ncol(data.b))
                knn.data.a.b.distance$aver.dist <- rowMeans(knn.data.a.b.distance)

                ## finding knn for data a ####
                printColoredMessage(
                    message = paste0(
                        '* Finding ',
                        min.sample.for.ps,
                        ' nearest neighbours for each sample of the "',
                        pairs.batch[1 , y],
                        '" using "' ,
                        pairs.batch[2 , y],
                        '" subgroups:'),
                    color = 'blue',
                    verbose = verbose
                )
                knn.data.b.a <- RANN::nn2(
                    data = t(data.b),
                    query = t(data.a),
                    k = min.sample.for.ps
                )
                printColoredMessage(
                    message = '** Obtaining the knn indexs.',
                    color = 'blue',
                    verbose = verbose
                )
                knn.data.b.a.index <- knn.data.b.a$nn.idx
                colnames(knn.data.b.a.index) <- paste0('knn', seq(min.sample.for.ps))
                row.names(knn.data.b.a.index) <- c(1:ncol(data.a))
                ### distance
                printColoredMessage(
                    message = '** Obtaining the distances.',
                    color = 'blue',
                    verbose = verbose
                )
                knn.data.b.a.distance <- as.data.frame(knn.data.b.a$nn.dists)
                colnames(knn.data.b.a.distance) <- paste0('knn', seq(min.sample.for.ps))
                row.names(knn.data.b.a.distance) <- c(1:ncol(data.a))
                knn.data.b.a.distance$aver.dist <- rowMeans(knn.data.b.a.distance)

                ## final mnn ####
                printColoredMessage(
                    message = paste0(
                        '* Finding MNN between the "',
                        pairs.batch[1 , y],
                        '" using "' ,
                        pairs.batch[2 , y],
                        '" subgroups:'),
                    color = 'blue',
                    verbose = verbose
                )
                all.mnn <- BiocNeighbors::findMutualNN(
                    data1 = t(data.a),
                    data2 = t(data.b),
                    k1 = min.mnn,
                    BNPARAM = mnn.nbparam,
                    BPPARAM = mnn.bpparam
                )
                printColoredMessage(
                    message = paste0(
                        '** ',
                        length(all.mnn$first),
                        ' MNN are found.'),
                    color = 'blue',
                    verbose = verbose
                )
                ## create prps index ####
                printColoredMessage(
                    message = paste0('* creating PRPS indexs and score.'),
                    color = 'blue',
                    verbose = verbose
                )
                all.prps.index <- lapply(
                    1:length(all.mnn$first),
                    function(z){
                        sample.to.ave.b <- knn.data.b.a.index[all.mnn$first[z] , ]
                        sample.to.ave.a <- knn.data.a.b.index[all.mnn$second[z] , ]
                        prps.index <- data.frame(
                            set.a = sample.to.ave.a,
                            set.b = sample.to.ave.b
                        )
                        samples.ids <- data.frame(
                            set.a = row.names(sample.annot.a)[sample.to.ave.a],
                            set.b = row.names(sample.annot.b)[sample.to.ave.b]
                        )
                        avre.dist.a <- knn.data.b.a.distance[all.mnn$first[z] , 'aver.dist' ]
                        avre.dist.b <- knn.data.a.b.distance[all.mnn$second[z] , 'aver.dist' ]
                        list(
                            prps.index = prps.index,
                            aver.dist = c(avre.dist.a + avre.dist.b)/2,
                            samples.ids = samples.ids
                        )
                    })
                names(all.prps.index) <- paste0('prps.set', 1:length(all.mnn$first))
                if (isTRUE(filter.prps.sets)){
                    printColoredMessage(
                        message = '* Filtering the number of PRPS sets.',
                        color = 'blue',
                        verbose = verbose
                    )
                    if (length(all.prps.index) >=  max.prps.sets){
                        printColoredMessage(
                            message = '* The number of the PRPS set is larger than the specified "max.prps.sets", filtering PRPS sets.',
                            color = 'blue',
                            verbose = verbose
                        )
                        aver.dists <- sapply(
                            1:length(all.prps.index),
                            function(p) all.prps.index[[p]]$aver.dist
                        )
                        names(aver.dists) <- 1:length(all.prps.index)
                        aver.dists <- aver.dists[order(aver.dists, decreasing = FALSE)]
                        select.prps <- names(aver.dists)[1:max.prps.sets]
                        all.prps.index <- all.prps.index[as.numeric(select.prps)]
                        printColoredMessage(
                            message = paste0(
                                '** ',
                                length(all.prps.index),
                                ' PRPS set are kept.'),
                            color = 'blue',
                            verbose = verbose
                        )
                    }
                }
                ## prps data
                # finding PRPS sets ####
                printColoredMessage(
                    message = '* Creating the PRPS data:',
                    color = 'blue',
                    verbose = verbose
                )
                data.a <- prps.data[ , colnames(data.a)]
                data.b <- prps.data[ , colnames(data.b)]
                if (isTRUE(apply.log)){
                    printColoredMessage(
                        message = '** Applying log transformation on the data before creating PRPS. ',
                        color = 'blue',
                        verbose = verbose
                    )
                    data.a <- log2(data.a + pseudo.count)
                    data.b <- log2(data.b + pseudo.count)
                } else {
                    printColoredMessage(
                        message = '** the data will be used without any transformation. ',
                        color = 'blue',
                        verbose = verbose
                    )
                }
                all.prps.data <- lapply(
                    1:length(all.prps.index),
                    function(a){
                        prps.set <- all.prps.index[[a]]
                        ps.a <- rowMeans(data.a[ , prps.set$prps.index$set.a, drop = FALSE])
                        ps.b <- rowMeans(data.b[ , prps.set$prps.index$set.b, drop = FALSE])
                        prps <- data.frame(ps.a, ps.b)
                        prps
                    })
                all.prps.data <- do.call(cbind, all.prps.data)
                colnames(all.prps.data) <- paste0(
                    pairs.batch[1 , y],
                    '||',
                    pairs.batch[2 , y],
                    '||',
                    rep(1:c(ncol(all.prps.data)/2), each = 2)
                )
                all.prps.sample.annot <- lapply(
                    1:length(all.prps.index),
                    function(f){
                        prps.set <- all.prps.index[[f]]$samples.ids
                    })
                names(all.prps.sample.annot) <- rep(
                    paste0(pairs.batch[1 , y], '_' , pairs.batch[2 , y]),
                    length(all.prps.sample.annot)
                )
                list(
                    all.prps.data = all.prps.data,
                    all.prps.sample.annot = all.prps.sample.annot
                )
            })

        all.prps.expr.data <- lapply(
            1:length(all.prps.data),
            function(x){
                all.prps.data[[x]]$all.prps.data
            })
        all.prps.expr.data <- do.call(cbind, all.prps.expr.data)
        se.obj[[main.uv.variable]] <- initial.variable
        sample.annotation <- as.data.frame(colData(x = se.obj))

        ### plot PRPS map ####
        printColoredMessage(
            message = '- Plotting the PRPS map:',
            color = 'magenta',
            verbose = verbose
        )
        all.prps.sample.annot <- lapply(
            1:length(all.prps.data),
            function(x){
                sample.annot <- all.prps.data[[x]]$all.prps.sample.annot
                sample.annot <- do.call(rbind, sample.annot )
                sample.annot$group <- paste0(
                    unique(names(all.prps.data[[x]]$all.prps.sample.annot)),
                    rep(paste0('||set', 1:min.sample.for.ps), each = min.sample.for.ps)
                    )
                sample.annot <- sample.annot %>%
                    pivot_longer(-c(group), values_to = 'var', names_to = 'rep')
                sample.annot$lib.size <- se.obj[[main.uv.variable]][ sample.annot$var]
                sample.annot
            })
        all.prps.sample.annot <- do.call(rbind, all.prps.sample.annot)
        prps.map <- ggplot(all.prps.sample.annot, aes(x = group, y = lib.size, color = rep)) +
            geom_boxplot() +
            geom_point() +
            facet_grid(~group, scales = 'free', space = 'free') +
            scale_x_discrete(expand = c(0, 0.5)) +
            xlab('') +
            ylab(main.uv.variable) +
            # xlim(c(
            #     min(se.obj[[main.uv.variable]]),
            #     max(se.obj[[main.uv.variable]])
            # )) +
            # geom_hline(yintercept = c(
            #     min(se.obj[[main.uv.variable]]),
            #     max(se.obj[[main.uv.variable]])), color = 'gray70') +
            theme_bw() +
            theme(
                legend.key = element_blank(),
                axis.line = element_line(colour = 'black', linewidth = 1),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.y = element_text(size = 14),
                axis.text.x = element_text(size = 14, angle = 90, vjust = 1, hjust = 1),
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 18),
                strip.text.y = element_text(size = 0),
                strip.text.x = element_text(size = 0)
            )
        if (isTRUE(plot.output)) print(prps.map)

        ## sample annotation ####
        all.prps.sample.annot <- lapply(
            1:length(all.prps.data),
            function(x){
                sample.annot <- all.prps.data[[x]]$all.prps.sample.annot
                sample.annot <- do.call(rbind, sample.annot)
                sample.annot$group <- all.prps.data[[x]]$groups
                sample.annot
            })
        all.prps.sample.annot <- do.call(rbind, all.prps.sample.annot)
    }

    # Save the results ####
    ## select output name ####
    if (is.null(output.name))
        output.name <- paste0(main.uv.variable, '|', 'mnn', '|', assay.name)
    if (is.null(prps.group))
        prps.group <- paste0('prps|mnn|', main.uv.variable)

    printColoredMessage(
        message = '-- Saving the PRPS data',
        color = 'magenta',
        verbose = verbose
        )
    ## save the PRPS data in the SummarizedExperiment object ####
    if (isTRUE(save.se.obj)) {
        printColoredMessage(
            message = 'Save all the PRPS data into the metadata of the SummarizedExperiment object.',
            color = 'blue',
            verbose = verbose
            )
        if (!'PRPS' %in% names(se.obj@metadata)) {
            se.obj@metadata[['PRPS']] <- list()
        }
        if (!'un.supervised' %in% names(se.obj@metadata[['PRPS']])) {
            se.obj@metadata[['PRPS']][['un.supervised']] <- list()
        }
        if (!prps.group %in% names(se.obj@metadata[['PRPS']][['un.supervised']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]] <- list()
        }
        if (!'prps.data' %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']] <- list()
        }
        if (!output.name %in% names(se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']])) {
            se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']][[output.name]] <- list()
        }
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']][[output.name]] <- all.prps.expr.data

        printColoredMessage(message = '------------The createPrPsByMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## output the PRPS data as matrix ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The createPrPsByMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(prps.data = all.prps.data)
    }
}

