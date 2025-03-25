#' Create PRPS sets using k and mutual nearest neighbors in RNA-seq data.

#' @author Ramyar Molania

#' @description
#' This function uses the k and mutual nearest neighbors approaches to create PRPS in the RNA-seq data. This function
#' can be used in situation that the biological variation are entirely unknown. The function applies the 'findKnn' function
#' to find similar samples per batch and then average them to create pseudo-samples. Then, function uses the 'findMnn' to
#' match up pseudo samples across batches to create pseudo-replicates.

#' @param se.obj A summarized experiment object.
#' @param assay.name Symbol. A symbol indicating the name of the assay in the SummarizedExperiment object. This assay will
#' be used to first find k nearest neighbors and them mutual nearest neighbors data. This data must the one that will be
#' used for the RUV-III normalization.
#' @param main.uv.variable Symbol. Indicates the name of a column in the sample annotation of the SummarizedExperiment object.
#' The 'uv.variable' can be either categorical and continuous. If 'uv.variable' is a continuous variable, this will be
#' divided into 'nb.clusters' groups using the 'clustering.method'.
#' @param clustering.method Symbol.A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param filter.prps.sets Logical. If 'TRUE', the number of PRPS sets across each pair of batches will be filtered if they are
#' higher than the 'max.prps.sets' value. The default is 'TRUE'. The high number of PRPS sets will just increase the
#' computational time for the RUV-III normalization.
#' @param max.prps.sets Numeric. A numeric value specifying the maximum number for PRPS set across each pair of batches.
#' The default is set to 10.
#' @param select.extereme.groups TTT
#' @param other.uv.variables Character. A character string or character vector representing the label of unwanted
#' variable(s), such as library size, etc., within the sample annotation (colData) of the SummarizedExperiment object.
#' This can comprise a vector containing either categorical, continuous, or a combination of both variables. These variables
#' will be considered when generating PRPS sets for the "main.uv.variable". This can help avoid potential contamination
#' from the "other.uv.variables" in the PRPS data of the "main.uv.variable". The default is set to 'NULL'.
#' @param other.uv.clustering.method Character. A character string indicating which clustering method should be used to
#' group each continuous unwanted variable, if specified in the "other.uv.variables". The options include 'kmeans', 'cut',
#' and 'quantile'. The default is 'kmeans'. We refer to the createHomogeneousUVGroups() function for more details.
#' @param nb.other.uv.clusters Numeric. A numeric value to specify the number of clusters/groups for each continuous unwanted
#' variable. The default is set to 3.
#' @param min.sample.for.ps Nummeric.
#' @param min.batch.to.cover Numeric. The maximum number of nearest neighbors to compute. The default is set 3.
#' @param data.input Character. A character string that indicates which data should be used as input for finding
#' the k nearest neighbors. Options include: 'expr' and 'pcs'. If 'pcs' is selected, the first 'nb.pcs' of PCs of
#' the data will be used as input. If 'expr' is selected, the expression data will be used as input. The default is set
#' to 'expr'.
#' @param nb.pcs Numeric. A numeric value that indicates the number of PCs to be used as data input for finding the k
#' nearest neighbors. The 'nb.pcs' must be set when "data.input = pcs". The default is set to 2.
#' @param center Logical. Indicates whether to scale the data or not. If center is TRUE, then centering is done by
#' subtracting the column means of the assay from their corresponding columns. The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD. If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is set to 'FALSE'.
#' @param svd.bsparam Character. A BiocParallelParam object specifying how parallelization should be performed. The default
#' is set to bsparam(). We refer to the 'runSVD' function from the BiocSingular R package for further details.
#' @param min.knn description
#' @param min.mnn description
#' @param hvg Vector. A vector of the names of the highly variable genes. These genes will be used to find the anchors
#' samples across the batches. The default is NULL.
#' @param normalization Symbol. Indicates which normalization methods should be applied before finding the knn. The default
#' is 'cpm'. If is set to NULL, no normalization will be applied.
#' @param apply.cosine.norm TTT
#' @param regress.out.variables Symbols. Indicates the columns names that contain biological variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by unwanted variation variable. The default is NULL, indicates the regression will not be applied.
#' @param check.prps.connectedness Logical. Indicates whether to assess the connectedness between the PRPS sets or not.
#' The default is set to TRUE. See the details for more information.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is TRUE.
#' Please, note, any RNA-seq data (assays) must be in log scale before computing RLE.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay(s) before applying
#' log transformation to avoid -Inf for measurements that are equal to 0. The default is 1.
#' @param mnn.bpparam Symbol. A BiocParallelParam object specifying how parallelization should be performed. The default
#' is SerialParam(). We refer to the 'findMutualNN' function from the BiocNeighbors R package for more details.
#' @param mnn.nbparam Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' . The default is KmknnParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment object
#' or to output the result as list. By default it is set to TRUE.
#' @param plot.output TTTT
#' @param output.name Symbol. A symbol specifying the name of output file. If is 'NULL', the function will select a name
#' based on "paste0(uv.variable, '|', 'anchor', '|', assay.name))".
#' @param prps.group Symbol. A symbol specifying the name of the output file. If is 'NULL', the function will select a name
#' based on "paste0('prps_mnn_', uv.variable)".
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @return The SummarizedExperiment object that contain all the PPRS data, knn, mnn and plot results in the metadata, or
#' a list of the results.

#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom BiocNeighbors findMutualNN KmknnParam
#' @importFrom SummarizedExperiment assay colData
#' @importFrom BiocParallel SerialParam
#' @importFrom stats dist
#' @importFrom RANN nn2
#' @export

createPrPsByKnnMnn <- function(
        se.obj,
        assay.name,
        main.uv.variable,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        select.extreme.groups = TRUE,
        other.uv.variables = NULL,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        min.sample.for.ps = 3,
        min.batch.to.cover = 'all',
        check.prps.connectedness = TRUE,
        data.input = 'expr',
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        nb.knn = 2,
        nb.mnn = 1,
        hvg = NULL,
        normalization = 'CPM',
        apply.cosine.norm = FALSE,
        regress.out.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        assess.se.obj = TRUE,
        plot.output = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        output.name = NULL,
        prps.group = NULL,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The createPrPsByKnnMnn function starts:',
                        color = 'white',
                        verbose = verbose)

    # Assessing the SummarizedExperiment object ####
    if (isTRUE(assess.se.obj)) {
        se.obj <- checkSeObj(
            se.obj = se.obj,
            assay.names = assay.name,
            variables = c(main.uv.variable, other.uv.variables),
            remove.na = remove.na,
            verbose = verbose
        )
    }
    # Assessing and grouping the main unwanted variable ####
    printColoredMessage(
        message = '- Assessing and grouping the unwanted variable:',
        color = 'magenta',
        verbose = verbose
        )
    if (isTRUE(select.extreme.groups)){
        se.obj.initial <- se.obj
    }

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
        if (isTRUE(select.extreme.groups)){
            printColoredMessage(
                message = paste0(
                    '- Selecting the two subgroups of the ',
                    main.uv.variable,
                    ' variable with highest and lowest values.'),
                color = 'blue',
                verbose = verbose
            )
            max.group <- se.obj[[main.uv.variable]][initial.variable == max(initial.variable)]
            min.group <- se.obj[[main.uv.variable]][initial.variable == min(initial.variable)]
            selected.samples <- se.obj[[main.uv.variable]] %in% c(max.group, min.group)
            se.obj <- se.obj[ , selected.samples]
            se.obj[[main.uv.variable]] <- droplevels(se.obj[[main.uv.variable]])
            initial.variable <- initial.variable[selected.samples]
        }
    }
    if (!is.numeric(initial.variable)){
        length.variable <- length(unique(initial.variable))
        if (length.variable == 1){
            stop('To create PRPS, the "main.uv.variable" must have at least two groups/levels.')
        } else if (length.variable > 1){
            printColoredMessage(
                message = paste0(
                    '- The "',
                    uv.variable,
                    '" is a categorical variable with ',
                    length(unique(se.obj[[uv.variable]])),
                    ' levels.'),
                color = 'blue',
                verbose = verbose
            )
            se.obj[[uv.variable]] <- factor(x = se.obj[[uv.variable]])
        }
    }

    # Creating PRPS data with KNN and MNN ####
    if (!is.null(other.uv.variables)){
        ## considering other unwanted variables ####
        printColoredMessage(
            message = '- Creating PRPS data by considering other specified variables',
            color = 'magenta',
            verbose = verbose
            )
        ## grouping the other unwanted variables ####
        printColoredMessage(
            message = '-- Assessing and grouping the other specified unwanted variable(s):',
            color = 'blue',
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
        min.sample.size <- max(
            min.sample.for.ps,
            nb.mnn,
            nb.knn
            )
        ## finding and plotting covered batches ####
        covered.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                subgroups.size <- findRepeatingPatterns(
                    vec = all.uv.groups$main.uv[all.uv.groups$other.uv == x],
                    n.repeat = min.sample.size
                )
            })
        names(covered.batches) <- unique(all.uv.groups$other.uv)
        covered.batches.plot <- table(all.uv.groups$main.uv, all.uv.groups$other.uv) %>%
            data.frame(.) %>%
            dplyr::mutate(selected = Freq >= min.sample.size) %>%
            ggplot(data = ., aes(x = Var2, y = Var1, color = selected)) +
            geom_point(size = 4) +
            geom_text(aes(label = Freq , hjust = 0.5, vjust = 0.5), color = 'black', size = 5) +
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
        if (isTRUE(plot.output)) print(covered.batches.plot)

        ## checking covered batches ####
        printColoredMessage(
            message = '-- Checking the distribution of the "main.uv.variable" across the "other.uv.variables":',
            color = 'blue',
            verbose = verbose
            )
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
                    min.sample.size,
                    ' samples across all sub-groups of the ',
                    main.uv.variable,
                    ' variable.'),
                color = 'blue',
                verbose = verbose
            )
        }
        ## finding K nearest neighbors ####
        printColoredMessage(
            message = '-- Finding k nearest neighbor by applying the findKnn function:',
            color = 'magenta',
            verbose = verbose
            )
        all.possible.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                possible.batch <- findRepeatingPatterns(
                    vec = all.uv.groups[all.uv.groups$other.uv == x, ]$main.uv,
                    n.repeat = max(min.sample.for.ps, nb.mnn)
                )
                if (length(possible.batch) > 1){
                    combn(x = possible.batch , m = 2)
                } else NA

            })
        names(all.possible.batches) <- unique(all.uv.groups$other.uv)
        all.possible.batches <- all.possible.batches[!is.na(all.possible.batches)]
        all.knn <- lapply(
            1:length(all.possible.batches),
            function(x){
                all.knn.samples <- findKnn(
                    se.obj = se.obj[ , homo.uv.groups == names(all.possible.batches)[x] ],
                    assay.name = assay.name,
                    uv.variable = main.uv.variable,
                    data.input = data.input,
                    nb.pcs = nb.pcs,
                    center = center,
                    scale = scale,
                    svd.bsparam = svd.bsparam,
                    clustering.method = clustering.method,
                    nb.clusters = nb.clusters,
                    nb.knn = nb.knn,
                    hvg = hvg,
                    normalization = normalization,
                    regress.out.variables = regress.out.variables,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    assess.se.obj = FALSE,
                    remove.na = remove.na,
                    output.name = output.name,
                    prps.group = prps.group,
                    save.se.obj = FALSE,
                    verbose = verbose
                )
                all.knn.samples$other.group <- names(all.possible.batches[x])
                all.knn.samples
            })
        names(all.knn) <- names(all.possible.batches)

        ## finding mutual nearest neighbors ####
        printColoredMessage(
            message = '-- Finding mutual nearest neighbors by applying the findMnn function:',
            color = 'magenta',
            verbose = verbose
            )
        all.mnn <- lapply(
            1:length(all.possible.batches),
            function(x){
                all.mnn.samples <- findMnn(
                    se.obj = se.obj[ , homo.uv.groups == names(all.possible.batches)[x] ],
                    assay.name = assay.name,
                    uv.variable = main.uv.variable,
                    nb.mnn = nb.mnn,
                    clustering.method = clustering.method,
                    nb.clusters = nb.clusters,
                    data.input = data.input,
                    nb.pcs = nb.pcs,
                    center = center,
                    scale = scale,
                    svd.bsparam = svd.bsparam,
                    normalization = normalization,
                    apply.cosine.norm = apply.cosine.norm,
                    regress.out.variables = regress.out.variables,
                    hvg = hvg,
                    apply.log = apply.log,
                    pseudo.count = pseudo.count,
                    mnn.bpparam = mnn.bpparam,
                    mnn.nbparam = mnn.nbparam,
                    assess.se.obj = assess.se.obj,
                    remove.na = remove.na,
                    plot.output = FALSE,
                    output.name = output.name,
                    prps.group = prps.group,
                    save.se.obj = FALSE,
                    verbose = verbose
                    )
                all.mnn.samples$mnn$other.group <- names(all.possible.batches[x])
                all.mnn.samples
            })
        ## plotting all mnn ####
        all.mnn.plots <- lapply(
            1:length(all.possible.batches),
            function(x){
                all.mnn[[x]]$mnn.plot + ggtitle(names(all.possible.batches[x]))
            })
        all.mnn.plots <- ggarrange(plotlist = all.mnn.plots)
        if (isTRUE(plot.output)) print(all.mnn.plots)
        ## obtaning mnn data ####
        all.mnn <- lapply(
            1:length(all.possible.batches),
            function(x){
                all.mnn[[x]]$mnn
            })
        names(all.mnn) <- names(all.possible.batches)

        ## matching mnn and knn data ####
        printColoredMessage(
            message = '- Matching the MNN sets with the corresponding KNN sets:',
            color = 'magenta',
            verbose = verbose
            )
        mnn.sets <- NULL
        all.prps.sets <- lapply(
            names(all.possible.batches),
            function(i){
                sub.all.knn <- all.knn[[i]]
                sub.all.mnn <- all.mnn[[i]]
                # sanity check ####
                if (length(unique(sub.all.knn$other.group)) > 1){
                    stop('There something wrong with knn and mnn.')
                }
                if (length(unique(sub.all.mnn$other.group)) > 1){
                    stop('There something wrong with knn and mnn.')
                }
                if (!unique(sub.all.knn$other.group) == unique(sub.all.mnn$other.group)){
                    stop('There something wrong with knn and mnn.')
                }
                sub.prps.sets <- lapply(
                    1:nrow(sub.all.mnn),
                    function(x) {
                        # ps set 1
                        ps.set.1 <- sub.all.knn[ ,   grep('sample.ids', colnames(sub.all.knn))] == sub.all.mnn$sample.ids.1[x]
                        ps.set.1 <- sub.all.knn[rowSums(ps.set.1) > 0 , ]
                        ps.set.1$mnn.sets <- paste0(paste0(
                            sort(c(sub.all.mnn[x , 3], sub.all.mnn[x , 4])),
                            collapse = '_'),
                            '_',
                            x,
                            '_',
                            i
                            )
                        ps.set.1$mnn.sets.data <- paste0(sort(c(sub.all.mnn[x , 1], sub.all.mnn[x , 2])), collapse = '_')
                        if (nrow(ps.set.1) > 1) {
                            ps.set.1 <- ps.set.1[ps.set.1$rank.aver.dist == min(ps.set.1$rank.aver.dist) , ]
                        }
                        # ps set 2
                        ps.set.2 <- sub.all.knn[ ,   grep('sample.ids', colnames(sub.all.knn))] == sub.all.mnn$sample.ids.2[x]
                        ps.set.2 <- sub.all.knn[rowSums(ps.set.2) > 0 , ]
                        ps.set.2$mnn.sets <- paste0(paste0(
                            sort(c(sub.all.mnn[x , 3], sub.all.mnn[x , 4])),
                            collapse = '_'),
                            '_',
                            x,
                            '_',
                            i
                            )
                        ps.set.2$mnn.sets.data <- paste0(sort(c(sub.all.mnn[x , 1], sub.all.mnn[x , 2])), collapse = '_')
                        if (nrow(ps.set.2) > 1) {
                            ps.set.2 <- ps.set.2[ps.set.2$rank.aver.dist == min(ps.set.2$rank.aver.dist) , ]
                        }
                        prps.set <- rbind(ps.set.1, ps.set.2)
                        # sanity check ####
                        if( sum(prps.set[1 , grep('sample.ids', colnames(sub.all.knn))] %in% prps.set[2 , grep('sample.ids', colnames(sub.all.knn))]) > 1 ){
                            stop('There something wrong with knn and mnn.')
                        }
                        prps.set
                    })
                sub.prps.sets <- do.call(rbind, sub.prps.sets)
            })
        names(all.prps.sets) <- names(all.possible.batches)

        if (is.null(all.prps.sets)) {
            stop('PRPS cannot be created. You may want to increase the value of the mnn.')
        }
        ## applying a sanity check ####
        printColoredMessage(
            message =  '- Applying a sanity check on the mnn and knn.',
            color = 'blue',
            verbose = verbose
        )
        sanity.check <- unlist(lapply(
            names(all.prps.sets),
            function(x){
                nrow(all.prps.sets[[x]]) != 2*nrow(all.mnn[[x]])
            }))
        if (sum(sanity.check) == 0){
            printColoredMessage(
                message =  '- The rows of the of the matched MNN and KNN is correct.',
                color = 'blue',
                verbose = verbose
            )
        } else {
            stop('For individual MNN set, the corresponding KNN sets cannot be found. Check the the input.')
        }

        ## adding the average distances of each knn sets for each PRPS set and then rank them ####
        printColoredMessage(
            message = '- Averaging the distances of each knn sets for each MNN set and then rank them:',
            color = 'blue',
            verbose = verbose
            )
        aver.mnn.sets <- NULL
        all.prps.sets <- do.call(rbind, all.prps.sets)

        all.prps.sets$aver.mnn.sets <- unlist(lapply(
            seq(1, nrow(all.prps.sets), 2),
            function(x)
                rep(mean(all.prps.sets$aver.dist[x:(x + 1)]), 2))
            )
        set.seed(2233)
        all.prps.sets$rank.aver.mnn.sets <- rank(
            x = all.prps.sets$aver.mnn.sets,
            ties.method = 'random'
            )
        ## filtering PRPS sets ####
        if (isTRUE(filter.prps.sets)) {
            printColoredMessage(
                message = '- Filtering the PRPS sets across each pair of batches:',
                color = 'orange',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    '- The maximum number of PRPS sets across each pairs of batches is ',
                    max.prps.sets,
                    '.'),
                color = 'blue',
                verbose = verbose
                )
            printColoredMessage(
                message = paste0(
                    '- The PRPS sets will be filtered based on the distances between each knn sets.'),
                color = 'blue',
                verbose = verbose
                )
            all.prps.sets <- lapply(
                unique(all.prps.sets$mnn.sets.data),
                function(x) {
                    temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                    if (length(unique(temp.prps.set$mnn.sets)) >= max.prps.sets) {
                        printColoredMessage(
                            message = paste0(
                                '* The number of PRPS sets across the batches "',
                                x,
                                '" is ',
                                length(unique(temp.prps.set$mnn.sets)), '.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- arrange(temp.prps.set, aver.mnn.sets, mnn.sets)
                        printColoredMessage(
                            message = paste0(
                                '* ',
                                length(unique(temp.prps.set$mnn.sets)) - max.prps.sets,
                                ' PRPS sets are removed.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- temp.prps.set[1:c(2 * max.prps.sets) , ]
                    } else if (length(unique(temp.prps.set$mnn.sets))  < max.prps.sets) {
                        printColoredMessage(
                            message = paste0(
                                '* The number of PRPS sets across the batches "',
                                x,
                                '" is ',
                                length(unique(temp.prps.set$mnn.sets)) , '.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                    }
                    return(temp.prps.set)
                })
            all.prps.sets <- do.call(rbind, all.prps.sets)
        }
        printColoredMessage(
            message = paste0(
                '- ',
                length(unique(all.prps.sets$mnn.sets)),
                ' PRPS stes are found in total.'),
            color = 'blue',
            verbose = verbose
            )
        ## creating PRPS data ####
        printColoredMessage(
            message = '-- Creating PRPS data matrix:',
            color = 'magenta',
            verbose = verbose
            )
        ## applying log ####
        printColoredMessage(
            message = '- Applying data log transformation before creating the PRPS expression data:',
            color = 'blue',
            verbose = verbose
            )
        if (isTRUE(apply.log) & !is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0(
                    '- Applying log2 on the "',
                    assay.name,
                    '" + ',
                    pseudo.count,
                    ' (pseudo.count)  data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        } else if (isTRUE(apply.log) & is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0(
                    'Applying log2 on the "',
                    assay.name,
                    '" data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log)) {
            printColoredMessage(
                message = paste0(
                    'The "',
                    assay.name,
                    '" data will be used without any log transformation.' ),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- assay(x = se.obj, i = assay.name)
        }
        printColoredMessage(
            message = '- Aeveraging samples to create pseudo samples:',
            color = 'blue',
            verbose = verbose
        )
        prps.data <- lapply(
            unique(all.prps.sets$mnn.sets),
            function(i){
                tep.mnn.sets <- all.prps.sets[all.prps.sets$mnn.sets == i , ]
                set.a <- tep.mnn.sets[1 , grep('sample.ids', colnames(all.prps.sets))]
                set.a <- rowMeans(expr.data[ , unlist(as.vector(set.a[1 , ])) ])
                set.b <- tep.mnn.sets[2 , grep('sample.ids', colnames(all.prps.sets))]
                set.b <- rowMeans(expr.data[ , unlist(as.vector(set.b[1 , ])) ])
                prps <- cbind(set.a, set.b)
                colnames(prps) <- rep(paste(
                    main.uv.variable,
                    i,
                    sep = '_'),  2)
                prps
            })
        prps.data <- do.call(cbind, prps.data)

        ## sanity check ####
        if (!sum(table(colnames(prps.data)) == 2) == ncol(prps.data) / 2) {
            stop('There someting wrong with PRPS sets.')
        }
        se.obj[[uv.variable]] <- initial.variable

        ## plotting the PRPS map ####
        prps.map <- lapply(
            unique(all.prps.sets$mnn.sets),
            function(i){
                tep.mnn.sets <- all.prps.sets[all.prps.sets$mnn.sets == i , ]
                set.a <- tep.mnn.sets[1 , grep('sample.ids', colnames(all.prps.sets))]
                set.a <- unlist(as.vector(set.a[1 , ]))
                set.b <- tep.mnn.sets[2 , grep('sample.ids', colnames(all.prps.sets))]
                set.b <- unlist(as.vector(set.b[1 , ]))
                initial.variable.set.a <- initial.variable[colnames(se.obj) %in% set.a]
                initial.variable.set.b <- initial.variable[colnames(se.obj) %in% set.b]
                data.frame(
                    group1 = initial.variable.set.a,
                    group2 = initial.variable.set.b,
                    set = rep(i, min.sample.for.ps)
                    )
            })
        prps.map <- do.call(rbind, prps.map)
        prps.map <- pivot_longer(prps.map, -set, names_to = 'group', values_to = 'var')
        prps.map$group2 <- 'PRPS sets'

        all.uv.variable <- data.frame(set = main.uv.variable, group = 'UV', var = se.obj.initial[[main.uv.variable]], group2 = 'UV')
        prps.map <- rbind(prps.map, all.uv.variable)
        prps.map.plot <- ggplot(prps.map, aes(x = set, y = var, color = group)) +
            geom_boxplot() +
            geom_point(size = 2) +
            xlab('Homogeneous groups') +
            ylab(main.uv.variable) +
            scale_color_manual(values = c('darkgreen', 'tomato', 'navy')) +
            facet_grid(.~group2, scales = 'free', space = 'free') +
            scale_x_discrete(expand = c(0, 0.5)) +
            theme_bw() +
            theme(
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 18),
                axis.line = element_line(colour = 'black', linewidth = .85),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5 , vjust = 0.5),
                axis.text.y = element_text(size = 12),
                legend.position = 'right') +
            guides(color = guide_legend(title = "PS"))
        if(isTRUE(verbose)) print(prps.map.plot)
    }

    ## considering only main unwanted variables ####
    if (is.null(other.uv.variables)){
        ### checking sample sizes of each sub group ####
        ### KNN ####
        sub.group.sample.size.knn <- findRepeatingPatterns(
            vec = se.obj[[main.uv.variable]],
            n.repeat = nb.knn + 1
            )
        if (length(sub.group.sample.size.knn) == 0){
            stop(paste0(
                'All subgroups of the unwanted variable have less than ',
                nb.knn + 1,
                ' (nb.knn + 1) samples. KNN cannot be found.')
                )
        } else if (length(sub.group.sample.size.knn) != length(unique(se.obj[[uv.variable]])) ){
            printColoredMessage(
                message = paste0(
                    'All or some subgroups of the unwanted variable have less than ',
                    nb.knn + 1,
                    ' (nb.knn + 1) samples. Then KNN for those sub-groups cannot be created.'),
                color = 'red',
                verbose = verbose
                )
        } else {
            printColoredMessage(
                message = paste0(
                    '- All the sub-groups of the unwanted variable have at least ',
                    nb.knn + 1,
                    ' nb.knn + 1 samples.'),
                color = 'blue',
                verbose = verbose
                )
        }
        ### MNN ####
        sub.group.sample.size.mnn <- findRepeatingPatterns(
            vec = se.obj[[uv.variable]],
            n.repeat = nb.mnn + 1
            )
        if (length(sub.group.sample.size.mnn) == 0){
            stop(paste0(
                'All subgroups of the unwanted variable have less than ',
                nb.mnn + 1,
                ' (nb.mnn + 1) samples. MNN cannot be found.')
                )
        } else if (length(sub.group.sample.size.mnn) != length(unique(se.obj[[uv.variable]])) ){
            printColoredMessage(
                message = paste0(
                    'All or some subgroups of the unwanted variable have less than ',
                    nb.mnn + 1,
                    ' (nb.mnn + 1) samples.'),
                color = 'red',
                verbose = verbose
            )
        } else {
            printColoredMessage(
                message = paste0(
                    '- All the subgroups of the unwanted variable have at least, ',
                    nb.mnn + 1,
                    ' (nb.mnn + 1) samples.'),
                color = 'blue',
                verbose = verbose
            )
        }
        ### finding k nearest neighbor ####
        printColoredMessage(
            message = '-- Finding k nearest neighbor by applying the findKnn function:',
            color = 'magenta',
            verbose = verbose
        )
        se.obj <- findKnn(
            se.obj = se.obj,
            assay.name = assay.name,
            uv.variable = main.uv.variable,
            data.input = data.input,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            svd.bsparam = svd.bsparam,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            nb.knn = nb.knn,
            hvg = hvg,
            normalization = normalization,
            regress.out.variables = regress.out.variables,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            output.name = output.name,
            prps.group = prps.group,
            save.se.obj = save.se.obj,
            verbose = verbose
            )
        if (is.null(output.name)) {
            output.name.knn <- paste0(uv.variable, '|' , assay.name)
        } else output.name.knn <- output.name
        if (is.null(prps.group)){
            prps.group.mnn <- paste0('prps|knnMnn|', uv.variable)
        } else prps.group.mnn <- prps.group
        all.knn <- se.obj@metadata$PRPS$un.supervised[[prps.group.mnn]]$KnnMnn$knn[[output.name.knn]]
        if (isFALSE(save.se.obj)) {
            se.obj@metadata$PRPS$un.supervised[[prps.group.mnn]]$KnnMnn$knn[[output.name.knn]] <- NULL
        }

        ### finding mutual nearest neighbor ####
        printColoredMessage(
            message = '-- Finding mutual nearest neighbors by applying the findMnn function:',
            color = 'magenta',
            verbose = verbose
            )
        se.obj <- findMnn(
            se.obj = se.obj,
            assay.name = assay.name,
            uv.variable = uv.variable,
            nb.mnn = nb.mnn,
            clustering.method = clustering.method,
            nb.clusters = nb.clusters,
            data.input = data.input,
            nb.pcs = nb.pcs,
            center = center,
            scale = scale,
            svd.bsparam = svd.bsparam,
            normalization = normalization,
            apply.cosine.norm = apply.cosine.norm,
            regress.out.variables = regress.out.variables,
            hvg = hvg,
            apply.log = apply.log,
            pseudo.count = pseudo.count,
            mnn.bpparam = mnn.bpparam,
            mnn.nbparam = mnn.nbparam,
            assess.se.obj = assess.se.obj,
            remove.na = remove.na,
            plot.output = plot.output,
            output.name = output.name,
            prps.group = prps.group,
            save.se.obj = save.se.obj,
            verbose = verbose
            )
        if (is.null(output.name)) {
            output.name.mnn <- paste0(uv.variable, '|' , assay.name)
            } else output.name.mnn <- output.name
        if (is.null(prps.group)){
            prps.group.mnn <- paste0('prps|knnMnn|', uv.variable)
            } else prps.group.mnn <- prps.group
        all.mnn <- se.obj@metadata$PRPS$un.supervised[[prps.group.mnn]]$KnnMnn$mnn[[output.name.mnn]]
        if (isFALSE(save.se.obj)) {
            se.obj@metadata$PRPS$un.supervised[[prps.group.mnn]]$KnnMnn$mnn[[output.name.mnn]] <- NULL
            }

        ### matching KNN and MNN  ####
        printColoredMessage(
            message = '-- Finding all possible similar samples across batches:',
            color = 'magenta',
            verbose = verbose
            )
        ### finding the knn for each mnn set ####
        printColoredMessage(
            message = '- Matching the results of the KNN and MNN data:',
            color = 'orange',
            verbose = verbose
            )
        printColoredMessage(
            message = '* Matching the MNN sets with the corresponding KNN sets:',
            color = 'blue',
            verbose = verbose
            )
        mnn.sets <- NULL
        ### finding all the prps sets ####
        all.prps.sets <- lapply(
            1:nrow(all.mnn),
            function(x) {
                # ps set 1
                ps.set.1 <- all.knn[ ,   grep('sample.ids', colnames(all.knn))] == all.mnn$sample.ids.1[x]
                ps.set.1 <- all.knn[rowSums(ps.set.1) > 0 , ]
                ps.set.1$mnn.sets <- paste0(paste0(
                    sort(c(all.mnn[x , 3], all.mnn[x , 4])),
                    collapse = '_'),
                    '_',
                    x
                    )
                ps.set.1$mnn.sets.data <- paste0(sort(c(all.mnn[x , 1], all.mnn[x , 2])), collapse = '_')
                if (nrow(ps.set.1) > 1) {
                    ps.set.1 <- ps.set.1[ps.set.1$rank.aver.dist == min(ps.set.1$rank.aver.dist) , ]
                }
                # ps set 2
                ps.set.2 <- all.knn[ ,   grep('sample.ids', colnames(all.knn))] == all.mnn$sample.ids.2[x]
                ps.set.2 <- all.knn[rowSums(ps.set.2) > 0 , ]
                ps.set.2$mnn.sets <- paste0(paste0(
                    sort(c(all.mnn[x , 3], all.mnn[x , 4])),
                    collapse = '_'),
                    '_',
                    x
                )
                ps.set.2$mnn.sets.data <- paste0(sort(c(all.mnn[x , 1], all.mnn[x , 2])), collapse = '_')
                if (nrow(ps.set.2) > 1) {
                    ps.set.2 <- ps.set.2[ps.set.2$rank.aver.dist == min(ps.set.2$rank.aver.dist) , ]
                }
                prps.set <- rbind(ps.set.1, ps.set.2)
                # sanity check ####
                if( sum(prps.set[1 , grep('sample.ids', colnames(all.knn))] %in% prps.set[2 , grep('sample.ids', colnames(all.knn))]) > 1 ){
                    stop('There something wrong with knn and mnn.')
                }
                prps.set
            })
        all.prps.sets <- do.call(rbind, all.prps.sets)

        if (is.null(all.prps.sets)) {
            stop('PRPS cannot be created. You may want to increase the value of the mnn.')
        }
        ### sanity check ####
        if(nrow(all.prps.sets) == 2*nrow(all.mnn)){
            printColoredMessage(
                message = paste0(
                    '* The nrow of the matched MNN and KNN is ',
                    nrow(all.prps.sets),
                    '.'),
                color = 'blue',
                verbose = verbose
            )
        } else {
            stop('For individual MNN set, the corresponding KNN sets cannot be found. Check the the input.')
        }

        ### adding the average of the knn sets for each PRPS set and then rank them ####
        printColoredMessage(
            message = '* Average the knn sets for each MNN set and then rank them:',
            color = 'blue',
            verbose = verbose
        )
        aver.mnn.sets <- NULL
        all.prps.sets$aver.mnn.sets <- unlist(lapply(
            seq(1, nrow(all.prps.sets), 2),
            function(x)
                rep(mean(all.prps.sets$aver.dist[x:(x + 1)]), 2))
        )
        set.seed(2233)
        all.prps.sets$rank.aver.mnn.sets <- rank(
            x = all.prps.sets$aver.mnn.sets,
            ties.method = 'random'
        )

        ### filtering PRPS sets ####
        if (isTRUE(filter.prps.sets)) {
            printColoredMessage(
                message = '- Filtering the PRPS sets across each pair of batches:',
                color = 'orange',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    '- The maximum number of PRPS sets across each pairs of batches is ',
                    max.prps.sets,
                    '.'),
                color = 'blue',
                verbose = verbose
            )
            printColoredMessage(
                message = paste0(
                    '- The PRPS sets will be filtered based on the distances between each knn sets.'),
                color = 'blue',
                verbose = verbose
            )
            all.prps.sets <- lapply(
                unique(all.prps.sets$mnn.sets.data),
                function(x) {
                    temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                    if (length(unique(temp.prps.set$mnn.sets)) >= max.prps.sets) {
                        printColoredMessage(
                            message = paste0(
                                '* The number of PRPS sets across the batches "',
                                x,
                                '" is ',
                                length(unique(temp.prps.set$mnn.sets)), '.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- arrange(temp.prps.set, aver.mnn.sets, mnn.sets)
                        printColoredMessage(
                            message = paste0(
                                '* ',
                                length(unique(temp.prps.set$mnn.sets)) - max.prps.sets,
                                ' PRPS sets are removed.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- temp.prps.set[1:c(2 * max.prps.sets) , ]
                    } else if (length(unique(temp.prps.set$mnn.sets))  < max.prps.sets) {
                        printColoredMessage(
                            message = paste0(
                                '* The number of PRPS sets across the batches "',
                                x,
                                '" is ',
                                length(unique(temp.prps.set$mnn.sets)) , '.' ),
                            color = 'blue',
                            verbose = verbose
                        )
                        temp.prps.set <- all.prps.sets[all.prps.sets$mnn.sets.data == x ,]
                    }
                    return(temp.prps.set)
                })
            all.prps.sets <- do.call(rbind, all.prps.sets)
        }
        printColoredMessage(
            message = paste0(
                '- ',
                length(unique(all.prps.sets$mnn.sets)),
                ' PRPS stes are found in total.'),
            color = 'blue',
            verbose = verbose
        )

        ### Ccreating PRPS data ####
        printColoredMessage(
            message = '-- Creating PRPS data:',
            color = 'magenta',
            verbose = verbose
        )
        ### applying log ####
        printColoredMessage(
            message = '- Applying data log transformation before creating the PRPS expression data:',
            color = 'blue',
            verbose = verbose
        )
        if (isTRUE(apply.log) & !is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0(
                    '- Applying log2 on the "',
                    assay.name,
                    '" + ',
                    pseudo.count,
                    ' (pseudo.count)  data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        } else if (isTRUE(apply.log) & is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0(
                    'Applying log2 on the "',
                    assay.name,
                    '" data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log)) {
            printColoredMessage(
                message = paste0(
                    'The "',
                    assay.name,
                    '" data will be used without any log transformation.' ),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- assay(x = se.obj, i = assay.name)
        }
        printColoredMessage(
            message = '- Aeveraging samples to create pseudo samples:',
            color = 'blue',
            verbose = verbose
        )
        prps.data <- lapply(
            unique(all.prps.sets$mnn.sets),
            function(x) {
                temp.prps <- all.prps.sets[all.prps.sets$mnn.sets == x, ]
                index.a <- unlist(unname(temp.prps[1, grep('overal', colnames(temp.prps))]))
                index.b <- unlist(unname(temp.prps[2, grep('overal', colnames(temp.prps))]))
                prps.a <- rowMeans(expr.data[, index.a])
                prps.b <- rowMeans(expr.data[, index.b])
                prps <- cbind(prps.a, prps.b)
                colnames(prps) <- paste(uv.variable, temp.prps$mnn.sets, sep = '_')
                return(prps)
            })
        prps.data <- do.call(cbind, prps.data)
        ### sanity check ####
        if (!sum(table(colnames(prps.data)) == 2) == ncol(prps.data) / 2) {
            stop('There someting wrong with PRPS sets.')
        }
        ## plotting the PRPS map ####
        prps.map <- lapply(
            unique(all.prps.sets$mnn.sets),
            function(i){
                tep.mnn.sets <- all.prps.sets[all.prps.sets$mnn.sets == i , ]
                set.a <- tep.mnn.sets[1 , grep('sample.ids', colnames(all.prps.sets))]
                set.a <- unlist(as.vector(set.a[1 , ]))
                set.b <- tep.mnn.sets[2 , grep('sample.ids', colnames(all.prps.sets))]
                set.b <- unlist(as.vector(set.b[1 , ]))
                initial.variable.set.a <- initial.variable[colnames(se.obj) %in% set.a]
                initial.variable.set.b <- initial.variable[colnames(se.obj) %in% set.b]
                data.frame(
                    group1 = initial.variable.set.a,
                    group2 = initial.variable.set.b,
                    set = rep(i, min.sample.for.ps)
                )
            })
        prps.map <- do.call(rbind, prps.map)
        prps.map <- pivot_longer(prps.map, -set, names_to = 'group', values_to = 'var')
        prps.map$group2 <- 'PRPS sets'

        all.uv.variable <- data.frame(set = main.uv.variable, group = 'UV', var = se.obj.initial[[main.uv.variable]], group2 = 'UV')
        prps.map <- rbind(prps.map, all.uv.variable)
        prps.map.plot <- ggplot(prps.map, aes(x = set, y = var, color = group)) +
            geom_boxplot() +
            geom_point(size = 2) +
            xlab('Homogeneous groups') +
            ylab(main.uv.variable) +
            scale_color_manual(values = c('darkgreen', 'tomato', 'navy')) +
            facet_grid(.~group2, scales = 'free', space = 'free') +
            scale_x_discrete(expand = c(0, 0.5)) +
            theme_bw() +
            theme(
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 18),
                axis.line = element_line(colour = 'black', linewidth = .85),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5 , vjust = 0.5),
                axis.text.y = element_text(size = 12),
                legend.position = 'right') +
            guides(color = guide_legend(title = "PS"))
        if(isTRUE(verbose)) print(prps.map.plot)
        se.obj[[uv.variable]] <- initial.variable
    }

    # Save the results ####
    ## select output name ####
    output.name <- paste0(uv.variable, '|', 'mnn', '|', assay.name)
    if (is.null(prps.group)) {
        prps.group <- paste0('prps|knnMnn|', uv.variable)
    }
    printColoredMessage(message = '-- Saving the PRPS data',
                        color = 'magenta',
                        verbose = verbose)
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
        se.obj@metadata[['PRPS']][['un.supervised']][[prps.group]][['prps.data']][[output.name]] <- prps.data

        printColoredMessage(message = '------------The createPrPsByKnnMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    ## output the PRPS data as matrix ####
    if (isFALSE(save.se.obj)) {
        printColoredMessage(message = '------------The createPrPsByKnnMnn function finished.',
                            color = 'white',
                            verbose = verbose)
        return(prps.data = prps.data)
    }
}
