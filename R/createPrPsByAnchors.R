#' Create PRPS data using integration anchors.

#' @author Ramyar Molania

#' @description
#' This functions employs the FindIntegrationAnchors function from the Seurat R package to create PRPS sets for the
#' RUV-III normalization of RNA-seq data. This function can be used in situations in which biological variation are unknown.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A string indicating the assay name within the SummarizedExperiment object for the
#' creation of PRPS data. The assay must be the one that will be used as data input for the RUV-III-PRPS normalization.
#' @param uv.variable Character. A string specifying the name of columns in the sample annotation of the SummarizedExperiment
#' object. This variable can be categorical or continuous. If a continuous variable is provided, it will be
#' divided into groups using the clustering method.
#' @param clustering.method Character. A string indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param hvg Vector. A vector containing the names of highly variable genes. These genes will be utilized to identify
#' anchor samples across different batches. The default value is set to 'NULL'.
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay before applying
#' log transformation to avoid -Inf for raw counts that are equal to 0. The default is 1.
#' @param anchor.features Numeric. A numeric value indicating the provided number of features to be used in anchor
#' finding. The default is 2000. We refer to the FindIntegrationAnchors R function for more details.
#' @param scale Logical. Whether or not to scale the features provided. Only set to FALSE if you have previously scaled
#' the features you want to use for each object in the object.list.
#' @param min.ps.samples Numeric. The minimum number of samples to be averaged to create a pseudo-sample. The default
#' is 3. The minimum value is 2.
#' @param max.ps.samples Numeric value indicating the maximum number of samples to be averaged for creating a pseudo-sample.
#' The default is 'inf'. Please note that averaging a high number of samples may lead to inaccurate PRPS estimates.
#' @param max.prps.sets Numeric. The maximum number of PRPS sets across batches. The default is 10.
#' @param normalization Character. Indicate which normalization methods should be used before finding the anchors.
#' The options are "LogNormalize" or "SCT". The default is "LogNormalize".
#' @param sct.clip.range Numeric. Numeric of length two specifying the min and max values the Pearson residual will be
#' clipped to. The default is 'NULL'. We refer to the FindIntegrationAnchors R function for more details.
#' @param reduction Character. Indicates which dimensional reduction to perform when finding anchors. The options are
#' "cca": canonical correlation analysis, "rpca": reciprocal PCA, and "rlsi": Reciprocal LSI. The default is "cca".
#' @param l2.norm Logical. Indicates whether to perform L2 normalization on the CCA sample embeddings after dimensional
#' reduction or not. The default is 'TRUE'.
#' @param dims Numeric. Indicates which dimensions to use from the CCA to specify the neighbor search space. The default
#' is 10.
#' @param k.anchor Numeric. How many neighbors (k) to use when picking anchors. The default is set to 2.
#' @param k.filter Numeric. How many neighbors (k) to use when filtering anchors. The default is 20.
#' @param k.score Numeric. How many neighbors (k) to use when scoring anchors. The default is 30.
#' @param max.features Numeric. The maximum number of features to use when specifying the neighborhood search space in
#' the anchor filtering. The default is 30.
#' @param nn.method Character. Method for nearest neighbor finding. Options include: "rann", "annoy". The default is
#' "annoy".
#' @param n.trees Numeric. More trees gives higher precision when using annoy approximate nearest neighbor search.
#' @param eps Numeric. Error bound on the neighbor finding algorithm (from RANN/Annoy).
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Character. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment
#' object or to output the result as a list. By default, it is set to TRUE.
#' @param plot.output. TTT
#' @param output.name Character. A string specifying the name of the output file. If it is 'NULL', the function will
#' select a name based on "paste0(uv.variable, '|', 'anchor', '|', assay.name)".
#' @param prps.group Character. A string specifying the name of the PRPS group. If it is 'NULL', the function will select
#' a name based on "paste0('prps|anchor|', uv.variable)".
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

#' @importFrom Seurat VariableFeatures FindIntegrationAnchors
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom dplyr group_by top_n desc
#' @importFrom Matrix rowMeans
#' @importFrom stats setNames
#' @export

createPrPsByAnchors <- function(
        se.obj,
        assay.name,
        main.uv.variable,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        min.sample.for.ps = 3,
        max.sample.for.ps = 10,
        select.extereme.groups = FALSE,
        other.uv.variables = NULL,
        other.uv.clustering.method = 'kmeans',
        check.prps.connectedness = TRUE,
        nb.other.uv.clusters = 2,
        hvg = NULL,
        min.batches.to.cover = 'all',
        apply.log = TRUE,
        pseudo.count = 1,
        anchor.features = 2000,
        scale = TRUE,
        normalization = "LogNormalize",
        sct.clip.range = NULL,
        reduction = "cca",
        l2.norm = TRUE,
        dims = 1:15,
        k.anchor = 2,
        k.filter = 10,
        k.score = 30,
        max.features = 200,
        nn.method = "annoy",
        n.trees = 50,
        eps = 0,
        assess.se.obj = TRUE,
        remove.na = 'both',
        save.se.obj = TRUE,
        plot.output = TRUE,
        output.name = NULL,
        prps.group = NULL,
        verbose = TRUE
        ) {
    printColoredMessage(message = '------------The createPrPsByAnchors function starts:',
                        color = 'white',
                        verbose = verbose)
    # Checking the inputs #####
    if (is.list(assay.name)) {
        stop('The "assay.name" must be the name of an assay in the SummarizedExperiment object.')
    }
    if (length(assay.name) > 1) {
        stop('The "assay.name" must be a single name of an assay in the SummarizedExperiment object..')
    }
    if (is.null(main.uv.variable)) {
        stop('The "main.uv.variable" cannot be empty. Note, unknown sources of UV can be found by the "identifyUnknownUV" function.')
    }
    if (length(main.uv.variable) > 1) {
        stop('The "main.uv.variable" must be a single variable name in the SummarizedExperiment object.')
    }
    if (!main.uv.variable %in% colnames(colData(se.obj))) {
        stop('The "main.uv.variable" cannot be found in the SummarizedExperiment object')
    }
    if (k.anchor == 0) {
        stop('The k.anchor cannot be 0.')
    }
    if (! normalization %in% c("LogNormalize", "SCT")) {
        stop('The " normalization" must be one of the "LogNormalize" or "SCT".')
    }
    if (!reduction %in% c("cca", "rpca", 'rlsi')) {
        stop('The "reduction" should be one of the "cca", "rpca" or "rlsi"')
    }
    if (sum(hvg %in% row.names(se.obj)) != length(hvg)) {
        stop('All the "hvg" genes are not found in the SummarizedExperiment object.')
    }

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

    # Assessing and grouping the unwanted variable ####
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
            perfix = '_group',
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
        length.variable <- length(unique(initial.variable))
        if (length.variable == 1){
            stop('To create PRPS, the "main.uv.variable" must have at least two groups/levels.')
        } else if (length.variable > 1){
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
    # Creating PRPS data ####
    if (!is.null(other.uv.variables)){
        ## other uv variable is TRUE ####
        ## assessing and grouping the other unwanted variable ####
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
            other.uv = homo.uv.groups,
            sample.ids = colnames(se.obj)
            )

        ## finding subgroups that has enough samples for PRPS ####
        printColoredMessage(
            message = '- Finding subgroups with respect to other unwanted variable(s) that have enough samples for PRPS:',
            color = 'magenta',
            verbose = verbose
            )
        covered.batches <- lapply(
            unique(all.uv.groups$other.uv),
            function(x){
                subgroups.size <- findRepeatingPatterns(
                    vec = all.uv.groups$main.uv[all.uv.groups$other.uv == x],
                    n.repeat = max(min.sample.for.ps, c(max(dims) + 1),  k.anchor, k.filter)
                )
            })
        names(covered.batches) <- unique(all.uv.groups$other.uv)
        printColoredMessage(
            message = paste0(
                '- There are ',
                length(covered.batches),
                ' subgroups that have enough samples for PRPS.'),
            color = 'blue',
            verbose = verbose
        )
        covered.batches.table <- as.data.frame(
            table(all.uv.groups$main.uv, all.uv.groups$other.uv)
        )
        covered.batches.table$selected <- covered.batches.table$Freq >= max(min.sample.for.ps)
        if(isTRUE(plot.output)){
            printColoredMessage(
                message = '- Plotting distribution of samples across subgroups with respect to other unwanted variable(s):',
                color = 'magenta',
                verbose = verbose
            )
        }
        covered.batches.table <- ggplot(covered.batches.table, aes(x = Var2, y = Var1, color = selected)) +
            geom_point(size = 4) +
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
                    min.samples = c(min.sample.for.ps),
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
                    max(min.sample.for.ps),
                    ' samples across all sub-groups of the ',
                    main.uv.variable,
                    ' variable.'),
                color = 'blue',
                verbose = verbose
            )
        }
        selected.covered.batches <- unlist(lapply(
            1:length(covered.batches),
            function(x) length(covered.batches[[x]]) > 1
            ))
        covered.batches <- covered.batches[selected.covered.batches]

        ## finding anchors across pairs of batches ####
        printColoredMessage(
            message = '-- Finding anchors between all the pairs of the subgroups :',
            color = 'magenta',
            verbose = verbose
            )
        ## splitting the data into groups and create SeuratObjects ####
        printColoredMessage(
            message = paste0(
                '- spliting the SummarizedExperiment object into ',
                length(covered.batches),
                ' groups and then select highly variable genes.'),
            color = 'blue',
            verbose = verbose
            )
        all.seurat.objects <- lapply(
            names(covered.batches),
            function(x) {
                all.uv.groups.sub <- all.uv.groups[all.uv.groups[['other.uv']] == x , ]
                all.seurat.objects.sub <- lapply(
                    covered.batches[[x]],
                    function(y){
                        samples.index <- all.uv.groups.sub$sample.ids[all.uv.groups.sub$main.uv == y]
                        seurat.obj <- SeuratObject::CreateSeuratObject(
                            counts = assay(x = se.obj[, samples.index], i = assay.name),
                            project = paste(x, y, '_')
                            )
                        seurat.obj <- Seurat::NormalizeData(
                            object = seurat.obj,
                            normalization.method = normalization,
                            verbose = verbose
                            )
                        if (!is.null(hvg)){
                            Seurat::VariableFeatures(seurat.obj) <- hvg
                            } else if (is.null(hvg)) {
                            seurat.obj <- Seurat::FindVariableFeatures(
                                object = seurat.obj,
                                nfeatures = anchor.features,
                                verbose = verbose
                            )
                        }
                        return(seurat.obj)
                    })
                names(all.seurat.objects.sub) <- covered.batches[[x]]
                all.seurat.objects.sub
            })
        names(all.seurat.objects) <- names(covered.batches)
        all.samples.index <- c(1:ncol(se.obj))
        ## finding anchors  ####
        printColoredMessage(
            message = '- Finding the anchors:',
            color = 'blue',
            verbose = verbose
            )
        all.anchors <- lapply(
            names(all.seurat.objects),
            function(x){
                all.anchors <- Seurat::FindIntegrationAnchors(
                    object.list = all.seurat.objects[[x]],
                    anchor.features = anchor.features,
                    reference = NULL,
                    scale = scale,
                    normalization = normalization,
                    sct.clip.range = sct.clip.range,
                    reduction = reduction,
                    l2.norm = l2.norm,
                    dims = dims,
                    k.anchor = k.anchor,
                    k.filter = k.filter,
                    k.score = k.score,
                    max.features = max.features,
                    nn.method = nn.method,
                    n.trees = n.trees,
                    eps = eps,
                    verbose = verbose
                    )
                score <- NULL
                all.anchors <- all.anchors@anchors
                # excluding repetitive scores ####
                half.data.1 <- all.anchors$score[1:(nrow(all.anchors)/2)]
                half.data.2 <- all.anchors$score[c(c(nrow(all.anchors)/2)+1):nrow(all.anchors)]
                if (all.equal(half.data.1, half.data.2)){
                    all.anchors <- all.anchors[ 1:(nrow(all.anchors)/2), ]
                }
                all.anchors$score <- round(x = all.anchors$score, digits = 3)
                all.anchors <- all.anchors[ , c(1,2,4,5, 3)]
                # renaming the columns ####
                colnames(all.anchors)[1:4] <- c(
                    'dataset.sample.1',
                    'dataset.sample.2',
                    'dataset.index.1',
                    'dataset.index.2'
                    )
                # applying a sanity check ####
                sanity.a <- lapply(
                    unique(all.anchors$dataset.name.1),
                    function(i){
                        max.anchor <- max(all.anchors$dataset.sample.1[all.anchors$dataset.name.1 == i])
                        max.sample <- ncol(all.seurat.objects[[x]][[i]])
                        if (max.anchor > max.sample){
                            stop('There something wrong with the anchors.')
                        }
                    })
                sanity.b <- lapply(
                    unique(all.anchors$dataset.name.2),
                    function(i){
                        max.anchor <- max(all.anchors$dataset.sample.2[all.anchors$dataset.name.2 == i])
                        max.sample <- ncol(all.seurat.objects[[x]][[i]])
                        if (max.anchor > max.sample){
                            stop('There something wrong with the anchors.')
                        }
                    })
                # adding data set names ####
                all.anchors$dataset.name.1 <- names(all.seurat.objects[[x]])[all.anchors$dataset.index.1]
                all.anchors$dataset.name.2 <- names(all.seurat.objects[[x]])[all.anchors$dataset.index.2]
                all.anchors$other.uv <- x
                # adding sample ids ####
                all.anchors$sample.ids.1 <- unlist(lapply(
                    1:nrow(all.anchors),
                    function(i){
                        data1 <- all.seurat.objects[[x]][[all.anchors$dataset.name.1[i]]]
                        colnames(data1)[all.anchors$dataset.sample.1[i]]
                    }))
                all.anchors$sample.ids.2 <- unlist(lapply(
                    1:nrow(all.anchors),
                    function(i){
                        data2 <- all.seurat.objects[[x]][[all.anchors$dataset.name.2[i]]]
                        colnames(data2)[all.anchors$dataset.sample.2[i]]
                    }))
                all.anchors$group <- paste(
                    all.anchors$dataset.name.1,
                    all.anchors$dataset.name.2,
                    all.anchors$other.uv,
                    sep = '.'
                    )
                printColoredMessage(
                    message = paste0(
                        nrow(all.anchors),
                        ' sample pairs (anchors) are found across all the sub-groups of the "',
                        main.uv.variable,
                        '" variable.'),
                    color = 'blue',
                    verbose = verbose
                    )
                all.anchors
            })
        names(all.anchors) <- names(all.seurat.objects)
        all.anchors <- do.call(rbind, all.anchors)
        ## applying a sanity check  ####
        printColoredMessage(
            message = '- Sanity check on the anchors.',
            color = 'blue',
            verbose = verbose
            )
        sanity.a <- lapply(
            all.anchors$sample.ids.1,
            function(x){
                a <- unique(all.anchors$other.uv[all.anchors$sample.ids.1 == x])
                b <- all.uv.groups$other.uv[all.uv.groups$sample.ids == x]
                if (a != b)  stop('There something wrong with the anchors.')
                a <- unique(all.anchors$dataset.name.1[all.anchors$sample.ids.1 == x])
                b <- all.uv.groups$main.uv[all.uv.groups$sample.ids == x]
                if (a != b)  stop('There something wrong with the anchors.')
            })
        sanity.b <- lapply(
            all.anchors$sample.ids.2,
            function(x){
                a <- unique(all.anchors$other.uv[all.anchors$sample.ids.2 == x])
                b <- all.uv.groups$other.uv[all.uv.groups$sample.ids == x]
                if (a != b) stop('ffff')
                a <- unique(all.anchors$dataset.name.2[all.anchors$sample.ids.2 == x])
                b <- all.uv.groups$main.uv[all.uv.groups$sample.ids == x]
                if (a != b) stop('ffff')
            })
        all.prps.sets <- split(
            x = all.anchors,
            f = all.anchors$group
            )
        groups.sets <- names(all.prps.sets)
        ## creating all possible PRPS sets ####
        all.prps.sets <- lapply(
            groups.sets,
            function(i){
                all.prps.sets.sub <- all.prps.sets[[i]]
                all.prps.sets.split <- split(
                    x = all.prps.sets.sub,
                    f = all.prps.sets.sub$sample.ids.1
                    )
                groups.sets.sub <- names(all.prps.sets.split)
                all.prps.sets.sub <- lapply(
                    groups.sets.sub,
                    function(x) {
                        temp.anchors <- all.prps.sets.split[[x]]
                        temp.anchors <- do.call(
                            rbind,
                            lapply(temp.anchors$sample.ids.2, function(j)
                                all.prps.sets.sub[all.prps.sets.sub$sample.ids.2 == j,])
                            )
                        temp.anchors <- rbind(
                            temp.anchors,
                            do.call(rbind,
                                    lapply(temp.anchors$sample.ids.1, function(j)
                                        all.prps.sets.sub[all.prps.sets.sub$sample.ids.1 == j,])
                            ))
                        n <- 3
                        counter <- 1
                        while (length(unique(temp.anchors$dataset.sample.1)) < min.sample.for.ps){
                            print(paste("Iteration:", counter))
                            temp.anchors <- rbind(
                                temp.anchors,
                                do.call(rbind,
                                        lapply(temp.anchors$sample.ids.2, function(j)
                                            all.prps.sets.sub[all.prps.sets.sub$sample.ids.2 == j,])
                                ))
                            temp.anchors
                            if(counter >= n)  break
                            counter <- counter + 1
                        }
                        n <- 3
                        counter <- 1
                        while (length(unique(temp.anchors$dataset.sample.2)) < min.sample.for.ps){
                            print(paste("Iteration:", counter))
                            temp.anchors <- rbind(
                                temp.anchors,
                                do.call(rbind,
                                        lapply(temp.anchors$sample.ids.1, function(j)
                                            all.prps.sets.sub[all.prps.sets.sub$sample.ids.1 == j,])
                                ))
                            temp.anchors
                            if(counter >= n) break
                            counter <- counter + 1
                        }

                        all.datasets <- sort(unique(c(temp.anchors$dataset.name.1, temp.anchors$dataset.name.2)))
                        temp.anchors <- temp.anchors[order(temp.anchors$score, decreasing = TRUE) , ]
                        anchor.sets <- lapply(
                            all.datasets,
                            function(x) {
                                unique(c(
                                    temp.anchors$sample.ids.1[temp.anchors$dataset.name.1 == x],
                                    temp.anchors$sample.ids.2[temp.anchors$dataset.name.2 == x]
                                ))
                            })
                        length.sets <- sapply(anchor.sets, length)
                        if (sum(length.sets >= min.sample.for.ps) != 2 ){
                            return(NULL)
                        } else {
                            # filtering PRPS sets based on PS sample size ####
                            printColoredMessage(
                                message = '- Filter the PRPS sets based on sample size.',
                                color = 'blue',
                                verbose = verbose
                                )
                            anchor.sets <- lapply(
                                1:2,
                                function(e){
                                    if (length(anchor.sets[[e]]) > max.sample.for.ps){
                                        anchor.sets[[e]][1:max.sample.for.ps]
                                    } else anchor.sets[[e]]
                                })
                            anchor.datasets <-
                                unlist(lapply(all.datasets, function(x) {
                                    unique(c(
                                        temp.anchors$dataset.name.1[temp.anchors$dataset.name.1  == x],
                                        temp.anchors$dataset.name.2[temp.anchors$dataset.name.2 == x]
                                    ))
                                }))
                            return(
                                list(
                                    anchor.sets = setNames(anchor.sets, anchor.datasets),
                                    length.sets = length.sets,
                                    anchor.datasets = anchor.datasets))
                        }
                    })
                names(all.prps.sets.sub) <- paste0('PRPSset_', groups.sets.sub)
                all.prps.sets.sub <- Filter(Negate(is.null), all.prps.sets.sub)
                if (length(all.prps.sets.sub) == 0) all.prps.sets.sub = NULL
                all.prps.sets.sub
            })
        names(all.prps.sets) <- groups.sets
        all.prps.sets <- Filter(Negate(is.null), all.prps.sets)
        printColoredMessage(
            message = paste0(
                '- ',  sum(sapply(all.prps.sets, length)),
                ' possible PRPS stes are found.'),
            color = 'blue',
            verbose = verbose
            )
        ## filtering number of PRPS sets ####
        groups.sets <- names(all.prps.sets)
        all.prps.sets <- lapply(
            names(all.prps.sets),
            function(x){
                if (length(all.prps.sets[[x]]) > max.prps.sets){
                    all.prps.sets[[x]][1:max.prps.sets]
                } else all.prps.sets[[x]]
            })
        names(all.prps.sets) <- groups.sets
        ## checking initial coverage ####
        printColoredMessage(
            message = '- Check the distribution of the PRPS sets across the batches.',
            color = 'orange',
            verbose = verbose
            )
        prps.coverage <- lapply(
            names(all.prps.sets),
            function(x){
                unique(unlist(lapply(
                    names(all.prps.sets[[x]]),
                    function(y){
                        all.prps.sets[[x]][[y]]$anchor.datasets
                    })))
            })
        all.batches <- length(unique(se.obj[[main.uv.variable]]))
        prps.coverage.length <- sapply(prps.coverage, length)
        if (max(prps.coverage.length) == all.batches) {
            printColoredMessage(
                message = '- all sub-groups of the main unwanted variable are covered by PRPS data.',
                color = 'blue',
                verbose = verbose
                )
        }
        if (max(prps.coverage.length) < all.batches) {
            printColoredMessage(
                message = '- all sub-groups of the main unwanted variable are not covered by PRPS data across each sub-groups of other unwanted variable.',
                color = 'blue',
                verbose = verbose
                )
            if (isTRUE(check.prps.connectedness)){
                printColoredMessage(
                    message = '- Assessing the connection between differet sets of PRPS:',
                    color = 'blue',
                    verbose = verbose
                )
            }
            if (isFALSE(check.prps.connectedness)){
                printColoredMessage(
                    message = '- Assessing the connection between differet sets of PRPS:',
                    color = 'blue',
                    verbose = verbose
                )
            }

        }
        ## creating PRPS data ####
        # data transformation and normalization ####
        printColoredMessage(
            message = '-- Data transformation and normalization:',
            color = 'magenta',
            verbose = verbose
        )
        ## apply log ####
        if (isTRUE(apply.log) & !is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0('- apply log2 the ', assay.name,' + ', pseudo.count, ' data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name) + pseudo.count)
        } else if (isTRUE(apply.log) & is.null(pseudo.count)) {
            printColoredMessage(
                message = paste0('- apply log2 on the ', assay.name,' data.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- log2(assay(x = se.obj, i = assay.name))
        } else if (isFALSE(apply.log)) {
            printColoredMessage(
                message = paste0('The ', assay.name, ' data will be used without any log transformation.'),
                color = 'blue',
                verbose = verbose
            )
            expr.data <- assay(x = se.obj, i = assay.name)
        }
        prps.data <- lapply(
            names(all.prps.sets),
            function(x) {
                prps.sets <- all.prps.sets[[x]]$anchor.sets
                temp.data <- lapply(
                    1:length(prps.sets),
                    function(y) {
                        if(length(prps.sets[[y]]) >= min.ps.samples){
                            rowMeans(expr.data[, prps.sets[[y]], drop = FALSE])
                        }
                    })
                temp.data <- do.call(cbind, temp.data)
                colnames(temp.data) <- rep(
                    x = paste0(uv.variable, '.unsu.', x),
                    ncol(temp.data)
                )
                return(temp.data)
            })
        prps.data <- do.call(cbind, prps.data)
        if (sum(table(colnames(prps.data)) == 1)) {
            stop( 'There are something wrong with the prps.data. All the column names of the prps.data are the same.')
        }
        if(sum(is.na(prps.data)) !=0){
            stop( 'There NA in the PRPS data, please check the data input and parameters.')
        }
        se.obj[[uv.variable]] <- initial.variable

    }


    # Saving the output ####
    ## select output name ####
    if(is.null(output.name))
        output.name <- paste0(uv.variable, '|', 'anchor', '|', assay.name)
    if (is.null(prps.group))
        prps.group <- paste0('prps|anchor|', uv.variable)
    if (isTRUE(save.se.obj)) {
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

        printColoredMessage(message = '------------The createPrPsByAnchors function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)
    }
    if (isFALSE(save.se.obj)){
        printColoredMessage(message = '------------The createPrPsByAnchors function finished.',
                            color = 'white',
                            verbose = verbose)
        return(list(prps.data = prps.data))
    }
}



