#' Creates unsupervised PRPS data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function uses the `createUnSupervisedPRPSbyAnchors()` or `createUnSupervisedPRPSbyMNN()` to create PRPS sets of all
#' the unwanted variables in situations when the biological variation are unknown.
#'
#' @details
#' Additional details...
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character indicating the assay name within the SummarizedExperiment object for
#' the creation of PRPS data. The assay must be the one that will be used as data input for the RUV-III-PRPS normalization.
#' @param uv.variables Character. A character or vector of strings specifying the name(s) of column(s) in the
#' sample annotation of the SummarizedExperiment object. This variable can be categorical or continuous. If a continuous
#' variable is provided, this will be divided into groups using the clustering method.
#' @param approach Character. Indicates which method to be used. The options are `anchor` and `mnn`.
#' @param data.input Character. Indicates which data should be used as input for finding the k nearest neighbors. Options
#' include: `expr` and `pcs`. If `pcs` is selected, the first PCs of the data will be used as input. If `expr` is selected,
#' the data will be selected as input.
#' @param clustering.method Character. A character indicating the choice of clustering method for grouping the `uv.variable`
#' if a continuous variable is provided. Options include `kmeans`, `cut`, and `quantile`. The default is set to `kmeans`.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the `uv.variable` is a
#' continuous variable. The default is 3.
#' @param other.uv.variables TTTT
#' @param other.uv.clustering.method TTT
#' @param nb.other.uv.clusters TTT
#' @param hvg Vector. A vector containing the names of highly variable genes. These genes will be utilized to identify
#' anchor samples across different batches. The default value is set to `NULL`
#' @param select.extreme.groups Logical. Indicates whether to select only the extreme groups e.g., highest and lowest
#' clusters, when the `main.uv.variable` is a continuous variable. Default is set to `TRUE`. This will increase the
#' variation between PR sets in order to better capture the unwanted variation.
#' @param min.batches.to.cover TTTT
#' @param filter.prps.sets TTTT
#' @param max.prps.sets Numeric. The maximum number of PRPS sets across batches. The default is 10.
#' @param min.sample.for.ps Numeric. The minimum number of samples to be averaged to create a pseudo-sample. The default
#' is 3. The minimum value is 2.
#' @param max.sample.for.ps Numeric value indicating the maximum number of samples to be averaged for creating a pseudo-sample.
#' The default is `inf`. Please note that averaging a high number of samples may lead to inaccurate PRPS estimates.
#' @param nb.knn Numeric. The maximum number of nearest neighbors to compute. The default is 3.
#' @param nb.mnn Numeric. The maximum number of mutual nearest neighbors to compute. The default is 3.
#' @param anchor.features Numeric. A numeric value indicating the number of features to be used in anchor finding.
#' The default is 2000. We refer to the `FindIntegrationAnchors` R function for more details.
#' @param scale Logical. Whether or not to scale the features provided. Only set to FALSE if you have previously scaled
#' the features you want to use for each object in the object list.
#' @param sct.clip.range Numeric. A numeric vector of length two specifying the min and max values the Pearson residual will be
#' clipped to. The default is `NULL`. See `FindIntegrationAnchors` for more.
#' @param reduction Character. Indicates which dimensional reduction to perform when finding anchors. Options are `cca`,
#' `rpca`, and `rlsi`. The default is `cca`.
#' @param l2.norm Logical. Indicates whether to perform L2 normalization on the CCA sample embeddings after dimensional
#' reduction. Default is `TRUE`.
#' @param dims Numeric. Indicates which dimensions to use from CCA to define the neighbor search space. Default is 10.
#' @param k.anchor Numeric. Number of neighbors (k) used to select anchors. Default is 3.
#' @param k.filter Numeric. Number of neighbors (k) used for anchor filtering. Default is 200.
#' @param k.score Numeric. Number of neighbors (k) used when scoring anchors. Default is 30.
#' @param max.features Numeric. Maximum number of features used for neighborhood search space in anchor filtering. Default is 30.
#' @param nn.method Character. Method for nearest neighbor finding. Options: `rann`, `annoy`. Default is `annoy`.
#' @param n.trees Numeric. Number of trees for `annoy` search (more trees = higher precision).
#' @param eps Numeric. Error bound on the neighbor finding algorithm (from RANN or Annoy).
#' @param nb.pcs Numeric. Number of principal components to use as input for nearest neighbor search. Must be set if
#' `data.input = pcs`. Default is 2.
#' @param center Logical. Whether to center the data (subtract column means). Default is `TRUE`.
#' @param scale Logical. Whether to scale the data before applying SVD. Default is `FALSE`.
#' @param svd.bsparam A `BiocParallelParam` object defining parallelization for SVD. Default: `bsparam()`. See `runSVD`
#' from the BiocSingular package.
#' @param normalization Character. Normalization method before computing nearest neighbors. Default: `cpm`. If `NULL`,
#' no normalization is applied.
#' @param apply.cosine.norm Logical. Indicates whether cosine normalization should be applied before finding MNN. Default
#' is set to `TRUE`.
#' @param check.prps.connectedness Logical. Indicates whether to assess the `connectedness` between the PRPS sets across
#' all batches. Default is set to `TRUE`, indicating if there is not connections between all PRPS sets across all batches,
#' the function will give error.We refer to the checkPRPSconnectedness() function for more details.
#' @param regress.out.variables Character. Column names containing biological variables to regress out before estimating
#'  unwanted variation. Default: `NULL`.
#' @param apply.log Logical. Whether to apply a log transformation. Default is `TRUE`.
#' @param pseudo.count Numeric. Value added to counts before log transform to avoid `-Inf`. Default is 1.
#' @param mnn.bpparam Character. A `BiocParallelParam` object for parallel MNN search. Default: `SerialParam()`. See `findMutualNN`.
#' @param mnn.nbparam Character. A `BiocParallelParam` object for nearest neighbor search. Default: `KmknnParam()`.
#' @param check.se.obj Logical. Whether to assess the SummarizedExperiment object using `checkSeObj`.
#' @param remove.na Character. Indicates whether to remove `NA` values. Options: `assays`, `none`. Default is `assays`.
#' @param knn.group.name Character. A character specifying the name of the knn  to which the current KNN belong.
#' If set to `NULL`, the function will automatically assign a name using  `main.uv.variable`.
#' @param knn.sets.name  Character. A character specifying the name of the knn set names to be saved in the metadata
#' of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0(uv.variable, '|', assay.name)`.
#' @param mnn.group.name Character. A character specifying the name of the mnn to which the current MNN belong.
#' If set to `NULL`, the function will automatically assign a name using  `main.uv.variable`.
#' @param mnn.sets.name A character specifying the name of the mnn set names to be saved in the metadata of the
#' SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' @param prps.group.name Character. A character specifying the name of the prps.group.name to which the current KNN belong.
#' If set to `NULL`, the function will automatically assign a name using  `main.uv.variable`.
#' @param prps.sets.name Character. A character specifying the name of the output file to be saved in the metadata
#' of the SummarizedExperiment object. If set to `NULL`, the function will select a name based on
#' `paste0(uv.variable, '|', assay.name)`.
#' @param plot.output TTTT
#' @param save.se.obj Logical. If `TRUE`, saves results in `metadata` of the SummarizedExperiment object. Otherwise,
#' returns results as a list. Default is `TRUE`.
#' @param verbose Logical. If `TRUE`, prints messages throughout function execution.
#' @param bio.dims TTT
#' @param sample.to.use TTT
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr count
#' @importFrom tidyr %>%
#' @export

createPrPsUnSupervised <- function(
        se.obj,
        assay.name,
        uv.variables,
        approach = 'mnn',
        data.input = 'expr',
        bio.dims = 50,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        other.uv.variables = NULL,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 3,
        hvg = NULL,
        select.extreme.groups = FALSE,
        min.batches.to.cover = 'all',
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        min.sample.for.ps = 3,
        max.sample.for.ps = 10,
        nb.mnn = 1,
        nb.knn = 2,
        anchor.features = 2000,
        sct.clip.range = NULL,
        reduction = "cca",
        l2.norm = TRUE,
        dims = 1:10,
        k.anchor = 3,
        k.filter = 5,
        k.score = 5,
        max.features = 200,
        nn.method = "annoy",
        n.trees = 50,
        eps = 0,
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        sample.to.use = 'all',
        svd.bsparam = bsparam(),
        normalization = 'CPM',
        apply.cosine.norm = FALSE,
        regress.out.variables = NULL,
        check.prps.connectedness = FALSE,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        check.se.obj = TRUE,
        remove.na = 'both',
        knn.group.name = NULL,
        knn.sets.name = NULL,
        mnn.group.name = NULL,
        mnn.sets.name = NULL,
        prps.group.name = NULL,
        prps.sets.name = NULL,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    if (approach == 'anchor') {
        for (i in uv.variables) {
            if (!is.null(other.uv.variables)) {
                other.uv.variables <- uv.variables[!uv.variables %in% i]
            } else other.uv.variables <- NULL
            se.obj <- createPrPsByAnchors(
                se.obj = se.obj,
                assay.name = assay.name,
                main.uv.variable = i,
                clustering.method = clustering.method,
                nb.clusters = nb.clusters,
                other.uv.variables = other.uv.variables,
                other.uv.clustering.method = other.uv.clustering.method,
                nb.other.uv.clusters = nb.other.uv.clusters,
                min.sample.for.ps = min.sample.for.ps,
                max.sample.for.ps = max.sample.for.ps,
                filter.prps.sets = filter.prps.sets,
                max.prps.sets = max.prps.sets,
                min.batches.to.cover = min.batches.to.cover,
                select.extreme.groups = select.extreme.groups,
                check.prps.connectedness = check.prps.connectedness,
                hvg = hvg,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                anchor.features = anchor.features,
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
                check.se.obj = check.se.obj,
                remove.na = remove.na,
                plot.output = plot.output,
                prps.group.name = prps.group.name,
                prps.sets.name = prps.sets.name,
                save.se.obj = save.se.obj,
                verbose = verbose
                )
        }
    }
    if (approach == 'knn.mnn') {
        for (i in uv.variables) {
            if (!is.null(other.uv.variables)) {
                other.uv.variables <- uv.variables[!uv.variables %in% i]
            } else other.uv.variables <- NULL
            se.obj <- createPrPsByKnnMnn (
                se.obj = se.obj,
                assay.name = assay.name,
                main.uv.variable = i,
                clustering.method = clustering.method,
                nb.clusters = nb.clusters,
                other.uv.variables = other.uv.variables,
                other.uv.clustering.method = other.uv.clustering.method,
                nb.other.uv.clusters = nb.other.uv.clusters,
                min.sample.for.ps = min.sample.for.ps,
                select.extreme.groups = select.extreme.groups,
                filter.prps.sets = filter.prps.sets,
                max.prps.sets = max.prps.sets,
                min.batches.to.cover = min.batches.to.cover,
                check.prps.connectedness = check.prps.connectedness,
                data.input = data.input,
                nb.pcs = nb.pcs,
                center = center,
                scale = scale,
                svd.bsparam = svd.bsparam,
                nb.knn = nb.knn,
                nb.mnn = nb.mnn,
                hvg = hvg,
                normalization = normalization,
                apply.cosine.norm = apply.cosine.norm,
                regress.out.variables = regress.out.variables,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                mnn.bpparam = mnn.bpparam,
                mnn.nbparam = mnn.nbparam,
                check.se.obj = check.se.obj,
                remove.na = remove.na,
                plot.output = plot.output,
                knn.group.name = knn.group.name,
                knn.sets.name = knn.sets.name,
                mnn.group.name = mnn.group.name,
                mnn.sets.name = mnn.sets.name,
                prps.group.name = prps.group.name,
                prps.sets.name = prps.sets.name,
                save.se.obj = save.se.obj
                )
        }
    }
    if (approach == 'mnn') {
        for (i in uv.variables) {
            if (!is.null(other.uv.variables)) {
                other.uv.variables <- uv.variables[!uv.variables %in% i]
            } else other.uv.variables <- NULL
            se.obj <- createPrPsByMnn(
                se.obj = se.obj,
                assay.name = assay.name,
                bio.dims = bio.dims,
                data.input = data.input,
                main.uv.variable = i,
                clustering.method = clustering.method,
                nb.clusters = nb.clusters,
                other.uv.variables = other.uv.variables,
                other.uv.clustering.method = other.uv.clustering.method,
                nb.other.uv.clusters = nb.other.uv.clusters,
                min.sample.for.ps = min.sample.for.ps,
                filter.prps.sets = filter.prps.sets,
                select.extreme.groups = select.extreme.groups,
                max.prps.sets = max.prps.sets,
                min.batches.to.cover = min.batches.to.cover,
                check.prps.connectedness = check.prps.connectedness,
                nb.mnn = nb.mnn,
                hvg = hvg,
                sample.to.use = sample.to.use,
                normalization = normalization,
                apply.cosine.norm = apply.cosine.norm,
                regress.out.variables = regress.out.variables,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                mnn.bpparam = mnn.bpparam,
                mnn.nbparam = mnn.nbparam,
                check.se.obj = check.se.obj,
                remove.na = remove.na,
                plot.output = plot.output,
                prps.group.name = prps.group.name,
                prps.sets.name = prps.sets.name,
                save.se.obj = save.se.obj,
                verbose = verbose
            )
        }
    }
    return(se.obj)
}
