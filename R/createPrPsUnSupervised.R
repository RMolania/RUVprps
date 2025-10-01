#' Creates unsupervised PRPS data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function uses the `createUnSupervisedPRPSbyAnchors()` or `createUnSupervisedPRPSbyMNN()` to create PRPS sets of all
#' the unwanted variables in situations when the biological variation is unknown.
#'
#' @details
#' Additional details...
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Character. A character indicating the assay name within the SummarizedExperiment object for
#' the creation of PRPS data. The assay must be the one that will be used as data input for the RUV-III-PRPS normalization.
#' @param uv.variables Character. A character or vector of strings specifying the name(s) of column(s) in the
#' sample annotation of the SummarizedExperiment object. This variable can be categorical or continuous. If a continuous
#' variable is provided, it will be divided into groups using the clustering method.
#' @param other.uv.variables Character. Additional unwanted variables to be included alongside `uv.variables`.
#' These may provide alternative axes of unwanted variation. Default is `NULL`.
#' @param approach Character. Indicates which method to be used. Options are `anchor` or `mnn`.
#' @param nb.pcs Numeric. Number of principal components to use as input for nearest neighbor search. Must be set if
#' `data.input = pcs`. Default is 2.
#' @param nb.cca Numeric. Number of canonical correlation components to compute when using the anchor approach. Default is 10.
#' @param data.input Character. Indicates which data should be used as input for finding nearest neighbors. Options:
#' `expr` or `pcs`. If `pcs` is selected, the first PCs of the data will be used as input. If `expr` is selected,
#' expression values will be used directly.
#' @param coordinates.to.use Character or integer. Specifies which dimensions (e.g., PCs or CCA components) to use
#' for neighbor search. Default is `NULL` (use all available).
#' @param samples.to.use Character, logical, or numeric. A subset of samples to be included in PRPS construction.
#' Default is `all`.
#' @param min.sample.for.ps Numeric. Minimum number of samples to average to create a pseudo-sample. Default is 3 (minimum 2).
#' @param max.sample.for.ps Numeric. Maximum number of samples to average for creating a pseudo-sample.
#' Default is `Inf`. Larger values may dilute unwanted variation signals.
#' @param select.extreme.groups Logical. If `TRUE`, selects only extreme groups (e.g., lowest and highest clusters) when
#' `uv.variables` is continuous. Default is `TRUE`.
#' @param filter.prps.sets Logical. Indicates whether to filter PRPS sets based on size or variability before saving.
#' Default is `FALSE`.
#' @param max.prps.sets Numeric. Maximum number of PRPS sets across batches. Default is 10.
#' @param min.batches.to.cover Numeric. Minimum number of batches each PRPS set must cover to be retained. Default is 2.
#' @param cover.all.batches Logical. If `TRUE`, requires each PRPS set to cover all batches. Default is `FALSE`.
#' @param check.prps.connectedness Logical. If `TRUE`, checks whether PRPS sets form a connected structure across batches.
#' Default is `TRUE`. See `checkPRPSconnectedness()`.
#' @param reference.group Character. Specifies a group within `uv.variables` to be used as the reference in pairwise
#' comparisons. Default is `NULL` (all groups compared).
#' @param hvg Character vector. Names of highly variable genes to use in anchor/MNN identification. Default is `NULL`
#' (all genes used).
#' @param scale.cca Logical. Indicates whether to scale features before applying CCA. Default is `TRUE`.
#' @param apply.ruviii.norm Logical. If `TRUE`, applies RUV-III normalization to pseudo-samples. Default is `FALSE`.
#' @param use.ruviii.norm.for.mnn Logical. If `TRUE`, applies RUV-III normalization prior to MNN search. Default is `FALSE`.
#' @param ncg Numeric. Number of negative control genes to include for RUV-III normalization. Default is `NULL`.
#' @param k Numeric. Number of nearest neighbors for kNN-based similarity. Default is 20.
#' @param center Logical. Whether to mean-center features before dimension reduction. Default is `TRUE`.
#' @param bio.dims Numeric. Number of biological dimensions to retain after removing unwanted variation. Default is `NULL`.
#' @param min.ps Numeric. Minimum number of pseudo-samples required per group. Default is 1.
#' @param nb.knn Numeric. Number of nearest neighbors to compute. Default is 3.
#' @param nb.mnn Numeric. Number of mutual nearest neighbors to compute. Default is 3.
#' @param min.nb.for.mnn Numeric. Minimum number of mutual nearest neighbors required to establish an anchor link. Default is 1.
#' @param similarity.approach Character. Method to compute similarity between samples/anchors. Options: `correlation`,
#' `euclidean`, or `cosine`. Default is `correlation`.
#' @param data.for.similarity Character. Specifies which data representation to use for similarity (e.g., `expr`, `pcs`).
#' Default is `expr`.
#' @param clustering.method Character. Clustering method for continuous `uv.variables`. Options: `kmeans`, `cut`, `quantile`.
#' Default is `kmeans`.
#' @param nb.clusters Numeric. Number of clusters for continuous `uv.variables`. Default is 3.
#' @param other.uv.clustering.method Character. Clustering method for `other.uv.variables`. Default matches `clustering.method`.
#' @param nb.other.uv.clusters Numeric. Number of clusters for `other.uv.variables`. Default matches `nb.clusters`.
#' @param nb.batches.to.cover Numeric. Minimum number of batches that must be represented in final PRPS sets. Default is 2.
#' @param svd.bsparam A `BiocParallelParam` object defining parallelization for SVD. Default: `bsparam()` (BiocSingular).
#' @param normalization Character. Normalization method prior to neighbor search. Default is `cpm`. If `NULL`, no normalization.
#' @param cosine.norm Logical. Whether to apply cosine normalization before MNN search. Default is `TRUE`.
#' @param regress.out.variables Character vector. Names of biological variables to regress out before PRPS estimation.
#' Default is `NULL`.
#' @param regress.out.rle.med Logical. If `TRUE`, regresses out median RLE values. Default is `FALSE`.
#' @param mnn.bpparam A `BiocParallelParam` object for parallel MNN search. Default: `SerialParam()`.
#' @param mnn.nbparam A nearest neighbor search parameter object, e.g., `KmknnParam()`. Default: `KmknnParam()`.
#' @param apply.log Logical. Whether to log-transform expression values. Default is `TRUE`.
#' @param apply.log.for.prps Logical. Whether to log-transform pseudo-samples after averaging. Default is `FALSE`.
#' @param pseudo.count Numeric. Value added prior to log-transform to avoid `-Inf`. Default is 1.
#' @param assess.variables.association Logical. Whether to assess associations between unwanted variables and PRPS sets.
#' Default is `FALSE`.
#' @param create.prps.map Logical. If `TRUE`, generates a PRPS similarity map across batches. Default is `FALSE`.
#' @param check.se.obj Logical. If `TRUE`, checks that `se.obj` is valid using `checkSeObj()`. Default is `TRUE`.
#' @param remove.na Character. Whether to remove `NA` values. Options: `assays`, `sample.annotation`, `both`, `none`.
#' Default is `assays`.
#' @param cca.set.name Character. Name to assign to CCA results in metadata. Default is constructed automatically.
#' @param knn.group.name Character. Name of the current KNN group. Default is auto-generated from `uv.variables`.
#' @param knn.sets.name Character. Name of the KNN set saved in metadata. Default is constructed automatically.
#' @param mnn.group.name Character. Name of the current MNN group. Default is auto-generated from `uv.variables`.
#' @param mnn.sets.name Character. Name of the MNN set saved in metadata. Default is constructed automatically.
#' @param prps.group.name Character. Name of the PRPS group. Default is auto-generated from `uv.variables`.
#' @param prps.sets.name Character. Name of the PRPS sets saved in metadata. Default is constructed automatically.
#' @param plot.output Logical. Whether to generate diagnostic plots of PRPS sets. Default is `FALSE`.
#' @param save.se.obj Logical. If `TRUE`, saves results in `metadata` of the SummarizedExperiment object. Otherwise,
#' returns results as a list. Default is `TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default is `TRUE`.

#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr count
#' @importFrom tidyr %>%
#' @export

createPrPsUnSupervised <- function(
        se.obj,
        assay.name,
        uv.variables,
        other.uv.variables = NULL,
        approach = 'cca',
        data.input = 'pcs',
        coordinates.to.use = 'both',
        nb.cca = 2,
        nb.pcs = 5,
        samples.to.use = 'all',
        min.sample.for.ps = 3,
        max.sample.for.ps = 10,
        select.extreme.groups = FALSE,
        filter.prps.sets = TRUE,
        max.prps.sets = 3,
        min.batches.to.cover = 'all',
        cover.all.batches = FALSE,
        check.prps.connectedness = FALSE,
        reference.group = NULL,
        hvg = NULL,
        scale.cca = TRUE,
        apply.ruviii.norm = TRUE,
        use.ruviii.norm.for.mnn = TRUE,
        ncg = NULL,
        k = 2,
        min.ps = 10,
        bio.dims = 3,
        nb.mnn = 3,
        nb.knn = 2,
        center = TRUE,
        min.nb.for.mnn = 1,
        similarity.approach = 'euclidean',
        data.for.similarity = 'ruv',
        clustering.method = 'cut',
        nb.clusters = 3,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        nb.batches.to.cover = 2,
        svd.bsparam = bsparam(),
        normalization = 'CPM',
        cosine.norm = FALSE,
        regress.out.variables = NULL,
        regress.out.rle.med = FALSE,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        apply.log = TRUE,
        apply.log.for.prps = TRUE,
        pseudo.count = 1,
        assess.variables.association = TRUE,
        create.prps.map = FALSE,
        check.se.obj = TRUE,
        remove.na = 'both',
        cca.set.name = NULL,
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
    if (approach == 'cca') {
        for (i in uv.variables) {
            if (!is.null(other.uv.variables)) {
                other.uv.variables <- uv.variables[!uv.variables %in% i]
            } else other.uv.variables <- NULL
            se.obj <- createPrPsUnSupervisedByCca(
                se.obj = se.obj,
                assay.name = assay.name,
                main.uv.variable = i,
                reference.group = reference.group,
                other.uv.variables = other.uv.variables,
                coordinates.to.use = coordinates.to.use,
                nb.cca = nb.cca,
                nb.pcs = nb.pcs,
                samples.to.use = samples.to.use,
                min.sample.for.ps = min.sample.for.ps,
                select.extreme.groups = select.extreme.groups,
                filter.prps.sets = filter.prps.sets,
                max.prps.sets = max.prps.sets,
                min.batches.to.cover = min.batches.to.cover,
                cover.all.batches = cover.all.batches,
                check.prps.connectedness = check.prps.connectedness,
                hvg = hvg,
                scale.cca = scale.cca,
                apply.ruviii.norm = apply.ruviii.norm,
                use.ruviii.norm.for.mnn = use.ruviii.norm.for.mnn,
                ncg = ncg,
                k = k ,
                nb.mnn = nb.mnn,
                min.ps = min.ps,
                min.nb.for.mnn = min.nb.for.mnn,
                similarity.approach = similarity.approach,
                data.for.similarity = data.for.similarity,
                clustering.method = clustering.method,
                nb.clusters = nb.clusters,
                other.uv.clustering.method = other.uv.clustering.method,
                nb.other.uv.clusters = nb.other.uv.clusters,
                nb.batches.to.cover = nb.batches.to.cover,
                normalization = normalization,
                cosine.norm = cosine.norm,
                regress.out.variables = regress.out.variables,
                regress.out.rle.med = regress.out.rle.med,
                apply.log = apply.log,
                apply.log.for.prps = apply.log.for.prps,
                pseudo.count = pseudo.count,
                assess.variables.association = assess.variables.association,
                create.prps.map = create.prps.map,
                plot.output = plot.output,
                mnn.bpparam = mnn.bpparam,
                mnn.nbparam = mnn.nbparam,
                check.se.obj = check.se.obj,
                remove.na = remove.na,
                cca.set.name = cca.set.name,
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
            se.obj <- createPrPsUnsupervisedByKnnMnn(
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
                apply.cosine.norm = cosine.norm,
                regress.out.variables = regress.out.variables,
                apply.log = apply.log,
                apply.log.for.prps = apply.log.for.prps,
                pseudo.count = pseudo.count,
                mnn.bpparam = mnn.bpparam,
                mnn.nbparam = mnn.nbparam,
                samples.to.use = samples.to.use,
                check.se.obj = check.se.obj,
                remove.na = remove.na,
                plot.output = plot.output,
                knn.group.name = knn.group.name,
                knn.sets.name = knn.sets.name,
                mnn.group.name = mnn.group.name,
                mnn.sets.name = mnn.sets.name,
                prps.group.name = prps.group.name,
                prps.sets.name = prps.sets.name,
                save.se.obj = save.se.obj,
                verbose = verbose
            )
        }
    }
    if (approach == 'mnn') {
        for (i in uv.variables) {
            if (!is.null(other.uv.variables)) {
                other.uv.variables <- uv.variables[!uv.variables %in% i]
            } else other.uv.variables <- NULL
            createPrPsUnsupervisedByMnn(
                se.obj = se.obj,
                assay.name = assay.name,
                main.uv.variable = i,
                other.uv.variables = other.uv.variables,
                data.input = data.input,
                bio.dims = bio.dims,
                clustering.method = clustering.method,
                nb.clusters = nb.clusters,
                other.uv.clustering.method = other.uv.clustering.method,
                nb.other.uv.clusters = nb.other.uv.clusters,
                select.extreme.groups = select.extreme.groups,
                min.sample.for.ps = min.sample.for.ps,
                filter.prps.sets = filter.prps.sets,
                max.prps.sets = max.prps.sets,
                min.batches.to.cover = min.batches.to.cover,
                check.prps.connectedness = check.prps.connectedness,
                nb.mnn = nb.mnn,
                hvg = hvg,
                samples.to.use = samples.to.use,
                normalization = normalization,
                apply.cosine.norm = cosine.norm,
                regress.out.variables = regress.out.variables,
                regress.out.rle.med = regress.out.rle.med,
                apply.log = apply.log,
                apply.log.for.prps = apply.log.for.prps,
                pseudo.count = pseudo.count,
                mnn.bpparam = mnn.bpparam,
                mnn.nbparam = mnn.nbparam,
                create.prps.map = create.prps.map,
                plot.output = plot.output,
                prps.group.name = prps.group.name,
                prps.sets.name = prps.sets.name,
                check.se.obj = check.se.obj,
                remove.na = remove.na,
                save.se.obj = save.se.obj,
                verbose = verbose
            )
        }
    }
    return(se.obj)
}
