#' Creates unsupervised pseudo-replicates of pseudo samples (PRPS).

#' @author Ramyar Molania

#' @description
#' This function uses the 'createUnSupervisedPRPSbyAnchors' or 'createUnSupervisedPRPSbyMNN' to create PRPS sets of all
#' the unwanted variables in situations when the biological variation are unknown.

#' @param se.obj A SummarizedExperiment object.
#' @param assay.name Symbol. A symbol indicating the assay name within the SummarizedExperiment object for the creation
#' of PRPS data. The assay must be the one that will be used as data input for the RUV-III-PRPS normalization.
#' @param uv.variables Symbol. A symbol or a vector of symbols specifying the name(s) of column(s) in the sample annotation
#' of the SummarizedExperiment object. This variable can to be categorical or continuous variable. If a continuous variable
#' is provide, this will be divided into groups using the clustering method.
#' @param approach Symbol. Indicates which method to be used. The options are 'anchor' and 'mnn'.
#' @param data.input Symbol. Indicates which data should be used an input for finding the k nearest neighbors data. Options
#' include: 'expr' and 'pcs'. If 'pcs' is selected, the first PCs of the data will be used as input. If 'expr' is selected,
#' the data will be selected as input.
#' @param clustering.method Symbol.A symbol indicating the choice of clustering method for grouping the 'uv.variable'
#' if a continuous variable is provided. Options include 'kmeans', 'cut', and 'quantile'. The default is set to 'kmeans'.
#' @param nb.clusters Numeric. A numeric value indicating how many clusters should be found if the 'uv.variable' is a
#' continuous variable. The default is 3.
#' @param other.uv.variables TTTT
#' @param other.uv.clustering.method TTT
#' @param nb.other.uv.clusters TTT
#' @param hvg Vector. A vector containing the names of highly variable genes. These genes will be utilized to identify
#' anchor samples across different batches. The default value is set to 'NULL'
#' @param select.extreme.groups TTTT
#' @param min.batches.to.cover TTTT
#' @param filter.prps.sets TTTT
#' @param max.prps.sets Numeric. The maximum number of PRPS sets across batches. The default in 10.
#' @param min.sample.for.ps Numeric. The minimum number of samples to be averaged to create a pseudo-sample. The default
#' is 3. The minimum value is 2.
#' @param max.sample.for.ps Numeric value indicating the maximum number of samples to be averaged for creating a pseudo-sample.
#' The default is 'inf'. Please note that averaging a high number of samples may lead to inaccurate PRPS estimates.
#' @param nb.knn Numeric.The maximum number of nearest neighbors to compute. The default is set 3.
#' @param nb.mnn Numeric.The maximum number of nearest neighbors to compute. The default is set 3.
#' @param anchor.features Numeric. A numeric value indicating the provided number of features to be used in anchor finding.
#' The default is 2000. We refer to the FindIntegrationAnchors R function for more details.
#' @param scale Logical. Whether or not to scale the features provided. Only set to FALSE if you have previously scaled
#' the features you want to use for each object in the object.list.
#' @param sct.clip.range Numeric. Numeric of length two specifying the min and max values the Pearson residual will be
#' clipped to. The default is 'NULL'. We refer to the FindIntegrationAnchors R function for more details.
#' @param reduction Symbol. Indicates which dimensional reduction to perform when finding anchors. The options are "cca":
#' canonical correlation analysis, "rpca": reciprocal PCA and "rlsi": Reciprocal LSI. The default is "cca".
#' @param l2.norm Logical. Indicates whether to perform L2 normalization on the CCA samples embeddings after dimensional
#' reduction or not. The default is 'TRUE'.
#' @param dims Numeric. Indicates which dimensions to use from the CCA to specify the neighbor search space. Th default
#' is 10.
#' @param k.anchor Numeric. How many neighbors (k) to use when picking anchors. Th default is 3.
#' @param k.filter Numeric. How many neighbors (k) to use when filtering anchors. Th default is 200.
#' @param k.score Numeric. How many neighbors (k) to use when scoring anchors. Th default is 30.
#' @param max.features Numeric. The maximum number of features to use when specifying the neighborhood search space in
#' the anchor filtering.Th default is 30.
#' @param nn.method Symbol.Method for nearest neighbor finding. Options include: "rann", "annoy". The defauly is "annoy".
#' @param n.trees Numeric. More trees gives higher precision when using annoy approximate nearest neighbor search
#' @param eps Numeric. Error bound on the neighbor finding algorithm (from RANN/Annoy)
#' @param nb.pcs Numeric. Indicates the number PCs should be used as data input for finding the k nearest neighbors. The
#' 'nb.pcs' must be set when the "data.input = PCs". The default is 2.
#' @param center Logical. Indicates whether to scale the data or not. If center is TRUE, then centering is done by
#' subtracting the column means of the assay from their corresponding columns. The default is TRUE.
#' @param scale Logical. Indicates whether to scale the data or not before applying SVD.  If scale is TRUE, then scaling
#' is done by dividing the (centered) columns of the assays by their standard deviations if center is TRUE, and the root
#' mean square otherwise. The default is FALSE
#' @param svd.bsparam A BiocParallelParam object specifying how parallelization should be performed. The default is bsparam().
#' We refer to the 'runSVD' function from the BiocSingular R package.
#' @param normalization Symbol. Indicates which normalization methods should be applied before finding the knn. The default
#' is 'cpm'. If is set to NULL, no normalization will be applied.
#' @param apply.cosine.norm TTTT
#' @param regress.out.variables Symbols. Indicates the columns names that contain biological variables in the
#' SummarizedExperiment object. These variables will be regressed out from the data before finding genes that are highly
#' affected by unwanted variation variable. The default is NULL, indicates the regression will not be applied.
#' @param check.prps.connectedness TTT
#' @param apply.log Logical. Indicates whether to apply a log-transformation to the data or not. The default is 'TRUE'.
#' @param pseudo.count Numeric. A value as a pseudo count to be added to all measurements of the assay before applying
#' log transformation to avoid -Inf for raw counts that are equal to 0. The default is 1.
#' @param mnn.bpparam Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' . The default is SerialParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.
#' @param mnn.nbparam Symbol. A BiocParallelParam object specifying how parallelization should be performed to find MNN.
#' . The default is KmknnParam(). We refer to the 'findMutualNN' function from the 'BiocNeighbors' R package.
#' @param assess.se.obj Logical. Indicates whether to assess the SummarizedExperiment object or not. See the checkSeObj
#' function for more details.
#' @param remove.na Symbol. To remove NA or missing values from the assays or not. The options are 'assays' and 'none'.
#' The default is "assays", so all the NA or missing values from the assay(s) will be removed before computing RLE. See
#' the checkSeObj function for more details.
#' @param output.name Symbol.
#' @param prps.group Symbol.
#' @param plot.output TTTT
#' @param save.se.obj Logical. Indicates whether to save the RLE results in the metadata of the SummarizedExperiment
#' object or to output the result as list. By default it is set to TRUE.
#' @param verbose Logical. If 'TRUE', shows the messages of different steps of the function.

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
        k.filter = 20,
        k.score = 5,
        max.features = 200,
        nn.method = "annoy",
        n.trees = 50,
        eps = 0,
        nb.pcs = 2,
        center = TRUE,
        scale = FALSE,
        svd.bsparam = bsparam(),
        normalization = 'CPM',
        apply.cosine.norm = FALSE,
        regress.out.variables = NULL,
        check.prps.connectedness = FALSE,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        assess.se.obj = TRUE,
        remove.na = 'both',
        output.name = NULL,
        prps.group = NULL,
        plot.output = TRUE,
        save.se.obj = TRUE,
        verbose = TRUE
        ){
    if (approach == 'anchor') {
        for (i in uv.variables) {
            se.obj <- createPrPsByAnchors(
                se.obj = se.obj,
                assay.name = assay.name,
                main.uv.variable = i,
                clustering.method = clustering.method,
                nb.clusters = nb.clusters,
                filter.prps.sets = filter.prps.sets,
                max.prps.sets = max.prps.sets,
                min.sample.for.ps = min.sample.for.ps,
                max.sample.for.ps = max.sample.for.ps,
                select.extreme.groups = select.extreme.groups,
                other.uv.variables = other.uv.variables,
                other.uv.clustering.method = other.uv.clustering.method,
                nb.other.uv.clusters = nb.other.uv.clusters,
                check.prps.connectedness = check.prps.connectedness,
                hvg = hvg,
                min.batches.to.cover = min.batches.to.cover,
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
                remove.na = remove.na,
                plot.output = plot.output,
                output.name = output.name,
                prps.group = prps.group,
                assess.se.obj = assess.se.obj,
                save.se.obj = save.se.obj,
                verbose = verbose
            )
        }
    }
    if (approach == 'knn.mnn') {
        for (i in uv.variables) {
            se.obj <- createPrPsByKnnMnn (
                se.obj = se.obj,
                assay.name = assay.name,
                main.uv.variable = i,
                clustering.method = clustering.method,
                nb.clusters = nb.clusters,
                filter.prps.sets = filter.prps.sets,
                max.prps.sets = max.prps.sets,
                select.extreme.groups = select.extreme.groups,
                other.uv.variables = other.uv.variables,
                other.uv.clustering.method = other.uv.clustering.method,
                nb.other.uv.clusters = nb.other.uv.clusters,
                min.sample.for.ps = min.sample.for.ps,
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
                assess.se.obj = assess.se.obj,
                remove.na = remove.na,
                plot.output = plot.output,
                output.name = output.name,
                prps.group = prps.group,
                save.se.obj = save.se.obj,
                verbose = verbose
            )
        }
    }
    if (approach == 'mnn') {
        for (i in uv.variables) {
            se.obj <- createPrPsByMnn(
                se.obj = se.obj,
                assay.name = assay.name,
                main.uv.variable = i,
                clustering.method = clustering.method,
                nb.clusters = nb.clusters,
                filter.prps.sets = filter.prps.sets,
                max.prps.sets = max.prps.sets,
                select.extreme.groups = select.extreme.groups,
                other.uv.variables = other.uv.variables,
                other.uv.clustering.method = other.uv.clustering.method,
                nb.other.uv.clusters = nb.other.uv.clusters,
                min.sample.for.ps = min.sample.for.ps,
                nb.mnn = nb.mnn,
                min.batches.to.cover = min.batches.to.cover,
                hvg = hvg,
                normalization = normalization,
                apply.cosine.norm = apply.cosine.norm,
                regress.out.variables = regress.out.variables,
                check.prps.connectedness = check.prps.connectedness,
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
        }
    }
    return(se.obj)
}
