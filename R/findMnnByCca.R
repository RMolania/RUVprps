





findMnnByCca <- function(
        se.obj,
        assay.name,
        main.uv.variable,
        nb.cca = 5,
        reference.group = NULL,
        hvg = NULL,
        scale = TRUE,
        regress.out.rle.med = FALSE,
        sample.to.use = 'all',
        nb.mnn = 3,
        clustering.method = 'kmeans',
        nb.clusters = 3,
        other.uv.variables = NULL,
        other.uv.clustering.method = 'kmeans',
        nb.other.uv.clusters = 2,
        normalization = 'CPM',
        cosine.norm = FALSE,
        regress.out.variables = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        mnn.bpparam = SerialParam(),
        mnn.nbparam = KmknnParam(),
        check.se.obj = TRUE,
        remove.na = 'both',
        plot.output = TRUE,
        cca.set.name = NULL,
        mnn.group.name = NULL,
        mnn.sets.name = NULL,
        save.se.obj = TRUE,
        verbose = TRUE
){

}
