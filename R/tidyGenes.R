#' Remove assays from a SummarizedExperiment object
#'
#' @param se.obj A SummarizedExperiment object.
#' @param keep.gene.type TTTT
#' @param gene.type.col.name TTTT
#' @param remove.duplicates.ids TTTT
#' @param ids.col.name TTTT
#' @param chaneg.row.names TTTT
#' @param new.row.names TTTT
#'
#' @importFrom SummarizedExperiment assays
#'
tidyGenes <- function(
        se.obj,
        keep.gene.type,
        gene.type.col.name,
        remove.duplicates.ids = TRUE,
        ids.col.name,
        change.row.names = TRUE,
        new.row.names = NULL
        ){
    # Checking the input ####
    gene.annot <- rowData(x = se.obj)
    # subset to protein coding genes and remove duplicated gene names
    to.remove <- gene.annot[[gene.type.col.name]] != keep.gene.type | duplicated(gene.annot[[ids.col.name]])
    se.obj <- se.obj[!to.remove, ]
    if (isTRUE(change.row.names)){
        row.names(se.obj) <- rowData(x = se.obj)[[new.row.names]]
    }
    return(se.obj)
}


