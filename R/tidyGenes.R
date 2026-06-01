#' Remove assays from a SummarizedExperiment object
#'
#' @description
#' This function removes assays or subsets rows from a `SummarizedExperiment` object based on gene type, identifiers,
#' or custom row names. It allows filtering while optionally handling duplicate IDs or renaming row identifiers.
#'
#' @param se.obj A `SummarizedExperiment` object.
#' @param keep.gene.type Character or NULL. Specifies which gene type(s) to retain (e.g., `"protein_coding"`). If `NULL`,
#' no filtering by gene type is applied.
#' @param gene.type.col.name Character. The column name in `rowData` that contains the gene type annotation used for filtering.
#' @param remove.duplicates.ids Logical. If `TRUE`, rows with duplicated identifiers in `ids.col.name` will be removed. The
#' default is `FALSE`.
#' @param ids.col.name Character. The column name in `rowData` containing the identifiers to check for duplicates.
#' @param change.row.names Logical. If `TRUE`, replaces the current row names of the `SummarizedExperiment` object with values
#' from `new.row.names`.
#' @param new.row.names Character. A column name in `rowData` containing the values to use as new row names. Required if
#' `change.row.names = TRUE`.
#'
#' @importFrom SummarizedExperiment assays
#'
#' @return A filtered `SummarizedExperiment` object with updated assays and/or row names.
#'
#' @export

tidyGenes <- function(
        se.obj,
        keep.gene.type,
        gene.type.col.name,
        remove.duplicates.ids = TRUE,
        ids.col.name,
        change.row.names      = TRUE,
        new.row.names         = NULL
        ){
    if (class(se.obj)[1] == 'SummarizedExperiment' | class(se.obj)[1] == 'RangedSummarizedExperiment'){
        # Checking the input ####
        gene.annot <- rowData(x = se.obj)
        # subset to protein coding genes and remove duplicated gene names
        to.remove <- gene.annot[[gene.type.col.name]] != keep.gene.type | duplicated(gene.annot[[ids.col.name]])
        se.obj <- se.obj[!to.remove, ]
        if (isTRUE(change.row.names)){
            row.names(se.obj) <- rowData(x = se.obj)[[new.row.names]]
            }
    } else stop("The 'se.obj' is not a class of 'SummarizedExperiment'.")
    return(se.obj)
}


