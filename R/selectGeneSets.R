#' Selects genes for partial correlation and gene scoring analysis.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function adds the TCGA RNA-seq batch information to the a SummarizedExperiment object.
#'
#' @param se.obj TTT
#' @param gene.set TTT
#' @param set.name TTT

selectGenesSets <- function(
        se.obj,
        gene.set = 'immune.stromal',
        set.name = NULL
        ){
    if (gene.set == 'immune.stromal'){
        purity.gene.set <- as.data.frame(rowData(se.obj)[c('immune.gene.signature', 'stromal.gene.signature')])
        purity.gene.set[is.na(purity.gene.set)] <- FALSE
        purity.gene.set <- c(
            row.names(purity.gene.set)[purity.gene.set$immune.gene.signature],
            row.names(purity.gene.set)[purity.gene.set$stromal.gene.signature]
        )
        purity.gene.set <- list(
            'Tumour.purity' = list(
                upset.genes = purity.gene.set,
                downset.genes = NULL)
        )
    }
    return(gene.sets = purity.gene.set)
}
