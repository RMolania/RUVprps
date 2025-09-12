#' Count the current NCG sets in SummarizedExperiment object.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function adds pre-selected sets of negative control genes (NCGs) to a SummarizedExperiment object.
#'
#' @details
#' A pre-selected set of negative control genes (NCGs) will be stored in the following location:
#' se.obj->metadata->NCG->pre.selected->subset.name. These genes can be used for various analyses, including identifying
#' unknown sources of variation, assessing variation, performing RUV normalization, and evaluating normalization steps. The
#' gene set will be stored in: metadata->NCG->pre.selected->subset.name->gene.set
#'
#' @param se.obj A SummarizedExperiment object.
#' @param ncg.selection Character. A character indicating which type of NCG should be selected. The options are: `supervised`
#' and `unsupervised`.
#' @param ncg.group.name Character. A character indicating which group of NCG should be selected.
#' @param create.venn.diagram Logical. Specifies whether to creaet a venn diagram using all NCG sets or not. The default
#' is ser to `TRUE`.
#' @param check.all
#'
#' @return A SummarizedExperiment object with a metadata that contains the NCGs.
#'
#' @importFrom ggVennDiagram ggVennDiagram
#' @export

countNCG <- function(
        se.obj,
        ncg.selection,
        ncg.group.name,
        create.venn.diagram = FALSE,
        check.all = FALSE){

    if (isFALSE(check.all)){
        ncgs.list <- names(se.obj@metadata$NCG[[ncg.selection]][[ncg.group.name]])
        ncgs.list <- ncgs.list[!ncgs.list %in% 'assessment.plot']
        count.ncgs <- sapply(
            ncgs.list,
            function(x){
                sum(se.obj@metadata$NCG[[ncg.selection]][[ncg.group.name]][[x]]$ncg.set)
            })
        count.ncgs
    }
    if (isTRUE(create.venn.diagram)){
        count.ncgs <- lapply(
            ncgs.list,
            function(x){
                row.names(se.obj)[se.obj@metadata$NCG[[ncg.selection]][[ncg.group.name]][[x]]$ncg.set]
            })
        names(count.ncgs) <- ncgs.list
        venn.plot <- ggVennDiagram(count.ncgs, label_alpha = 0) +
            scale_fill_gradient(low = "white", high = "darkgreen") +
            xlab('') +
            ylab('') +
            theme(
                panel.background = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank()
            )
    }
    if (isTRUE(create.venn.diagram)){
        return(list(count.ncgs = count.ncgs, venn.plot = venn.plot))
    } else return(list(count.ncgs = count.ncgs))
}








