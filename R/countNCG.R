#' Count and check the overlap the NCG sets in SummarizedExperiment object.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function counts and evaluates the overlap of the NGCs stored in the metadata of the SummarizedExperiment object.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param ncg.selection Character. A character string indicating which type of NCG should be selected. The options are:
#' `supervised` and `unsupervised`.
#' @param ncg.group.name Character. A character string indicating which group of NCGs should be selected. Group names are
#' specified in the functions that identifies NCGs. If set to `NULL`, all current NCG group will be selected.
#' @param create.venn.diagram Logical. Specifies whether to create a venn diagram using all NCG sets or not. The default
#' is set to `TRUE`.
#' @param verbose Logical. Indicates whether to display output messages during function execution. The default is set to
#' `TRUE`.
#'
#' @export

countNCG <- function(
        se.obj,
        ncg.selection,
        ncg.group.name = NULL,
        create.venn.diagram = FALSE,
        verbose = TRUE
        ){
    if (!is.character(ncg.selection) | length(ncg.selection) > 1){
        stop('The "ncg.selection" must be one of the "supervised" or "unsupervised".')
    }
    if (!ncg.selection %in% c('supervised', 'unsupervised')){
        stop('The "ncg.selection" must be one of the "supervised" or "unsupervised".')
    }
    if (!is.null(ncg.group.name)){
        ncgs.list <- names(se.obj@metadata$NCG[[ncg.selection]][[ncg.group.name]])
        ncgs.list <- ncgs.list[!ncgs.list %in% 'assessment.plot']
        count.ncgs <- sapply(
            ncgs.list,
            function(x){
                sum(se.obj@metadata$NCG[[ncg.selection]][[ncg.group.name]][[x]]$ncg.set)
            })
        count.ncgs
    }
    if (is.null(ncg.group.name)){
        ncg.group.name <- names(se.obj@metadata$NCG[[ncg.selection]])
        ncgs.set.names <- lapply(
            ncg.group.name,
            function(x){
                ncgs.list <- names(se.obj@metadata$NCG[[ncg.selection]][[x]])
                ncgs.list <- ncgs.list[!ncgs.list %in% 'assessment.plot']
                ncgs.list
            })
        count.ncgs <- sapply(
            ncg.group.name,
            function(x){
                ncgs.list <- names(se.obj@metadata$NCG[[ncg.selection]][[x]])
                ncgs.list <- ncgs.list[!ncgs.list %in% 'assessment.plot']
                ncgs.list <- sapply(
                    ncgs.list,
                    function(y){
                        sum(se.obj@metadata$NCG[[ncg.selection]][[x]][[y]]$ncg.set)
                    })
                ncgs.list
            })
        count.ncgs
    }
    if (isTRUE(create.venn.diagram)){
        if (!requireNamespace("ggVennDiagram", quietly = TRUE)){
            stop("Package 'ggVennDiagram' is required for creating a Venn diagram. Please install it.")
        }
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








