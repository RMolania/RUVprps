#' Selects genes for partial correlation and gene scoring analysis.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function adds the TCGA RNA-seq batch information to the a SummarizedExperiment object.
#'
#' @param se.obj TTT
#' @param assay.names TTT
#' @param variables TTTT
#' @param groups TTTTsss
#' @param cor.cutoff TTTT
#' @param abs.cor TTTT
#' @param top.genes TTTTT
#' @param method TTTTT
#' @param a TTTTT
#' @param rho TTTTT
#' @param apply.log TTTTT
#' @param pseudo.count TTTTT
#' @param check.se.object TTTTT
#' @param remove.na TTTTT
#' @param verbose TTTTT

selectGenesForPPcorr <- function(
        se.obj,
        assay.names,
        variables,
        groups = NULL,
        cor.cutoff = 0.5,
        abs.cor = TRUE,
        top.genes = NULL,
        method = "spearman",
        a = 0.05,
        rho = 0,
        apply.log = TRUE,
        pseudo.count = 1,
        check.se.object = TRUE,
        remove.na = 'both',
        verbose = TRUE
        ){
    names(variables) <- assay.names
    groups <- se.obj[[groups]]
    if (!is.null(groups)){
        all.corr <- list()
        for(i in 1:length(variables)){
            corr.genes <- sapply(
                levels(groups),
                function(x){
                    gene.var.corr <- computeGenesVariableCorrelation(
                        se.obj = se.obj[ , groups == x],
                        assay.names = names(variables)[i],
                        variable = variables[i],
                        method = method,
                        a = a, rho = rho,
                        apply.log = apply.log,
                        pseudo.count = pseudo.count,
                        check.se.obj = FALSE,
                        save.se.obj = FALSE
                        )[[names(variables)[i]]]
                    if (isTRUE(abs.cor)){
                        abs(gene.var.corr[ , 'correlation']) > cor.cutoff
                    } else gene.var.corr[ , 'correlation'] > cor.cutoff
                })
            all.corr[[variables[i]]] <- row.names(corr.genes)[rowSums(corr.genes)!= 0]
        }
    }
    return(all.corr)
}



