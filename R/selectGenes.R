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
#' @param select.tp.genes TTTTT
#' @param method TTTTT
#' @param a TTTTT
#' @param rho TTTTT
#' @param apply.log TTTTT
#' @param pseudo.count TTTTT
#' @param check.se.object TTTTT
#' @param remove.na TTTTT
#' @param verbose TTTTT

selectGenes <- function(
        se.obj,
        assay.names,
        variables,
        groups = NULL,
        cor.cutoff = 0.5,
        abs.cor = TRUE,
        top.genes = NULL,
        select.tp.genes = FALSE,
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
    if (isTRUE(select.tp.genes)){
        purity.gene.set <- as.data.frame(rowData(read.se.obj)[c('immune.gene.signature', 'stromal.gene.signature')])
        purity.gene.set[is.na(purity.gene.set)] <- FALSE
        purity.gene.set <- c(
            row.names(purity.gene.set)[purity.gene.set$immune.gene.signature],
            row.names(purity.gene.set)[purity.gene.set$stromal.gene.signature]
        )
        purity.gene.set <- list(
            'Tumour purity' = list(
                upset.genes = purity.gene.set,
                downset.genes = NULL)
        )
    }
    all.gene.list <- list(
        non.tp.genes = all.corr,
        tp.genes = purity.gene.set
    )
    return(all.gene.list)
}

