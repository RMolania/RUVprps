#' Estimate PAM50 subtypes
#'
#' @author Ramyar Molania
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character or character vector. The name(s) of the assay(s) in the `SummarizedExperiment` object.
#' The default is set to 'all'.
#' @param entrez.gene.id description
#' @param probe description
#' @param tissue.type description
#' @param groups description
#' @param apply.log description
#' @param pseudo.count description
#' @param normalization description
#' @param regress.out.variables description
#' @param regress.out.rle.med description
#' @param out.put.name description
#' @param check.se.obj description
#' @param remove.na description
#' @param verbose description
#'
#' @importFrom SummarizedExperiment colData
#' @import genefu

estimatePAM50 <- function(
        se.obj,
        assay.names = 'all',
        entrez.gene.id,
        probe,
        tissue.type = NULL,
        groups = NULL,
        apply.log = TRUE,
        pseudo.count = 1,
        normalization = 'CPM',
        regress.out.variables = NULL,
        regress.out.rle.med = FALSE,
        out.put.name = NULL,
        check.se.obj = TRUE,
        remove.na = 'none',
        verbose = TRUE
        ){
    gene.annot <- as.data.frame(rowData(se.obj))
    egi <- which(colnames(gene.annot) == entrez.gene.id)
    colnames(gene.annot)[egi] <- 'EntrezGene.ID'
    probe.col <- which(colnames(gene.annot) == probe)
    colnames(gene.annot)[probe.col] <- 'probe'
    pam50.genfu <- lapply(
        names(assays(se.obj)),
        function(x){
            expr.data <- preProcessData(
                se.obj = se.obj,
                assay.name = x,
                normalization = normalization,
                regress.out.variables = regress.out.variables,
                regress.out.rle.med = regress.out.rle.med,
                apply.log = apply.log,
                pseudo.count = pseudo.count,
                check.se.obj = FALSE,
                remove.na = 'none',
                verbose = verbose
                )
            row.names(expr.data)[row.names(expr.data) == 'NDC80'] <-'KNTC2'
            row.names(expr.data)[row.names(expr.data) == 'NUF2'] <- 'CDCA1'
            row.names(expr.data)[row.names(expr.data) == 'ORC6'] <- 'ORC6L'
            gene.annot$probe[gene.annot$probe == 'NDC80'] <- 'KNTC2'
            gene.annot$probe[gene.annot$probe == 'NUF2'] <- 'CDCA1'
            gene.annot$probe[gene.annot$probe == 'ORC6'] <- 'ORC6L'
            genefu::molecular.subtyping(
                sbt.model = "pam50",
                data = t(expr.data),
                annot = gene.annot,
                do.mapping = TRUE)
        })
    names(pam50.genfu) <- paste0(names(assays(se.obj)), '.PAM50')
    for(i in names(pam50.genfu)){
        SummarizedExperiment::colData(se.obj)[i] <- as.character(pam50.genfu[[i]]$subtype)
        index <- SummarizedExperiment::colData(se.obj)[[i]] == 'Normal'
        SummarizedExperiment::colData(se.obj)[[i]][index] <- 'Normal.like'
        if (!is.null(tissue.type)){
            index.normal <- SummarizedExperiment::colData(se.obj)[[tissue.type]] == 'Normal'
            SummarizedExperiment::colData(se.obj)[[i]][index.normal] <- 'Normal'
            SummarizedExperiment::colData(se.obj)[[i]] <- factor(
                x = SummarizedExperiment::colData(se.obj)[ , i],
                levels = c('Basal','Her2','LumA','LumB', 'Normal.like', 'Normal'))
        } else {
            SummarizedExperiment::colData(se.obj)[[i]] <- factor(
                x = SummarizedExperiment::colData(se.obj)[ , i],
                levels = c('Basal','Her2','LumA','LumB', 'Normal.like'))
        }
    }
    return(se.obj)
}


# usethis::use_data(pam50.robust, overwrite = TRUE)
# usethis::use_r("pam50.robust.R")
# data(pam50.robust)
