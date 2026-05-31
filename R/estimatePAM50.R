#' Estimate PAM50 subtypes
#'
#' @author Ramyar Molania
#'
#' @description
#' This function estimates PAM50 molecular subtypes from RNA-seq or microarray data stored in a `SummarizedExperiment` object
#' using the `genefu` package. It supports preprocessing steps such as normalization, log-transformation, and regression
#' of unwanted variables.
#'
#' @param se.obj A SummarizedExperiment object.
#' @param assay.names Character or character vector. The name(s) of the assay(s) in `se.obj` to use. Default is `"all"`.
#' @param entrez.gene.id Character. Column name in `rowData(se.obj)` containing Entrez gene IDs for mapping to PAM50 genes.
#' @param probe Character. Column name in `rowData(se.obj)` containing probe or gene symbols if Entrez IDs are not used.
#' @param tissue.type Character. Tissue type for samples (e.g., `"breast"`) to select the appropriate PAM50 reference.
#' @param groups Character or NULL. Optional grouping variable(s) in `colData(se.obj)` to perform subtype estimation within groups.
#' @param apply.log Logical. If `TRUE`, log2-transforms the assay data prior to subtype estimation. Default is `TRUE`.
#' @param pseudo.count Numeric. Pseudo-count added to assay values before log2 transformation to avoid `-Inf`. Default is 1.
#' @param normalization Character. Normalization method applied before subtype estimation (e.g., `"CPM"`, `"TMM"`). Default is `"CPM"`.
#' @param regress.out.variables Character vector. Column names in `colData(se.obj)` representing unwanted variation variables to regress out. Default is `NULL`.
#' @param regress.out.rle.med Logical. If `TRUE`, regress out the median relative log expression (RLE) per sample. Default is `FALSE`.
#' @param out.put.name Character. Name to assign to the results stored in the `SummarizedExperiment` metadata.
#' @param check.se.obj Logical. If `TRUE`, validates the structure of `se.obj` before analysis. Default is `TRUE`.
#' @param remove.na Character. Determines how to handle missing values. Options: `"assays"`, `"sample.annotation"`, `"both"`, or `"none"`. Default is `"both"`.
#' @param verbose Logical. If `TRUE`, prints progress messages during subtype estimation. Default is `TRUE`.
#'
#' @return A `SummarizedExperiment` object with PAM50 subtype assignments stored in the `colData`.
#'
#' @importFrom SummarizedExperiment colData
#' @import genefu
#'
#' @export

estimatePAM50 <- function(
        se.obj,
        assay.names           = 'all',
        entrez.gene.id,
        probe,
        tissue.type           = NULL,
        groups                = NULL,
        apply.log             = TRUE,
        pseudo.count          = 1,
        normalization         = 'CPM',
        regress.out.variables = NULL,
        regress.out.rle.med   = FALSE,
        out.put.name          = NULL,
        check.se.obj          = TRUE,
        remove.na             = 'none',
        verbose               = TRUE
        ){
    gene.annot                      <- as.data.frame(rowData(se.obj))
    egi                             <- which(colnames(gene.annot) == entrez.gene.id)
    colnames(gene.annot)[egi]       <- 'EntrezGene.ID'
    probe.col                       <- which(colnames(gene.annot) == probe)
    colnames(gene.annot)[probe.col] <- 'probe'

    # Checking the assays ####
    if (length(assay.names) == 1 && assay.names == 'all') {
        assay.names <- factor(x = names(assays(se.obj)), levels = names(assays(se.obj)))
    } else  assay.names <- factor(x = assay.names , levels = assay.names)
    if(!sum(assay.names %in% names(assays(se.obj))) == length(assay.names)){
        stop('The "assay.names" cannot be found in the SummarizedExperiment object.')
    }
    pam50.genfu <- lapply(
        levels(assay.names),
        function(x){
            expr.data <- preProcessData(
                se.obj                = se.obj,
                assay.name            = x,
                normalization         = normalization,
                regress.out.variables = regress.out.variables,
                regress.out.rle.med   = regress.out.rle.med,
                apply.log             = apply.log,
                pseudo.count          = pseudo.count,
                check.se.obj          = FALSE,
                remove.na             = 'none',
                verbose               = verbose
                )
            row.names(expr.data)[row.names(expr.data) == 'NDC80'] <-'KNTC2'
            row.names(expr.data)[row.names(expr.data) == 'NUF2']  <- 'CDCA1'
            row.names(expr.data)[row.names(expr.data) == 'ORC6']  <- 'ORC6L'
            gene.annot$probe[gene.annot$probe == 'NDC80']         <- 'KNTC2'
            gene.annot$probe[gene.annot$probe == 'NUF2']          <- 'CDCA1'
            gene.annot$probe[gene.annot$probe == 'ORC6']          <- 'ORC6L'
            genefu::molecular.subtyping(
                sbt.model  = "pam50",
                data       = t(expr.data),
                annot      = gene.annot,
                do.mapping = TRUE)
        })
    names(pam50.genfu) <- paste0(assay.names, '.PAM50')
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
        rm(i)
    }
    return(se.obj)
}
