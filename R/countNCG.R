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
#' @param ncg.selection TTTT
#' @param ncg.group.name TTTT
#' @param create.venn.diagram TTTT
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


# m <- read_fst(path = '../gc_txi.counts.fst')
# dim(m)
# m[1:3,1:3]
# plot(log2(colSums(m)) , col = factor(annot$Library.type))
# colSums(m)[1:3]
# annot <- read.csv('../filtered_samples.csv.csv')
# dim(annot)
# plot(annot$totalcounts, colSums(m))
# dim(m)
# dim(annot)
# table(annot$Library.type)
#
#
#
#
#
#
# set.seed(42)
#
# # Parameters
# n_genes <- 1000
# n_samples_per_group <- 10
# groups <- c("A", "B", "C", "D")
# total_samples <- n_samples_per_group * length(groups)
#
# # Baseline expression
# baseline <- matrix(rnbinom(n_genes * total_samples, mu = 100, size = 1),
#                    nrow = n_genes, ncol = total_samples)
#
# # Add group-specific signals
# expr_matrix <- baseline
# colnames(expr_matrix) <- paste0("Sample_", 1:total_samples)
# rownames(expr_matrix) <- paste0("Gene_", 1:n_genes)
#
# # Index samples by group
# sample_labels <- rep(groups, each = n_samples_per_group)
#
# # Add subtype signals
# # 50 DE genes for each subtype
# genes_A <- 1:50
# genes_B <- 51:100
# genes_C <- 101:150
#
# # Increase expression for subtypes
# expr_matrix[genes_A, sample_labels == "A"] <- expr_matrix[genes_A, sample_labels == "A"] + 100
# expr_matrix[genes_B, sample_labels == "B"] <- expr_matrix[genes_B, sample_labels == "B"] + 100
# expr_matrix[genes_C, sample_labels == "C"] <- expr_matrix[genes_C, sample_labels == "C"] + 100
#
# # Create D as weighted mix of A and B
# A_samples <- which(sample_labels == "A")
# B_samples <- which(sample_labels == "B")
# D_samples <- which(sample_labels == "D")
#
# for (i in seq_along(D_samples)) {
#     a <- expr_matrix[, sample(A_samples, 1)]
#     b <- expr_matrix[, sample(B_samples, 1)]
#     w <- runif(1, 0.3, 0.7)  # random mix weight
#     expr_matrix[, D_samples[i]] <- round(w * a + (1 - w) * b)
# }
#
# # Create metadata
# metadata <- data.frame(
#     sample_id = colnames(expr_matrix),
#     group = sample_labels
# )
#
# # Done: expr_matrix contains simulated counts
#
# library(DESeq2)
# dds <- DESeqDataSetFromMatrix(expr_matrix, colData = metadata, design = ~ group)
# vsd <- vst(dds, blind = TRUE)
#
# plotPCA(vsd, intgroup = "group")
#
#
#
#
#
# set.seed(42)
#
# # Parameters
# n_genes <- 1000
# n_samples_per_group <- 10
# groups <- c("A", "B", "C", "D")
# total_samples <- n_samples_per_group * length(groups)
#
# # Assign genes
# housekeeping_genes <- 201:n_genes
# genes_A <- 1:50
# genes_B <- 51:100
# genes_C <- 101:150
#
# # Simulate baseline expression with gene-specific dispersion
# gene_means <- runif(n_genes, 50, 500)
# gene_dispersion <- runif(n_genes, 0.5, 2)
#
# expr_matrix <- sapply(1:total_samples, function(i) {
#     rnbinom(n_genes, mu = gene_means, size = 1 / gene_dispersion)
# })
# expr_matrix <- round(expr_matrix)
# rownames(expr_matrix) <- paste0("Gene_", 1:n_genes)
# colnames(expr_matrix) <- paste0("Sample_", 1:total_samples)
# sample_labels <- rep(groups, each = n_samples_per_group)
#
# # Batch effect: 2 batches
# batch_labels <- rep(c("batch1", "batch2"), length.out = total_samples)
# batch_effect <- rnorm(n_genes, 1, 0.2)
#
# # Add subtype-specific expression
# expr_matrix[genes_A, sample_labels == "A"] <- expr_matrix[genes_A, sample_labels == "A"] + 150
# expr_matrix[genes_B, sample_labels == "B"] <- expr_matrix[genes_B, sample_labels == "B"] + 150
# expr_matrix[genes_C, sample_labels == "C"] <- expr_matrix[genes_C, sample_labels == "C"] + 150
#
# # Add batch effects (e.g., batch2 samples get small increase)
# expr_matrix[, batch_labels == "batch2"] <- t(t(expr_matrix[, batch_labels == "batch2"]) * batch_effect)
#
# # Mix A+B+C to generate group D
# group_A <- which(sample_labels == "A")
# group_B <- which(sample_labels == "B")
# group_C <- which(sample_labels == "C")
# group_D <- which(sample_labels == "D")
#
# for (i in group_D) {
#     w <- runif(3)  # random weights
#     w <- w / sum(w)  # normalize to 1
#     a <- expr_matrix[, sample(group_A, 1)]
#     b <- expr_matrix[, sample(group_B, 1)]
#     c <- expr_matrix[, sample(group_C, 1)]
#     expr_matrix[, i] <- round(w[1] * a + w[2] * b + w[3] * c)
# }
#
# # Metadata
# metadata <- data.frame(
#     sample_id = colnames(expr_matrix),
#     group = sample_labels,
#     batch = batch_labels
# )
# rownames(metadata) <- metadata$sample_id
#
# # Done
# library(DESeq2)
#
# dds <- DESeqDataSetFromMatrix(round(expr_matrix, digits = 0), colData = metadata, design = ~ batch + group)
# vsd <- vst(dds, blind = TRUE)
#
# plotPCA(vsd, intgroup = c("group", "batch"))
#
#
# library(SummarizedExperiment)
# sce <- zellkonverter::readH5AD("../prostate_portal_300921.h5ad")
# table(sce$batch, sce$group)
# count.data <- assay(x = sce, i = 'soupx_counts')
# sample.annot <- as.data.frame(colData(sce))
# all.equal(row.names(sample.annot), colnames(count.data))
# sample.annot <- dplyr::arrange(sample.annot, 'batch', 'group', 'patient', 'celltype')
# count.data <- count.data[ , row.names(sample.annot)]
# all.equal(row.names(sample.annot), colnames(count.data))
# plot(log2(colSums(count.data)))
#
#
# gg <- rowSums(count.data)
# count.data <- count.data[rowSums(count.data) != 0 , ]
# sum(gg == 0)
#
# library(Seurat)
# index <- sample.annot$group == 'tumor'
# count.data <- count.data[ , index]
# sample.annot <- sample.annot[index , ]
#
# ps <- sapply(
#     unique(sample.annot$patient),
#     function(x){
#         temp.annot <- sample.annot[sample.annot$patient == x , ]
#         rowSums(count.data[ , row.names(temp.annot[temp.annot$celltype == 'LE-KLK4' , ]), drop = F])
#     })
# colnames(ps) <- unique(sample.annot$patient)
# ps <- edgeR::cpm(y = ps, log = T)
# p <- prcomp(t(ps))
# plot(p$x)
#
#
# sum(index)
# s1 <- CreateSeuratObject(counts = count.data[ , index], meta.data = sample.annot[index , ])
# s1 <- NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = 10000)
# s1 <- FindVariableFeatures(s1, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(s1)
# s1 <- ScaleData(s1, features = all.genes)
# s1 <- RunPCA(s1, features = VariableFeatures(object = s1))
# s1 <- FindNeighbors(s1, dims = 1:10)
# s1 <- FindClusters(s1, resolution = 0.5)
# s1 <- RunUMAP(s1, dims = 1:10)
# DimPlot(s1, reduction = "umap", group.by = 'batch', label = TRUE)
# DimPlot(s1, reduction = "umap", group.by = 'group', label = TRUE)
# DimPlot(s1, reduction = "umap", group.by = 'celltype', label = TRUE)
# DimPlot(s1, reduction = "umap", group.by = 'Biopsy.Location.BRIEF', label = TRUE)
# DimPlot(s1, reduction = "umap", group.by = 'name', label = TRUE)
# table(s1$group, s1$celltype)
#










