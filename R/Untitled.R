# all.data <- qs::qread(file = '/Users/molania.r/Documents/Current_Projects/Project_DFCI_2_RUVprpsApplication/Application_RUVprps/ProstateCancerAtlas_BulkRNASeq_RawCount_SampleAnnot_12012025.qs')
# expr.data <- all.data$expr.data
# sample.annot <- all.data$sample.annot
#
# library(edgeR)
# expressed.genes <- sapply(
#     unique(sample.annot$dataset),
#     function(x){
#         rowSums(expr.data[ ,sample.annot$dataset == x ])
#     })
# expr.data <- expr.data[rowSums(expressed.genes == 0) == 0 , ] # 33937  1365
#
# ## lowly expressed genes
# keep.genes <- filterByExpr(
#     y = expr.data,
#     group = sample.annot$Molecular.Type
# )
# sum(keep.genes) # 32300
# expr.data <- expr.data[keep.genes , ] # 21744542
#
#
#
# sample.annot$Molecular.Type[sample.annot$Molecular.Type == 'NORMAL'] <- 'Normal'
# sample.annot$Molecular.Type[sample.annot$Molecular.Type == 'PRIMARY'] <- 'Primary'
# sample.annot$Molecular.Type[sample.annot$Molecular.Type == 'Normal' & sample.annot$dataset == 'GTEx'] <- 'Healthy'
#
# sample.annot$Tissues <- 'Cancer'
# sample.annot$Tissues[sample.annot$Molecular.Type == 'PRIMARY'] <- 'Cancer'
# sample.annot$Tissues[sample.annot$Molecular.Type == 'Healthy'] <- 'Normal'
# sample.annot$Tissues[sample.annot$Molecular.Type == 'Normal'] <- 'Normal'
#
#
#
#
# ncg <- singscore::getStableGenes(n_stable = 1000)
# ncg <- intersect(ncg, row.names(expr.data)) # 994
# length(ncg) == 994
#
# ### creating SE obj
# row.names(sample.annot) <- sample.annot$sample.ids
# row.names(expr.data)[is.na(row.names(expr.data))] <- 'na.genes'
# all.equal(colnames(expr.data), row.names(sample.annot))
# gene.annot <- data.frame(hgnc_symbol = row.names(expr.data))
# pad.se <- SummarizedExperiment::SummarizedExperiment(
#     assays = list(RawCount = expr.data),
#     colData = sample.annot,
#     rowData = gene.annot
# )
# pad.se <- prepareSeObj(
#     data = pad.se,
#     raw.count.assay.name = 'RawCount',
#     remove.lowly.expressed.genes = FALSE,
#     estimate.tumor.purity = 'both',
#     assay.name.to.estimate.purity = 'RawCount',
#     scale.singscore.values = FALSE,
#     calculate.library.size = TRUE,
#     create.gene.annotation = TRUE,
#     add.immun.stroma.genes = TRUE,
#     add.gene.details = TRUE,
#     add.housekeeping.genes = TRUE,
#     column.name = 'gene_name',
#     gene.group = 'hgnc_symbol'
# ) # dim: 32300 1365
# pad.se <- identifyUnknownUV(
#     se.obj = pad.se,
#     assay.name = 'RawCount',
#     approach = 'rle',
#     nb.pcs = 3,
#     apply.log = TRUE,
#     ncg = ncg,
#     chronological.detection = TRUE,
#     clustering.methods = 'mclust',
#     assess.bio.association = FALSE,
#     assess.uv.association = FALSE,
#     generate.association.plot = FALSE,
#     add.to.sample.annotation = TRUE,
#     col.name = 'Estimated.batches',
#     output.name = 'iteration1',
#     cpt.minseglen = 3
# )
#
# pccat.fig <- 'PCCAT/'
# pdf(
#     paste0(pccat.fig, 'RUVprps_PCCAT_EstimateUnWantedVariation_RLEapproach_NCG.pdf'),
#     height = 5,
#     width = 10
# )
# pad.se@metadata$UnKnownUV$RawCount$iteration1$plots$batch.plot + ylab('RLE medians') +
#     ylim(-4,4)
# dev.off()
#
#
# ## confusion matrix
# conf.mat <- table(
#     pad.se@metadata$UnKnownUV$RawCount$iteration1$batches ,
#     pad.se$dataset
# )
# conf.mat <- as.matrix(conf.mat)
# max.col <- apply(conf.mat, 1, which.max)
# df <- data.frame(
#     batch = rownames(conf.mat),
#     matched.col = colnames(conf.mat)[max.col],
#     value = mapply(function(r, c) conf.mat[r, c], rownames(conf.mat), max.col)
# )
# df <- df[order(df$value, decreasing = TRUE), ]
# conf.mat.ordered <- conf.mat[df$batch, unique(df$matched.col), drop = FALSE]
#
# library(grid)
# pdf(
#     paste0(pccat.fig, 'RUVprps_PCCAT_EstimateUnWantedVariation_ConfusionMatrix_DatasetsAndEstimated.pdf'),
#     height = 5,
#     width = 10
# )
# ComplexHeatmap::Heatmap(
#     matrix = log2(conf.mat.ordered + 1),
#     name = 'Frequency',
#     cluster_columns = F,
#     cluster_rows = F,
#     show_row_dend = F,
#     show_column_dend = F,
#     col = c('white', 'tomato'),
#     heatmap_legend_param = list(
#         title = "Frequency",
#         title_gp = gpar(fontsize = 14, fontface = "bold"),
#         labels_gp = gpar(fontsize = 14),
#         legend_height = unit(5, "cm"),
#         legend_width = unit(1.2, "cm")
#     )
# )
# dev.off()
#
#
# ## confusion matrix
# conf.mat <- table(pad.se$Molecular.Type , pad.se$dataset)
# conf.mat <- as.matrix(conf.mat)
#
# pdf(
#     paste0(pccat.fig, 'RUVprps_PCCAT_ConfusionMatrix_DatasetsAndMolecularTypes.pdf'),
#     height = 5,
#     width = 10
# )
# ComplexHeatmap::Heatmap(
#     matrix = log2(conf.mat + 1) ,
#     name = 'Frequency',
#     cluster_columns = T,
#     cluster_rows = T,
#     show_row_dend = F,
#     show_column_dend = F,
#     col = c('white', 'tomato'),
#     heatmap_legend_param = list(
#         title = "Frequency",
#         title_gp = gpar(fontsize = 14, fontface = "bold"),
#         labels_gp = gpar(fontsize = 14),
#         legend_height = unit(5, "cm"),
#         legend_width = unit(1.2, "cm")                       # make legend wider (if horizontal)
#     )
# )
# dev.off()
#
#
#
#
#
# library(SummarizedExperiment)
# pad.se$library.size
# gene.annot <- data.frame(rowData(pad.se))
# keep.genes <- gene.annot$gene_biotype == 'protein_coding'
# keep.genes[is.na(keep.genes)] <- FALSE
# pad.se <- pad.se[keep.genes , ]
#
#
#
#
#
#
#
#
#
# gene.annot <- rowData(pad.se)
# hk.normal <- row.names(gene.annot)[gene.annot$single.cell.rnaseq.hk.genes == TRUE]
# hk.tumours <- singscore::getStableGenes(n_stable = 2000)
# hk <- intersect(hk.normal, hk.tumours)
# hk <- intersect(hk, row.names(pad.se)) # 466
#
# data.sets <- unique(pad.se$dataset)
# pca.un.esti <- lapply(
#     data.sets,
#     function(x){
#         identifyUnknownUV(
#             se.obj = pad.se[ , pad.se$dataset == x],
#             assay.name = 'RawCount',
#             ncg = hk,
#             clustering.methods = "mclust",
#             approach = 'pca',
#             nb.pcs = 2,
#             save.se.obj = FALSE)
#     })
# names(pca.un.esti) <- data.sets
# pca.un.esti.data <- lapply(
#     data.sets,
#     function(x){
#         data.frame(
#             batch = pca.un.esti[[x]]$batches,
#             samples = row.names(pca.un.esti[[x]]$input.data)
#         )
#     })
# pca.un.esti.data <- do.call(
#     rbind,
#     pca.un.esti.data
# )
# row.names(pca.un.esti.data) <- pca.un.esti.data$samples
# pca.un.esti.data <- pca.un.esti.data[colnames(pad.se) , ]
# pad.se$pca.batch <- pca.un.esti.data$batch
# pad.se$pca.batch <- paste0(pad.se$pca.batch,'_' ,pad.se$dataset)
# pad.se <- orderSeObj(
#     se.obj = pad.se,
#     factors.to.order = c('dataset', 'pca.batch')
# )
# plot(pad.se$library.size, col = factor(pad.se$pca.batch))
# plot(pad.se$library.size, col = factor(pad.se$dataset))
#
#
#
#
#
#
# pad.se <- identifyUnknownUV(
#     se.obj = pad.se,
#     assay.name = 'RawCount',
#     approach = 'rle',
#     chronological.detection = TRUE,
#     add.to.sample.annotation = TRUE,
#     col.name = 'Estimated.batches'
# )
# pad.se$Tumour.purity <- pad.se$tumour.purity.estimate
#
# pdf(
#     paste0(
#         pccat.fig,
#         'RUVprps_PCCAT_EstimatedUnwantedVar_TwoStepsPcaAndRLE.pdf'
#     ),
#     width = 10,
#     height = 4
# )
# pad.se@metadata$UnKnownUV$RawCount$`23batches|rle.median|AllGenes|nbClust.kmeansClustering`$plots$batch.plot
# dev.off()
#
#
#
#
#
#
#
#
# pad.se$Biological.groups <- paste0(
#     pad.se$Molecular.Type,
#     '_' ,
#     pad.se$Body.site
# )
# pad.se <- RUVprps:::renameVariables(
#     se.obj = pad.se,
#     current.names = c('Molecular.Type', 'library.size'),
#     new.names = c('Molecular.type', 'Library.size')
# )
# ncg.options <- c(
#     'AnovaCorr.AcrossAllSamples',
#     'AnovaCorr.PerBatchPerBiology',
#     'LinearMixedModel',
#     'TwoWayAnova'
# )
# pad.se@metadata$NCG <- NULL
# t1 <- Sys.time()
# for(i in 1:4){
#     # pdf(
#     #     paste0(pccat.fig, "RUVprps_PCCAT_NcgIdentification_Approach_", ncg.options[i], ".pdf"),
#     #     height = 5,
#     #     width = 6)
#     pad.se <- findNcgSupervised(
#         se.obj = pad.se,
#         assay.name = 'RawCount',
#         bio.variables = c('Tumour.purity', 'Biological.groups'),
#         uv.variables = c('Estimated.batches', 'Library.size', 'Library.type'),
#         samples.to.use = 'all',
#         approach = ncg.options[i],
#         form = ~ Tumour.purity + (1|Biological.groups) + (1|Estimated.batches) + Library.size + (1|Library.type),
#         use.rank = FALSE,
#         bio.percentile = NULL,
#         uv.percentile = 0.8,
#         pseudo.count = 0.5,
#         adjust.data = FALSE,
#         adjustment.variables = "uv",
#         ncg.selection.method = 'quantile',
#         nb.ncg = 0.05,
#         ncg.group.name = 'pa.scenario1',
#         ncg.set.name = ncg.options[i],
#         nb.cores = 15,
#         check.se.obj = FALSE,
#         plot.ncg.assessment = TRUE,
#         assess.ncg = TRUE
#     )
#     # dev.off()
# }
# t2 <- Sys.time()
# t2 - t1
#
#
#
#
#
# pad.se <- createPrPsSupervised(
#     se.obj = pad.se,
#     assay.name = 'RawCount',
#     bio.variables = c('Tumour.purity', 'Biological.groups'),
#     samples.to.use = 'all',
#     uv.variables = c('Estimated.batches', 'Library.type', 'Library.size'),
#     bio.clustering.method = "kmeans",
#     nb.bio.clusters = 3,
#     apply.log = TRUE,
#     nb.other.uv.clusters = 4,
#     check.prps.connectedness = FALSE,
#     apply.other.uv.variables = FALSE,
#     prps.group.name = 'pa.scenario1',
#     verbose = TRUE
# )
#
#
#
#
# t1 <- Sys.time()
# pad.se <- RUVIIIprps(
#     se.obj = pad.se,
#     assay.name = 'RawCount',
#     control.sample.types = 'prps',
#     prps.type = 'supervised',
#     prps.group  = 'pa.scenario1',
#     prps.set.names = 'all',
#     ncg.type = 'supervised',
#     ncg.group.names = 'pa.scenario1',
#     ncg.set.names = 'TwoWayAnova',
#     k = c(1:10, 15, 20, 30, 40),
#     apply.log = TRUE,
#     eta = NULL,
#     data.to.log = 'assay',
#     pseudo.count = 0.5,
#     return.wa = TRUE,
#     check.se.obj = FALSE
# )
# t2 <- Sys.time()
# t2 - t1
#
#
#
#
# pad.se <- qs::qread('pad.se.tes.qs')
# pad.se$dataset <- as.factor(pad.se$dataset)
# ppcorr.gene.sets <- selectGenesForPPcorr(
#     se.obj = pad.se,
#     assay.names = c('RawCount', 'RawCount'),
#     variables = c('Library.size', 'Tumour.purity'),
#     cor.cutoff = 0.98,
#     groups = 'dataset',
#     apply.log = TRUE
# )
# lapply(ppcorr.gene.sets, length)
# tumour.purity.genes.sets <- RUVprps:::selectGenesSets(
#     se.obj = pad.se,
#     gene.set = "immune.stromal"
# )
#
# t1 <- Sys.time()
# pad.se <- assessVariation(
#     se.obj = pad.se,
#     assay.names = 'all',
#     assessment.level = 'L1',
#     bio.variables =  c('Tumour.purity', 'Biological.groups'),
#     uv.variables =  c('Library.size','Library.type', 'Estimated.batches'),
#     pcorr.genes = ppcorr.gene.sets,
#     gene.set.score.list = tumour.purity.genes.sets,
#     general.points.size = 3,
#     apply.log = FALSE,
#     override.check = FALSE,
#     check.se.obj = FALSE
# )
# t2 <- Sys.time()
# t2-t1 # 1.16
#
# pad.se <- assessNormalization(
#     se.obj = pad.se,
#     assay.names = 'all',
#     assessment.level = 'L1',
#     select.top.ruv = FALSE,
#     bio.variables = c("Tumour.purity", "Biological.groups"),
#     uv.variables = c("Library.size", "Library.type", "Estimated.batches"),
#     bio.weight = 0.6,
#     uv.weight = 0.4,
#     corr.cutoff = list(Tumour.purity = .5, Library.size = .3),
#     output.name = 'gtex.assessNormalization'
# )
