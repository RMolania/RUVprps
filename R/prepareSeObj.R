#' Prepares a SummarizedExperiment object from tabular data.
#'
#' @author Ramyar Molania
#'
#' @description
#' This function creates a SummarizedExperiment object from tabular expression data set(s) and sample annotation (if
#' available). The function can identify and remove lowly expressed genes if a raw count data is provided, add a range
#' of annotations for individual genes e.g. biotype, chromosome names, ... . , estimate tumor purity, and provide several
#' sets of housekeeping genes and a immune and stromal gene signatures. The housekeeping gene sets could be suitable
#' negative control genes for identification of sources of unknown unwanted variation and the implementation of the RUV
#' methods. The immune and stromal gene signature could be used to estimate tumor purity variation in cancer RNA-seq data.
#'
#' @details
#' The SummarizedExperiment object is a data structure used in the R for representing and manipulating high-dimensional
#' experimental data. Here are some key features and components of the SummarizedExperiment object:
#' Assays:
#' The SummarizedExperiment allows for the incorporation of multiple assays (data). Each assay is a separate matrix of
#' data associated with the same features and samples. Rows typically represent features (e.g., genes, transcripts), and
#' columns represent samples or experimental conditions.
#' Row and Column Metadata:
#' The SummarizedExperiment object includes metadata associated with both rows and columns. Row metadata can contain
#' information about the features, such as gene annotations or genomic coordinates. Column metadata may include sample
#' information, experimental conditions, or other relevant details.
#' Metadata
#' The Metadata of the SummarizedExperiment allows for flexibility in terms of data types and structures. This makes it
#' suitable for saving metrics and plots in the RUVIIIPRPS R package. We refer to the SummarizedExperiment R package for
#' more details.
#'
#' @param data List or a SummarizedExperiment object. Either a list containing data sets (assays) or expression data. For
#' individual datasets, genes should be arranged in rows and samples in columns. If multiple datasets are provided, ensure
#' that the row names of the assays are in the same order or a SummarizedExperiment object.
#' @param sample.annotation Data frame. Contains sample-specific information for the assays. The order of rows should
#' match the column names of the assay(s).
#' @param raw.count.assay.name Character. The name of the raw RNA-seq count data within the list of assays or expression
#' data. Raw counts data is required if `remove.lowly.expressed.genes` or `calculate.library.size` is set to `TRUE`.
#' @param remove.lowly.expressed.genes Logical. If `TRUE`, lowly expressed genes will be identified and removed from
#' the assay specified in the `raw.count.assay.name` argument. The default is set to `FALSE`.
#' @param count.cutoff Numeric. Minimum count required for at least some sample groups. If the `biological.group` argument
#' is set to `NULL`, all samples will be considered as a single group. Otherwise, the smallest subgroups of the
#' `biological.group` will be considered for filtering lowly expressed genes. This follows the `filterByExpr` function
#' from the edgeR package.
#' @param biological.group Character. The column name in the sample annotation specifying the biological groups.
#' The smallest biological groups will be considered for filtering lowly expressed genes. If `NULL`, all samples will
#' be considered as a single group.
#' @param minimum.proportion Numeric. The minimum proportion of samples within a group in which a gene must be expressed.
#' The default is set to 0.5.
#' @param calculate.library.size Logical. If `TRUE`, the library size (total counts) will be calculated using the data
#' provided in `raw.count.assay.name`. Library size will be computed after removing lowly expressed genes.
#' @param estimate.tumor.purity Character. A method to estimate tumor purity. Options include: `estimate`, `singscore`,
#' `both`, and `NULL`. If `estimate`, the ESTIMATE method will be used. If `singscore`, the `singscore` method will be
#' applied. If `both`, both methods will be applied. The default is set to `NULL`.
#' @param gene.ids.to.estimate.tumor.purity TTT
#' @param scale.singscore.values Logical. Indicating to change the range of the tumor purity estimated by the singscore
#' method based on the  the one estimated by the ESTIMATE method. The default is set to `TRUE`.
#' @param assay.name.to.estimate.purity Character. The name of an assay within the list of assays or expression data
#' to be used for estimating tumor purity.
#' @param create.sample.annotation Logical. If `TRUE`, a sample annotation will be generated, initially containing
#' the column names of the assay(s). Additional variables, such as library size, will be included in the annotation.
#' @param gene.annotation Data frame. Contains gene-specific details (e.g., chromosome names, GC content). The list of
#' housekeeping genes and immune/stromal genes will be included in the gene annotation.
#' @param create.gene.annotation Logical. If `TRUE`, a gene annotation will be generated, initially containing row names
#' of the assay(s). Additional gene details and lists will be added to the annotation.
#' @param add.gene.details Logical. If `TRUE`, preset or provided gene details in the `gene.details` argument will be
#' added to the gene annotation.
#' @param gene.group Character. The name of a gene class from the row names of the assay(s). This must be one of the
#' following: `entrezgene_id`, `hgnc_symbol`, or `ensembl_gene_id`.
#' @param column.name TTTT
#' @param gene.details Character or character vector. Indicates the gene details to be included in the gene annotation.
#' @param add.housekeeping.genes Logical. If `TRUE`, publicly available ``housekeeping`` gene sets will be added to the gene
#' annotation. These genes may be used as negative controls for RUV normalization.
#' @param add.immun.stroma.genes Logical. If `TRUE`, immune and stromal gene signatures from Kosuke Yoshihara et al. will
#' be added to the gene annotation. These signatures can be used to estimate tumor purity in cancer RNA-seq data.
#' @param metaData Any. Metadata in any format and dimensions.
#' @param verbose Logical. If `TRUE`, the function will print messages at various steps of the process.
#'
#' @return A SummarizedExperiment object containing assay(s) and also samples annotation, gene annotation and metadata
#' if there are specified.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assay rowData colData rowData<-
#' @importFrom tidyestimate filter_common_genes estimate_score
#' @importFrom biomaRt getBM useMart useDataset
#' @importFrom singscore rankGenes simpleScore
#' @importFrom Matrix colSums rowSums
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr left_join
#' @importFrom knitr kable
#' @importFrom edgeR cpm
#' @export

# raw.count <- data.table::fread(
#     '/Users/molania.r/Documents/Current_Projects/Project_DFCI_2_RUVprpsApplication/Application_RUVprps/data/Tumor-25.01-Polya_ensembl_counts_60499genes_2025-03-03.tsv')
# raw.count <- as.data.frame(raw.count)
# gene.names <- raw.count$Gene
# raw.count <- raw.count[ , -1]
# raw.count <- data.matrix(raw.count)
# row.names(raw.count) <- gene.names
#
#
# sample.annot <- read.delim('/Users/molania.r/Documents/Current_Projects/Project_DFCI_2_RUVprpsApplication/Application_RUVprps/data/GSE294351_clinical_Treehouse-Tumor-Compendium-25.01-PolyA_20250131v1.tsv')
# row.names(sample.annot) <- sample.annot$th_dataset_id
# dim(sample.annot)
# common.samples <- intersect(sample.annot$th_dataset_id, colnames(raw.count))
# sample.annot <- sample.annot[common.samples , ]
# raw.count <- raw.count[ , common.samples ]
# all.equal(colnames(raw.count), sample.annot$th_dataset_id)
#
# row.names(raw.count) <- sub("\\..*$", "", row.names(raw.count))
#
# data = raw.count
# sample.annotation = sample.annot
# raw.count.assay.name = 'RawCount'
# remove.lowly.expressed.genes = FALSE
# count.cutoff = 10
# biological.group = NULL
# minimum.proportion = 0.5
# calculate.library.size = TRUE
# estimate.tumor.purity = 'both'
# scale.singscore.values = TRUE
# assay.name.to.estimate.purity = 'RawCount'
# create.sample.annotation = FALSE
# gene.annotation = NULL
# create.gene.annotation = TRUE
# add.gene.details = TRUE
# gene.group = 'ensembl_gene_id'
# column.name = NULL
# gene.details = NULL
# add.housekeeping.genes = TRUE
# add.immun.stroma.genes = TRUE
# metaData = NULL
# verbose = TRUE
# gene.ids.to.estimate.tumor.purity = 'hgnc_symbol'
# gene.group = 'ensembl_gene_id'
#
# se.obj <- prepareSeObj(
#     data = raw.count,
#     sample.annotation = sample.annot,
#     raw.count.assay.name = 'RawCount',
#     remove.lowly.expressed.genes = FALSE,
#     calculate.library.size = TRUE,
#     estimate.tumor.purity = 'both',
#     gene.ids.to.estimate.tumor.purity = 'ensembl_gene_id',
#     scale.singscore.values = TRUE,
#     gene.annotation = NULL,
#     create.gene.annotation = TRUE,
#     add.gene.details = TRUE,
#     add.housekeeping.genes = TRUE,
#     assay.name.to.estimate.purity = 'RawCount',
#     gene.group = 'ensembl_gene_id')


prepareSeObj <- function(
        data = NULL,
        sample.annotation = NULL,
        raw.count.assay.name = NULL,
        remove.lowly.expressed.genes = FALSE,
        count.cutoff = 10,
        biological.group = NULL,
        minimum.proportion = 0.5,
        calculate.library.size = FALSE,
        estimate.tumor.purity = NULL,
        gene.ids.to.estimate.tumor.purity = 'hgnc_symbol',
        scale.singscore.values = TRUE,
        assay.name.to.estimate.purity = NULL,
        create.sample.annotation = FALSE,
        gene.annotation = NULL,
        create.gene.annotation = FALSE,
        add.gene.details = FALSE,
        gene.group = NULL,
        column.name = NULL,
        gene.details = NULL,
        add.housekeeping.genes = FALSE,
        add.immun.stroma.genes = FALSE,
        metaData = NULL,
        verbose = TRUE
        ){
    printColoredMessage(message = '------------The prepareSeObj function starts:',
                        color = 'white',
                        verbose = verbose
                        )
    if (inherits(data, what = 'list')){
        # data input is list ####
        ## Checking the function inputs ####
        ### dimensions and orders of the assays ####
        if (is.null(names(data))){
            stop('The data must have a name for individual element.')
        }
        if (length(names(data)) != length(unique(names(data))) ){
            stop('The name(s) of gene expression matrix in the "data" object should be unique.' )
        }
        if (length(data) > 1){
            dim.data <- unlist(lapply(
                names(data),
                function(x) c(nrow(data[[x]]), ncol(data[[x]])))
            )
            if (length(unique(dim.data)) !=2 ){
                stop('All the expression matrix in the "data" object must have the same dimensions.')
            }
            all.assays <- combn(x = 1:length(data), m = 2)
            m.out <- lapply(
                1:ncol(all.assays),
                function(x){
                    if (!all.equal(row.names(data[[all.assays[ 1, x]]]), row.names(data[[all.assays[ 2, x]]])))
                        stop('The row names of the "data" must be in the same order.')
                })
            m.out <- lapply(
                1:ncol(all.assays),
                function(x){
                    if (!all.equal(colnames(data[[all.assays[ 1, x]]]),colnames(data[[all.assays[ 2, x]]])))
                        stop('The column names of the "data" must be in the same order.')
                })
        }
        ### raw count assay name ####
        if (!is.null(raw.count.assay.name)){
            if (is.logical(raw.count.assay.name)){
                stop('The "raw.count.assay.name" must be a name of gene expression matrix in the "data" object.')
            }
            if (length(raw.count.assay.name) > 1){
                stop('The "raw.count.assay.name" must contain a single assay name.')
            }
            if (!raw.count.assay.name %in% names(data)){
                stop('The "raw.count.assay.name" must be a name of gene expression matrix in the "data" object.')
            }
        }

        ### library size calculation ####
        if (!is.logical(calculate.library.size)){
            stop('The "calculate.library.size" must be logical: "TRUE" or "FALSE".')
        }
        if (isTRUE(calculate.library.size)){
            if (is.null(raw.count.assay.name)){
                stop('To calculate library size, "raw.count.assay.name" must be provided.')
            }
            if (is.null(sample.annotation) & isFALSE(create.sample.annotation)){
                stop('To calculate library size, either a "sample.annotation" must be provided or the "create.sample.annotation" must be set to "TRUE".')
            }
            if (!is.null(sample.annotation)){
                if (nrow(sample.annotation) != ncol(data[[raw.count.assay.name]]) ){
                    stop('The number of rows in the "sample.annotation" must be the same as the number fo the columns in the data.')
                }
                if (!all.equal(colnames(data[[1]]), row.names(sample.annotation))){
                    stop('The colum names of the dataset(s) must be the same as row names of the sample.annotation.')
                }
            }
            if (isTRUE(remove.lowly.expressed.genes)){
                printColoredMessage(
                    message = '- Note that, the library size will be calculated after removing lowly expressed genes.',
                    color = 'blue',
                    verbose = verbose
                )
            }
        }
        ### tumor purity estimation ####
        if (is.logical(estimate.tumor.purity)){
            stop('The "estimate.tumor.purity" must be one of the "estimate", "singscore", "both" or "NULL"')
        }
        if (!is.null(estimate.tumor.purity)){
            if (is.null(sample.annotation) & isFALSE(create.sample.annotation)){
                stop('To add the calculated tumour purity, either a "sample.annotation" should be provided or set the create.sample.annotation=TRUE.')
            } else if (is.null(assay.name.to.estimate.purity)){
                stop('To estimate the tumour purity, the "assay.name.to.estimate.purity" must be provided.')
            } else if (!estimate.tumor.purity %in% c("estimate", "singscore", "both")){
                stop('The "estimate.tumor.purity" must be one of the "estimate", "singscore", "both" or "NULL"')
            } else if (is.null(gene.ids.to.estimate.tumor.purity)){
                stop('To estimate purity, the "gene.ids.to.estimate.tumor.purity" must be specified to one of the "entrezgene_id", "hgnc_symbol" or "ensembl_gene_id".')
            }
        }
        if (isTRUE(add.gene.details)){
            if (is.null(gene.annotation) & isFALSE(create.gene.annotation)){
                stop('To add gene details, either a "gene.annotation" should be provided or set the "create.gene.annotation=TRUE".')
            }
        }
        if (isTRUE(add.housekeeping.genes)){
            if (is.null(gene.annotation) & isFALSE(create.gene.annotation)){
                stop('To add housekeeping genes, either a "gene.annotation" should be provided or set the "create.gene.annotation=TRUE".')
            }
        }
        if (isTRUE(add.immun.stroma.genes)){
            if (is.null(gene.annotation) & isFALSE(create.gene.annotation)){
                stop('To add immun stroma genes, either a "gene.annotation" should be provided or set the "create.gene.annotation=TRUE".')
            }
        }

        ## remove lowly expressed genes ####
        if (isTRUE(remove.lowly.expressed.genes)){
            if (is.null(raw.count.assay.name)){
                stop('To find and remove lowly expressed genes, the "raw.count.assay.name" must be provided.')
            }
            if (!is.null(biological.group)){
                if (is.null(sample.annotation)){
                    stop('To find "biological.group", please provide a sample annotation that contains the "biological.group" variable.')
                }
                if (!biological.group %in% colnames(sample.annotation)){
                    stop('The "biological.group" cannot be found in the sample annotation.')
                }
            }
            if (minimum.proportion > 1 | minimum.proportion < 0){
                stop('The "minimum.proportion" must between 0 to 1.')
            }
            if (is.null(count.cutoff)){
                stop('The count.cutoff cannot be empty.')
            }
            if (count.cutoff < 0){
                stop('The value of the "count.cutoff" cannot be negative.')
            }
        }
        ## sample annotation ####
        if (!is.null(sample.annotation)){
            if (nrow(sample.annotation) != ncol(data[[1]])){
                stop('The number of rows in "sample.annotation" and the number of columns in the "data" must be the same.')
            } else if (!all.equal(row.names(sample.annotation), colnames(data[[1]])) ){
                stop('The order and lables of the row names of "sample.annotation" and colummn names of the "data" should be the same.')
            }
            if (create.sample.annotation){
                stop('A "sample.annotation" is provided, then "create.sample.annotation" must be set to "FALSE".')
            }
        }
        ## gene annotation ####
        if (!is.null(gene.annotation) & isTRUE(create.gene.annotation)){
            stop('A gene annotation is provided, then "create.gene.annotation" must be "FALSE".')
        }
        if (isTRUE(add.gene.details) & is.null(gene.annotation) & !create.gene.annotation){
            stop('To add "add.gene.details" , the "create.gene.annotation" must be set to "TRUE".')
        }
        if (isTRUE(add.housekeeping.genes) & is.null(gene.annotation) & !create.gene.annotation){
            stop('To add "add.housekeeping.genes", gene annotation must be provided or "create.gene.annotation" must be set to "TRUE" .')
        }
        if (isTRUE(add.immun.stroma.genes) & is.null(gene.annotation) & !create.gene.annotation){
            stop('To add "add.immun.stroma.genes", gene annotation must be provided or "create.gene.annotation" must set to "TRUE".')
        }
        if (isTRUE(add.gene.details) & is.null(gene.group)){
            stop('To add gene details, the "gene.group" must be specified ("entrezgene_id", "hgnc_symbol" and "ensembl_gene_id").')
        }
        if (isTRUE(add.gene.details) &!is.null(gene.group)){
            if (!gene.group %in% c('entrezgene_id', 'hgnc_symbol', 'ensembl_gene_id')){
                stop('The "gene.group" must be one of the "entrezgene_id", "hgnc_symbol" or "ensembl_gene_id".')
            }
        }
        if (!is.null(gene.annotation)) {
            if (nrow(gene.annotation) != nrow(data[[1]])) {
                stop('The number of rows in "gene.annotation" and the number of rows in the "data" must be the same.')
            } else if (!sum(row.names(gene.annotation) == row.names(data[[1]])) == nrow(gene.annotation))  {
                stop('The row names in "gene.annotation" and the row names datastes should be identical.')
            } else if (!gene.group %in% colnames(gene.annotation)) {
                stop('The "gene.group" must be a column name in the "gene.annotation".')
            }
        }
        if (isTRUE(add.housekeeping.genes) & is.null(gene.group)){
            stop('To add housekeeping genes, the "gene.group" must be specified.')
        }
        if (isTRUE(add.housekeeping.genes) & is.null(gene.annotation) & !create.gene.annotation){
            stop('To add housekeeping genes, the "create.gene.annotation" must be "TRUE".')
        }
        # grammar
        if (length(data) == 1){
            assay.n <- 'assay'
        } else assay.n <- 'assays'
        # remove lowly expressed genes ####
        if (isTRUE(remove.lowly.expressed.genes)){
            printColoredMessage(
                message = paste0('-- Removing lowly expressed genes from the ', raw.count.assay.name, ' assay.'),
                color = 'magenta',
                verbose = verbose
                )
            library.size <- Matrix::colSums(data[[raw.count.assay.name]])
            cpm.data <- edgeR::cpm(y = data[[raw.count.assay.name]], lib.size = NULL)
            cpm.cutoff <- round(count.cutoff/median(library.size) * 1e6, digits = 2)
            if (!is.null(biological.group)){
                sample.size <- min(table(sample.annotation[[biological.group]]))
            } else if (!is.null(minimum.proportion)){
                sample.size <- round(ncol(cpm.data) * minimum.proportion, digits = 0)
            }
            if (is.null(minimum.proportion) & is.null(biological.group)){
                sample.size <- ncol(cpm.data)
            }
            keep.genes <- Matrix::rowSums(cpm.data >= cpm.cutoff) >= sample.size
            printColoredMessage(
                message = paste0(
                    sum(keep.genes),
                    ' of ',
                    nrow(cpm.data),
                    ' genes with expression cpm cutoff => ',
                    cpm.cutoff,
                    ' in at least ',
                    sample.size,
                    ' samples are kept as highly expressed genes.'),
                color = 'blue',
                verbose = verbose
            )
            names.assays <-  names(data)
            data <- lapply(
                names.assays,
                function(x) data[[x]][keep.genes ,])
            names(data) <- names.assays
            rm(cpm.data, library.size)
            if (!is.null(gene.annotation)){
                gene.annotation <- gene.annotation[keep.genes , ]
            }
        }
        ### calculate library size
        if (isTRUE(calculate.library.size)){
            printColoredMessage(
                message = ' -- Calculating library size.',
                color = 'magenta',
                verbose = verbose
            )
            printColoredMessage(
                message = 'Note, to calculate library size, a raw count data without any transformation should be provided.',
                color = 'red',
                verbose = verbose
            )
            library.size <- Matrix::colSums(data[[raw.count.assay.name]])
            printColoredMessage(
                message = 'The library size is calculated, with summaries (in millions):',
                color = 'blue',
                verbose = verbose
            )
            if (isTRUE(verbose)) print(summary(library.size/10^6), color = 'blue')
        }
        # sample annotation ####
        printColoredMessage(
            message = '-- Sample annotation:',
            color = 'magenta',
            verbose = verbose
        )
        if (!is.null(sample.annotation)){
            printColoredMessage(
                message = 'The "sample.annotation" will be added to the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            if (calculate.library.size){
                sample.annotation[['library.size']] <- library.size
            }
        }
        if (is.null(sample.annotation) & create.sample.annotation){
            printColoredMessage(
                message = 'A sample.annotation will be created and added to the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose
            )
            if (calculate.library.size){
                sample.annotation <- data.frame(
                    sample.ids = colnames(data[[1]]),
                    library.size = library.size
                )
            } else sample.annotation <- data.frame(sample.ids = colnames(data[[1]]))
        }
        if (is.null(sample.annotation) & !create.sample.annotation){
            printColoredMessage(
                message = 'The SummarizedExperiment object will not contain sample annotation.',
                color = 'blue',
                verbose = verbose)
        }
        # gene annotation ####
        printColoredMessage(
            message = '-- Gene annotation:',
            color = 'magenta',
            verbose = verbose
        )
        if (!is.null(gene.annotation)){
            printColoredMessage(
                message = 'The gene.annotation will be added to the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose)
        } else if (is.null(gene.annotation) & isTRUE(create.gene.annotation)){
            printColoredMessage(
                message = 'A gene.annotation (rowData) will be created and added to the SummarizedExperiment object.',
                color = 'blue',
                verbose = verbose)
            gene.annotation <- data.frame(gene.ids = row.names(data[[1]]))
            colnames(gene.annotation) <- gene.group
        } else if (is.null(gene.annotation)){
            printColoredMessage(
                message = 'The SummarizedExperiment object will not contain gene annotation.',
                color = 'blue',
                verbose = verbose)
        }
        # add gene details ####
        printColoredMessage(
            message = '-- Gene details:',
            color = 'magenta',
            verbose = verbose
        )
        if (isTRUE(add.gene.details)){
            if (is.null(gene.details)){
                printColoredMessage(
                    message = 'The gene.details is not specified, some pre-set details will be added to the gene annotation.',
                    color = 'blue',
                    verbose = verbose
                    )
                printColoredMessage(
                    message = 'Obtain the pre-set gene details from the bioMart R package, this may take a few minutes.',
                    color = 'blue',
                    verbose = verbose
                    )
                ensembl <- useMart('ensembl') # 24 March 2022
                ensembl <- useDataset(
                    mart = ensembl,
                    'hsapiens_gene_ensembl'
                    )
                bioMart.geneAnnot <- biomaRt::getBM(
                    attributes = c(
                        'entrezgene_id',
                        'hgnc_symbol',
                        'gene_biotype',
                        'ensembl_gene_id',
                        'description',
                        'chromosome_name'),
                    mart = ensembl
                    )
                bioMart.geneAnnot <- bioMart.geneAnnot[!duplicated(bioMart.geneAnnot[[gene.group]]), ]
                gene.annotation <- dplyr::left_join(
                    x = gene.annotation,
                    y = bioMart.geneAnnot,
                    by = gene.group,
                    multiple = 'first')
            } else if (!is.null(gene.details)){
                printColoredMessage(
                    message = 'Obtain gene details from the bioMart R package.',
                    color = 'blue',
                    verbose = verbose)
                ensembl <- useMart('ensembl') # 24 March 2022
                ensembl <- useDataset(
                    mart = ensembl,
                    'hsapiens_gene_ensembl'
                )
                attributes.list <- biomaRt::listAttributes(mart = ensembl)
                if (sum(gene.details %in% attributes.list$name) == 0){
                    stop('Non of the provided gene.details are found in the attributes list (biomaRt::listAttributes) in the biomaRt.')
                } else {
                    printColoredMessage(
                        message = paste0(
                            sum(gene.details %in% attributes.list$name),
                            ' of ',
                            length(gene.details),
                            'gene.details are found.'),
                        color = 'blue',
                        verbose = verbose)
                }
                gene.details <- unique(gene.details, gene.group)
                bioMart.geneAnnot <- biomaRt::getBM(
                    attributes = gene.details,
                    mart = ensembl
                )
                bioMart.geneAnnot <- bioMart.geneAnnot[!duplicated(bioMart.geneAnnot[[gene.group]]), ]
                gene.annotation <- dplyr::left_join(
                    x = gene.annotation,
                    y = bioMart.geneAnnot,
                    multiple = 'first',
                    by = gene.group)
            }
        } else{
            printColoredMessage(
                message = 'Any extra gene details are not specified.',
                color = 'blue',
                verbose = verbose
            )
        }
        # add housekeeping genes list ####
        if (isTRUE(add.housekeeping.genes)){
            printColoredMessage(
                message = '-- Add several lists of housekeeping genes to the gene annotation:',
                color = 'magenta',
                verbose = verbose
            )
            hk.im.genes <- hk_immunStroma
            keep.cols <- c(which(colnames(hk.im.genes) %in% gene.group), 4:9)
            gene.annotation <- dplyr::left_join(
                x = gene.annotation,
                y = hk.im.genes[ , keep.cols],
                by = gene.group,
                multiple = 'first'
            )
            printColoredMessage(
                message = 'Seven different lists of housekeeping genes are added to the gene annotation.',
                color = 'blue',
                verbose = verbose
            )
            for(g in colnames(hk.im.genes)[4:9])
                gene.annotation[[g]][is.na(gene.annotation[[g]])] <- FALSE
            nb.hk.genes <- lapply(
                colnames(hk.im.genes)[4:9],
                function(x) sum(gene.annotation[[x]]))
            names(nb.hk.genes) <- colnames(hk.im.genes)[4:9]
            if (verbose) print(kable(unlist(nb.hk.genes),
                                    caption = 'Number of genes in each list of housekeeping genes:',
                                    col.names = 'nb.genes'))
        }
        # add immune and stroma genes signatures ####
        if (isTRUE(add.immun.stroma.genes)){
            printColoredMessage(
                message = '-- Add immune and stromal genes signature to the gene annotation:',
                color = 'magenta',
                verbose = verbose
            )
            hk.im.genes <- hk_immunStroma
            keep.cols <- c(
                which(colnames(hk.im.genes) %in% gene.group),
                10:ncol(hk.im.genes))
            gene.annotation <- as.data.frame(dplyr::left_join(
                x = gene.annotation,
                y = hk.im.genes[ , keep.cols],
                by = gene.group,
                multiple = 'first'
            ))
            printColoredMessage(
                message = 'The immune and stromal genes signature from Kosuke Yoshihara et.al are added.',
                color = 'blue',
                verbose = verbose
            )
            for(g in colnames(hk.im.genes)[10:11])
                gene.annotation[[g]][is.na(gene.annotation[[g]])] <- FALSE
            nb.genes <- lapply(
                colnames(hk.im.genes)[10:11],
                function(x) sum(gene.annotation[x]))
            names(nb.genes) <- colnames(hk.im.genes)[10:11]
            if (verbose) print(
                kable(unlist(nb.genes),
                      caption = 'Number of genes in the immune and stromal gene signatures:',
                      col.names = 'nb.genes'))
        }
        # estimate tumor purity ####
        if (!is.null(estimate.tumor.purity)){
            printColoredMessage(
                message = '-- Estimating tumour purity:',
                color = 'magenta',
                verbose = verbose
                )
            current.row.names <- NULL
            if (!gene.group %in% c("entrezgene_id", "hgnc_symbol")){
                if (!is.null(gene.annotation)){
                    if (sum(colnames(gene.annotation) %in% c("entrezgene_id", "hgnc_symbol")) > 0){
                        gene.ids <- gene.annotation[[gene.ids.to.estimate.tumor.purity]]
                        dup.ids <- duplicated(gene.ids)
                        gene.ids[dup.ids] <- paste('gene', 1:sum(dup.ids))
                    }
                    current.row.names <- row.names(data[[assay.name.to.estimate.purity]])
                    row.names(data[[assay.name.to.estimate.purity]]) <- gene.ids
                    if (sum(colnames(gene.annotation) %in% c("entrezgene_id", "hgnc_symbol")) == 0){
                        stop('To estimate tumour purity, the "gene.ids.to.estimate.tumor.purity" must be either "entrezgene_id" or "hgnc_symbol".')
                    }
                }
                if (is.null(gene.annotation)){
                    stop('To estimate tumour purity, the "gene.ids.to.estimate.tumor.purity" must be either "entrezgene_id" or "hgnc_symbol".')
                }
            }
            if (estimate.tumor.purity == 'estimate'){
                printColoredMessage(
                    message = '-- Estimating tumour purity using the ESTIMATE method:',
                    color = 'blue',
                    verbose = verbose
                    )
                tumour.purity <- tidyestimate::filter_common_genes(
                    df = data[[assay.name.to.estimate.purity]],
                    id = gene.ids.to.estimate.tumor.purity,
                    tidy = FALSE,
                    tell_missing = verbose,
                    find_alias = TRUE
                    )
                tumour.purity <- tidyestimate::estimate_score(
                    df = tumour.purity,
                    is_affymetrix = TRUE
                    )
                tumour.purity <- tumour.purity$purity
                sample.annotation[['tumour.purity']] <- tumour.purity
                if (!is.null(current.row.names)){
                    row.names(data[[assay.name.to.estimate.purity]]) <- current.row.names
                }
            } else if (estimate.tumor.purity == 'singscore'){
                printColoredMessage(
                    message = '-- Estimating tumour purity using the singscore method:',
                    color = 'blue',
                    verbose = verbose
                    )
                im.str.gene.sig <- hk_immunStroma$immune.gene.signature == 'TRUE' |
                    hk_immunStroma$stromal.gene.signature == 'TRUE'
                if (gene.ids.to.estimate.tumor.purity == "entrezgene_id"){
                    im.str.gene.sig <- hk_immunStroma$entrezgene_id[im.str.gene.sig]
                } else if (gene.ids.to.estimate.tumor.purity == 'hgnc_symbol'){
                    im.str.gene.sig <- hk_immunStroma$hgnc_symbol[im.str.gene.sig]
                } else if (gene.ids.to.estimate.tumor.purity == 'ensembl_gene_id')
                    im.str.gene.sig <- hk_immunStroma$ensembl_gene_id[im.str.gene.sig]
                tumour.purity <- singscore::rankGenes(data[[assay.name.to.estimate.purity]])
                tumour.purity <- singscore::simpleScore(
                    rankData = tumour.purity,
                    upSet = im.str.gene.sig)
                tumour.purity <- tumour.purity$TotalScore
                sample.annotation[['tumour.purity']] <- 1 - tumour.purity
                if (!is.null(current.row.names)){
                    row.names(data[[assay.name.to.estimate.purity]]) <- current.row.names
                }
            } else if (estimate.tumor.purity == 'both'){
                printColoredMessage(
                    message = '- Estimating tumour purity using both ESTIMATE and singscore methods:',
                    color = 'blue',
                    verbose = verbose
                    )
                printColoredMessage(
                    message = '- Applying the Estimate method:',
                    color = 'blue',
                    verbose = verbose
                    )
                tumour.purity <- tidyestimate::filter_common_genes(
                    df = data[[assay.name.to.estimate.purity]],
                    id = gene.ids.to.estimate.tumor.purity,
                    tidy = FALSE,
                    tell_missing = verbose,
                    find_alias = TRUE
                    )
                tumour.purity <- tidyestimate::estimate_score(
                    df = tumour.purity,
                    is_affymetrix = TRUE)
                tumour.purity.estimate <- tumour.purity$purity
                printColoredMessage(
                    message = '- Applying the singscore method:',
                    color = 'blue',
                    verbose = verbose
                    )
                im.str.gene.sig <- hk_immunStroma$immune.gene.signature == 'TRUE' |
                    hk_immunStroma$stromal.gene.signature == 'TRUE'
                if (gene.ids.to.estimate.tumor.purity == "entrezgene_id"){
                    im.str.gene.sig <- hk_immunStroma$entrezgene_id[im.str.gene.sig]
                } else if (gene.ids.to.estimate.tumor.purity == 'hgnc_symbol'){
                    im.str.gene.sig <- hk_immunStroma$hgnc_symbol[im.str.gene.sig]
                } else if (gene.ids.to.estimate.tumor.purity == 'ensembl_gene_id')
                    im.str.gene.sig <- hk_immunStroma$ensembl_gene_id[im.str.gene.sig]
                tumour.purity <- singscore::rankGenes(data[[assay.name.to.estimate.purity]])
                tumour.purity <- singscore::simpleScore(
                    rankData = tumour.purity,
                    upSet = im.str.gene.sig
                    )
                tumour.purity.singscore <- tumour.purity$TotalScore
                sample.annotation[['tumour.purity.estimate']] <- tumour.purity.estimate
                sample.annotation[['tumour.purity.singscore']] <- 1 - tumour.purity.singscore
                if (isTRUE(scale.singscore.values)){
                    tps <- sample.annotation[['tumour.purity.singscore']]
                    tpe <- sample.annotation[['tumour.purity.estimate']]
                    rtps <- range(tps)
                    rtpe <- range(tpe)
                    stps <- (tps - min(rtps)) / (max(rtps) - min(rtps)) * (max(rtpe) - min(rtpe)) + min(rtpe)
                    sample.annotation[['tumour.purity.singscore.scaled']] <- stps
                }
            }
            if (!is.null(current.row.names)){
                row.names(data[[assay.name.to.estimate.purity]]) <- current.row.names
            }
        }
        # outputs ####
        printColoredMessage(
            message = '-- Creating a SummarizedExperiment object:',
            color = 'magenta',
            verbose = verbose
        )
        if (is.null(gene.annotation) & is.null(sample.annotation) & is.null(metaData)) {
            se.obj <- SummarizedExperiment::SummarizedExperiment(assays = data)
        } else if (!is.null(gene.annotation) & is.null(sample.annotation) & is.null(metaData)) {
            se.obj <- SummarizedExperiment::SummarizedExperiment(
                assays = data,
                rowData = gene.annotation)
        } else if (is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
            se.obj <- SummarizedExperiment::SummarizedExperiment(
                assays = data,
                colData = DataFrame(sample.annotation),
                metadata = metaData
            )
        } else if (is.null(gene.annotation) & is.null(sample.annotation) & !is.null(metaData)) {
            se.obj <- SummarizedExperiment::SummarizedExperiment(
                assays = data,
                metadata = metaData)
        } else if (!is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
            se.obj <- SummarizedExperiment::SummarizedExperiment(
                assays = data,
                rowData = gene.annotation,
                colData = sample.annotation)
        } else if (is.null(gene.annotation) & !is.null(sample.annotation) & is.null(metaData)) {
            se.obj <- SummarizedExperiment::SummarizedExperiment(
                assays = data,
                rowData = gene.annotation,
                colData = DataFrame(sample.annotation)
            )
        } else if (!is.null(gene.annotation) & !is.null(sample.annotation) & !is.null(metaData)) {
            se.obj <- SummarizedExperiment::SummarizedExperiment(
                assays = data,
                rowData = gene.annotation,
                colData = sample.annotation,
                metadata = metaData
            )
        }
        printColoredMessage(
            message = paste0('A SummarizedExperiment object is created with:'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '-',
                nrow(se.obj),
                ' measurements (e.g. genes) and ',
                ncol(se.obj),
                ' assays (e.g.samples)'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '-',
                length(assays(se.obj)),
                ' data sets (assays)'),
            color = 'blue',
            verbose = verbose
        )
        if (!is.null(sample.annotation))
            printColoredMessage(
                message = paste0(
                    '-',
                    ncol(SummarizedExperiment::colData(se.obj)),
                    ' annotations for the samples'),
                color = 'blue',
                verbose = verbose
            )
        if (!is.null(gene.annotation))
            printColoredMessage(
                message = paste0(
                    '-',
                    ncol(SummarizedExperiment::rowData(se.obj)),
                    ' annotations for the genes'),
                color = 'blue',
                verbose = verbose
            )
        if (!is.null(metaData))
            printColoredMessage(message = '- a metadata',
                                color = 'blue',
                                verbose = verbose)
        printColoredMessage(message = '------------The createSeObj function finished.',
                            color = 'white',
                            verbose = verbose)
        return(se.obj)

    }
    # Data input is SummarizedExperiment ####
    if (class(data)[1] == 'SummarizedExperiment' | class(data)[1] == 'RangedSummarizedExperiment'){
        if (!is.null(gene.group)){
            gene.annotation.a <- rowData(data)
            gene.annotation.a[[gene.group]] <- gene.annotation.a[[column.name]]
            rowData(data) <- gene.annotation.a
        }
        ## Checking and adding sample annotation ####
        if (is.logical(sample.annotation)){
            stop('The "sample.annotation" cannot be logical.')
        }
        if (!is.null(sample.annotation)){
            if (ncol(SummarizedExperiment::colData(x = data)) == 0){
                if (nrow(sample.annotation) != ncol(data)){
                    stop(paste('The number of rows in the provided "sample.annotation is not equal',
                               'to the number of columns in the the SummarizedExperiment object.'))
                }
                if (!all.equal(row.names(sample.annotation), colnames(data))){
                    stop(paste0(
                        'The row names of the provided sample annotation are not matched with',
                        'the column names of the SummarizedExperiment.'))
                } else if (all.equal(row.names(sample.annotation), colnames(data))){
                    SummarizedExperiment::colData(data) <- DataFrame(sample.annotation)
                    if (isTRUE(calculate.library.size) & !is.null(biological.group)){
                        if (!is.character(biological.group)){
                            stop('The "biological group" must be a name of a column that conatins biological variable.')
                        }
                        if (length(biological.group) > 1){
                            stop('The "biological group" must be a single name of a column that conatins biological variable.')
                        }
                        if (!biological.group %in% colnames(SummarizedExperiment::colData(x = data))){
                            stop('The "biological group" cannot be found in the SummarizedExperiment object.')
                        }
                    }
                }
            } else {
                printColoredMessage(
                    message = paste0(
                        '** The SummarizedExperiment already contains a sample annotation file.',
                        'The provided sample annotation will not be added.'
                        ),
                    color = 'red',
                    verbose = verbose
                )
                if (isTRUE(calculate.library.size) & !is.null(biological.group)){
                    if (!is.character(biological.group)){
                        stop('The "biological group" must be a name of a column that conatins biological variable.')
                    }
                    if (length(biological.group) > 1){
                        stop('The "biological group" must be a single name of a column that conatins biological variable.')
                    }
                    if (!biological.group %in% colnames(SummarizedExperiment::colData(x = data))){
                        stop('The "biological group" cannot be found in the SummarizedExperiment object.')
                    }
                }
                sample.annotation <- as.data.frame(SummarizedExperiment::colData(data))
            }
        }
        if (is.null(sample.annotation)){
            if (isTRUE(calculate.library.size) & !is.null(biological.group)){
                stop('A sample annotation wiht at least a column contating "biological.group" must be provided.')
            }
            if (isTRUE(calculate.library.size) | isTRUE(estimate.tumor.purity)){
                if (ncol(SummarizedExperiment::colData(x = data)) == 0){
                    printColoredMessage(
                        message = paste0(
                            '- The SummarizedExperiment does no contain sample annotation.',
                            'A sample annotation will be created and added to the object.'),
                        color = 'blue',
                        verbose = verbose
                    )
                    SummarizedExperiment::colData(data) <- DataFrame(sample.ids = colnames(data))
                }
            }
        }
        # Checking and adding gene annotation ####
        if (is.logical(gene.annotation)){
            stop('The "gene.annotation" canot be logical.')
        }
        if (!is.null(gene.annotation) & isTRUE(create.gene.annotation)){
            stop('Both the "gene.annotation" and "create.gene.annotation" canot be specified.')
        }
        if (!is.null(gene.annotation)){
            if (ncol(SummarizedExperiment::rowData(x = data)) == 0){
                if (nrow(gene.annotation) != nrow(data)){
                    stop(paste('The number of rows in the provided "gene.annotation is not equal',
                               'to the number of rows in the the SummarizedExperiment object.'))
                }
                if (!all.equal(row.names(gene.annotation), row.names(data))){
                    stop(paste0(
                        'The row names of the provided gene annotation are not matched with',
                        'the row names of the SummarizedExperiment')
                        )
                } else if (all.equal(row.names(gene.annotation), row.names(data))){
                    SummarizedExperiment::rowData(data) <- DataFrame(gene.annotation)
                }
            } else {
                printColoredMessage(
                    message = paste0(
                        'The SummarizedExperiment contains a gene annotation file.',
                        'The provided gene annotation file will not be added.'),
                    color = 'red',
                    verbose = verbose
                )
            }
        }
        if (is.null(gene.annotation) & isFALSE(create.gene.annotation)){
            if (isTRUE(add.gene.details) | isTRUE(add.housekeeping.genes) | isTRUE(add.immun.stroma.genes)){
                if (ncol(SummarizedExperiment::rowData(x = data)) == 0){
                    stop(paste0(
                        'The SummarizedExperiment object dose not contain a gene annotation, then',
                        'The "gene.annotation" or "create.gene.annotation" must be specified if any of',
                        '"add.gene.details" or "add.housekeeping.genes" or "add.immun.stroma.genes" are specified.'))
                }
                if (is.null(gene.group)){
                    stop(paste0(
                        'The "gene.group" must be provided if any of ',
                        '"add.gene.details" or "add.housekeeping.genes" or "add.immun.stroma.genes" are specified.')
                        )
                }
                if (!gene.group %in% c('entrezgene_id', 'hgnc_symbol', 'ensembl_gene_id')){
                    stop('The "gene.group" must be one of the "entrezgene_id", "hgnc_symbol" or "ensembl_gene_id".')
                }
                if (!gene.group %in% colnames(SummarizedExperiment::rowData(x = data))){
                    stop('The "gene.group" cannot be found in the gene annotation file of the SummarizedExperiment object.')
                }
            }
        }
        if (is.null(gene.annotation) & isTRUE(create.gene.annotation)){
            if (ncol(SummarizedExperiment::rowData(data)) != 0){
                stop(paste0('The SummarizedExperiment contains a gene annotation, so',
                            'the "create.gene.annotation" must be set to "FALSE".'))
            }
            if (ncol(SummarizedExperiment::rowData(data)) == 0){
                printColoredMessage(
                    message = paste0(
                        'The SummarizedExperiment does no contain gene annotation.',
                        'A gene annotation will be created and added to the object.'),
                    color = 'blue',
                    verbose = verbose
                )
                gene.annotation <- data.frame(gene.ids = row.names(data))
                if (isTRUE(add.gene.details) | isTRUE(add.housekeeping.genes) | isTRUE(add.immun.stroma.genes)){
                    if (is.null(gene.group)){
                        stop(paste0(
                            'The "gene.group" must be provided if any of the ', '"add.gene.details" or ',
                            '"add.housekeeping.genes" or "add.immun.stroma.genes" are specified.'))
                    }
                    colnames(gene.annotation) <- gene.group
                }
                row.names(gene.annotation) <- row.names(data)
                SummarizedExperiment::rowData(data) <- DataFrame(gene.annotation)
            }
        }

        # Checking estimating of tumor purity ####
        if (is.logical(estimate.tumor.purity)){
            stop('The "estimate.tumor.purity" must be one of the "estimate", "singscore", "both" or "NULL"')
        }
        if (!is.null(estimate.tumor.purity)){
            if (!estimate.tumor.purity %in% c("estimate", "singscore", "both")){
                stop('The "estimate.tumor.purity" must be one of the "estimate", "singscore", "both" or "NULL"')
            }
            if (is.logical(assay.name.to.estimate.purity)){
                stop('The "assay.name.to.estimate.purity" must be a assay name in the SummarizedExperiment object.')
            }
            if (length(assay.name.to.estimate.purity) > 1){
                stop('The "assay.name.to.estimate.purity" must be a single assay name in the SummarizedExperiment object.')
            }
            if (is.null(assay.name.to.estimate.purity)){
                stop('To estimate the tumour purity, the "assay.name.to.estimate.purity" must be provided.')
            }
            if (!assay.name.to.estimate.purity %in% names(assays(data))){
                stop('The "assay.name.to.estimate.purity" cannot be found in the SummarizedExperiment object.')
            }
        }
        # Checking raw count data ####
        if (!is.null(raw.count.assay.name)){
            if (is.logical(raw.count.assay.name)){
                stop('The "raw.count.assay.name" cannot be logical.')
            }
            if (length(raw.count.assay.name) > 1){
                stop('The "raw.count.assay.name" must be a single assay name in the SummarizedExperiment object.')
            }
            if (!raw.count.assay.name %in% names(assays(data))){
                stop('The "raw.count.assay.name" canot be found in the SummarizedExperiment object.')
            }
        }
        # Removing lowly expressed genes ####
        if (isTRUE(remove.lowly.expressed.genes)){
            if (is.null(raw.count.assay.name)){
                stop('To find and remove lowly expressed genes, the "raw.count.assay.name" must be provided.')
            }
            if (!is.null(biological.group)){
                if (!is.character(biological.group)){
                    stop('The "biological.group" must be a name of a column that conatins a biological variable.')
                }
                if (length(biological.group) > 1){
                    stop('The "biological.group" must be single name of a column that conatins a biological variable.')
                }
                if (!biological.group %in% colnames(SummarizedExperiment::colData(data))){
                    stop('The "biological.group" cannot be found in the sample annotation of the SummarizedExperiment object.')
                }
            }
            if (!is.numeric(minimum.proportion)){
                stop('The "minimum.proportion" must be a numeric value between 0 to 1. ')
            }
            if (minimum.proportion > 1 | minimum.proportion < 0){
                stop('The "minimum.proportion" must between 0 to 1.')
            }
            if (is.null(count.cutoff)){
                stop('The count.cutoff cannot be empty.')
            }
            if (!is.numeric(count.cutoff)){
                stop('The "count.cutoff" must be a postive numeric value. ')
            }
            if (count.cutoff < 0){
                stop('The value of the "count.cutoff" cannot be negative.')
            }
        }

        # Adding gene details ####
        if (isTRUE(add.gene.details) & is.null(gene.group)){
            stop('To add gene details, the "gene.group" must be specified as one of the "entrezgene_id", "hgnc_symbol" and "ensembl_gene_id".')
        }
        if (isTRUE(add.gene.details) & !is.null(gene.group)){
            if (!gene.group %in% c('entrezgene_id', 'hgnc_symbol', 'ensembl_gene_id')){
                stop('The "gene.group" must be one of the "entrezgene_id", "hgnc_symbol" or "ensembl_gene_id".')
            }
        }
        if (isTRUE(add.housekeeping.genes) & is.null(gene.group)){
            stop('To add housekeeping genes, the "gene.group" must be specified.')
        }
        # Removing lowly expressed genes ####
        if (isTRUE(remove.lowly.expressed.genes)){
            printColoredMessage(
                message = paste0(
                    '-- Removing lowly expressed genes from the ',
                    raw.count.assay.name,
                    ' data'),
                color = 'magenta',
                verbose = verbose
                )
            if (sum(is.na(assay(x = data , i = raw.count.assay.name))) != 0){
                stop(paste0(
                    'There are NAs in the ',
                    raw.count.assay.name,
                    'data. Lowly expressed genes cannot be found.'))
            }
            library.size <- Matrix::colSums(assay(x = data , i = raw.count.assay.name))
            cpm.data <- edgeR::cpm(
                y = assay(x = data , i = raw.count.assay.name),
                lib.size = NULL
                )
            cpm.cutoff <- round(
                x = count.cutoff/median(library.size) * 1e6,
                digits = 2
                )
            if (!is.null(biological.group)){
                sample.size <- min(table(sample.annotation[[biological.group]]))
            } else if (!is.null(minimum.proportion)){
                sample.size <- round(ncol(cpm.data) * minimum.proportion, digits = 0)
            }
            if (is.null(minimum.proportion) & is.null(biological.group)){
                sample.size <- ncol(cpm.data)
            }
            keep.genes <- Matrix::rowSums(cpm.data >= cpm.cutoff) >= sample.size
            printColoredMessage(
                message = paste0(
                    '- ',
                    sum(keep.genes),
                    ' of ',
                    nrow(cpm.data),
                    ' genes with expression cpm cutoff => ',
                    cpm.cutoff,
                    ' in at least ',
                    sample.size,
                    ' samples are kept as highly expressed genes.'),
                color = 'blue',
                verbose = verbose
            )
            data <- data[keep.genes , ]
            gene.annotation <- as.data.frame(SummarizedExperiment::rowData(data))
        }
        # Calculating library size ####
        if (isTRUE(calculate.library.size)){
            printColoredMessage(
                message = '-- Calculating library size (totall reads).',
                color = 'magenta',
                verbose = verbose
                )
            printColoredMessage(
                message = '**Note, to calculate library size, a raw count data without any transformation should be provided.',
                color = 'red',
                verbose = verbose
                )
            library.size <- Matrix::colSums(assay(x = data, i = raw.count.assay.name))
            printColoredMessage(
                message = '* The library size is calculated, with summaries (in millions):',
                color = 'blue',
                verbose = verbose
            )
            if (isTRUE(verbose)) print(summary(library.size/10^6), color = 'blue')
            if ('library.size' %in% colnames(SummarizedExperiment::colData(data))){
                stop('There is a column named "library.size" in the SummarizedExperiment object.')
            } else data[['library.size']] <- log2(library.size)

        }
        # Adding gene details ####
        if (isTRUE(add.gene.details)){
            printColoredMessage(
                message = '-- Adding gene details:',
                color = 'magenta',
                verbose = verbose
            )
            if (is.null(gene.details)){
                printColoredMessage(
                    message = 'The gene.details is not specified, some pre-set details will be added to the gene annotation.',
                    color = 'blue',
                    verbose = verbose
                    )
                printColoredMessage(
                    message = 'Obtaining the pre-set gene details from the bioMart R package, this may take a few minutes.',
                    color = 'blue',
                    verbose = verbose
                    )
                ensembl <- useMart('ensembl')
                ensembl <- useDataset(
                    mart = ensembl,
                    'hsapiens_gene_ensembl'
                )
                bioMart.geneAnnot <- biomaRt::getBM(
                    attributes = c(
                        'entrezgene_id',
                        'hgnc_symbol',
                        'gene_biotype',
                        'ensembl_gene_id',
                        'description',
                        'chromosome_name'),
                    mart = ensembl
                )
                bioMart.geneAnnot <- bioMart.geneAnnot[!duplicated(bioMart.geneAnnot[[gene.group]]), ]
                gene.annotation <- as.data.frame(SummarizedExperiment::rowData(data))
                gene.annotation$gene.order <- c(1:nrow(gene.annotation))
                gene.annotation <- dplyr::left_join(
                    x = gene.annotation,
                    y = bioMart.geneAnnot,
                    by = gene.group,
                    multiple = 'first'
                    )
                if (isTRUE(all.equal(gene.annotation$gene.order, 1:nrow(gene.annotation)))) {
                    SummarizedExperiment::rowData(data) <- as.data.frame(gene.annotation)
                } else {
                    gene.annotation <- gene.annotation[order(gene.annotation$gene.order) , ]
                    SummarizedExperiment::rowData(data) <- as.data.frame(gene.annotation)
                }

            } else if (!is.null(gene.details)){
                printColoredMessage(
                    message = 'Obtaining the specified gene details from the bioMart R package.',
                    color = 'blue',
                    verbose = verbose
                    )
                ensembl <- useMart('ensembl')
                ensembl <- useDataset(
                    mart = ensembl,
                    'hsapiens_gene_ensembl'
                    )
                attributes.list <- biomaRt::listAttributes(mart = ensembl)
                if (sum(gene.details %in% attributes.list$name) == 0){
                    stop('Non of the provided "gene.details" are found in the attributes list (biomaRt::listAttributes) in the biomaRt.')
                } else {
                    printColoredMessage(
                        message = paste0(
                            sum(gene.details %in% attributes.list$name),
                            ' of ',
                            length(gene.details),
                            'gene.details are found.'),
                        color = 'blue',
                        verbose = verbose)
                }
                gene.details <- unique(gene.details, gene.group)
                bioMart.geneAnnot <- biomaRt::getBM(
                    attributes = gene.details,
                    mart = ensembl
                    )
                bioMart.geneAnnot <- bioMart.geneAnnot[!duplicated(bioMart.geneAnnot[[gene.group]]), ]
                gene.annotation <- as.data.frame(SummarizedExperiment::rowData(data))
                gene.annotation$gene.order <- c(1:nrow(gene.annotation))
                gene.annotation <- dplyr::left_join(
                    x = gene.annotation,
                    y = bioMart.geneAnnot,
                    multiple = 'first',
                    by = gene.group
                    )
                if (isTRUE(all.equal(gene.annotation$gene.order, 1:nrow(gene.annotation)))) {
                    SummarizedExperiment::rowData(data) <- as.data.frame(gene.annotation)
                } else {
                    gene.annotation <- gene.annotation[order(gene.annotation$gene.order) , ]
                    SummarizedExperiment::rowData(data) <- as.data.frame(gene.annotation)
                }
            }
        }
        # Adding housekeeping genes list ####
        if (isTRUE(add.housekeeping.genes)){
            printColoredMessage(
                message = '-- Adding several lists of publicly avaiable housekeeping genes to the gene annotation:',
                color = 'magenta',
                verbose = verbose
            )
            hk.im.genes <- hk_immunStroma
            keep.cols <- c(which(colnames(hk.im.genes) %in% gene.group), 4:9)
            gene.annotation <- as.data.frame(SummarizedExperiment::rowData(data))
            gene.annotation$gene.order <- c(1:nrow(gene.annotation))
            gene.annotation <- dplyr::left_join(
                x = gene.annotation,
                y = hk.im.genes[ , keep.cols],
                by = gene.group,
                multiple = 'first'
            )
            if (isTRUE(all.equal(gene.annotation$gene.order, 1:nrow(gene.annotation)))) {
                SummarizedExperiment::rowData(data) <- as.data.frame(gene.annotation)
            } else {
                gene.annotation <- gene.annotation[order(gene.annotation$gene.order) , ]
                SummarizedExperiment::rowData(data) <- as.data.frame(gene.annotation)
            }
            printColoredMessage(
                message = '* Seven different lists of housekeeping genes are added to the gene annotation.',
                color = 'blue',
                verbose = verbose
            )
            for(g in colnames(hk.im.genes)[4:9])
                gene.annotation[[g]][is.na(gene.annotation[[g]])] <- FALSE
            nb.hk.genes <- lapply(
                colnames(hk.im.genes)[4:9],
                function(x) sum(gene.annotation[[x]])
                )
            names(nb.hk.genes) <- colnames(hk.im.genes)[4:9]
            if (isTRUE(verbose)) print(
                kable(unlist(nb.hk.genes),
                      caption = 'Number of genes in each list of housekeeping genes:',
                      ol.names = 'nb.genes'))
        }
        # Adding immune and stroma genes signatures ####
        if (isTRUE(add.immun.stroma.genes)){
            printColoredMessage(
                message = '-- Adding immune and stromal genes signature to the gene annotation:',
                color = 'magenta',
                verbose = verbose
            )
            hk.im.genes <- hk_immunStroma
            keep.cols <- c(
                which(colnames(hk.im.genes) %in% gene.group),
                10:ncol(hk.im.genes))
            gene.annotation <- as.data.frame(SummarizedExperiment::rowData(data))
            gene.annotation$gene.order <- c(1:nrow(gene.annotation))
            gene.annotation <- as.data.frame(dplyr::left_join(
                x = gene.annotation,
                y = hk.im.genes[ , keep.cols],
                by = gene.group,
                multiple = 'first'
            ))
            if (isTRUE(all.equal(gene.annotation$gene.order, 1:nrow(gene.annotation)))) {
                SummarizedExperiment::rowData(data) <- as.data.frame(gene.annotation)
            } else {
                gene.annotation <- gene.annotation[order(gene.annotation$gene.order) , ]
                SummarizedExperiment::rowData(data) <- as.data.frame(gene.annotation)
            }
            printColoredMessage(
                message = '* The immune and stromal genes signature from Kosuke Yoshihara et.al are added.',
                color = 'blue',
                verbose = verbose
            )
            for (g in colnames(hk.im.genes)[10:11])
                gene.annotation[[g]][is.na(gene.annotation[[g]])] <- FALSE
            nb.genes <- lapply(
                colnames(hk.im.genes)[10:11],
                function(x) sum(gene.annotation[x])
                )
            names(nb.genes) <- colnames(hk.im.genes)[10:11]
            if (isTRUE(verbose)) print(
                kable(unlist(nb.genes),
                      caption = 'Number of genes in the immune and stromal gene signatures:',
                      col.names = 'nb.genes')
                )
        }
        # Estimating tumor purity ####
        if (!is.null(estimate.tumor.purity)){
            printColoredMessage(
                message = '-- Estimating tumour purity (this is just for cancer RNAseq data):',
                color = 'magenta',
                verbose = verbose
                )
            if (estimate.tumor.purity == 'estimate'){
                printColoredMessage(
                    message = '- Estimating tumour purity using the ESTIMATE method:',
                    color = 'blue',
                    verbose = verbose
                    )
                tumour.purity <- tidyestimate::filter_common_genes(
                    df = assay(x = data, i = assay.name.to.estimate.purity),
                    id = gene.ids.to.estimate.tumor.purity,
                    tidy = FALSE,
                    tell_missing = verbose,
                    find_alias = TRUE
                    )
                tumour.purity <- tidyestimate::estimate_score(
                    df = tumour.purity,
                    is_affymetrix = TRUE
                    )
                tumour.purity <- tumour.purity$purity
                sample.annotation[['tumour.purity.estimate']] <- tumour.purity
            }
            if (estimate.tumor.purity == 'singscore'){
                printColoredMessage(
                    message = '* Estimating tumour purity using the singscore method:',
                    color = 'blue',
                    verbose = verbose
                    )
                im.str.gene.sig <- hk_immunStroma$immune.gene.signature == 'TRUE' |
                    hk_immunStroma$stromal.gene.signature == 'TRUE'
                if (gene.ids.to.estimate.tumor.purity == "entrezgene_id"){
                    im.str.gene.sig <- hk_immunStroma$entrezgene_id[im.str.gene.sig]
                } else if (gene.ids.to.estimate.tumor.purity == 'hgnc_symbol'){
                    im.str.gene.sig <- hk_immunStroma$hgnc_symbol[im.str.gene.sig]
                } else if (gene.ids.to.estimate.tumor.purity == 'ensembl_gene_id')
                    im.str.gene.sig <- hk_immunStroma$ensembl_gene_id[im.str.gene.sig]
                tumour.purity <- singscore::rankGenes(
                    assay(x = data, i = assay.name.to.estimate.purity)
                     )
                tumour.purity <- singscore::simpleScore(
                    rankData = tumour.purity,
                    upSet = im.str.gene.sig
                    )
                tumour.purity <- tumour.purity$TotalScore
                sample.annotation[['tumour.purity']] <- 1 - tumour.purity
            } else if (estimate.tumor.purity == 'both'){
                printColoredMessage(
                    message = '- Estimating tumour purity using both ESTIMATE and singscore methods:',
                    color = 'blue',
                    verbose = verbose
                    )
                printColoredMessage(
                    message = '-- Estimating tumour purity using the ESTIMATE method:',
                    color = 'blue',
                    verbose = verbose
                    )
                tumour.purity <- tidyestimate::filter_common_genes(
                    df = assay(x = data, i = assay.name.to.estimate.purity),
                    id = gene.ids.to.estimate.tumor.purity,
                    tidy = FALSE,
                    tell_missing = verbose,
                    find_alias = TRUE
                    )
                tumour.purity <- tidyestimate::estimate_score(
                    df = tumour.purity,
                    is_affymetrix = TRUE
                    )
                tumour.purity.estimate <- tumour.purity$purity
                printColoredMessage(
                    message = '- Estimating tumour purity using the singscore method:',
                    color = 'blue',
                    verbose = verbose
                    )
                im.str.gene.sig <- hk_immunStroma$immune.gene.signature == 'TRUE' |
                    hk_immunStroma$stromal.gene.signature == 'TRUE'
                if (gene.group == "entrezgene_id"){
                    im.str.gene.sig <- hk_immunStroma$entrezgene_id[im.str.gene.sig]
                } else if (gene.group == 'hgnc_symbol'){
                    im.str.gene.sig <- hk_immunStroma$hgnc_symbol[im.str.gene.sig]
                } else if (gene.group == 'ensembl_gene_id')
                    im.str.gene.sig <- hk_immunStroma$ensembl_gene_id[im.str.gene.sig]
                tumour.purity <- singscore::rankGenes(assay(x = data, i = assay.name.to.estimate.purity))
                tumour.purity <- singscore::simpleScore(
                    rankData = tumour.purity,
                    upSet = im.str.gene.sig
                    )
                tumour.purity.singscore <- tumour.purity$TotalScore
                data[['tumour.purity.estimate']] <- tumour.purity.estimate
                data[['tumour.purity.singscore']] <- 1 - tumour.purity.singscore
                if (isTRUE(scale.singscore.values)){
                    tps <- data[['tumour.purity.singscore']]
                    tpe <- data[['tumour.purity.estimate']]
                    rtps <- range(tps)
                    rtpe <- range(tpe)
                    stps <- (tps - min(rtps)) / (max(rtps) - min(rtps)) * (max(rtpe) - min(rtpe)) + min(rtpe)
                    data[['tumour.purity.singscore.scaled']]<- stps
                }
            }
        }
        # Saving the outputs ####
        printColoredMessage(
            message = paste0('The SummarizedExperiment object:'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '-',
                nrow(data),
                ' measurements (e.g. genes) and ',
                ncol(data),
                ' assays (e.g.samples)'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '-',
                length(assays(data)),
                ' data sets (assays)'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '-',
                ncol(SummarizedExperiment::colData(data)),
                ' annotations for the samples'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(
            message = paste0(
                '-',
                ncol(SummarizedExperiment::rowData(data)),
                ' annotations for the genes'),
            color = 'blue',
            verbose = verbose
        )
        printColoredMessage(message = '------------The prepareSeObj function finished.',
                            color = 'white',
                            verbose = verbose)
        return(data)

    }
    # else {
    #     stop('The "data" must be either a SummarizedExperiment object or a list of data set (s).')
    # }
}


