#' Adds the TCGA RNA-seq batch information
#'
#' @author Ramyar Molania
#'
#' @description
#' This function adds the TCGA RNA-seq batch information to a `SummarizedExperiment` object, using metadata
#' provided by TCGA sample barcodes or external batch information files.
#'
#' @param se.obj A `SummarizedExperiment` object containing TCGA RNA-seq data.
#' @param batch.info A data frame or character string specifying batch information. This can be a metadata table
#' with batch IDs for each sample or a file path to such a table.
#' @param cancer.type Character. A TCGA cancer type code (e.g., `"BRCA"`, `"PRAD"`, `"COAD"`) used to match batch
#' annotations with the corresponding dataset.
#' @param keep.all.samples Logical. If `TRUE`, all samples in the `SummarizedExperiment` object are retained even if
#' batch information is missing. If `FALSE`, samples without batch annotations are removed. The default is `TRUE`.
#' @param use.lables Logical. Indicates whether to use TCGA barcode labels to infer batch information directly.
#' If `TRUE`, batch identifiers are parsed from the TCGA sample barcodes.
#' @param missing.samples.name Character. A label to assign to samples without batch information. By default,
#' `"Unknown"`.
#'
#' @importFrom stringr str_split
#' @importFrom tidyr separate
#' @importFrom purrr map_int
#'
#' @return A `SummarizedExperiment` object with TCGA RNA-seq batch information added to the `colData`.
#'
#' @export

addTcgaBatchInfo <- function(
        se.obj,
        batch.info = 'all',
        cancer.type = NULL,
        keep.all.samples = TRUE,
        use.lables = FALSE,
        missing.samples.name = 'DF'
        ){
    colnames(tcga.batch.info)[11:13] <- c('Years', 'Months', 'Days')
    tcga.batchs <- c('Years','Months','Days','Plates','TSS','Center')
    # Checking batch information ####
    if (length(batch.info) > 1){
        if (sum(batch.info %in% tcga.batchs) != length(batch.info)){
            stop('All or some of the "batch.info" cannot be found in the TCGA batch information.')
        }
    }
    if (identical(x = batch.info, 'all')){
        batch.info = tcga.batchs
    }
    # Selecting cancer type ####
    if (!is.null(cancer.type)){
        tcga.batch.info <- tcga.batch.info[tcga.batch.info$Cancer.type == cancer.type , , drop = FALSE]
    }
    # Checking sample bar codes ####
    ## TCGA id
    samples.barcode <- colnames(se.obj)
    tcga.id <- substr(x = samples.barcode, start = 1, stop = 4)
    if (isFALSE(sum(tcga.id == 'TCGA') == ncol(se.obj))){
        stop('The colum names of the SummarizedExperimet object must start with "TCGA"')
    }
    ## Checking the length of the barcodes
    samples.barcode <- gsub("[_.]", "-", samples.barcode)
    length.barcodes <- map_int(samples.barcode, ~ length(str_split(.x, "-")[[1]]))
    if (length(table(length.barcodes)) > 1){
        stop('The sample ids do not have the consistent TCGA information, please check.')
    }
    ## plot an example
    # df <- tibble(barcode = tcga.batch.info$Sample.ids.full.a[1]) %>%
    #     separate(
    #         barcode,
    #         into = c("Project", "TSS", "Participant", "SampleType", "Portion", "Plates", "Center"), sep = "-") %>%
    #     mutate(ID = row_number())
    # df.long <- df %>%
    #     pivot_longer(cols = -c(ID), names_to = "Component", values_to = "Value") %>%
    #     group_by(ID) %>%
    #     mutate(Position = row_number())
    # separator.positions <- data.frame(
    #     Position = seq(1.5, max(df.long$Position) - 0.5, by = 1),
    #     ID = 1,
    #     sep = "-"
    # )
    # p1 <- ggplot(df.long, aes(x = Position, y = factor(ID))) +
    #     geom_text(aes(label = Value), size = 6, color = "black") +
    #     geom_text(data = separator.positions, aes(x = Position, y = factor(ID), label = sep),
    #               size = 8, inherit.aes = FALSE) +
    #     scale_x_continuous(breaks = df.long$Position, labels = df.long$Component) +
    #     scale_y_discrete(expand = expansion(mult = c(0., 0.))) +  # reduces vertical padding
    #     guides(fill = "none") +
    #     xlab('-------------------Common TCGA barcode -------------------') +
    #     theme_minimal() +
    #     theme(
    #         axis.title.y = element_blank(),
    #         axis.text.y = element_blank(),
    #         axis.ticks = element_blank(),
    #         panel.grid = element_blank(),
    #         plot.title = element_text(size = 20),
    #         axis.title.x = element_text(size = 20),
    #         axis.text.x = element_text(margin = margin(t = 1, b = 25), size = 18)
    #     )
    # ggarrange(p1, p1, ncol = 1)

    ## Checking the provided information
    a.barcode <- length.barcodes[1]
    if (a.barcode < 2){
        stop('There are not enough information in the sample ids to find batch information.')
    }
    # Extract the 4th part of each barcode
    if (length(a.barcode) == 3){
        index <- intersect(tcga.batch.info$Sample.ids.3, colnames(se.obj))
        tcga.batch.info <- tcga.batch.info[tcga.batch.info$Samples %in% index, ]
        read.se.obj <- read.se.obj[ , tcga.batch.info$Samples]
        ## add batch information to the SummarizedExperiment Object
        read.se.obj$Years <- as.factor(tcga.batch.info$Year)
        read.se.obj$Plates <- tcga.batch.info$Plates
        read.se.obj$TSS <- tcga.batch.info$TSS
    }
    if (length(a.barcode) == 5){
        stop('ffff')
    }

    if (length(a.barcode) == 6){
        stop('fffff')
    }

    if (a.barcode == 7){
        common.samples <- intersect(tcga.batch.info$Sample.ids.full.a, colnames(se.obj))
        if (length(common.samples) == ncol(se.obj)){
            tcga.batch.info <- tcga.batch.info[tcga.batch.info$Sample.ids.full.a %in% common.samples, ]
            match.samples <- match(colnames(se.obj), tcga.batch.info$Sample.ids.full.a)
            tcga.batch.info <- tcga.batch.info[match.samples , ]
            # if (isTRUE(include.sample.info)){
            #     for(i in sample.info) tcga.batch.info[i] <- se.obj[[i]]
            # }
            for(i in batch.info) se.obj[[i]] <- tcga.batch.info[[i]]
        }
        if (length(common.samples) < ncol(se.obj)){
            tcga.batch.info <- tcga.batch.info[tcga.batch.info$Sample.ids.full.a %in% common.samples, ]
            match.samples <- match(colnames(se.obj), tcga.batch.info$Sample.ids.full.a)
            tcga.batch.info <- tcga.batch.info[match.samples , ]
            tcga.batch.info$New.ids <- colnames(se.obj)
            for(i in batch.info) se.obj[[i]] <- tcga.batch.info[[i]]
        }
    }
    return(se.obj)
}
