#' Adds the TCGA RNA-seq batch information
#'
#' @author Ramyar Molania
#'
#' @description
#' This function adds the TCGA RNA-seq batch information to a `SummarizedExperiment` object, using metadata
#' provided by TCGA sample barcodes or external batch information files.
#'
#' @param se.obj A `SummarizedExperiment` object containing TCGA RNA-seq data.
#' @param tissue.type A data frame or character string specifying batch information. This can be a metadata table
#' with batch IDs for each sample, or a file path to such a table.
#' @param factors Character. A TCGA cancer type code (e.g., `"BRCA"`, `"PRAD"`, `"COAD"`) used to match batch
#' annotations with the corresponding dataset.
#' @param cancer.type Logical. If `TRUE`, all samples in the `SummarizedExperiment` object are retained, even if
#' batch information is missing. If `FALSE`, samples without batch annotations are removed. The default is `TRUE`.
#' @param missing.samples.name Character. A label to assign to samples without batch information. The default is
#' `"Unknown"`.
#' @param verbose UUU
#'
#' @importFrom stringr str_split
#' @importFrom tidyr separate
#' @importFrom purrr map_int
#'
#' @return A `SummarizedExperiment` object with TCGA RNA-seq batch information added to the `colData`.
#'
#' @export


addTcgaClinicalInfo <- function(
        se.obj,
        tissue.type,
        factors = c('OS', 'OS.time'),
        cancer.type = NULL,
        missing.samples.name = 'DF',
        verbose = TRUE
        ){
    # tcga.clinic.info <- read.csv('../RequiredData_Temp/TCGA.Liuetal.ClinicalData.PanCancer.csv')
    # usethis::use_data(tcga.clinic.info)
    all.tcga.clinic <- colnames(tcga.clinic.info)
    if (sum(factors %in% colnames(tcga.clinic.info)) != length(factors)){
        stop('All or some of the "factors" cannot be found in the TCGA clinical information data.')
    }
    if (!is.null(cancer.type)){
        if (!is.character(cancer.type) | length(cancer.type) > 1){
            stop('The cancer type must be a character string of a TCGA cacncer types.')
        }
        if (!cancer.type %in% unique(tcga.clinic.info$type)){
            stop('The "cancer.type" cannot be found in the TCGA clinical information data')
        }
        tcga.clinic.info <- tcga.clinic.info[tcga.clinic.info$type == cancer.type , ]
        printColoredMessage(
            message = paste0(
                '- The specifed cancer type has clinical information for ',
                nrow(tcga.clinic.info),
                ' samples.'),
            color = 'blue',
            verbose = verbose)
    }

    # Checking sample barcodes
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
    # df <- tibble(barcode = tcga.clinic.info$bcr_patient_barcode[10]) %>%
    #     separate(
    #         barcode,
    #         into = c("Project", "TSS", "Participant"), sep = "-") %>%
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
    if (a.barcode > 4){
        samples.barcode <- sub("^(([^-]+-){2}[^-]+)-.*", "\\1", samples.barcode)
    }
    common.samples <- intersect(samples.barcode, tcga.clinic.info$bcr_patient_barcode)
    printColoredMessage(
        message = paste0(
            'There are ',
            length(common.samples),
            ' samples in common.'),
        color = 'blue',
        verbose = verbose
        )
    se.obj[['OS']] <- -1
    se.obj[['OS.time']] <- -1
    for (i in common.samples){
        index <- tcga.clinic.info$bcr_patient_barcode == i
        se.obj$OS[samples.barcode == i & se.obj[[tissue.type]] == 'Tumor'] <-  tcga.clinic.info$OS[index]
        se.obj$OS.time[samples.barcode == i & se.obj[[tissue.type]] == 'Tumor'] <- tcga.clinic.info$OS.time[index]
    }
    return(se.obj)
}





