#' Prepare variables for PRPS
#'
#' @author Ramyar Molania
#'
#' @description
#' This function prepares biological and unwanted variation variables for Partial Residualization and Partial Scoring (PRPS)
#' analysis. It includes optional clustering, extreme group selection, batch coverage, and variable association assessment.
#'
#' @param se.obj A `SummarizedExperiment` object containing gene expression data and sample annotations.
#' @param main.variable Character. The primary variable of interest to be used in PRPS.
#' @param other.variables Character vector. Other covariates or unwanted variation variables to include in the analysis.
#' @param nb.mnn Numeric. Number of nearest neighbors to consider when assessing variable relationships. Default is 10.
#' @param min.sample.for.ps Numeric. Minimum number of samples required in a group to perform partial scoring. Default is 5.
#' @param clustering.method Character. Clustering method to apply to continuous variables (e.g., `"kmeans"`, `"hierarchical"`). Default is `"kmeans"`.
#' @param nb.clusters Numeric. Number of clusters for the main variable if clustering is applied. Default is 3.
#' @param other.uv.clustering.method Character. Clustering method for other unwanted variation variables. Default is `"kmeans"`.
#' @param nb.other.uv.clusters Numeric. Number of clusters for other unwanted variation variables. Default is 3.
#' @param select.extreme.groups Logical. If `TRUE`, selects groups with extreme values for the main variable. Default is `FALSE`.
#' @param cover.all.batches Logical. If `TRUE`, ensures that selected groups cover all batches present in the dataset. Default is `TRUE`.
#' @param nb.batches.to.cover Numeric. Minimum number of batches to cover when `cover.all.batches = TRUE`. Default is 1.
#' @param assess.variables.association Logical. If `TRUE`, computes associations between variables using `canCorPairs`. Default is `TRUE`.
#' @param plot.output Logical. If `TRUE`, generates heatmaps or plots summarizing variable clusters or associations. Default is `TRUE`.
#' @param verbose Logical. If `TRUE`, displays messages during function execution. Default is `TRUE`.
#'
#' @return A list containing processed main and other variables, clusters, and optionally plots and variable associations,
#' ready for PRPS analysis.
#'
#' @importFrom variancePartition canCorPairs
#' @importFrom ComplexHeatmap Heatmap
#' @export


prepareVariableForPrPs <- function(
        se.obj,
        main.variable,
        other.variables = NULL,
        nb.mnn = 1,
        min.sample.for.ps = 3,
        clustering.method = 'cut',
        nb.clusters = 3,
        other.uv.clustering.method = 'cut',
        nb.other.uv.clusters = 3,
        select.extreme.groups = FALSE,
        cover.all.batches = FALSE,
        nb.batches.to.cover = 2,
        assess.variables.association = TRUE,
        plot.output = TRUE,
        verbose = TRUE
        ){
    # Assessing association between variables ####
    variables <- c(main.variable, other.variables)
    if (isTRUE(assess.variables.association) & length(variables) == 1){
        printColoredMessage(
            message = '- Only one variable is provided, then the "assess.variables.association" cannot be applied.',
            color = 'green',
            verbose = verbose
            )
    }
    if (isTRUE(assess.variables.association) & length(variables) > 1){
        sample.annotation <- colData(se.obj)
        variables <- gsub(' ', '.', variables)
        form <- as.formula(paste0(
            '~',
            paste0(variables, collapse = '+')
            ))
        var.association <- canCorPairs(
            formula = form,
            data = sample.annotation,
            showWarnings = FALSE
            )
        ht.vars <- ComplexHeatmap::Heatmap(
            matrix = round(x = var.association, digits = 3),
            col = c('grey', 'darkgreen'),
            heatmap_legend_param = list(title = 'Correlation')
            )
        if (isTRUE(plot.output)) print(ht.vars)
        cor.df <- as.data.frame(as.table(var.association))
        colnames(cor.df) <- c('Variable1', 'Variable2', 'Corr')
        cor.df <- cor.df[cor.df$Variable1 != cor.df$Variable2, ]
        cor.df <- cor.df[!duplicated(t(apply(cor.df[,1:2], 1, sort))), ]
        if (sum(cor.df$Corr > 0.8) > 0 ){
            high.corr.vars <- cor.df[abs(cor.df$Corr) > .8, ]
            high.corr.vars <- lapply(
                1:nrow(high.corr.vars),
                function(x){
                    printColoredMessage(
                        message = paste0(
                            'The "',
                            high.corr.vars$Variable1[x],
                            '" and "'
                            ,
                            high.corr.vars$Variable2[x],
                            '" variables are highly correlated, this can be a problem when "apply.other.variables" is set to TRUE.'
                        ),
                        color = 'red',
                        verbose = verbose)
                })
        }
    }
    # Considering only on variable at the time ####
    if (is.null(other.variables)){
        ## Continuous variable ####
        if (is.numeric(se.obj[[main.variable]])){
            ## Keeping original values of the main unwanted variable ####
            initial.variable <- se.obj[[main.variable]]
            ## Grouping main the main unwanted variable ####
            grouped.variable <- groupContinuousVariable(
                se.obj = se.obj,
                variable = main.variable,
                nb.clusters = nb.clusters,
                clustering.method = clustering.method,
                perfix = '.',
                plot.output = plot.output,
                verbose = verbose
            )
            ## Selecting subgroups of the main unwanted variable with highest and lowest values ####
            if (isTRUE(select.extreme.groups)){
                printColoredMessage(
                    message = paste0(
                        '- Selecting the two subgroups of the "',
                        main.variable,
                        '" variable with the highest and lowest values.'),
                    color = 'blue',
                    verbose = verbose
                    )
                max.group <- grouped.variable[initial.variable == max(initial.variable)]
                if (sum(grouped.variable == max.group) < max(min.sample.for.ps, nb.mnn)){
                    stop(paste0(
                        'The subgroup with highest values of the ',
                        main.variable,
                        'variable does not have enough samples for finding MNN or creatinf PS.'))
                }
                min.group <- grouped.variable[initial.variable == min(initial.variable)]
                if (sum(grouped.variable == min.group) < max(min.sample.for.ps, nb.mnn)){
                    stop(paste0(
                        'The subgroup with highest values of the ',
                        main.variable,
                        'variable does not have enough samples for finding MNN or creatinf PS.'))
                }
                selected.subgroups <- c(min.group, max.group)
                 all.variables <- data.frame(
                    sampl.ids = colnames(se.obj),
                    variable = initial.variable,
                    groups = grouped.variable,
                    selected = ifelse(grouped.variable %in% c(selected.subgroups), 'TRUE', 'FALSE')
                    )
            }
            if (isFALSE(select.extreme.groups)) {
                 all.variables <- data.frame(
                    sampl.ids = colnames(se.obj),
                    variable = initial.variable,
                    groups = grouped.variable,
                    selected = 'TRUE'
                    )
            }
        }
        ## Categorical variable ####
        if (!is.numeric(se.obj[[main.variable]])){
            if (length(unique(se.obj[[main.variable]])) == 1){
                stop(paste0(
                    'To create MNN, the "',
                    main.variable,
                    '" variable must have at least two groups/levels.')
                )
            }
            covered.batches <- findRepeatingPatterns(
                vec = se.obj[[main.variable]],
                n.repeat = max(min.sample.for.ps, nb.mnn)
                )
            if (isTRUE(cover.all.batches)){
                if (length(covered.batches) != length(unique(se.obj[[main.variable]]))){
                    stop(paste0(
                        'All or some subgroups/batches of the "',
                        main.variable,
                        '" variable does not have enough samples for finding MNN or creatinf PS')
                        )
                } else {
                    all.variables <- data.frame(
                        sampl.ids = colnames(se.obj),
                        variable =  se.obj[[main.variable]],
                        groups = se.obj[[main.variable]],
                        selected = ifelse(se.obj[[main.variable]] %in% c(covered.batches), 'TRUE', 'FALSE')
                        )
                }
            }
            if (isFALSE(cover.all.batches)){
                if ( nb.batches.to.cover > length(unique(se.obj[[main.variable]]))){
                    printColoredMessage(
                        message = paste0(
                            '- The " nb.batches.to.cover" is larger than the number of subgroups of the',
                            main.variable,
                            'variable, then " nb.batches.to.cover" is capped.'),
                        color = 'blue',
                        verbose = verbose
                        )
                    nb.batches.to.cover <- length(unique(se.obj[[main.variable]]))
                }
                if (length(covered.batches) <  nb.batches.to.cover){
                    stop(paste0(
                        nb.batches.to.cover,
                        'subgroups/batches of the "',
                        main.variable,
                        '" variable does not have enough samples for finding MNN or creatinf PS')
                    )
                }
                if (length(covered.batches) >=  nb.batches.to.cover){
                    all.batches <- sort(table(se.obj[[main.variable]]), decreasing = TRUE)
                    selected.subgroups <- names(all.batches)[1: nb.batches.to.cover]
                     all.variables <- data.frame(
                        sampl.ids = colnames(se.obj),
                        variable =  se.obj[[main.variable]],
                        groups = se.obj[[main.variable]],
                        selected = ifelse(se.obj[[main.variable]] %in% c(covered.batches), 'TRUE', 'FALSE')
                        )
                     all.variables
                }
            }
        }
    }
    # Considering all variables ####
    if (!is.null(other.variables)){
        if (is.numeric(se.obj[[main.variable]])){
            initial.variable <- se.obj[[main.variable]]
            grouped.variable <- groupContinuousVariable(
                se.obj = se.obj,
                variable = main.variable,
                nb.clusters = nb.clusters,
                clustering.method = clustering.method,
                perfix = '.',
                plot.output = plot.output,
                verbose = verbose
                )
        } else grouped.variable <- se.obj[[main.variable]]

        if (isTRUE(select.extreme.groups) & is.numeric(se.obj[[main.variable]])){
            printColoredMessage(
                message = paste0(
                    '- Selecting the two subgroups of the "',
                    main.variable,
                    '" variable with the highest and lowest values.'),
                color = 'blue',
                verbose = verbose
                )
            max.group <- grouped.variable[initial.variable == max(initial.variable)]
            if (sum(grouped.variable == max.group) < max(min.sample.for.ps, nb.mnn)){
                stop(paste0(
                    'The subgroup with highest values of the ',
                    main.variable,
                    'variable does not have enough samples for finding MNN or creatinf PS.'))
                }
            min.group <- grouped.variable[initial.variable == min(initial.variable)]
            if (sum(grouped.variable == min.group) < max(min.sample.for.ps, nb.mnn)){
                stop(paste0(
                    'The subgroup with highest values of the ',
                    main.variable,
                    'variable does not have enough samples for finding MNN or creatinf PS.'))
                }
            selected.subgroups <- c(min.group, max.group)
        }
        ## grouping other variables ####
        homo.uv.groups <- createHomogeneousUVGroups(
            se.obj = se.obj,
            uv.variables = other.variables,
            nb.clusters = nb.other.uv.clusters,
            clustering.method = other.uv.clustering.method,
            check.se.obj = FALSE,
            save.se.obj = FALSE,
            remove.na = 'none',
            verbose = verbose
            )
        all.uv.groups <- data.frame(
            sample.ids = colnames(se.obj),
            main.uv = grouped.variable,
            other.uv = homo.uv.groups
            )
        if (isTRUE(select.extreme.groups) & is.numeric(se.obj[[main.variable]])){
            all.uv.groups <- all.uv.groups[all.uv.groups$main.uv %in% selected.subgroups , ]
            covered.batches <- lapply(
                unique(all.uv.groups$other.uv),
                function(s){
                    subgroups.size <- findRepeatingPatterns(
                        vec = all.uv.groups$main.uv[all.uv.groups$other.uv == s],
                        n.repeat = max(min.sample.for.ps, nb.mnn))
                })
            names(covered.batches) <- unique(all.uv.groups$other.uv)
            length.covered.batches <- sapply(covered.batches, length)
            if (sum(length.covered.batches == 2) == 0){
                stop(paste0(
                    'There is no at least two subgroups of the "',
                    main.variable,
                    '" variable that have enough samples within the other variable group to find MNN or creating PRPS')
                )
            }
            selected.subgroups <- covered.batches[length.covered.batches == 2]
            all.uv.groups <- all.uv.groups[all.uv.groups$other.uv %in% names(selected.subgroups) , ]
             all.variables <- data.frame(
                sampl.ids = colnames(se.obj),
                variable =  se.obj[[main.variable]],
                main.variable = grouped.variable,
                other.variables = homo.uv.groups,
                selected = ifelse(colnames(se.obj) %in% all.uv.groups$sample.ids, 'TRUE', 'FALSE')
            )
        } else{
            covered.batches <- lapply(
                unique(all.uv.groups$other.uv),
                function(s){
                    subgroups.size <- findRepeatingPatterns(
                        vec = all.uv.groups$main.uv[all.uv.groups$other.uv == s],
                        n.repeat = max(min.sample.for.ps, nb.mnn))
                })
            names(covered.batches) <- unique(all.uv.groups$other.uv)
            length.covered.batches <- sapply(covered.batches, length)
            if (sum(length.covered.batches == 1) == length(unique(all.uv.groups$other.uv))){
                stop(paste0(
                    'There is no at least two subgroups of the "',
                    main.variable,
                    '" variable that have enough samples within the other variable group to find MNN or creating PRPS')
                )
            }
            if (isTRUE(cover.all.batches) & isFALSE(select.extreme.groups)){
                if (sum(length.covered.batches == length(unique(all.uv.groups$main.uv))) == 0 ){
                    stop(paste0(
                        'There is no any subgroups of other variables that have enough samples for all the subgropus of "',
                        main.variable,
                        '" variable,to find MNN or creating PRPS.')
                    )
                }
            }
            if (isFALSE(cover.all.batches)){
                if(is.null(nb.batches.to.cover)){
                    selected.subgroups <- names(covered.batches[length.covered.batches > 1])
                    all.uv.groups <- all.uv.groups[all.uv.groups$other.uv %in% selected.subgroups, ]
                }
                if (is.numeric( nb.batches.to.cover)){
                    if (sum(length.covered.batches >=  nb.batches.to.cover) == 0){
                        stop(paste0(
                            'There is no any subgroups of other variables that have enough samples for all the subgropus of "',
                            main.variable,
                            '" variable,to find MNN or creating PRPS.')
                        )
                    }
                    covered.batches <- covered.batches[length.covered.batches >=  nb.batches.to.cover]
                    selected.subgroups <- lapply(
                        names(covered.batches),
                        function(x){
                            all.uv.groups.sub <- droplevels(all.uv.groups[all.uv.groups$other.uv == x ,  ])
                            if (length(unique(all.uv.groups.sub$other.uv)) <  nb.batches.to.cover){
                                new.nb.batches.to.cover = length(unique(all.uv.groups.sub$main.uv))
                            } else new.nb.batches.to.cover <-  nb.batches.to.cover
                            names(sort(table(all.uv.groups.sub$main.uv), decreasing = TRUE))[1:new.nb.batches.to.cover]
                        })
                    selected.subgroups <- names(covered.batches)
                    all.variables <- data.frame(
                        sampl.ids = colnames(se.obj),
                        variable =  se.obj[[main.variable]],
                        main.variable = grouped.variable,
                        other.variables = homo.uv.groups,
                        selected = ifelse(homo.uv.groups %in% c(selected.subgroups), 'TRUE', 'FALSE'))
                    gg <- unique(unlist(unname(covered.batches)))
                    all.variables$selected[!all.variables$main.variable %in% gg] <- 'FALSE'
                }
            }
        }
         all.variables
    }
    return(all.variables)
}


# printColoredMessage(
#     message = '- Assessing and grouping the other specified unwanted variable(s):',
#     color = 'magenta',
#     verbose = verbose
# )
#
#
# ## Finding subgroups that has enough samples for generating PS ####
# printColoredMessage(
#     message = '- Finding subgroups with respect to other unwanted variable(s) that have enoug samples for PRPS:',
#     color = 'magenta',
#     verbose = verbose
# )
#
#
# printColoredMessage(
#     message = paste0(
#         '- There are ',
#         length(covered.batches),
#         ' subgroups with respect to other unwanted variables that have enough samples for generating PS.'),
#     color = 'blue',
#     verbose = verbose
# )
# printColoredMessage(
#     message = '- Checking the distribution on main unwanted variable across subgroups with respect to other unwanted variables.',
#     color = 'blue',
#     verbose = verbose
# )
# covered.batches.table <- as.data.frame(
#     table(all.uv.groups$main.uv, all.uv.groups$other.uv)
# )
# covered.batches.table$groups <- covered.batches.table$Freq >= max(min.sample.for.ps, nb.mnn)
#
# ## Plotting the distribution of th main unwanted variable across the subgroups of other unwanted variable ####
# if (isTRUE(plot.output)){
#     printColoredMessage(
#         message = '- Plotting distribution of samples across subgroups with respect to other unwanted variable(s):',
#         color = 'magenta',
#         verbose = verbose
#     )
# }
# covered.batches.table$groups <- ifelse(covered.batches.table$groups == 'TRUE', 'selected', 'unselected')
# covered.batches.plot <- ggplot(covered.batches.table, aes(x = Var1, y = Var2, color = groups)) +
#     geom_point(size = 6) +
#     geom_text(aes(label = Freq , hjust = 0.5, vjust = 0.5), color = 'black') +
#     xlab(main.uv.variable) +
#     ylab('Homogeneous groups\n(other unwanted variables)') +
#     theme_bw() +
#     theme(
#         legend.key = element_blank(),
#         axis.line = element_line(colour = 'black', linewidth = 1),
#         axis.title.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 14),
#         axis.text.x = element_text(size = 14, angle = 90, vjust = 1, hjust = 1),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 18),
#         strip.text.y = element_text(size = 0)
#     )
# if (isTRUE(plot.output)) print(covered.batches.plot)
#
# ## Finding subgroups of other unwanted variable that have enough samples for PRPS ####
# printColoredMessage(
#     message = '- Finding subgroups of other unwanted variable that have enough samples for generating PRPS:',
#     color = 'magenta',
#     verbose = verbose
# )
# selected.covered.batches <- lapply( covered.batches, length)
# if (sum(selected.covered.batches == 1) == length(selected.covered.batches)){
#     stop(paste0(
#         ' Non of the sample groups with respect to the other unwanted variables that have at least ',
#         max(min.sample.for.ps, nb.mnn),
#         ' samples across at least two sub-groups of the ',
#         main.uv.variable,
#         ' variable.'))
# }
# if (sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) == 0){
#     printColoredMessage(
#         message = paste0(
#             '- Non of the sample groups with respect to the other unwanted variables have at least ',
#             max(min.sample.for.ps, nb.mnn),
#             ' samples across all the sub-groups of the main unwanted variable: "',
#             main.uv.variable,
#             '" variable.'),
#         color = 'blue',
#         verbose = verbose
#     )
#     if (isFALSE(check.prps.connectedness)){
#         printColoredMessage(
#             message = '-- We recommend specifiying the "check.prps.connectedness"',
#             color = 'red',
#             verbose = verbose
#         )
#     } else if (isTRUE(check.prps.connectedness)){
#         checkPRPSconnectedness(
#             data.input = table(all.uv.groups$main.uv, all.uv.groups$other.uv),
#             min.samples = c(nb.mnn, min.sample.for.ps),
#             batch.name = main.uv.variable,
#             verbose = verbose
#         )
#     }
# }
# if (sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) > 0  ){
#     printColoredMessage(
#         message = paste0(
#             '-- There are ',
#             sum(selected.covered.batches == length(unique(all.uv.groups$main.uv))) ,
#             ' groups with respect to the other unwanted variables that have at least ',
#             max(min.sample.for.ps, nb.mnn),
#             ' samples across all sub-groups of the main unwanted variable: ',
#             main.uv.variable,
#             ' variable.'),
#         color = 'blue',
#         verbose = verbose
#     )
# }

