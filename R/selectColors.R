#' Print colored messages.
#'
#' @param nb.color Message to be printed.
#' @param group The color of the message.
#'
#' @importFrom RColorBrewer brewer.pal
#'
#' @export


selectColors <- function(
        nb.color,
        group = 'pan.selection.a'
        ){
    if (group == 'nrc'){
        current.colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF")
    }
    if (group == 'pan.selection.a'){
        current.colors <- c(
            RColorBrewer::brewer.pal(8, "Dark2")[-5],
            RColorBrewer::brewer.pal(10, "Paired"),
            RColorBrewer::brewer.pal(12, "Set3"),
            RColorBrewer::brewer.pal(9, "Blues")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Oranges")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Greens")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Purples")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Reds")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Greys")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "BuGn")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "PuRd")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "BuPu")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "YlGn")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(10, "Paired")
        )
    }
    if (group == 'pan.selection.b'){
        current.colors <- c(
            "#E7B800", "#2E9FDF", 'red4',
            RColorBrewer::brewer.pal(8, "Dark2")[-5],
            RColorBrewer::brewer.pal(10, "Paired"),
            RColorBrewer::brewer.pal(12, "Set3"),
            RColorBrewer::brewer.pal(9, "Blues")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Oranges")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Greens")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Purples")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Reds")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "Greys")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "BuGn")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "PuRd")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "BuPu")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(9, "YlGn")[c(8, 3, 7, 4, 6, 9, 5)],
            RColorBrewer::brewer.pal(10, "Paired")
        )
    }
    if (group == 'pan.selection.c'){
        current.colors <- rev(nc_colors <- c(
            "#564547", "#E0E0E0", "#FE002A", "#1C47FF", "#0DFB26",
            "#FC0DDD", "#FDA200", "#AE164F", "#16C0FC", "#D9E75C",
            "#1CFCC9", "#AF75D2", "#006C1C", "#C416FF", "#FEA1CE",
            "#A7833B", "#009093", "#DE4F00", "#9B0081", "#224B80",
            "#FD0DA3", "#77DD6D", "#FD9587", "#5F7DE4", "#BDCBFE",
            "#16F9FD", "#ADC487", "#F12268", "#8E3D55", "#FDD138",
            "#7300C7", "#FF86EE", "#F2BBAF", "#A784AB", "#AAF90D",
            "#A0F0DE", "#B100BD", "#FA66B6", "#95380D", "#607E63",
            "#F400FF", "#F6DF92", "#65C58C", "#FF6C79", "#6C9B00",
            "#FAB9FC", "#FBB675", "#91E0F5", "#774B2E", "#723B83",
            "#C46CF8", "#636616", "#8B9FAF", "#D56C1C", "#3591B8",
            "#00FE9F", "#56329B", "#FACAE5", "#E1708E", "#264F5A",
            "#A2A126", "#8C16FB", "#A4917A", "#7684CA", "#AE6900",
            "#BB69AD", "#DC35B6", "#EBDBB8", "#AB7885", "#DB1C0D",
            "#896DFB", "#FF66FD", "#BE7566", "#0DB9A0", "#C3EF90",
            "#F34787", "#AD3273", "#1632CD", "#AF90FE", "#820099",
            "#B2221C", "#009816", "#9CD12E", "#5CDC16", "#F816BC",
            "#9EC2AD", "#0DD8F8", "#0DA7FF", "#F4E500", "#1C79FC",
            "#C6AFFC", "#68698C", "#812E6C", "#D1C6E5", "#B6EEBE",
            "#E987FB", "#C49500", "#C6007C", "#88BDFD", "#F7A1AD"))

    }
    return(current.colors[nb.color])
}





