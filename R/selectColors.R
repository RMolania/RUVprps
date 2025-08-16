#' Print colored messages.
#'
#' @param nb.color Message to be printed.
#' @param group The color of the message.
#'
#' @importFrom Polychrome createPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggsci pal_npg
#'
#' @export


selectColors <- function(
        nb.color,
        group = 'pan.selection.a'
        ){
    if (group == 'nrc'){
        current.colors <- pal_npg("nrc")(10)
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
        current.colors <- rev(createPalette(100, seedcolors = c("#000000", "#FFFFFF")))

    }
    return(current.colors[nb.color])
}
