#' plotVenn
#'
#' Plot a venn diagram with VennDiagram package using some nice defaults
#'
#' @param x matrix of values (similar to what is fed to VennDiagram in limma package).
#' @param ... additional parameters passed down to venn.diagrama()
#' @param add.universe logical; whether to add a group containing all elements.
#'
#' @return NULL
#' @export
#'
#' @examples
#' NULL
plot_venn <- function(x, ...) {
  UseMethod("plot_venn")
}

#' @rdname plot_venn
#' @export
plot_venn.matrix <- function(x, add.universe = FALSE, euler = FALSE, scaled = FALSE, filename = NULL, fontfamily = "sans", cat.fontfamily = "sans", main.fontfamily = "sans", ...) {
  m2l <- function(x, add.universe = FALSE) {
    l <- lapply(seq_len(ncol(x)), function(i) {
      if (is.null(rownames(x)))
        (1:nrow(x))[x[,i] != 0]
      else
        rownames(x)[x[,i] != 0]
    })
    if (is.null(colnames(x)))
      names(l) <- paste0("group-", 1:length(l))
    else
      names(l) <- colnames(x)
    if (add.universe)
      l$total <- rownames(x)
    l
  }
  plot_venn(
    m2l(x, add.universe = add.universe),
    euler = euler,
    scaled = scaled,
    filename = filename,
    fontfamily = fontfamily,
    cat.fontfamily = cat.fontfamily,
    main.fontfamily = main.fontfamily,
    ...
  )
}

#' @rdname plot_venn
#' @export
plot_venn.data.frame <- function(x, euler = FALSE, scaled = FALSE, filename = NULL, fontfamily = "sans", cat.fontfamily = "sans", main.fontfamily = "sans", ...) {
  plot_venn(
    data.matrix(x),
    euler = euler,
    scaled = scaled,
    filename = filename,
    fontfamily = fontfamily,
    cat.fontfamily = cat.fontfamily,
    main.fontfamily = main.fontfamily,
    ...
  )
}

#' @rdname plot_venn
#' @export
plot_venn.list <- function(x, euler = FALSE, scaled = FALSE, filename = NULL, fontfamily = "sans", cat.fontfamily = "sans", main.fontfamily = "sans", ...) {
  flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  grid::grid.newpage()
  grid::grid.draw(
    VennDiagram::venn.diagram(
      x,
      euler = euler,
      scaled = scaled,
      filename = filename,
      fontfamily = fontfamily,
      cat.fontfamily = cat.fontfamily,
      main.fontfamily = main.fontfamily,
      ...
    )
  )
}
