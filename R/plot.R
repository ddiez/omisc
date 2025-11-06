#' plot_enrichr
#' 
#' Plots the results of enrichR with some different choice of colors.
#'
#' @param x enrichR result.
plot_enrichr <- function(x, db) {
  x$Term <- factor(x$Term, rev(x$Term))
  x$Count <- as.integer(sub("\\/.*", "", x$Overlap))
  ggplot(x, aes(.data[["Count"]], .data[["Term"]], fill=-log10(.data[["Adjusted.P.value"]]))) +
    geom_col() + 
    scale_fill_distiller(palette="Spectral") +
    labs(y=NULL)
}

#' plot_heatmap
#'
#' @param x object to plot.
#' @param log whether to use log CPMs.
#' @param batch batch variable for removeBatchEffect.
#' @param design design matrix for removeBatchEffect.
#' @param symbol.col column name with gene symbols.
#' @param keep.ids whether to keep original feature ids.
#' @param top_ann names of columns to be used as top annotations.
#' @param top_ann_col color definition for the categories in the top annotations.
#' @param show_column_names whether to show column names (default: TRUE).
#' @param show_ann_legend whether to show the annotation legends.
#' @param ...
#'
#' @export
plot_heatmap <- function(x, ...) {
  UseMethod("plot_heatmap")
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.DGEList <- function(x, log = TRUE, batch = NULL, design = NULL, symbol.col = "SYMBOL", keep.ids = FALSE, top_ann = NULL, top_ann_col = NULL, show_ann_legend=TRUE, ...) {

  if (!is.null(top_ann)) {
    df <- x$sample |> select(top_ann)
    top_ann <- ComplexHeatmap::columnAnnotation(df = df, col = top_ann_col, show_legend=show_ann_legend)
  }


  m <- edgeR::cpm(x, log = log)

  if (!is.null(batch)) {
    if (is.null(design))
      design <- matrix(1, ncol(m),1)
    m <- removeBatchEffect(m, batch = batch, design = design)
  }

  symbol <- x$genes[[symbol.col]]

  if (keep.ids) {
    symbol <- paste0(rownames(m), ":", symbol)
  } else {
    sel.na <- is.na(symbol)
    symbol[sel.na] <- rownames(m)[sel.na]
  }

  rownames(m) <- symbol

  plot_heatmap(m, top_ann = top_ann, ...)
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.matrix <- function(x, scale = TRUE, show_column_names = TRUE, ...) {
  if (scale)
    x <- t(scale(t(x)))

  ComplexHeatmap::Heatmap(x, name = "logCPM", show_column_names = show_column_names, ...)
}


#' plot_hist
#'
#' Plots histogram from various types of objects.
#'
#' @param x an R object.
#' @param coef name of coefficient.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_hist <- function(x, ...) {
  UseMethod("plot_hist")
}

#' @rdname plot_hist
#' @export
plot_hist.MArrayLM <- function(x, coef = NULL, remove_intercept=TRUE, ...) {
  d <- to_tidy(x$p.value, "gene", "group", "p.value")

  if (remove_intercept) {
    d <- d |> filter(! grepl("intercept", .data[["group"]], ignore.case=TRUE))
  }

  if (!is.null(coef)) {
    d <- d %>% filter(.data[["group"]] %in% coef)
  }

  #groups <- unique(d[["group"]])
  ggplot(d, aes(.data[["p.value"]])) +
    geom_histogram(binwidth=0.01) +
    facet_wrap("group", ...)
  #lapply(groups, function(g) {
  #  ggplot(d |> filter(.data[["group"]] == g), aes(.data[["p.value"]])) +
  #    geom_histogram(binwidth = .01) + labs(title=g)
  #}) |> patchwork::wrap_plots()
}


#' plot_ma
#'
#' MA plot from various types of objects.
#'
#' @param x an R object.
#' @param coef name of the coefficient.
#' @param color.by column to color points.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_ma <- function(x, ...) {
  UseMethod("plot_ma")
}

#' @rdname plot_ma
#' @export
plot_ma.MArrayLM <- function(x, coef = 1, color.by = NULL) {
  d <- limma::topTable(x, coef = coef, n = Inf)

  if (!is.null(color.by)) {
    d <- d %>% arrange(.data[[color.by]])
  }

  p <- ggplot(d, aes(AveExpr, logFC))

  if (!is.null(color.by)) {
    p <- p + geom_point(aes(color = .data[[color.by]]))
  } else {
    p <- p + geom_point(color = "lightgrey")
  }

  p + geom_hline(yintercept = c(-1, 1), lty = "dotted") +
    geom_hline(yintercept = 0, color = "limegreen") +
    labs(title = coef)
}


#' plot_volcano
#'
#' Volcano plot from various types of objects.
#'
#' @param x an R object.
#' @param coef name of the coefficient.
#' @param top_genes number of genes to highlight.
#' @param lfc logFC cutoff for top genes.
#' @param fdr FDR cutoff for top genes.
#' @param use.column column containing gene names (default: symbol).
#' @param color.by column to color points.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_volcano <- function(x, ...) {
  UseMethod("plot_volcano")
}

#' @rdname plot_volcano
#' @export
plot_volcano.MArrayLM <- function(x, coef = 1, top_genes = NULL, lfc = 1, fdr = 0.01, use.column = "symbol", color.by = NULL) {
  d <- limma::topTable(x, coef, number = Inf)

  if (!is.null(color.by)) {
    d <- d %>% arrange(.data[[color.by]])
  }

  p <- ggplot(d, aes(logFC, -log10(P.Value)))

  if (!is.null(color.by)) {
    p <- p + geom_point(aes(color = .data[[color.by]]))
  } else {
    p <- p + geom_point(color = "black", size = .1)
  }

  p <- p + geom_hline(yintercept = -log10(1e-3), lty = "dotted") +
    geom_vline(xintercept = c(-1, 1), color = "black", lty = "dashed")

  if (!is.null(top_genes)) {
    top.up <- d %>% filter(logFC >= lfc, adj.P.Val < fdr) %>% head(top_genes)
    p <- p + ggrepel::geom_text_repel(aes(label = .data[[use.column]]), color = "red", data = top.up, min.segment.length = 0, max.overlaps = Inf)
    top.down <- d %>% filter(logFC <= -lfc, adj.P.Val < fdr) %>% head(top_genes)
    p <- p + ggrepel::geom_text_repel(aes(label = .data[[use.column]]), color = "blue", data = top.down, min.segment.length = 0, max.overlaps = Inf)
  }
  p + labs(title = coef)
}

#' plot_correlation
#'
#' Correlation heatmap from samples.
#'
#' @param x an R object.
#' @param log whether to use log CPMs.
#' @param ... arguments passed down to methods.
#'
#' @export
plot_correlation <- function(x, ...) {
  UseMethod("plot_correlation")
}

#' @rdname plot_volcano
#' @export
plot_correlation.DGEList <- function(x, log=TRUE, ...) {
  col_fun <- circlize::colorRamp2(c(-1,0,1), c("blue", "white", "red"))
  ComplexHeatmap::Heatmap(cor(cpm(x, log=log)), name="correlation", col=col_fun, ...)
}
