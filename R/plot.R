#' plot_heatmap
#'
#' @param x object to plot.
#' @param log whether to use log CPMs.
#' @param batch batch variable for removeBatchEffect.
#' @param design design matrix for removeBatchEffect.
#' @param symbol.col column name with gene symbols.
#' @param top_ann names of columns to be used as top annotations.
#' @param top_ann_col color definition for the categories in the top annotations.
#' @param show_column_names whether to show column names (default: TRUE).
#' @param ...
#'
#' @export
plot_heatmap <- function(x, ...) {
  UseMethod("plot_heatmap")
}

#' @rdname plot_heatmap
#' @export
plot_heatmap.DGEList <- function(x, log = TRUE, batch = NULL, design = NULL, symbol.col = "SYMBOL", top_ann = NULL, top_ann_col = NULL, ...) {

  if (!is.null(top_ann)) {
    df <- x$sample[, top_ann]
    top_ann <- ComplexHeatmap::columnAnnotation(df = df, col = top_ann_col)
  }


  m <- edgeR::cpm(x, log = log)

  if (!is.null(batch)) {
    if (is.null(design))
      design <- matrix(1, ncol(m),1)
    m <- removeBatchEffect(m, batch = batch, design = design)
  }

  ids <- x$genes[[symbol.col]]
  sel.na <- is.na(ids)
  ids[sel.na] <- rownames(x)[sel.na]
  rownames(m) <- ids

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
#' @return
#' @export
plot_hist <- function(x, ...) {
  UseMethod("plot_hist")
}

#' @rdname plot_hist
#' @export
plot_hist.MArrayLM <- function(x, coef = NULL) {
  d <- to_tidy(fit2$p.value, "gene", "group", "p.value")

  if (!is.null(coef)) {
    d <- d %>% filter(group %in% coef)
  }

  ggplot(d, aes(p.value)) +
    geom_histogram(binwidth = .01) +
    facet_wrap(~group)
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
#' @return
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
#' @return
#' @export
plot_volcano <- function(x, ...) {
  UseMethod("plot_vocano")
}

#' @rdname plot_volcano
#' @export
plot_volcano <- function(x, coef = 1, top_genes = NULL, lfc = 1, fdr = 0.01, use.column = "symbol", color.by = NULL) {
  d <- limma::topTable(x, coef, number = Inf)

  if (!is.null(color.by)) {
    d <- d %>% arrange(.data[[color.by]])
  }

  p <- ggplot(d, aes(logFC, -log10(P.Value)))

  if (!is.null(color.by)) {
    p <- p + geom_point(aes(color = .data[[color.by]]))
  } else {
    p <- p + geom_point(color = "lightgrey")
  }

  p <- p + geom_hline(yintercept = -log10(1e-3), lty = "dotted") +
    geom_vline(xintercept = c(-1, 1), color = "black", lty = "dashed")

  if (!is.null(top_genes)) {
    top.up <- d %>% filter(logFC >= lfc, adj.P.Val < fdr) %>% head(top_genes)
    p <- p + ggrepel::geom_text_repel(aes(label = .data[[use.column]]), color = "red", data = top.up, min.segment.length = 0, max.overlaps = Inf)
    top.down <- d %>% filter(logFC <= -lfc, adj.P.Val < fdr) %>% head(top_genes)
    p <- p + ggrepel::geom_text_repel(aes(label = .data[[use.column]]), color = "blue", data = top.down, min.segment.length = 0, max.overlaps = Inf)
  }
  p
}

