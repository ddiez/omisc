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
