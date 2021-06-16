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
