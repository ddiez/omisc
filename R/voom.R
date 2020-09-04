#' Plots the mean-variance trend.
#'
#' This function uses the information stored in the object returned by voom()
#' when called with save.plot = TRUE.
#'
#' @param x Output of voom() function.
#'
#' @return NULL
#' @export
#'
plot_voom <- function(x) {
  UseMethod("plot_voom")
}

#' @rdname plot_voom
#' @export
plot_voom.EList <- function(x) {
  d <- data.frame(x = x$voom.xy$x, y = x$voom.xy$y)
  l <- x$voom.line

  ggplot(d, aes(x, y)) +
    geom_point(size = .1) +
    annotate("line", x = l$x, y = l$y, color = "red", size = 1) +
    labs(x = d$xlab, y = d$ylab, title = "voom: Mean variance trend")
}
