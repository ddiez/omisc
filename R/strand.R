#' plot_strand_qc
#'
#' Plots to facilitate identifying the strandness of the sequencing data.
#'
#' @param x a data.frame like object obtained with read_star().
#'
#' @return
#' @export
#'
plot_strand_qc <- function(x) {
  p0 <- x %>% gather(type, count, -id) %>%
    ggplot(aes(log10(count + 1))) +
    geom_histogram(binwidth = .1) +
    scale_y_log10() +
    facet_wrap(~type)

  p1 <- ggplot(x, aes(unstranded, fr_strand)) +
    geom_point(size = .1) +
    scale_x_log10() +
    scale_y_log10()

  p2 <- ggplot(x, aes(unstranded, sr_strand)) +
    geom_point(size = .1) +
    scale_x_log10() +
    scale_y_log10()

  p3 <- ggplot(x, aes(fr_strand, sr_strand)) +
    geom_point(size = .1) +
    scale_x_log10() +
    scale_y_log10()

  p0 / (p1 + p2 + p3)
}
