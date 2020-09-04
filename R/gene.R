#' Basic gene ploting function
#'
#' Plots a gene from a specified assembly using Gviz package. Optionaly plots some data (e.g. read alignments). Enables zooming.
#'
#' @param symbol symbol of the gene to plot. Will be used to search using BiomartGeneRegionTrack.
#' @param genome genome assembly (e.g. mm10 or hg38).
#' @param tracks additional tracks to add.
#' @param add.data list with data to be added as a DataTrack.
#' @param add.ann annotation information.
#' @param add.ideogram whether to add an IdeogramTrack (FALSE by default to spead up testing).
#' @param from start genomic coordinates.
#' @param to end genomic coordinates.
#' @param biomart mart object (optional).
#' @param ylim sets the limits for y-axis (default: NULL).
#' @param ... further arguments passed to plotTracks.
#'
#' @return nothing but produces a plot as side effect.
#' @export
#'
#' @examples
#' NULL
plot_gene <- function(symbol, genome, tracks = NULL, add.data = NULL, add.ann = NULL, add.ideogram = TRUE, from = NULL, to = NULL, biomart = NULL, ylim = NULL, ...) {
  itrack <- NULL
  atrack <- GenomeAxisTrack(add35 = TRUE, add53 = TRUE)
  if (!is.null(biomart)) {
    bmtrack <- BiomartGeneRegionTrack(symbol = symbol,
                                      biomart = biomart,
                                      name = "EnsEMBL",
                                      geneSymbols = TRUE,
                                      showTranscriptID = TRUE,
                                      col.line = NULL,
                                      col = NULL,
                                      col.title = "black",
                                      background.title = "white")
  } else {
    bmtrack <- BiomartGeneRegionTrack(symbol = symbol,
                                      genome = genome,
                                      name = "EnsEMBL",
                                      geneSymbols = TRUE,
                                      showTranscriptID = TRUE,
                                      col.line = NULL,
                                      col = NULL,
                                      col.title = "black",
                                      background.title = "white")
  }

  if (add.ideogram)
    itrack <- IdeogramTrack(genome = genome, chromosome = chromosome(bmtrack))

  dtrack <- NULL
  if (!is.null(add.data)) {

    if (!is.list(add.data))
      add.data <- as.list(add.data)

    if (is.null(names(add.data)))
      names(add.data) <- paste0("DataTrack-", seq_len(length(add.data)))

    dtrack <- lapply(names(add.data), function(n) {
      x <- add.data[[n]]
      DataTrack(x,
                name = n,
                type = "histogram",
                ylim = ylim,
                col.histogram = "cornflowerblue",
                fill.histogram = "cornflowerblue",
                col.title = "black",
                col.axis = "black",
                background.title = "white")
    })


    if (!is.null(from) || !is.null(to)) {
      if (is.null(from))
        from <- min(start(bmtrack), na.rm = TRUE)

      if (is.null(to))
        to <- max(end(bmtrack), na.rm = TRUE)
    }
  }

  anntrack <- NULL
  if (!is.null(add.ann)) {
    anntrack <- get_ann_track(add.ann, symbol)
  }

  tl <- c(
    itrack,
    atrack,
    anntrack,
    tracks,
    bmtrack,
    dtrack
  )
  plotTracks(tl, chromosome = chromosome(bmtrack), from = from, to = to, ...)
  invisible(bmtrack)
}

get_ann_track <- function(x, name) {
  x <- x %>% distinct(seqnames, start, end, strand, .keep_all = TRUE) %>% filter(symbol == !!name)
  AnnotationTrack(
    GRanges(x),
    name = "csaw",
    col.line = NULL,
    col = NULL,
    fill = "black",
    col.title = "black",
    background.title = "white"
  )
}
