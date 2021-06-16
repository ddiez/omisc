#' @export
as_immarchset <- function(x, chain = NULL, meta = NULL) {
  data <- lapply(x, as_immarch, chain = chain)
  if (is.null(meta))
    meta <- data.frame(Sample = names(x))

  list(
    data = data,
    meta = meta
  )
}

#' @export
as_immarch <- function(x, chain = NULL) {
  if (is.null(chain)) {
    stop("specify one of [", paste(unique(x$chain), collapse = " "), "]")
  }
  x %>% filter(chain == !!chain) %>%
    count(cdr3, cdr3_seq, v_gene, d_gene, j_gene, sort = TRUE) %>%
    add_tally(n, name = "total") %>%
    mutate(proportion = n / total) %>%
    select(Clones = n, Proportion = proportion, CDR3.nt = cdr3_seq, CDR3.aa = cdr3, V.name = v_gene, D.name = d_gene, J.name = j_gene)
}

#' @export
get_pathway_genes <- function(pathway) {
  require(KEGGREST)

  pathway <- sub(".*:", "", pathway)
  d <- KEGGREST::keggGet(pathway)
  genes <- d[[1]]$GENE
  genes <- grep(";", genes, value = TRUE)
  genes <- sub(";.*", "", genes)
  genes <- sort(unique(genes))
  genes
}

#' @export
plot_enrichment <- function(x, n = 10, cutoff = 0.05, ontology = "BP", title = NULL) {
  if (colnames(x)[1] == "Term") {
    x <- x %>% filter(Ont == !!ontology)
  }

  top.up <- x %>% arrange(P.Up) %>% head(n)
  top.down <- x %>% arrange(P.Down) %>% head(n)

  if (colnames(x)[1] == "Pathway") {
    d <- bind_rows(top.up, top.down) %>%
      select(term = Pathway, up = P.Up, down = P.Down)
  } else {
    d <- bind_rows(top.up, top.down) %>%
      select(term = Term, up = P.Up, down = P.Down)
  }

  d <- d %>%
    gather(direction, p.value, up, down) %>%
    mutate(score = ifelse(direction == "up", -1 * log10(p.value), log10(p.value))) %>%
    mutate(term = fct_reorder(term, score))

  ggplot(d, aes(term, score, fill = direction)) +
    geom_hline(yintercept = 0, lty = "dotted") +
    geom_col() +
    geom_hline(yintercept = c(-log10(cutoff), log10(cutoff)), lty = "dotted") +
    scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
    coord_flip() +
    labs(x = NULL, y = NULL, title = title)
}

# #' @export
# plot_gene <- function(x, name = NULL) {
#   d <- cpm(x, log = TRUE) %>% as_tibble(rownames = "entrezgene")
#   d <- d %>% gather(samplename, logcpm, -entrezgene)
#   d <- d %>% left_join(x$genes, by = "entrezgene") %>% left_join(x$samples, by = "samplename")
#
#   d <- d %>% filter(symbol %in% !!name)
#   ggplot(d, aes(background, logcpm, color = background)) +
#     geom_boxplot() +
#     geom_point() +
#     facet_wrap(~region) +
#     labs(y = paste("logCPM (", name, ")"))
# }

# plot_heatmap <- function(x, features = NULL, cluster_rows = TRUE, ...) {
#   m <- cpm(x, log = TRUE)
#
#   sel.ok <- x$genes$symbol %in% features
#   m <- m[sel.ok, , drop = FALSE]
#
#   m <- m[! is.na(x[rownames(m), ]$genes$symbol), , drop = FALSE]
#
#   rownames(m) <- x[rownames(m), ]$genes$symbol
#   m <- m[features, , drop = FALSE]
#   m <- t(scale(t(m)))
#
#
#   background <- HeatmapAnnotation(background = x$samples$background, col = list(background = c("WT" = "green", "Tg" = "orange", "KO" = "red")))
#   region <- HeatmapAnnotation(region = x$samples$region, col = list(region = c("LZ" = "limegreen", "DZ" = "steelblue")))
#
#   ComplexHeatmap::Heatmap(m, name = "logCPM",  top_annotation = c(background, region), cluster_column_slices = FALSE, ...)
# }

#' @export
plot_volcano <- function(x, coef = 1, top_genes = NULL, lfc = 1, fdr = 0.01, use.column = "symbol") {
  d <- limma::topTable(x, coef, number = Inf)

  p <- ggplot(d, aes(logFC, -log10(P.Value))) +
    geom_point(size = .1) +
    geom_hline(yintercept = -log10(1e-3), lty = "dotted") +
    geom_vline(xintercept = c(-1, 1), color = "black", lty = "dashed")

  if (!is.null(top_genes)) {
    top.up <- d %>% filter(logFC >= lfc, adj.P.Val < fdr) %>% head(top_genes)
    p <- p + ggrepel::geom_text_repel(aes(label = .data[[use.column]]), color = "red", data = top.up, min.segment.length = 0, max.overlaps = Inf)
    top.down <- d %>% filter(logFC <= -lfc, adj.P.Val < fdr) %>% head(top_genes)
    p <- p + ggrepel::geom_text_repel(aes(label = .data[[use.column]]), color = "blue", data = top.down, min.segment.length = 0, max.overlaps = Inf)
  }
  p
}


# plot_volcano <- function(x, coef = NULL, cutoff = 0.05, logfc = 1) {
#   d <- topTable(x, coef = coef, n = Inf)
#
#   ggplot(d, aes(logFC, -log10(P.Value), color = AveExpr)) +
#     geom_point(size = .1) +
#     geom_vline(xintercept = c(-1, 1), lty = "dotted") +
#     geom_hline(yintercept = -log10(1e-3), lty = "dotted") +
#     scale_color_viridis_c() +
#     geom_point(pch = 21, color = "red", data = d %>% filter(adj.P.Val < cutoff, abs(logFC) > 1)) +
#     labs(title = coef)
# }

#' @export
plot_ma <- function(x, coef = NULL, cutoff = 0.05, logfc = 1) {
  d <- topTable(x, coef = coef, n = Inf)

  ggplot(d, aes(AveExpr, logFC)) +
    geom_point(size = .1) +
    geom_hline(yintercept = c(-1, 1), lty = "dotted") +
    geom_hline(yintercept = 0, color = "limegreen") +
    scale_color_viridis_c() +
    geom_point(pch = 21, color = "red", data = d %>% filter(adj.P.Val < cutoff, abs(logFC) > 1)) +
    labs(title = coef)
}

#' @export
plot_result <- function(x) {
  UseMethod("plot_result")
}

#' @export
plot_result.TestResults <- function(x) {
  plot_result(unclass(x))
}

#' @export
plot_result.matrix <- function(x) {
  ord <- do.call(order, as.list(as.data.frame(x)))
  x <- x[ord, ]
  y <- to_tidy(x)
  y <- y %>% mutate(row = factor(row, rownames(x))) %>%
    mutate(value = factor(value, c("-1", "0", "1")))
  ggplot(y, aes(col, row, fill = value)) +
    geom_tile() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = c(`-1` = "blue3", `0` = "white", `1` = "red3")) +
    labs(x = "", y = "")
}

#' @export
compute_mds <- function(x, ...) {
  UseMethod("compute_mds")
}

#' @export
compute_mds.DGEList <- function(x, ...) {
  compute_mds(t(cpm(x, ...)))
}

#' @export
compute_mds.matrix <- function(x) {
  x <- cmdscale(dist(x))
  colnames(x) <- c("MDS_1", "MDS_2")
  x
}

#' @export
compute_pca<- function(x, ...) {
  UseMethod("compute_pca")
}

#' @export
compute_pca.DGEList <- function(x, ...) {
  compute_pca(t(cpm(x, ...)))
}

#' @export
compute_pca.matrix <- function(x) {
  prcomp(dist(x), scale = TRUE)
}

#' to_tidy
#'
#' Convert a matrix, data.frame or tibble into a tidy tibble.
#'
#' @param x a matrix, data.frame or tibble object.
#' @param row.name name for row data.
#' @param col.name name for column data.
#' @param value.name name use for value column.
#' @param stringsAsFactors logical; whether to convert col/row names to factors (preserving original ordering).
#' @param ... arguments passed to methods.
#' @export
to_tidy <- function(x, ...) {
  UseMethod("to_tidy")
}

#' @rdname to_tidy
#' @export
to_tidy.matrix <- function(x, row.name = "row", ...) {
  if (is.null(rownames(x)))
    x <- as_tibble(x) %>% rownames_to_column(var = row.name)
  else
    x <- as_tibble(x, rownames = row.name)
  to_tidy(x, row.name = row.name, ...)
}

#' @rdname to_tidy
#' @export
to_tidy.data.frame <- function(x, row.name = "row", ...) {
  y <- as_tibble(x) %>% rownames_to_column(var = row.name)
  to_tidy(y, row.name = row.name, ...)
}


#' @rdname to_tidy
#' @export
to_tidy.tbl_df <- function(x, row.name = "row", col.name = "col", value.name = "value", stringsAsFactors = FALSE, ...) {
  y <- x %>% gather(!!col.name, !!value.name, -!!row.name, ...)

  if (stringsAsFactors) {
    rn <- x[[row.name]]
    cn <- colnames(x)[-1]

    y <- y %>%
      mutate(!!row.name := factor(.data[[row.name]], levels = rn)) %>%
      mutate(!!col.name := factor(.data[[col.name]], levels = cn))
  }

  y
}

