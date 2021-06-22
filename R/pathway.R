#' get_pathway_genes
#'
#' @param pathway
#'
#' @export
get_pathway_genes <- function(pathway) {
  genes <- KEGGREST::keggGet(pathway)[[1]][["GENE"]]
  sel <- grepl(";", genes)
  genes <- data.frame(entrezgene = genes[!sel], description = genes[sel])
  genes |>
    mutate(symbol = sub(";.*", "", description)) |>
    mutate(description = sub(".*; ", "", description)) |>
    select(entrezgene, symbol, description)
}
