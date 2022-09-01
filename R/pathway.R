#' get_pathway_genes
#'
#' @param pathway
#'
#' @export
get_pathway_genes <- function(pathway) {
  genes <- KEGGREST::keggGet(pathway)[[1]][["GENE"]]
  if (is.null(genes)) return(NULL)


  sel <- grepl(";", genes)
  genes <- data.frame(entrezgene = genes[!sel], description = genes[sel])
  genes <- genes |>
    mutate(symbol = sub(";.*", "", description)) |>
    mutate(description = sub(".*; ", "", description)) |>
    select(entrezgene, symbol, description)
  return(genes)
}
