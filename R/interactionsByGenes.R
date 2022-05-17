#' Filters GenomicInteractions object by baits
#'
#' his function filters a GenomicInteractions object from load_interactions by a set of given bait(s)
#'
#' @param interactions a GenomicInteractions object from \code{\link{load_interactions}}
#' @param baits character vector containing bait names
#'
#' @return GenomicInteractions object, subset of the original based on the give baits
#'
#' @importFrom stringr str_split
#'
#' @export
interactionsByBaits <- function(interactions, baits)
{
  one <- sapply(stringr::str_split(as.character(interactions$gene_I),","), function(x) any(x %in% baits))
  two <- sapply(stringr::str_split(as.character(interactions$gene_II),","), function(x) any(x %in% baits))
  interactions_baits <- interactions[one | two]

  if (length(interactions_baits) == 0)
  {
    message("There is no interaction with the given bait(s)")
    message("Could Interactions not be annotated with the same bait nomenclature?")

  }
  else
  {
    message(paste("There are",length(interactions_baits), "interaction with the given bait(s)"))
  }
  return(gi)

}
