#' Load genome
#'
#' This function load the given genome if installed
#'
#' @param genome character with the name of the genome version
#'
#' @return BSgenome object
#'
#'
#' @importFrom BSgenome installed.genomes
#'
#'
#' @export
load_genome <- function(genome) {
  genomes <- BSgenome::installed.genomes(splitNameParts = TRUE)$genome
  if (length(genomes)) {
    names(genomes) <- gsub("^BSgenome.", "", BSgenome::installed.genomes())
  }
  if (!length(genomes)) {
    stop("No genomes installed!")
  }
  if (any(genomes %in% genome)) {
    pkg <- paste0("BSgenome.", names(genomes[genomes %in%
      genome]))[[1]]
    suppressPackageStartupMessages(library(pkg,
      character.only = TRUE,
      quietly = TRUE
    ))
    genome_output <- get(pkg)
    return(genome_output)
  } else {
    stop(paste("Genome name doesn't match any genome installed \nGenomes installed:\n", paste(genomes, collapse = ", ")))
  }
}
