#' Computes the intersect between a vector of genes and interactions
#'
#' This function computes the intersect between a given vector of genes and a GenomicInteractions objects using an annotaion file
#'
#' @param interactions a GenomicInteractions object from \code{\link{load_interactions}}
#' @param genes character vector containing genes
#' @param annotation full path to annotation file which MUST include header as follows: chr, start, end, and at least a colum with genes_annot
#' @param genes_annot column name of annation file with the annotation in the same format as the genes
#'
#' @return GenomicInteractions object, subset of the original based on the give genes
#'
#' @importFrom stringr str_split
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom data.table fread
#'
#' @export
interactionsByGenes <- function(interactions, genes, annotation, genes_annot)
{
  message("Remember: Annotation must be, at least, chr start end genes_annot\n MUST INCLUDE HEADER")
  annot <- data.table::fread(annotation, header=T)
  if(any(colnames(annot) %in% "V1"))
  {
    stop("Annotation file MUST include header")
  }
  else
  {
    mat <- sapply(stringr::str_split(as.character(annot[[genes_annot]]),","), function(x) any(x %in% genes))
    matGR <- GenomicRanges::makeGRangesFromDataFrame(annot[mat,], seqnames.field = "chr",start.field = "start",end.field = "end")
    message(paste("There are",length(matGR),"HindIII fragments with genes"))
    gi <- IRanges::subsetByOverlaps(interactions,matGR)
    if (length(gi) == 0)
    {
      message("There is no interaction with the given gene(s)")
    }
    else
    {
      message(paste("There are",length(gi), "interaction with the given gene(s)"))
    }
    return(gi)
  }
}
