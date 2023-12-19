#' Creates a list of HiCaptuRe objects based on a peakmatrix
#'
#' This function creates a list of HiCaptuRe objects separating all the samples included on a peakmatrix.
#'
#' @param peakmatrix HiCaptuRe object of a peakmatrix
#' @param cutoff Chicago score cutoff to export interactions
#'
#' @return list with HiCaptuRe objects
#'
#' @importFrom GenomicRanges mcols
#'
#' @export
peakmatrix2list <- function(peakmatrix,cutoff = 5)
{
  type <- interactions@parameters$load["type"]
  if (type != "peakmatrix")
  {
    stop("Input is not a peakmatrix")
  }

  m <- GenomicRanges::mcols(peakmatrix)
  initial <- grep("CS_",names(m))[1]-1
  CS <- grep("CS_",names(m))
  final <- which(names(m) %in% c("counts","int","distance"))

  int_list <- lapply(CS, function(x)
  {
    sub_int <- peakmatrix[m[,x] >= cutoff, c(1:initial,x,final)]
    names(GenomicRanges::mcols(sub_int))[initial+1] <- "CS"
    sub_int
  })

  names(int_list) <- gsub("CS_","",names(m)[CS])
  return(int_list)
}

