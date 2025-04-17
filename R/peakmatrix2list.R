#' Creates a list of HiCaptuRe objects based on a peakmatrix
#'
#' This function creates a list of HiCaptuRe objects separating all the samples included on a peakmatrix.
#'
#' @param peakmatrix HiCaptuRe object of a peakmatrix
#' @param cutoff Chicago score cut-off to filter interactions
#'
#' @return A named list of \code{HiCaptuRe} objects, each corresponding to a sample (inferred from CS_* columns in the peakmatrix).
#'
#' @importFrom GenomicRanges mcols
#'
#' @export
peakmatrix2list <- function(peakmatrix, cutoff = 5) {
  parameters <- getParameters(peakmatrix)
  format <-  parameters$load["format"]
  if (tolower(format) != "peakmatrix") {
    stop(paste0("Input is not a peakmatrix, it is a ",format))
  } else if (!is.null(parameters$peakmatrix2list)){
    stop("The input object already appears to have been processed by 'peakmatrix2list()', which means it contains interactions for a single sample only. Please provide an unprocessed 'peakmatrix' object.")
  }

  m <- GenomicRanges::mcols(peakmatrix)
  CS <- grep("CS_", names(m))
  if (is.na(CS))
  {
    stop("No CS_* columns found in the provided HiCaptuRe object.")
  }

  int_list <- lapply(CS, function(x) {
    sub_int <- peakmatrix[m[, x] >= cutoff, which(!1:ncol(m) %in% CS[CS != x])]
    names(GenomicRanges::mcols(sub_int))[6] <- "CS"
    parameters$peakmatrix2list <- c(split_by=names(m)[x])
    setParameters(sub_int) <- parameters
    sub_int
  })

  names(int_list) <- gsub("CS_", "", names(m)[CS])
  return(int_list)
}
