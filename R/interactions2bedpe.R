#' Generate bedpe file
#'
#' This function creates data.frame and files in ibed format
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param file full path to desired output file
#' @param over.write T/F to over write the output file
#'
#' @return tibble object with the bedpe table and save it in the desired output file
#'
#' @importFrom dplyr as_tibble
#' @importFrom data.table fwrite
#'
#' @export
interactions2bedpe <- function(interactions, file, over.write=F)
{
  if (file.exists(file) & !over.write)
  {
    stop(paste(basename(file),"already exist, if you want to overwrite it select over.write=T"))
  }
  else
  {
    interactions$name <- 1:length(interactions)
    int_df <- dplyr::as_tibble(interactions)[,c("seqnames1","start1","end1",
                                                "seqnames2","start2","end2",
                                                "name","CS","strand1","strand2")]

    data.table::fwrite(int_df,file = file, col.names = F, row.names = F, quote = F, sep = "\t")
    return(int_df)
  }
}
