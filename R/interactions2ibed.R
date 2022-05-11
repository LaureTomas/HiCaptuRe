#' Generate ibed file
#'
#' This function creates data.frame and files in ibed format
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param file full path to desired output file
#' @param over.write T/F to over write the output file
#'
#' @return tibble object with the ibed table and save it in the desired output file
#'
#' @importFrom dplyr as_tibble
#' @importFrom data.table fwrite
#'
#' @export
interactions2ibed <- function(interactions, file, over.write=F)
{
  if (file.exists(file) & !over.write)
  {
    stop(paste(basename(file),"already exist, if you want to overwrite it select over.write=T"))
  }
  else
  {
    int_df <- dplyr::as_tibble(interactions)[,c("seqnames1","start1","end1","gene_I",
                                             "seqnames2","start2","end2","gene_II",
                                             "reads","CS")]
    colnames(int_df) <- c("bait_chr", "bait_start", "bait_end", "bait_name",
                          "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
                          "N_reads", "score")
    data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")
    return(int_df)
  }
}
