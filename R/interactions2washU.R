#' Generate file to upload to WashU Browser
#'
#' This function creates data.frame and files to load as track in the WashU Browser
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param file full path to desired output file
#'
#' @return tibble object with the necessary table for the WashU Browser, and save it in the desired output file
#'
#' @importFrom dplyr as_tibble arrange
#' @importFrom Rsamtools bgzip indexTabix
#' @importFrom data.table fwrite
#'
#' @export
interactions2washUold <- function(interactions,file)
{
  int_df <- as.data.frame(interactions)
  washU <- data.frame(paste(paste("chr",int_df$seqnames1,sep = ""),":",int_df$start1,",",int_df$end1,sep = ""),
                      paste(paste("chr",int_df$seqnames2,sep = ""),":",int_df$start2,",",int_df$end2, sep = ""),
                      int_df$CS)
  colnames(washU) <- c("regionI","regionII","CS")
  data.table::fwrite(washU, file = file, col.names = F, row.names = F, quote = F, sep = "\t")

  message(paste0("washU track file: ",file,"\n"))
  message(paste0('browseURL("http://epigenomegateway.wustl.edu/legacy/")'))
  message("Click on Apps > File Upload > Choose File > Setup > Pairwise Interaction > Add as Track")

  return(dplyr::as_tibble(washU))
}
#' @export
interactions2washU <- function(interactions, file)
{
  int_df <- dplyr::as_tibble(interactions)
  int_df <- dplyr::arrange(int_df, seqnames1, start1, end1, seqnames2, start2, end2)
  washU <- data.frame(paste(paste("chr",int_df$seqnames1,sep = ""),int_df$start1,int_df$end1,sep = "\t"),
                      paste(paste("chr",int_df$seqnames2,sep = ""),":",int_df$start2,"-",int_df$end2,",",int_df$CS, sep = ""))
  colnames(washU) <- c("regionI","regionIICS")
  data.table::fwrite(washU, file = file, col.names = F, row.names = F, quote = F, sep = "\t")

  # a <- Rsamtools::bgzip(file)
  # b <- Rsamtools::indexTabix(a, format = "bed")
  # file.remove(file)
  # message(paste0("washU track file: ",a,"\nindex file: ",b,"\n"))
  message(paste0("washU track file: ",file,"\n"))
  message(paste0('browseURL("http://epigenomegateway.wustl.edu/browser/")'))
  message("Click on Tracks > Local Text Tracks > Choose long-range text > Load washU file")

  return()
}

