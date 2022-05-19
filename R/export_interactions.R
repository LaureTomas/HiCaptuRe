#' Export interactions in the desired output format
#'
#' This function exports interactions in different formats
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param file full path to desired output file (ibed, peakmatrix, washU, washUold, cytoscape, bedpe)
#' @param type type of output format
#' @param over.write T/F to over write the output file
#'
#' @return tibble object with the ibed table and save it in the desired output file
#'
#' @importFrom GenomicInteractions export.igraph
#' @importFrom dplyr as_tibble arrange bind_rows
#' @importFrom data.table fwrite
#' @importFrom igraph simplify as_edgelist
#' @importFrom stats setNames
#'
#' @export
export_interactions <- function(interactions, file, type = "ibed", over.write=F)
{
  if (file.exists(file) & !over.write)
  {
    user_input <- readline(paste(basename(file),"already exists. Do you want to overwrite it? (y/n)   "))
    if (user_input != "y")
    {
      stop()
    }
    if (type == "peakmatrix")
    {
      int_df <- dplyr::as_tibble(interactions)

      int_df <- int_df[,c("seqnames1","start1","end1","ID_I","gene_I",
                          "seqnames2","start2","end2","ID_II","gene_II",
                          "dist",grep("CS_",colnames(int_df),value = T))]
      colnames(int_df) <- c("baitChr", "baitStart", "baitEnd","baitID", "baitName",
                            "oeChr", "oeStart", "oeEnd","oeID", "oeName",
                            "dist", gsub("CS_","",colnames(int_df)[12:ncol(int_df)]))
      data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")

    }else
    {
      if (any(grepl("ID_I",colnames(interactions@elementMetadata))))
      {
        message("Interactions are peakmatrix. Exporting as peakmatrix...")
        int_df <- dplyr::as_tibble(interactions)

        int_df <- int_df[,c("seqnames1","start1","end1","ID_I","gene_I",
                            "seqnames2","start2","end2","ID_II","gene_II",
                            "dist",grep("CS_",colnames(int_df),value = T))]
        colnames(int_df) <- c("baitChr", "baitStart", "baitEnd","baitID", "baitName",
                              "oeChr", "oeStart", "oeEnd","oeID", "oeName",
                              "dist", gsub("CS_","",colnames(int_df)[12:ncol(int_df)]))
        data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")
      }else
      {
        if (type == "ibed")
        {
          int_df <- dplyr::as_tibble(interactions)[,c("seqnames1","start1","end1","gene_I",
                                                      "seqnames2","start2","end2","gene_II",
                                                      "reads","CS")]
          colnames(int_df) <- c("bait_chr", "bait_start", "bait_end", "bait_name",
                                "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
                                "N_reads", "score")
          data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")
        }
        if (type == "washU")
        {
          int_df <- dplyr::as_tibble(interactions)
          int_df <- dplyr::arrange(int_df, seqnames1, start1, end1, seqnames2, start2, end2)
          washU <- data.frame(paste(paste("chr",int_df$seqnames1,sep = ""),int_df$start1,int_df$end1,sep = "\t"),
                              paste(paste("chr",int_df$seqnames2,sep = ""),":",int_df$start2,"-",int_df$end2,",",int_df$CS, sep = ""))
          colnames(washU) <- c("regionI","regionIICS")
          data.table::fwrite(washU, file = file, col.names = F, row.names = F, quote = F, sep = "\t")

          message(paste0("washU track file: ",file,"\n"))
          message(paste0('browseURL("http://epigenomegateway.wustl.edu/browser/")'))
          message("Click on Tracks > Local Text Tracks > Choose long-range text > Load washU file")
        }
        if (type == "washUold")
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
        }
        if (type == "cytoscape")
        {
          net <- igraph::simplify(export.igraph(interactions))
          nodes_edges <- igraph::as_edgelist(net)
        }
        if (type == "bedpe")
        {
          interactions$name <- 1:length(interactions)
          int_df <- dplyr::as_tibble(interactions)[,c("seqnames1","start1","end1",
                                                      "seqnames2","start2","end2",
                                                      "name","CS","strand1","strand2")]

          data.table::fwrite(int_df,file = file, col.names = F, row.names = F, quote = F, sep = "\t")
        }
        if (type == "seqmonk")
        {
          int_df <- dplyr::as_tibble(interactions)[,c("seqnames1","start1","end1","gene_I",
                                                      "seqnames2","start2","end2","gene_II",
                                                      "reads","CS")]

          int_df$ID <- 1:nrow(int_df)
          df1 <- int_df[,c("seqnames2","start2","end2","gene_II",
                           "reads","CS","ID")]

          df2 <- int_df[,c("seqnames1","start1","end1","gene_I",
                           "reads","CS","ID")]

          seqmonk <- dplyr::bind_rows(df1, stats::setNames(df2, names(df1))) %>% dplyr::arrange(ID)
          data.table::fwrite(seqmonk[,-ncol(seqmonk)],file = file, col.names = F, row.names = F, quote = F, sep = "\t")
        }
      }
    }

  }else
  {
    if (type == "peakmatrix")
    {
      int_df <- dplyr::as_tibble(interactions)

      int_df <- int_df[,c("seqnames1","start1","end1","ID_I","gene_I",
                          "seqnames2","start2","end2","ID_II","gene_II",
                          "dist",grep("CS_",colnames(int_df),value = T))]
      colnames(int_df) <- c("baitChr", "baitStart", "baitEnd","baitID", "baitName",
                            "oeChr", "oeStart", "oeEnd","oeID", "oeName",
                            "dist", gsub("CS_","",colnames(int_df)[12:ncol(int_df)]))
      data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")
    }else
    {
      if (any(grepl("ID_I",colnames(interactions@elementMetadata))))
      {
        message("Interactions are peakmatrix. Exporting as peakmatrix...")
        int_df <- dplyr::as_tibble(interactions)

        int_df <- int_df[,c("seqnames1","start1","end1","ID_I","gene_I",
                            "seqnames2","start2","end2","ID_II","gene_II",
                            "dist",grep("CS_",colnames(int_df),value = T))]
        colnames(int_df) <- c("baitChr", "baitStart", "baitEnd","baitID", "baitName",
                              "oeChr", "oeStart", "oeEnd","oeID", "oeName",
                              "dist", gsub("CS_","",colnames(int_df)[12:ncol(int_df)]))
        data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")
      }else
      {
        if (type == "ibed")
        {
          int_df <- dplyr::as_tibble(interactions)[,c("seqnames1","start1","end1","gene_I",
                                                      "seqnames2","start2","end2","gene_II",
                                                      "reads","CS")]
          colnames(int_df) <- c("bait_chr", "bait_start", "bait_end", "bait_name",
                                "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
                                "N_reads", "score")
          data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")
        }
        if (type == "washU")
        {
          int_df <- dplyr::as_tibble(interactions)
          int_df <- dplyr::arrange(int_df, seqnames1, start1, end1, seqnames2, start2, end2)
          washU <- data.frame(paste(paste("chr",int_df$seqnames1,sep = ""),int_df$start1,int_df$end1,sep = "\t"),
                              paste(paste("chr",int_df$seqnames2,sep = ""),":",int_df$start2,"-",int_df$end2,",",int_df$CS, sep = ""))
          colnames(washU) <- c("regionI","regionIICS")
          data.table::fwrite(washU, file = file, col.names = F, row.names = F, quote = F, sep = "\t")

          message(paste0("washU track file: ",file,"\n"))
          message(paste0('browseURL("http://epigenomegateway.wustl.edu/browser/")'))
          message("Click on Tracks > Local Text Tracks > Choose long-range text > Load washU file")
        }
        if (type == "washUold")
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
        }
        if (type == "cytoscape")
        {
          net <- igraph::simplify(export.igraph(interactions))
          nodes_edges <- igraph::as_edgelist(net)
        }
        if (type == "bedpe")
        {
          interactions$name <- 1:length(interactions)
          int_df <- dplyr::as_tibble(interactions)[,c("seqnames1","start1","end1",
                                                      "seqnames2","start2","end2",
                                                      "name","CS","strand1","strand2")]

          data.table::fwrite(int_df,file = file, col.names = F, row.names = F, quote = F, sep = "\t")
        }
        if (type == "seqmonk")
        {
          int_df <- dplyr::as_tibble(interactions)[,c("seqnames1","start1","end1","gene_I",
                                                      "seqnames2","start2","end2","gene_II",
                                                      "reads","CS")]

          int_df$ID <- 1:nrow(int_df)
          df1 <- int_df[,c("seqnames2","start2","end2","gene_II",
                           "reads","CS","ID")]

          df2 <- int_df[,c("seqnames1","start1","end1","gene_I",
                           "reads","CS","ID")]

          seqmonk <- dplyr::bind_rows(df1, stats::setNames(df2, names(df1))) %>% dplyr::arrange(ID)
          data.table::fwrite(seqmonk[,-ncol(seqmonk)],file = file, col.names = F, row.names = F, quote = F, sep = "\t")
        }
      }
    }

  }


}
