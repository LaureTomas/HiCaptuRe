#' Export interactions in the desired output format
#'
#' This function exports interactions in different formats
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param file full path to desired output file (ibed, peakmatrix, washU, washUold, cytoscape, bedpe)
#' @param format type of output format (ibed, peakmatrix, washU, washUold, cytoscape, bedpe)
#' @param over.write T/F to over write the output file
#' @param washU_seqname string to add to the seqnames to export to washU format
#' @param cutoff Chicago score cutoff to export interactions
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
export_interactions <- function(interactions, file, format = "ibed", over.write=F,washU_seqname="chr",cutoff=5)
{
  if (file.exists(file) & !over.write)
  {
    user_input <- readline(paste(basename(file),"already exists. Do you want to overwrite it? (y/n)   "))
    if (user_input != "y")
    {
      stop()
    }
  }

  export_parameters(interactions,file)
  type <- interactions@parameters$load["type"]

  if (format == "peakmatrix" &  type == "peakmatrix")
  {
    m <- GenomicRanges::mcols(interactions)
    CS_m <- m[,grep("CS_",names(m))]
    interactions <- interactions[apply(CS_m, 1, function(x) any(x >= cutoff))]

    int_df <- dplyr::as_tibble(interactions)

    int_df <- int_df[,c("seqnames1","start1","end1","ID_1","bait_1",
                        "seqnames2","start2","end2","ID_2","bait_2",
                        "distance",grep("CS_",colnames(int_df),value = T))]
    colnames(int_df) <- c("baitChr", "baitStart", "baitEnd","baitID", "baitName",
                          "oeChr", "oeStart", "oeEnd","oeID", "oeName",
                          "dist", gsub("CS_","",colnames(int_df)[12:ncol(int_df)]))
    data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")

  }
  if (format != "peakmatrix" &  type == "peakmatrix")
  {
    message("Interactions is a peakmatrix. Exporting in individual files")

    int_list <- peakmatrix2list(peakmatrix = interactions,cutoff = cutoff)
    files <- paste0("~/Downloads/prueba_ibed_",names(int_list),".ibed")

    if (format == "ibed")
    {
      invisible(mapply(HiCaptuRe:::export_ibed,int_list,files))
    }
    if (format == "washU")
    {
      invisible(mapply(HiCaptuRe:::export_washU,int_list,files,washU_seqname))
    }
    if (format == "washUold")
    {
      invisible(mapply(HiCaptuRe:::export_washUold,int_list,files,washU_seqname))
    }
    if (format == "cytoscape")
    {
      invisible(mapply(HiCaptuRe:::export_citoscape,int_list,files))
    }
    if (format == "bedpe")
    {
      invisible(mapply(HiCaptuRe:::export_bedpe,int_list,files))
    }
    if (format == "seqmonk")
    {
      invisible(mapply(HiCaptuRe:::export_seqmonk,int_list,files))
    }

  }else
  {
    interactions <- interactions[interactions$CS >= cutoff]
    if (format == "ibed")
    {
      export_ibed(interactions,file)
    }
    if (format == "washU")
    {
      export_washU(interactions,file,washU_seqname)
    }
    if (format == "washUold")
    {
      export_washUold(interactions,file,washU_seqname)
    }
    if (format == "cytoscape")
    {
      export_citoscape(interactions,file)
    }
    if (format == "bedpe")
    {
      export_bedpe(interactions,file)
    }
    if (format == "seqmonk")
    {
      export_seqmonk(interactions,file)
    }
  }
}

export_parameters <- function(interactions,file)
{
  file <- paste0(file,".parameters")
  param <- getParameters(interactions)

  for (i in 1:length(param))
  {
    df <- as.data.frame(param[[i]])
    colnames(df) <- paste("#",toupper(names(param)[i]),"#")
    df <- rbind(df,"")
    rownames(df)[nrow(df)] <- ""
    if (i == 1)
    {
      data.table::fwrite(df,file = file,quote = F,sep = "\t",col.names = T)
    }else{
      suppressWarnings(data.table::fwrite(df,file = file,quote = F,sep = "\t",col.names = T,append = T))
    }
  }
}


export_ibed <- function(ints,file)
{
  int_df <- dplyr::as_tibble(ints)[,c("seqnames1","start1","end1","bait_1",
                                              "seqnames2","start2","end2","bait_2",
                                              "reads","CS")]
  colnames(int_df) <- c("bait_chr", "bait_start", "bait_end", "bait_name",
                        "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
                        "N_reads", "score")
  data.table::fwrite(int_df,file = file, col.names = T, row.names = F, quote = F, sep = "\t")
}

export_washU <- function(ints,file,washU_seqname = "chr")
{

  int_df <- dplyr::as_tibble(ints)
  int_df <- dplyr::arrange(int_df, seqnames1, start1, end1, seqnames2, start2, end2)
  washU <- data.frame(paste(paste(washU_seqname,int_df$seqnames1,sep = ""),int_df$start1,int_df$end1,sep = "\t"),
                      paste(paste(washU_seqname,int_df$seqnames2,sep = ""),":",int_df$start2,"-",int_df$end2,",",int_df$CS, sep = ""))
  colnames(washU) <- c("regionI","regionIICS")
  data.table::fwrite(washU, file = file, col.names = F, row.names = F, quote = F, sep = "\t")

  message(paste0("washU track file: ",file,"\n"))
  message(paste0('browseURL("http://epigenomegateway.wustl.edu/browser/")'))
  message("Click on Tracks > Local Text Tracks > Choose long-range text > Load washU file")
}

export_washUold <- function(ints,file,washU_seqname = "chr")
{
  int_df <- as.data.frame(ints)
  washU <- data.frame(paste(paste(washU_seqname,int_df$seqnames1,sep = ""),":",int_df$start1,",",int_df$end1,sep = ""),
                      paste(paste(washU_seqname,int_df$seqnames2,sep = ""),":",int_df$start2,",",int_df$end2, sep = ""),
                      int_df$CS)
  colnames(washU) <- c("regionI","regionII","CS")
  data.table::fwrite(washU, file = file, col.names = F, row.names = F, quote = F, sep = "\t")

  message(paste0("washU track file: ",file,"\n"))
  message(paste0('browseURL("http://epigenomegateway.wustl.edu/legacy/")'))
  message("Click on Apps > File Upload > Choose File > Setup > Pairwise Interaction > Add as Track")
}

export_bedpe <- function(ints,file)
{
  ints$name <- 1:length(ints)
  int_df <- dplyr::as_tibble(ints)[,c("seqnames1","start1","end1",
                                              "seqnames2","start2","end2",
                                              "name","CS","strand1","strand2")]

  data.table::fwrite(int_df,file = file, col.names = F, row.names = F, quote = F, sep = "\t")
}

export_citoscape <- function(ints,file)
{
  net <- igraph::simplify(export.igraph(ints))
  nodes_edges <- igraph::as_edgelist(net)
  data.table::fwrite(nodes_edges,file = file, col.names = F, row.names = F, quote = F, sep = "\t")
}

export_seqmonk <- function(ints,file)
{

  int_df$ID <- 1:nrow(int_df)
  df1 <- int_df[,c("seqnames2","start2","end2","bait_2",
                   "reads","CS","ID")]

  df2 <- int_df[,c("seqnames1","start1","end1","bait_1",
                   "reads","CS","ID")]

  seqmonk <- dplyr::bind_rows(df1, stats::setNames(df2, names(df1))) %>% dplyr::arrange(ID)
  data.table::fwrite(seqmonk[,-ncol(seqmonk)],file = file, col.names = F, row.names = F, quote = F, sep = "\t")
}
