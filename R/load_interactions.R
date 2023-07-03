#' Loads interaction file into GenomicInteractions Object
#'
#' This function loads interaction files from Chicago R package into a GenomicInteractions Object, and remove possible duplicated interactions
#'
#' @param file full path to the interaction file (seqmonk, ibed, washU)
#'
#' @return GenomicInteractions object
#'
#' @importFrom magrittr `%>%`
#' @importFrom tidyr separate
#' @importFrom dplyr as_tibble group_by filter slice n
#' @importFrom stringr str_replace_all
#' @importFrom GenomicInteractions GenomicInteractions
#' @importFrom GenomicRanges makeGRangesFromDataFrame split mcols
#' @importFrom data.table fread
#' @importFrom progressr progressor
#' @importFrom S4Vectors elementMetadata
#'
#' @export
load_interactions <- function(file)
{
  if (!file.exists(file))
  {
    stop(paste(basename(file), "does not exist"))
  }else
  {
    ## Setting pipe operator from magrittr package
    `%>%` <- magrittr::`%>%`

    ## Reading file and detecting file format depending of the number of columns
    ## Tranforming all file formats to seqmonk to proceed with the cleaning

    p <- progressr::progressor(steps = 10)

    p(sprintf("Reading File"))
    data <- data.table::fread(file = file, header = T, stringsAsFactors = F, na.strings = "")

    if (ncol(data) > 10)
    {
      type <- "peakmatrix"
      message(paste(basename(file), "is a peakmatrix"))
      a <- data.table::fread(file = file, header = T, stringsAsFactors = F, na.strings = "")
      a1 <- a[, c(1:5, 11:ncol(a)),with=F]
      a2 <- a[, c(6:ncol(a)),with=F]
      colnames(a2) <- colnames(a1)
      df <- rbind(data.frame(a1, index = 1:nrow(a1)), data.frame(a2,index = 1:nrow(a2)))
      df <- df[order(df$index), ]
      data <- df[, -c(which(colnames(df)=="index"))]
      data$rownames <- 1:nrow(data)
      p(sprintf("Preparing Data"))


      ## Extracting name of all cell types in this peakmatrix
      cell_types <- paste0("CS_",colnames(data)[7:(ncol(data)-1)])

      ## Putting together in one line each interactions and duplicating them
      new_data <- rbind(cbind(data[seq(1, nrow(data), 2), ], data[seq(2, nrow(data), 2), ]),
                        cbind(data[seq(2, nrow(data), 2), ], data[seq(1, nrow(data), 2), ]))
      p(sprintf("Sorting Data"))

      ## Ordering by the original line that came
      new_data <- dplyr::as_tibble(new_data[order(as.numeric(rownames(new_data))), ], .name_repair = "minimal")
      colnames(new_data) <- c("chr_I", "start_I", "end_I","ID_I","gene_I","dist_I",  paste0("CS_I_ct",1:length(cell_types)),"rownames1",
                              "chr_II", "start_II", "end_II","ID_II","gene_II","dist_II",paste0("CS_II_ct",1:length(cell_types)),"rownames2")
    }else
    {
      type <- "non_peakmatrix"
      if (ncol(data)==6)
      {
        data <- data.table::fread(file = file, header = F, stringsAsFactors = F, na.strings = "")
        message(paste(basename(file), "is in seqmonk format"))
        data$rownames <- 1:nrow(data)
        p(sprintf("Preparing Data"))

      }
      if (ncol(data)==10)
      {
        message(paste(basename(file), "is in ibed format"))
        a <- data
        a1 <- a[,c(1:4,9:10),with=F]
        a2 <- a[,c(5:10),with=F]
        colnames(a2) <- colnames(a1)
        df <- rbind(data.frame(a1, index = 1:nrow(a1)), data.frame(a2, index = 1:nrow(a2)))
        df <- df[order(df$index),]
        data <- df[,1:6]
        data$rownames <- 1:nrow(data)
        p(sprintf("Preparing Data"))

      }
      if (ncol(data)==3)
      {
        if(!grepl(":",data[1,1]))
        {
          message(paste(basename(file), "is in washU new format"))
          warning("We do not recommend to use washU format from Chicago \n The ibed output must be annotated, see ??annotate_interactions")

          data <- data.table::fread(file = file, header = F, stringsAsFactors = F, sep = "\t") #reading file

          data <- tidyr::separate(data, 4, into=c("a","b"), sep = ":", remove = TRUE,
                                  convert = FALSE, extra = "warn", fill = "warn") %>%
            tidyr::separate(5, into=c("b","c"), sep = "-", remove = TRUE,
                            convert = FALSE, extra = "warn", fill = "warn") %>%
            tidyr::separate(6, into=c("c","d"), sep = ",", remove = TRUE,
                            convert = FALSE, extra = "warn", fill = "warn")
          p(sprintf("Preparing Data"))

        }
        if(grepl(":",data[1,1]))
        {
          message(paste(basename(file), "is in washU old format"))
          warning("We do not recommend to use washU format from Chicago \n The ibed output must be annotated, see ??annotate_interactions")

          data <- data.table::fread(file = file, header = F, stringsAsFactors = F, sep = "\t") #reading file

          data <- tidyr::separate(data, 1, into=c("a","b"), sep = ":", remove = TRUE,
                                  convert = FALSE, extra = "warn", fill = "warn") %>%
            tidyr::separate(2, into=c("b","c"), sep = ",", remove = TRUE,
                            convert = FALSE, extra = "warn", fill = "warn") %>%
            tidyr::separate(4, into=c("d","e"), sep = ":", remove = TRUE,
                            convert = FALSE, extra = "warn", fill = "warn") %>%
            tidyr::separate(5, into=c("e","f"), sep = ",", remove = TRUE,
                            convert = FALSE, extra = "warn", fill = "warn")
          p(sprintf("Preparing Data"))

        }

        data[,c(1,4)] <- lapply(data[,c(1,4)], function(x) as.character(gsub("chr", "", x)))
        data[,c(2:3,5:7)] <- as.data.frame(apply(data[,c(2:3,5:7)], 2, function(x) as.numeric(x)))

        annotations <- rep("non-annotated", nrow(data))
        reads <- rep(0, nrow(data))

        df1 <- cbind(data[,1:3],annotations,reads,data[,7])
        df2 <- cbind(data[,4:6],annotations,reads,data[,7])
        colnames(df2) <- colnames(df1)

        df <- rbind(data.frame(df1, index = 1:nrow(df1)), data.frame(df2, index = 1:nrow(df2)))
        df <- df[order(df$index),]
        data <- df[,1:6]
        data$rownames <- 1:nrow(data)
      }
      ## Putting together in one line each interactions and duplicating them
      new_data <- rbind(cbind(data[seq(1,nrow(data),2),],data[seq(2,nrow(data),2),]),
                        cbind(data[seq(2,nrow(data),2),],data[seq(1,nrow(data),2),]))
      p(sprintf("Sorting Data"))

      ## Ordering by the original line that came
      new_data <- dplyr::as_tibble(new_data[order(as.numeric(rownames(new_data))),],
                                   .name_repair = "minimal")

      colnames(new_data) <- c("chr_I","start_I","end_I","gene_I", "R_I","CS_I","rownames1",
                              "chr_II","start_II","end_II","gene_II","R_II","CS_II","rownames2")
    }

    p(sprintf("Real Duplicates"))

    ## Here, could be real duplicated interactions, check with unique
    new_data$gene_I <- stringr::str_replace_all(new_data$gene_I, "\\|",",")
    new_data$gene_II <- stringr::str_replace_all(new_data$gene_II, "\\|",",")
    ## After correcting the annotation the number of real duplicates increase

    ## Removing real duplicates, if exist, in the file
    new_data <- data.table::as.data.table(new_data)
    new_data <- unique(new_data)

    dup_real <- (nrow(data)-nrow(new_data))/2
    p(sprintf("CS Duplicates"))

    ## Filtering those duplicated interactions with different CS, by the higher one
    new_data <- new_data[, lapply(.SD, max),by=list(chr_I, start_I, end_I, chr_II, start_II, end_II)]
    new_data <- new_data[order(new_data$rownames1),]



    dup_CS <- ((nrow(data) - nrow(new_data))/2) - dup_real


    if(all.equal(new_data[,grep("CS_I(?!I).*",colnames(new_data),perl = T),with=F],new_data[,grep("CS_II.*",colnames(new_data)),with=F],check.attributes = F))
    {
      ## Keeping only one of the artificial inverted duplications
      new_data <- dplyr::slice(new_data, seq(1,dplyr::n(),2))
      p(sprintf("Cleaning"))

      ## Removing interactions that involve the MT chromosome

      message(paste0("\n",(nrow(data)/2) - nrow(new_data))," interactions removed\n\t- ",dup_real," interactions were real duplicates\n\t- ",dup_CS," interactions duplicated with different Chicago Score")

      if (type == "peakmatrix")
      {
        new_data <- new_data[, !colnames(new_data) %in% c("rownames1","rownames2","dist_I",paste0("CS_I_ct",1:length(cell_types))), with=F]
        colnames(new_data)[11:ncol(new_data)] <- c("dist", cell_types)
        new_data <- new_data[,c("chr_I","start_I","end_I","gene_I","ID_I",
                                "chr_II","start_II","end_II","gene_II","ID_II",
                                "dist",cell_types),with=F]
        new_data[new_data$gene_I == ".", ] <- c(new_data[new_data$gene_I == ".", c(6:10, 1:5, 11:ncol(new_data))])

      }else
      {
        new_data <- new_data[,!colnames(new_data) %in% c("rownames1","rownames2","R_I","CS_I"), with=F]
        colnames(new_data)[9:10] <- c("reads","CS")
        new_data <- new_data[,c("chr_I","start_I","end_I","gene_I",
                                "chr_II","start_II","end_II","gene_II",
                                "reads","CS"),with=F]
        new_data[new_data$gene_I == ".",] <- c(new_data[new_data$gene_I == ".",c(5:8,1:4,9,10)])
      }

      p(sprintf("Sorting"))



      ## Creating the genomic interactions object

      region1 <- GenomicRanges::makeGRangesFromDataFrame(new_data[,1:(grep("chr_II",colnames(new_data))-1)],seqnames.field = "chr_I", start.field = "start_I", end.field = "end_I", keep.extra.columns = T)
      region2 <- GenomicRanges::makeGRangesFromDataFrame(new_data[,grep("chr_II",colnames(new_data)):ncol(new_data)],seqnames.field = "chr_II", start.field = "start_II", end.field = "end_II", keep.extra.columns = T)

      # region1$ID_I <- subjectHits(findOverlaps(region1,hindGR))
      # region2$ID_II <- subjectHits(findOverlaps(region2,hindGR))
      #
      # gi <- GenomicInteractions::GenomicInteractions(region1, region2[,c(1,4,2,3)])
      gi <- GenomicInteractions::GenomicInteractions(region1, region2)
      p(sprintf("GenomicInteractions"))

      names(GenomicRanges::mcols(gi)) <- gsub(x = names(GenomicRanges::mcols(gi)), pattern = "anchor[1-2]\\.", "")

      ## Annotating regions with P, OE or uce

      gi <- annotate_POEuce(gi)


      p(sprintf("Type of Interactions"))

      ## Sorting interactions P_P
      cond <- (gi@anchor1 > gi@anchor2) & (!grepl("OE",gi$int))

      a1 <- gi@anchor1[cond]
      a2 <- gi@anchor2[cond]

      gi@anchor1[cond] <- a2
      gi@anchor2[cond] <- a1


      cols <- sort(grep("_I",colnames(S4Vectors::elementMetadata(gi[cond]))[1:5],value = T))
      S4Vectors::elementMetadata(gi[cond])[cols] <- S4Vectors::elementMetadata(gi[cond])[cols[c(rbind(seq(2,length(cols),2),seq(1,length(cols),2)))]]



      gi@elementMetadata <- gi@elementMetadata[,-which(colnames(gi@elementMetadata) %in% c("counts"))]

      p(sprintf("Finishing"))

      return(gi)
    }
    else(
      stop("Loading Interactions file error")
    )
  }
}

