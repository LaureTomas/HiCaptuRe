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
#' @importFrom dplyr as_tibble group_by filter slice ungroup n
#' @importFrom stringr str_replace_all
#' @importFrom GenomicInteractions GenomicInteractions annotateInteractions
#' @importFrom GenomicRanges makeGRangesFromDataFrame split mcols
#' @importFrom ddpcr quiet
#' @importFrom data.table fread
#'
#' @export
load_interactions <- function(file)
{
  if (!file.exists(file))
  {
    stop(paste(basename(file), "does not exist"))
  }
  else
  {
    ## Setting pipe operator from magrittr package
    `%>%` <- magrittr::`%>%`

    ## Reading file and detecting file format depending of the number of columns
    ## Tranforming all file formats to seqmonk to proceed with the cleaning

    progress_bar = txtProgressBar(min=1, max=11, style = 3, char="=")

    data <- data.table::fread(file = file, header = F, stringsAsFactors = F, na.strings = "")
    if (ncol(data)==6)
    {
      message(paste(basename(file), "is in seqmonk format"))
      setTxtProgressBar(progress_bar, value = 2)

    }
    if (ncol(data)==10)
    {
      message(paste(basename(file), "is in ibed format"))
      a <- data.table::fread(file = file, header = T, stringsAsFactors = F, na.strings = "")
      a1 <- a[,c(1:4,9:10)]
      a2 <- a[,c(5:10)]
      colnames(a2) <- colnames(a1)
      df <- rbind(data.frame(a1, index = 1:nrow(a1)), data.frame(a2, index = 1:nrow(a2)))
      df <- df[order(df$index),]
      data <- df[,1:6]
      rownames(data) <- 1:nrow(data)
      setTxtProgressBar(progress_bar, value = 2)

    }
    if (ncol(data)==3)
    {
      if(!grepl(":",data[1,1]))
      {
        message(paste(basename(file), "is in washU new format"))
        warning("We do not recommend to use washU format from Chicago \n The ibed output must be annotated, see annote_ibed function")

        data <- data.table::fread(file = file, header = F, stringsAsFactors = F, sep = "\t") #reading file

        data <- tidyr::separate(data, 4, into=c("a","b"), sep = ":", remove = TRUE,
                                convert = FALSE, extra = "warn", fill = "warn") %>%
          tidyr::separate(5, into=c("b","c"), sep = "-", remove = TRUE,
                          convert = FALSE, extra = "warn", fill = "warn") %>%
          tidyr::separate(6, into=c("c","d"), sep = ",", remove = TRUE,
                          convert = FALSE, extra = "warn", fill = "warn")
        setTxtProgressBar(progress_bar, value = 2)

      }
      if(grepl(":",data[1,1]))
      {
        message(paste(basename(file), "is in washU old format"))
        warning("We do not recommend to use washU format from Chicago \n The ibed output must be annotated, see annote_ibed function")

        data <- data.table::fread(file = file, header = F, stringsAsFactors = F, sep = "\t") #reading file

        data <- tidyr::separate(data, 1, into=c("a","b"), sep = ":", remove = TRUE,
                                convert = FALSE, extra = "warn", fill = "warn") %>%
          tidyr::separate(2, into=c("b","c"), sep = ",", remove = TRUE,
                          convert = FALSE, extra = "warn", fill = "warn") %>%
          tidyr::separate(4, into=c("d","e"), sep = ":", remove = TRUE,
                          convert = FALSE, extra = "warn", fill = "warn") %>%
          tidyr::separate(5, into=c("e","f"), sep = ",", remove = TRUE,
                          convert = FALSE, extra = "warn", fill = "warn")
        setTxtProgressBar(progress_bar, value = 2)

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
      rownames(data) <- 1:nrow(data)
    }

    ## Putting together in one line each interactions and duplicating them
    data <- as.data.frame(data)
    new_data <- rbind(cbind(data[seq(1,nrow(data),2),],data[seq(2,nrow(data),2),]),
                      cbind(data[seq(2,nrow(data),2),],data[seq(1,nrow(data),2),]))
    setTxtProgressBar(progress_bar, value = 3)

    ## Ordering by the original line that came
    new_data <- dplyr::as_tibble(new_data[order(as.numeric(rownames(new_data))),],
                                 .name_repair = "minimal")

    colnames(new_data) <- c("chr_I","start_I","end_I","gene_I", "R_I","CS_I",
                            "chr_II","start_II","end_II","gene_II","R_II","CS_II")
    setTxtProgressBar(progress_bar, value = 4)

    ## Here, could be real duplicated interactions, check with unique
    new_data$gene_I <- stringr::str_replace_all(new_data$gene_I, "\\|",",")
    new_data$gene_II <- stringr::str_replace_all(new_data$gene_II, "\\|",",")
    ## After correcting the annotation the number of real duplicates increase


    ## Removing real duplicates, if exist, in the file
    new_data <- unique(new_data)

    dup_real <- (nrow(data)-nrow(new_data))/2

    ## Filtering those duplicated interactions with different CS, by the higher one
    new_data <- new_data %>% dplyr::group_by(chr_I, start_I, end_I, chr_II, start_II, end_II) %>% dplyr::filter(CS_I == max(CS_I))
    setTxtProgressBar(progress_bar, value = 5)

    dup_CS <- ((nrow(data) - nrow(new_data))/2) - dup_real

    if(identical(new_data$CS_I, new_data$CS_II))
    {
      message("\nCleaning interactions file correct")

      ## Keeping only one of the artificial inverted duplications
      new_data <- dplyr::slice(dplyr::ungroup(new_data), seq(1,dplyr::n(),2))
      setTxtProgressBar(progress_bar, value = 6)

      ## Removing interactions that involve the MT chromosome
      new_data <- new_data[!(grepl("MT",new_data$chr_I)|grepl("MT",new_data$chr_II)),]

      dup_MT <- (nrow(data)/2) - nrow(new_data) - dup_real - dup_CS

      message(paste0("\n",(nrow(data)/2) - nrow(new_data))," interactions removed\n\t- ",dup_real," interactions were real duplicates\n\t- ",dup_CS," interactions duplicated with different Chicago Score\n\t- ",dup_MT," interactions involve MT Chromosome")
      new_data <- new_data[,-c(5:6)]
      colnames(new_data)[9:10] <- c("reads","CS")
      setTxtProgressBar(progress_bar, value = 7)


      ## reordering to have always captured baits as first end

      new_data[new_data$gene_I == ".",] <- c(new_data[new_data$gene_I == ".",c(5:8,1:4,9,10)])

      ## Creating the genomic interactions object

      region1 <- GenomicRanges::makeGRangesFromDataFrame(new_data[,1:4],seqnames.field = "chr_I", start.field = "start_I", end.field = "end_I", keep.extra.columns = T)
      region2 <- GenomicRanges::makeGRangesFromDataFrame(new_data[,5:10],seqnames.field = "chr_II", start.field = "start_II", end.field = "end_II", keep.extra.columns = T)

      gi <- GenomicInteractions::GenomicInteractions(region1, region2)
      setTxtProgressBar(progress_bar, value = 8)

      ## Annotating regions with P, OE or uce

      colnames(region2@elementMetadata)[1] <- "gene_I"
      regions <- unique(c(region1, region2[,1]))
      regions$gene_I[is.na(regions$gene_I)] <- "non-annotated"
      prom <- regions[(regions$gene_I != ".") & (regions$gene_I != "uce")]
      oe <- regions[regions$gene_I == "."]
      uce <- regions[regions$gene_I == "uce"]

      proml <- GenomicRanges::split(prom[,-1], as.factor(prom$gene_I))
      oel <- GenomicRanges::split(oe[,-1], as.factor(oe$gene_I))
      ucel <- GenomicRanges::split(uce[,-1], as.factor(uce$gene_I))

      annotation.features = list(P=proml,OE=oel,uce=ucel)

      ddpcr::quiet(GenomicInteractions::annotateInteractions(gi, annotation.features))

      names(GenomicRanges::mcols(gi)) <- gsub(x = names(GenomicRanges::mcols(gi)), pattern = "anchor[1-2]\\.","")
      setTxtProgressBar(progress_bar, value = 9)

      ## Defining the types of interactions
      gi$int <- ifelse((gi$gene_I == "." | gi$gene_II == ".") &
                         (gi$gene_I != "uce" & gi$gene_II != "uce"), "P_OE",
                       ifelse((gi$gene_I == "." | gi$gene_II == ".") &
                                (gi$gene_I == "uce" | gi$gene_II == "uce"), "uce_OE",
                              ifelse(!gi$gene_I %in% c(".","uce") & !gi$gene_II %in% c(".","uce"),"P_P",
                                     ifelse((gi$gene_I == "uce" | gi$gene_II == "uce") &
                                              (!gi$gene_I %in% c(".","uce") | !gi$gene_II %in% c(".","uce")),"P_uce",
                                            ifelse(gi$gene_I=="uce" & gi$gene_II=="uce","uce_uce",NA)))))
      setTxtProgressBar(progress_bar, value = 10)

      ## Sorting interactions P_P
      cond <- (gi@anchor1 > gi@anchor2) & (!grepl("OE",gi$int))

      a1 <- gi@anchor1[cond]
      a2 <- gi@anchor2[cond]

      gi@anchor1[cond] <- a2
      gi@anchor2[cond] <- a1

      elementMetadata(gi[cond])[c("gene_I","gene_II")] <- elementMetadata(gi[cond])[c("gene_II","gene_I")]

      gi@elementMetadata <- gi@elementMetadata[,-which(colnames(gi@elementMetadata) %in% c("counts"))]
      setTxtProgressBar(progress_bar, value = 11)
      close(progress_bar)
      message("\n")
      return(gi)
    }
    else(
      stop("Loading Interactions file error")
    )
  }
}

