#' Loads peakmatrix file into GenomicInteractions Object
#'
#' This function loads peakmatrix files from Chicago R package into a GenomicInteractions Object, and remove possible duplicated interactions
#'
#' @param file full path to the peakmatrix file
#'
#' @return GenomicInteractions object
#'
#' @importFrom magrittr `%>%`
#' @importFrom tidyr separate
#' @importFrom dplyr as_tibble group_by filter slice ungroup n summarise_all
#' @importFrom stringr str_replace_all
#' @importFrom GenomicInteractions GenomicInteractions
#' @importFrom GenomicRanges makeGRangesFromDataFrame split mcols
#' @importFrom data.table fread
#'
#' @export
load_peakmatrix <- function (file)
{
  if (!file.exists(file)) {
    stop(paste(basename(file), "does not exist"))
  }
  else {
    ## Setting pipe operator from magrittr package
    `%>%` <- magrittr::`%>%`

    message(paste(basename(file), "is a peakmatrix"))
    a <- data.table::fread(file = file, header = T, stringsAsFactors = F, na.strings = "")
    a1 <- a[, c(1:5, 11:ncol(a)),with=F]
    a2 <- a[, c(6:ncol(a)),with=F]
    colnames(a2) <- colnames(a1)
    df <- rbind(data.frame(a1, index = 1:nrow(a1)), data.frame(a2,index = 1:nrow(a2)))
    df <- df[order(df$index), ]
    data <- df[, -c(which(colnames(df)=="index"))]
    rownames(data) <- 1:nrow(data)
    data$rownames <- 1:nrow(data)


    ## Extracting name of all cell types in this peakmatrix
    cell_types <- colnames(data)[7:(ncol(data)-1)]

    ## Putting together in one line each interactions and duplicating them
    new_data2 <- rbind(cbind(data[seq(1, nrow(data), 2), ], data[seq(2, nrow(data), 2), ]),
                       cbind(data[seq(2, nrow(data), 2), ], data[seq(1, nrow(data), 2), ]))

    ## Ordering by the original line that came
    new_data2 <- dplyr::as_tibble(new_data2[order(as.numeric(rownames(new_data2))), ], .name_repair = "minimal")
    colnames(new_data2) <- c("chr_I", "start_I", "end_I","ID_I","gene_I","dist_I",  paste0("CS_I_ct",1:length(cell_types)),"rownames1",
                             "chr_II", "start_II", "end_II","ID_II","gene_II","dist_II",paste0("CS_II_ct",1:length(cell_types)),"rownames2")

    new_data2$gene_I <- stringr::str_replace_all(new_data2$gene_I, "\\|", ",")
    new_data2$gene_II <- stringr::str_replace_all(new_data2$gene_II, "\\|", ",")

    ## Removing real duplicates, if exist, in the file
    new_data2 <- data.table::as.data.table(new_data2)
    new_data2 <- unique(new_data2)
    dup_real <- (nrow(data) - nrow(new_data2))/2

    ## Filtering those duplicated interactions with different CS, by the higher one in all samples
    new_data2 <- new_data2[, lapply(.SD, max),by=list(ID_I,ID_II)]
    new_data2 <- new_data2[order(new_data2$rownames1),]

    dup_CS <- ((nrow(data) - nrow(new_data2))/2) - dup_real
    if (identical(new_data2$CS_I_ct1, new_data2$CS_II_ct1)) {

      ## Keeping only one of the artificial inverted duplications
      new_data2 <- dplyr::slice(dplyr::ungroup(new_data2),
                                seq(1, dplyr::n(), 2))


      message(paste0((nrow(data)/2) - nrow(new_data2)),
              " interactions removed\n\t- ", dup_real, " interactions were real duplicates\n\t- ",
              dup_CS, " interactions duplicated with different Chicago Score")

      new_data2 <- new_data2[, !colnames(new_data2) %in% c("rownames1","rownames2","dist_I",paste0("CS_I_ct",1:length(cell_types))), with=F]
      colnames(new_data2)[11:ncol(new_data2)] <- c("dist", cell_types)
      new_data2 <- new_data2[,c("chr_I","start_I","end_I","gene_I","ID_I",
                                "chr_II","start_II","end_II","gene_II","ID_II",
                                "dist",cell_types),with=F]

      ## reordering to have always captured baits as first end
      new_data2[new_data2$gene_I == ".", ] <- c(new_data2[new_data2$gene_I == ".", c(6:10, 1:5, 11:ncol(new_data2)),with=F])

      ## Creating the genomic interactions object
      region1 <- GenomicRanges::makeGRangesFromDataFrame(new_data2[,1:5], seqnames.field = "chr_I", start.field = "start_I",
                                                         end.field = "end_I", keep.extra.columns = T)
      region2 <- GenomicRanges::makeGRangesFromDataFrame(new_data2[,6:ncol(new_data2)], seqnames.field = "chr_II", start.field = "start_II",
                                                         end.field = "end_II", keep.extra.columns = T)
      gi <- GenomicInteractions::GenomicInteractions(region1, region2)


      names(GenomicRanges::mcols(gi)) <- gsub(x = names(GenomicRanges::mcols(gi)), pattern = "anchor[1-2]\\.", "")

      ## Annotating regions with P, OE or uce

      gi <- annotate_POEuce(gi)

      ## Sorting interactions P_P
      cond <- (gi@anchor1 > gi@anchor2) & (!grepl("OE", gi$int))

      a1 <- gi@anchor1[cond]
      a2 <- gi@anchor2[cond]

      gi@anchor1[cond] <- a2
      gi@anchor2[cond] <- a1

      elementMetadata(gi[cond])[grep("_I",colnames(elementMetadata(gi[cond]))[1:5],value = T)[c(1,3,2,4)]] <- elementMetadata(gi[cond])[grep("_I",colnames(elementMetadata(gi[cond]))[1:5],value = T)[c(3,1,4,2)]]

      gi@elementMetadata <- gi@elementMetadata[, -which(colnames(gi@elementMetadata) %in% c("counts"))]

      return(gi)
    }
    else (stop("Loading Peakmatrix file error"))
  }
}
