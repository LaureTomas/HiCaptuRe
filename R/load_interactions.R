#' Loads interaction file into GenomicInteractions Object
#'
#' This function loads interaction files from Chicago R package into a GenomicInteractions Object, and remove possible duplicated interactions
#'
#' @param file full path to the interaction file (seqmonk, ibed, washU)
#' @param washU_seqname character to remove from seqnames column when using washU format
#' @param ... arguments to pass to \link[HiCaptuRe]{digest_genome}
#'
#' @return HiCaptuRe object
#'
#' @importFrom tidyr separate
#' @importFrom dplyr as_tibble slice n
#' @importFrom stringr str_replace_all
#' @importFrom GenomicInteractions GenomicInteractions calculateDistances
#' @importFrom GenomicRanges makeGRangesFromDataFrame split mcols findOverlaps seqnames
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom data.table fread
#' @importFrom progressr progressor
#' @importFrom S4Vectors elementMetadata subjectHits
#'
#' @export
load_interactions <- function(file, washU_seqname = "chr", ...) {
  if (!file.exists(file)) {
    stop(paste(basename(file), "does not exist"))
  } else {
    ## Setting pipe operator from magrittr package
    `%>%` <- magrittr::`%>%`

    ## Reading file and detecting file format depending of the number of columns
    ## Transforming all file formats to seqmonk to proceed with the cleaning

    p <- progressr::progressor(steps = 10)

    p(sprintf("Reading File"))
    data <- data.table::fread(file = file, header = T, stringsAsFactors = F, na.strings = "")

    if (ncol(data) > 10) {
      type <- "peakmatrix"

      data <- data[, !grepl("dist", colnames(data)), with = F]
      reads <- rep(0, nrow(data))
      df1 <- cbind(data[, 1:5], reads, data[, 11:ncol(data), with = F])
      df2 <- cbind(data[, 6:10], reads, data[, 11:ncol(data), with = F])
      colnames(df2) <- colnames(df1)

      df <- rbind(data.frame(df1, index = 1:nrow(df1)), data.frame(df2, index = 1:nrow(df2)))

      df <- df[order(df$index), ]
      data <- df[, -c(which(colnames(df) == "index"))]
      data$rownames <- 1:nrow(data)
      p(sprintf("Preparing Data"))


      ## Extracting name of all cell types in this peakmatrix
      cell_types <- paste0("CS_", colnames(data)[7:(ncol(data) - 1)])

      ## Putting together in one line each interactions and duplicating them
      new_data <- rbind(
        cbind(data[seq(1, nrow(data), 2), ], data[seq(2, nrow(data), 2), ]),
        cbind(data[seq(2, nrow(data), 2), ], data[seq(1, nrow(data), 2), ])
      )
      p(sprintf("Sorting Data"))

      ## Ordering by the original line that came
      new_data <- dplyr::as_tibble(new_data[order(as.numeric(rownames(new_data))), ], .name_repair = "minimal")
      colnames(new_data) <- c(
        "chr_1", "start_1", "end_1", "ID_1", "bait_1", "read_1", paste0("CS_1_ct", 1:length(cell_types)), "rownames1",
        "chr_2", "start_2", "end_2", "ID_2", "bait_2", "read_2", paste0("CS_2_ct", 1:length(cell_types)), "rownames2"
      )
    } else {
      if (ncol(data) == 6) {
        type <- "seqmonk"

        data <- data.table::fread(file = file, header = F, stringsAsFactors = F, na.strings = "")

        message(paste(basename(file), "is in seqmonk format"))
        data$rownames <- 1:nrow(data)
        p(sprintf("Preparing Data"))
      }
      if (ncol(data) == 10 & any(grepl("bait", colnames(data)))) {
        type <- "ibed"

        message(paste(basename(file), "is in ibed format"))
        a <- data
        a1 <- a[, c(1:4, 9:10), with = F]
        a2 <- a[, c(5:10), with = F]
        colnames(a2) <- colnames(a1)
        df <- rbind(data.frame(a1, index = 1:nrow(a1)), data.frame(a2, index = 1:nrow(a2)))
        df <- df[order(df$index), ]
        data <- df[, 1:6]
        data$rownames <- 1:nrow(data)
        p(sprintf("Preparing Data"))
      }
      if (ncol(data) == 10 & !any(grepl("bait", colnames(data)))) {
        type <- "bedpe"

        message(paste(basename(file), "is in bedpe format"))
        warning("We do not recommend to use bedpe format \n The HiCaptuRe output must be annotated, see ??annotate_interactions")

        data <- data.table::fread(file = file, header = F, stringsAsFactors = F, sep = "\t") # reading file
        annotations <- rep("non-annotated", nrow(data))
        reads <- rep(0, nrow(data))

        df1 <- cbind(data[, 1:3], annotations, reads, data[, 8])
        df2 <- cbind(data[, 4:6], annotations, reads, data[, 8])
        colnames(df2) <- colnames(df1)

        df <- rbind(data.frame(df1, index = 1:nrow(df1)), data.frame(df2, index = 1:nrow(df2)))
        df <- df[order(df$index), ]
        data <- df[, 1:6]
        data$rownames <- 1:nrow(data)
        p(sprintf("Preparing Data"))
      }


      if (ncol(data) %in% 3:4) ## if washU
        {
          if (!grepl(":", data[1, 1])) {
            type <- "washU"

            message(paste(basename(file), "is in washU new format"))
            warning("We do not recommend to use washU format from Chicago \n The HiCaptuRe output must be annotated, see ??annotate_interactions")

            data <- data.table::fread(file = file, header = F, stringsAsFactors = F, sep = "\t") # reading file

            data <- tidyr::separate(data, 4,
              into = c("a", "b"), sep = ":", remove = TRUE,
              convert = FALSE, extra = "warn", fill = "warn"
            ) %>%
              tidyr::separate(5,
                into = c("b", "c"), sep = "-", remove = TRUE,
                convert = FALSE, extra = "warn", fill = "warn"
              ) %>%
              tidyr::separate(6,
                into = c("c", "d"), sep = ",", remove = TRUE,
                convert = FALSE, extra = "warn", fill = "warn"
              )
            p(sprintf("Preparing Data"))
          }
          if (grepl(":", data[1, 1])) {
            type <- "washU_old"

            message(paste(basename(file), "is in washU old format"))
            warning("We do not recommend to use washU format from Chicago \n The HiCaptuRe output must be annotated, see ??annotate_interactions")

            data <- data.table::fread(file = file, header = F, stringsAsFactors = F, sep = "\t") # reading file

            data <- tidyr::separate(data, 1,
              into = c("a", "b"), sep = ":", remove = TRUE,
              convert = FALSE, extra = "warn", fill = "warn"
            ) %>%
              tidyr::separate(2,
                into = c("b", "c"), sep = ",", remove = TRUE,
                convert = FALSE, extra = "warn", fill = "warn"
              ) %>%
              tidyr::separate(4,
                into = c("d", "e"), sep = ":", remove = TRUE,
                convert = FALSE, extra = "warn", fill = "warn"
              ) %>%
              tidyr::separate(5,
                into = c("e", "f"), sep = ",", remove = TRUE,
                convert = FALSE, extra = "warn", fill = "warn"
              )
            p(sprintf("Preparing Data"))
          }

          data[, c(1, 4)] <- lapply(data[, c(1, 4)], function(x) as.character(gsub(washU_seqname, "", x)))
          data[, c(2:3, 5:7)] <- as.data.frame(apply(data[, c(2:3, 5:7)], 2, function(x) as.numeric(x)))

          annotations <- rep("non-annotated", nrow(data))
          reads <- rep(0, nrow(data))

          df1 <- cbind(data[, 1:3], annotations, reads, data[, 7])
          df2 <- cbind(data[, 4:6], annotations, reads, data[, 7])
          colnames(df2) <- colnames(df1)

          df <- rbind(data.frame(df1, index = 1:nrow(df1)), data.frame(df2, index = 1:nrow(df2)))
          df <- df[order(df$index), ]
          data <- df[, 1:6]
          data$rownames <- 1:nrow(data)
        } ## end if washU

      ## Putting together in one line each interactions and duplicating them
      new_data <- rbind(
        cbind(data[seq(1, nrow(data), 2), ], data[seq(2, nrow(data), 2), ]),
        cbind(data[seq(2, nrow(data), 2), ], data[seq(1, nrow(data), 2), ])
      )
      p(sprintf("Sorting Data"))

      ## Ordering by the original line that came
      new_data <- dplyr::as_tibble(new_data[order(as.numeric(rownames(new_data))), ],
        .name_repair = "minimal"
      )

      colnames(new_data) <- c(
        "chr_1", "start_1", "end_1", "bait_1", "R_1", "CS_1", "rownames1",
        "chr_2", "start_2", "end_2", "bait_2", "R_2", "CS_2", "rownames2"
      )
    } # end non-peakmatrix

    p(sprintf("Real Duplicates"))

    ## Here, could be real duplicated interactions, check with unique
    new_data$bait_1 <- stringr::str_replace_all(new_data$bait_1, "\\|", ",")
    new_data$bait_2 <- stringr::str_replace_all(new_data$bait_2, "\\|", ",")

    ## After correcting the annotation the number of real duplicates increase

    ## Removing real duplicates, if exist, in the file
    new_data <- data.table::as.data.table(new_data)
    new_data <- unique(new_data)

    dup_real <- (nrow(data) - nrow(new_data)) / 2
    p(sprintf("CS Duplicates"))

    ## Filtering those duplicated interactions with different CS, by the higher one
    new_data <- new_data[, lapply(.SD, max), by = list(chr_1, start_1, end_1, chr_2, start_2, end_2)]
    new_data <- new_data[order(new_data$rownames1), ]

    dup_CS <- ((nrow(data) - nrow(new_data)) / 2) - dup_real

    if (all.equal(new_data[, grep("CS_1.*", colnames(new_data), perl = T), with = F], new_data[, grep("CS_2.*", colnames(new_data)), with = F], check.attributes = F)) {
      ## Keeping only one of the artificial inverted duplications
      new_data <- dplyr::slice(new_data, seq(1, dplyr::n(), 2))
      p(sprintf("Cleaning"))

      ## Removing interactions that involve the MT chromosome

      message(paste0("\n", (nrow(data) / 2) - nrow(new_data)), " interactions removed\n\t- ", dup_real, " interactions were real duplicates\n\t- ", dup_CS, " interactions duplicated with different Chicago Score")

      if (type == "peakmatrix") {
        new_data <- new_data[, !colnames(new_data) %in% c("rownames1", "rownames2", "read_1", paste0("CS_1_ct", 1:length(cell_types))), with = F]
        colnames(new_data)[11:ncol(new_data)] <- c("reads", cell_types)
        new_data <- new_data[, c(
          "chr_1", "start_1", "end_1", "bait_1", "ID_1",
          "chr_2", "start_2", "end_2", "bait_2", "ID_2",
          "reads", cell_types
        ), with = F]
        flip <- c(6:10, 1:5, 11:ncol(new_data))
        new_data[new_data$bait_1 == ".", ] <- new_data[new_data$bait_1 == ".", ..flip]
      } else {
        new_data <- new_data[, !colnames(new_data) %in% c("rownames1", "rownames2", "R_1", "CS_1"), with = F]
        colnames(new_data)[9:10] <- c("reads", "CS")
        new_data <- new_data[, c(
          "chr_1", "start_1", "end_1", "bait_1",
          "chr_2", "start_2", "end_2", "bait_2",
          "reads", "CS"
        ), with = F]
        flip <- c(5:8, 1:4, 9, 10)
        new_data[new_data$bait_1 == ".", ] <- new_data[new_data$bait_1 == ".", ..flip]
      }

      p(sprintf("Digesting Genome"))

      digest <- tryCatch(
        {
          digest_genome(...)
        },
        error = function(e) {
          e$call[1] <- call("digest_genome")
          stop(e)
        }
      )

      digestGR <- GenomicRanges::makeGRangesFromDataFrame(digest$digest, keep.extra.columns = T)

      ## Creating the genomic interactions object

      region1 <- GenomicRanges::makeGRangesFromDataFrame(new_data[, 1:(grep("chr_2", colnames(new_data)) - 1)], seqnames.field = "chr_1", start.field = "start_1", end.field = "end_1", keep.extra.columns = T, seqinfo = digest$seqinfo[GenomicRanges::seqnames(GenomeInfoDb::seqinfo(digestGR))])
      region2 <- GenomicRanges::makeGRangesFromDataFrame(new_data[, grep("chr_2", colnames(new_data)):ncol(new_data)], seqnames.field = "chr_2", start.field = "start_2", end.field = "end_2", keep.extra.columns = T, seqinfo = digest$seqinfo[GenomicRanges::seqnames(GenomeInfoDb::seqinfo(digestGR))])

      ID1 <- GenomicRanges::findOverlaps(region1, digestGR)
      ID2 <- GenomicRanges::findOverlaps(region2, digestGR)

      if (length(ID1) == 0 | length(ID2) == 0) {
        stop("No fragment found in digest.\nMaybe the genome version is not correct")
      }
      if (length(unique(region1)) != length(unique(S4Vectors::subjectHits(ID1))) | length(unique(region2)) != length(unique(S4Vectors::subjectHits(ID2)))) {
        stop("Digest does not perfectly match with fragments in data.\n Some fragments from digest overlap more than one fragment of your data, or viceversa")
      }

      region1$ID_1 <- digestGR$fragment_ID[S4Vectors::subjectHits(ID1)]
      region2$ID_2 <- digestGR$fragment_ID[S4Vectors::subjectHits(ID2)]

      cols <- names(GenomicRanges::mcols(region2))
      order <- c(grep("bait_2", cols), grep("ID_2", cols), grep("bait_2|ID_2", cols, invert = T))
      S4Vectors::elementMetadata(region2) <- S4Vectors::elementMetadata(region2)[, order]
      # GenomicRanges::mcols(region2) <- GenomicRanges::mcols(region2)

      p(sprintf("Type of Interactions"))

      gi <- GenomicInteractions::GenomicInteractions(region1, region2)

      names(GenomicRanges::mcols(gi)) <- gsub(x = names(GenomicRanges::mcols(gi)), pattern = "anchor[1-2]\\.", "")

      ## Annotating regions with P, OE or uce
      gi <- annotate_BOE(gi)

      p(sprintf("Sorting"))

      ## Sorting interactions P_P
      cond <- ((gi$ID_1 > gi$ID_2) & gi$int == "B_B") | ((gi$ID_1 < gi$ID_2) & gi$int == "OE_B")

      a1 <- gi@anchor1[cond]
      a2 <- gi@anchor2[cond]

      gi@anchor1[cond] <- a2
      gi@anchor2[cond] <- a1


      cols <- sort(grep("_", colnames(S4Vectors::elementMetadata(gi[cond]))[1:4], value = T))
      S4Vectors::elementMetadata(gi[cond])[cols] <- S4Vectors::elementMetadata(gi[cond])[cols[c(rbind(seq(2, length(cols), 2), seq(1, length(cols), 2)))]]

      # gi@elementMetadata <- gi@elementMetadata[,-which(colnames(gi@elementMetadata) %in% c("counts"))]

      p(sprintf("Finishing"))
      final <- HiCaptuRe(genomicInteractions = gi, parameters = list(digest = digest$parameters, load = c(file = normalizePath(file), type = type)), ByBaits = list(), ByRegions = list())
      final$distance <- GenomicInteractions::calculateDistances(final)

      final <- final[order(final$ID_1, final$ID_2)]

      return(final)
    } else {
      (
        stop("Loading Interactions file error")
      )
    }
  }
}
