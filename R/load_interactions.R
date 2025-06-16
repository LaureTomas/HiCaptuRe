#' Loads interaction file into GenomicInteractions Object
#'
#' This function loads interaction files from Chicago R package into a GenomicInteractions Object, and remove possible duplicated interactions
#'
#' @param file full path to the interaction file (seqmonk, ibed, washU)
#' @param sep separator to read the file
#' @param ... arguments to pass to \link[HiCaptuRe]{digest_genome}
#'
#' @return HiCaptuRe object
#'
#' @importFrom tidyr separate
#' @importFrom dplyr as_tibble slice n
#' @importFrom stringr str_replace_all
#' @importFrom GenomicInteractions GenomicInteractions calculateDistances
#' @importFrom GenomicRanges makeGRangesFromDataFrame split mcols findOverlaps seqnames
#' @importFrom data.table fread
#' @importFrom progressr progressor
#' @importFrom S4Vectors elementMetadata subjectHits
#'
#' @examples
#' ibed1 <- system.file("extdata", "ibed1_example.zip", package="HiCaptuRe")
#' interactions <- load_interactions(ibed1, select_chr = "19")
#'
#' @export
load_interactions <- function(file, sep = "\t", ...) {
  if (!file.exists(file)) {
    stop(paste(basename(file), "does not exist"))
  }

  data <- data.table::fread(file = file, sep = sep, stringsAsFactors = F, na.strings = "")
  format <- .detect_format(data)
  process_function <- switch(format,
    ibed = .process_ibed,
    washU = .process_washU,
    washUold = .process_washUold,
    peakmatrix = .process_peakmatrix,
    bedpe = .process_bedpe,
    seqmonk = .process_seqmonk
  )

  new_datadd <- process_function(data)

  digest <- tryCatch(
    {
      digest_genome(...)
    },
    error = function(e) {
      e$call[1] <- call("digest_genome")
      stop(e)
    }
  )

  gi <- .generate_GInteractions(new_datadd, digest)

  final <- HiCaptuRe(genomicInteractions = gi, parameters = list(digest = digest$parameters, load = c(file = normalizePath(file), format = format)), ByBaits = list(), ByRegions = list())
  final$distance <- GenomicInteractions::calculateDistances(final)

  final <- final[order(final$ID_1, final$ID_2)]

  return(final)
}


.detect_format <- function(data) {
  if (ncol(data) > 11) {
    peakmatrix_columns <- c(
      "baitChr", "baitStart", "baitEnd", "baitID", "baitName",
      "oeChr", "oeStart", "oeEnd", "oeID", "oeName",
      "dist"
    )
    if (all(peakmatrix_columns == (colnames(data)[1:11]))) {
      format <- "peakmatrix"
    } else {
      stop(paste("File has the number of columns of a peakmatrix file but not the correct columns names. Peakmatrix files col.names should be:", paste(peakmatrix_columns, collapse = " ")))
    }
  } else if (ncol(data) == 10) {
    ibed_columns <- c(
      "bait_chr", "bait_start", "bait_end", "bait_name",
      "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
      "N_reads", "score"
    )
    if (all(colnames(data) == ibed_columns)) {
      format <- "ibed"
    } else if (all(colnames(data) == paste0("V", 1:10))) {
      format <- "bedpe"
    } else {
      stop(paste("File has the number of columns of an ibed file or a bedpe file but not the correct columns names. Bedpe files should not have header, and ibed files col.names should be:", paste(ibed_columns, collapse = " ")))
    }
  } else if (ncol(data) == 6) {
    if (all(colnames(data) == paste0("V", 1:6))) {
      format <- "seqmonk"
    } else {
      stop("File has the number of columns of a seqmonk file but should not have header")
    }
  } else if (ncol(data) %in% 3:4) {
    if (!grepl(":", data[1, 1])) {
      format <- "washU"
    } else if (grepl(":", data[1, 1])) {
      format <- "washU_old"
    } else {
      stop("File has the number of columns of a washU file but not the proper format")
    }
  }
  return(format)
}

.deduplicate_interactions <- function(new_data) {
  check1 <- unique(new_data[, c("chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2")])
  check2 <- unique(new_data[, c("chr_1", "start_1", "end_1", "bait_1", "chr_2", "start_2", "end_2", "bait_2")])

  if (nrow(check1) != nrow(check2)) {
    warning("There are fragments with the same coordinates but different annotations. Use annotate_interactions() to properly annotate")
  }
  ## Removing real duplicates, if exist, in the file
  new_data <- unique(new_data)

  # dup_real <- (nrow(data) - nrow(new_data)) / 2

  ## Filtering those duplicated interactions with different CS, by the higher one
  new_data <- new_data[, lapply(.SD, max), by = list(chr_1, start_1, end_1, chr_2, start_2, end_2)]
  new_data <- new_data[order(new_data$rownames1), ]

  if (all.equal(new_data[, grep("CS_1.*", colnames(new_data), perl = T), with = F], new_data[, grep("CS_2.*", colnames(new_data)), with = F], check.attributes = F)) {
    ## Keeping only one of the artificial inverted duplications
    new_data <- dplyr::slice(new_data, seq(1, dplyr::n(), 2))
    return(new_data)
  } else {
    stop("Problem deduplicating interactions")
  }
}


.process_peakmatrix <- function(data) {
  data <- data[, !grepl("dist|baitID|oeID", colnames(data)), with = F]
  score_names <- paste0("CS_", colnames(data)[9:(ncol(data))])

  reads <- rep(0, nrow(data))
  warning("reads column set to 0 because peakmatrix format does not contain this info")
  data <- cbind(data[, 1:8], reads, data[, 9:ncol(data), with = F])

  data <- .formating_data(data)
  new_datadd <- .process_data(data, score_names)
  return(new_datadd)
}

.process_bedpe <- function(data) {
  annotations <- rep("non-annotated", nrow(data))
  reads <- rep(0, nrow(data))
  warning("reads column set to 0 and annotation set to 'non-annotated' because bedpe format does not contain this info")
  data <- cbind(data[, 1:3], annotations, data[, 4:6], annotations, reads, data[, 8])

  score_names <- "CS"
  data <- .formating_data(data)
  new_datadd <- .process_data(data, score_names)
  return(new_datadd)
}

.process_ibed <- function(data) {
  score_names <- "CS"
  data <- .formating_data(data)
  new_datadd <- .process_data(data, score_names)
  return(new_datadd)
}

.process_seqmonk <- function(data) {
  new_datadd <- .process_data(data, score_names)

  data$rownames <- 1:nrow(data)
}

.process_washU <- function(data) {
  data <- tidyr::separate(data, 4,
    into = c("a", "b"), sep = ":", remove = TRUE,
    convert = FALSE, extra = "warn", fill = "warn"
  ) |>
    tidyr::separate(5,
      into = c("b", "c"), sep = "-", remove = TRUE,
      convert = FALSE, extra = "warn", fill = "warn"
    ) |>
    tidyr::separate(6,
      into = c("c", "d"), sep = ",", remove = TRUE,
      convert = FALSE, extra = "warn", fill = "warn"
    )

  annotations <- rep("non-annotated", nrow(data))
  reads <- rep(0, nrow(data))

  warning("reads column set to 0 and annotation set to 'non-annotated' because washU format does not contain this info")
  data <- cbind(data[, 1:3], annotations, data[, 4:6], annotations, reads, data[, 7])

  score_names <- "CS"
  data <- .formating_data(data.table::as.data.table(data))
  new_datadd <- .process_data(data, score_names)
  return(new_datadd)
}


.process_washUold <- function(data) {
  data <- tidyr::separate(data, 1,
    into = c("a", "b"), sep = ":", remove = TRUE,
    convert = FALSE, extra = "warn", fill = "warn"
  ) |>
    tidyr::separate(2,
      into = c("b", "c"), sep = ",", remove = TRUE,
      convert = FALSE, extra = "warn", fill = "warn"
    ) |>
    tidyr::separate(4,
      into = c("d", "e"), sep = ":", remove = TRUE,
      convert = FALSE, extra = "warn", fill = "warn"
    ) |>
    tidyr::separate(5,
      into = c("e", "f"), sep = ",", remove = TRUE,
      convert = FALSE, extra = "warn", fill = "warn"
    )

  annotations <- rep("non-annotated", nrow(data))
  reads <- rep(0, nrow(data))

  warning("reads column set to 0 and annotation set to 'non-annotated' because washU format does not contain this info")
  data <- cbind(data[, 1:3], annotations, data[, 4:6], annotations, reads, data[, 7])

  score_names <- "CS"
  data <- .formating_data(data.table::as.data.table(data))
  new_datadd <- .process_data(data, score_names)
  return(new_datadd)
}

.formating_data <- function(data) {
  df1 <- data[, c(1:4, 9:ncol(data)), with = F]
  df2 <- data[, c(5:ncol(data)), with = F]
  colnames(df2) <- colnames(df1)

  df <- rbind(data.table::data.table(df1, index = 1:nrow(df1)), data.table::data.table(df2, index = 1:nrow(df2)))

  df <- df[order(df$index), ]
  data <- df[, !grepl("index", colnames(df)), with = F]
  return(data)
}

.process_data <- function(data, score_names) {
  data$rownames <- 1:nrow(data)

  ## Putting together in one line each interactions and duplicating them
  new_data <- rbind(
    cbind(data[seq(1, nrow(data), 2), ], data[seq(2, nrow(data), 2), ]),
    cbind(data[seq(2, nrow(data), 2), ], data[seq(1, nrow(data), 2), ])
  )

  ## Ordering by the original line that came
  new_data <- new_data[order(new_data$rownames), ]
  colnames(new_data) <- c(
    "chr_1", "start_1", "end_1", "bait_1", "read_1", paste0("CS_1_ct", 1:length(score_names)), "rownames1",
    "chr_2", "start_2", "end_2", "bait_2", "read_2", paste0("CS_2_ct", 1:length(score_names)), "rownames2"
  )

  new_datadd <- .deduplicate_interactions(new_data)

  new_datadd <- new_datadd[, !colnames(new_datadd) %in% c("rownames1", "rownames2", "read_1", paste0("CS_1_ct", 1:length(score_names))), with = F]
  colnames(new_datadd)[9:ncol(new_datadd)] <- c("reads", score_names)
  new_datadd <- new_datadd[, c(
    "chr_1", "start_1", "end_1", "bait_1",
    "chr_2", "start_2", "end_2", "bait_2",
    "reads", score_names
  ), with = F]
  flip <- c(5:8, 1:4, 9:ncol(new_datadd))
  new_datadd[new_datadd$bait_1 == ".", ] <- new_datadd[new_datadd$bait_1 == ".", ..flip]
  return(new_datadd)
}


.generate_GInteractions <- function(new_datadd, digest) {
  seqnames_data <- unique(c(new_datadd$chr_1, new_datadd$chr_2))
  seqnames_digest <- GenomicRanges::seqnames(digest$seqinfo)

  if (!all(seqnames_data %in% seqnames_digest)) {
    stop(paste("Some chromosomes from data are not present in the digest genome. Missing chromosomes:", paste(seqnames_data[!seqnames_data %in% seqnames_digest], collapse = ", "), ". Check digest_genome() arguments."))
  }

  digestGR <- GenomicRanges::makeGRangesFromDataFrame(digest$digest, keep.extra.columns = T)

  region1 <- GenomicRanges::makeGRangesFromDataFrame(new_datadd[, 1:4], seqnames.field = "chr_1", start.field = "start_1", end.field = "end_1", keep.extra.columns = T, seqinfo = digest$seqinfo)
  region2 <- GenomicRanges::makeGRangesFromDataFrame(new_datadd[, 5:ncol(new_datadd)], seqnames.field = "chr_2", start.field = "start_2", end.field = "end_2", keep.extra.columns = T, seqinfo = digest$seqinfo)

  ID1 <- GenomicRanges::findOverlaps(region1, digestGR)
  ID2 <- GenomicRanges::findOverlaps(region2, digestGR)

  if (length(ID1) == 0 | length(ID2) == 0) {
    stop("No fragment found in digest. Maybe the genome version is not correct")
  }
  if (length(unique(region1)) != length(unique(S4Vectors::subjectHits(ID1))) | length(unique(region2)) != length(unique(S4Vectors::subjectHits(ID2)))) {
    stop("Digest does not perfectly match with fragments in data. Some fragments from digest overlap more than one fragment in data, or viceversa. Check digest_genome() arguments.")
  }

  region1$ID_1 <- digestGR$fragment_ID[S4Vectors::subjectHits(ID1)]
  region2$ID_2 <- digestGR$fragment_ID[S4Vectors::subjectHits(ID2)]

  cols <- names(GenomicRanges::mcols(region2))
  order <- c(grep("bait_2", cols), grep("ID_2", cols), grep("bait_2|ID_2", cols, invert = T))
  S4Vectors::elementMetadata(region2) <- S4Vectors::elementMetadata(region2)[, order]

  gi <- GenomicInteractions::GenomicInteractions(region1, region2)

  names(GenomicRanges::mcols(gi)) <- gsub(x = names(GenomicRanges::mcols(gi)), pattern = "anchor[1-2]\\.", "")

  ## Annotating regions with B or OE
  gi <- .annotate_BOE(gi)

  ## Sorting interactions B_B
  cond <- ((gi$ID_1 > gi$ID_2) & gi$int == "B_B") | ((gi$ID_1 < gi$ID_2) & gi$int == "OE_B")

  a1 <- gi@anchor1[cond]
  a2 <- gi@anchor2[cond]

  gi@anchor1[cond] <- a2
  gi@anchor2[cond] <- a1

  cols <- sort(grep("_", colnames(S4Vectors::elementMetadata(gi[cond]))[1:4], value = T))
  S4Vectors::elementMetadata(gi[cond])[cols] <- S4Vectors::elementMetadata(gi[cond])[cols[c(rbind(seq(2, length(cols), 2), seq(1, length(cols), 2)))]]

  return(gi)
}
