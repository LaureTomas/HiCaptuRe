#' Filters HiCaptuRe object by overlaping regions
#'
#' This function filters a HiCaptuRe object from load_interactions by overlaping regions
#'
#' @param interactions HiCaptuRe object
#' @param regions full path to regions file (bed format) or a GRanges object
#' @param chr column name of chromosome values
#' @param start column name of start positions
#' @param end column name of end positions
#' @param invert T/F if need those interactions that do NOT overlaps with any regions
#'
#' @return HiCaptuRe object filtered by regions, by default with additional columns regarding overlap on each node. If invert=T no additional columns. And an additional slot ByRegions with region-centric statistics
#'
#' @importFrom GenomicInteractions anchorOne anchorTwo
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps overlapsAny countOverlaps pintersect mergeByOverlaps
#' @importFrom data.table fread
#' @importFrom dplyr as_tibble group_by mutate summarise n rename filter
#' @importFrom tidyr separate_rows
#' @importFrom S4Vectors elementMetadata width
#'
#' @export
interactionsByRegions <- function(interactions, regions, chr = NULL, start = NULL, end = NULL, invert = F) {
  ## Setting pipe operator from magrittr package
  `%>%` <- magrittr::`%>%`

  if (is(regions, "GRanges")) {
    regions_name <- deparse(substitute(regions))
    regionsGR <- regions
  } else if (is(regions, "character")) {
    regions_name <- regions

    ## Reading regions and transforming to Genomic Ranges
    if (!is.null(chr) & !is.null(start) & !is.null(end)) ## If the file has header
      {
        regions_df <- data.table::fread(regions, header = T)
        regionsGR <- GenomicRanges::makeGRangesFromDataFrame(regions_df, seqnames.field = chr, start.field = start, end.field = end, keep.extra.columns = T)
      } else if (is.null(chr) & is.null(start) & is.null(end)) ## If the file has not header
      {
        chr <- "V1"
        start <- "V2"
        end <- "V3"
        regions_df <- data.table::fread(regions)
        regionsGR <- GenomicRanges::makeGRangesFromDataFrame(regions_df, seqnames.field = chr, start.field = start, end.field = end, keep.extra.columns = T)
      } else ## Message if the arguments are no correctly filled
    {
      stop("chr, start, end arguments must be all filled if regions file has a header")
    }
  } else {
    stop("`regions` must be either a GRanges object or a character path to a BED-like file.")
  }
  regionsGR$regionID <- 1:length(regionsGR)

  if (invert) {
    interactions_regions <- unique(IRanges::subsetByOverlaps(interactions, regionsGR, invert = invert))

    if (length(interactions_regions) != length(interactions)) {
      anchor1 <- suppressWarnings(unique(IRanges::mergeByOverlaps(GenomicInteractions::anchorOne(interactions), regionsGR)))
      anchor1$intersect <- suppressWarnings(S4Vectors::width(IRanges::pintersect(anchor1[, 1], anchor1$regionsGR)))
      if (nrow(anchor1) != 0) {
        anchor1$B.id <- unlist(anchor1$B.id)
      }
      anchor2 <- suppressWarnings(unique(IRanges::mergeByOverlaps(GenomicInteractions::anchorTwo(interactions), regionsGR)))
      anchor2$intersect <- suppressWarnings(S4Vectors::width(IRanges::pintersect(anchor2[, 1], anchor2$regionsGR)))
      if (nrow(anchor2) != 0) {
        anchor2$B.id <- unlist(anchor2$B.id)
      }
      df <- dplyr::as_tibble(unique(rbind(anchor1[, c("regionID", "fragmentID", "B.id")], anchor2[, c("regionID", "fragmentID", "B.id")]))) %>%
        dplyr::group_by(regionID) %>%
        dplyr::mutate(annot = ifelse(is.na(B.id), ".", B.id)) %>%
        dplyr::summarise(
          Nfragment = dplyr::n(),
          NOE = sum(annot == "."),
          fragmentID = paste(fragmentID, collapse = ","),
          fragmentAnnot = paste(unique(B.id), collapse = ",")
        )

      byregions <- makeGRangesFromDataFrame(merge(regionsGR, df, by = "regionID", all = T), keep.extra.columns = T)
    } else {
      byregions <- regionsGR
      S4Vectors::elementMetadata(byregions) <- cbind(S4Vectors::elementMetadata(byregions), data.frame(Nfragment = NA, NOE = NA, fragmentID = NA, fragmentAnnot = NA))
    } ## end if/else length == interactions
  } ## end if invert T
  else {
    ## Subseting ibed by overlap with regions
    interactions_regions <- unique(IRanges::subsetByOverlaps(interactions, regionsGR, invert = invert))

    if (length(interactions_regions) != 0) {
      anchor1 <- suppressWarnings(unique(IRanges::mergeByOverlaps(GenomicInteractions::anchorOne(interactions), regionsGR)))
      anchor1$intersect <- suppressWarnings(S4Vectors::width(IRanges::pintersect(anchor1[, 1], anchor1[, "regionsGR"])))
      if (nrow(anchor1) != 0) {
        anchor1$B.id <- unlist(anchor1$B.id)
      }
      df1 <- dplyr::as_tibble(anchor1[, c("fragmentID", "regionID", "intersect")]) %>%
        dplyr::rename(ID_1 = fragmentID) %>%
        dplyr::group_by(ID_1) %>%
        dplyr::summarise(
          region_1 = T,
          Nregion_1 = dplyr::n(),
          regionID_1 = paste(regionID, collapse = ","),
          regionCov_1 = sum(intersect)
        )

      m1 <- merge(S4Vectors::elementMetadata(interactions_regions), df1, all = T)


      anchor2 <- suppressWarnings(unique(IRanges::mergeByOverlaps(GenomicInteractions::anchorTwo(interactions), regionsGR)))
      anchor2$intersect <- suppressWarnings(S4Vectors::width(IRanges::pintersect(anchor2[, 1], anchor2[, "regionsGR"])))
      if (nrow(anchor2) != 0) {
        anchor2$B.id <- unlist(anchor2$B.id)
      }
      df2 <- dplyr::as_tibble(anchor2[, c("fragmentID", "regionID", "intersect")]) %>%
        dplyr::rename(ID_2 = fragmentID) %>%
        dplyr::group_by(ID_2) %>%
        dplyr::summarise(
          region_2 = T,
          Nregion_2 = dplyr::n(),
          regionID_2 = paste(regionID, collapse = ","),
          regionCov_2 = sum(intersect)
        )

      m2 <- merge(S4Vectors::elementMetadata(interactions_regions), df2, all = T)

      m <- merge(m1, m2, by = colnames(S4Vectors::elementMetadata(interactions_regions)))
      m <- m[order(m$ID_1, m$ID_2), ]
      m$region_1[is.na(m$region_1)] <- F
      m$region_2[is.na(m$region_2)] <- F
      m$Nregion_1[is.na(m$Nregion_1)] <- 0
      m$Nregion_2[is.na(m$Nregion_2)] <- 0
      m$regionCov_1[is.na(m$regionCov_1)] <- 0
      m$regionCov_2[is.na(m$regionCov_2)] <- 0


      if (all(as.data.frame(m[, 1:ncol(S4Vectors::elementMetadata(interactions_regions))]) == as.data.frame(S4Vectors::elementMetadata(interactions_regions)), na.rm = T)) {
        S4Vectors::elementMetadata(interactions_regions) <- m
      } else {
        stop("Data order is not correct. Maybe some bait was initially clasified as Other-End")
      }

      df <- dplyr::as_tibble(unique(rbind(anchor1[, c("regionID", "fragmentID", "B.id")], anchor2[, c("regionID", "fragmentID", "B.id")]))) %>%
        dplyr::group_by(regionID) %>%
        dplyr::mutate(annot = ifelse(is.na(B.id), ".", B.id)) %>%
        dplyr::summarise(
          Nfragment = dplyr::n(),
          NOE = sum(annot == "."),
          fragmentID = paste(fragmentID, collapse = ","),
          fragmentAnnot = paste(unique(B.id), collapse = ",")
        )

      m1 <- dplyr::as_tibble(m[, c("regionID_1", "ID_1", "ID_2")]) %>%
        dplyr::mutate(ID = paste(ID_1, ID_2, sep = "_")) %>%
        dplyr::rename(regionID = regionID_1)
      m2 <- dplyr::as_tibble(m[, c("regionID_2", "ID_1", "ID_2")]) %>%
        dplyr::mutate(ID = paste(ID_1, ID_2, sep = "_")) %>%
        dplyr::rename(regionID = regionID_2)
      mfinal <- rbind(m1[, c(1, 4)], m2[, c(1, 4)]) %>%
        dplyr::filter(!is.na(regionID)) %>%
        tidyr::separate_rows(regionID, sep = ",") %>%
        dplyr::mutate(regionID = as.numeric(regionID)) %>%
        dplyr::group_by(regionID) %>%
        dplyr::summarise(N_int = length(unique(ID)))

      df <- merge(mfinal, df, by = "regionID")

      byregions <- GenomicRanges::makeGRangesFromDataFrame(merge(regionsGR, df, by = "regionID", all = T), keep.extra.columns = T)
    } else {
      byregions <- regionsGR
      S4Vectors::elementMetadata(byregions) <- cbind(S4Vectors::elementMetadata(byregions), data.frame(Nfragment = NA, NOE = NA, fragmentID = NA, fragmentAnnot = NA))
    } ## end if/else length == != 0
  } ## end else invert F


  # updating slots

  ByRegions_list <- getByRegions(interactions_regions)

  if (length(ByRegions_list) == 0) {
    ByRegions_list[[1]] <- byregions
  } else {
    ByRegions_list[[length(ByRegions_list) + 1]] <- byregions
  }

  interactions_regions <- setByRegions(interactions_regions, ByRegions_list)

  param <- getParameters(interactions_regions)
  param[[paste0("ByRegions_", length(ByRegions_list))]] <- c(interactions = deparse(substitute(interactions)), regions = regions_name, chr = chr, start = start, end = end, invert = invert)
  interactions_regions <- setParameters(interactions_regions, param)

  return(interactions_regions)
}
