#' Filters GenomicInteractions object by integrating a regions file
#'
#' This function filters a GenomicInteractions object from load_interactions by integrating regions and classifies the interactions by type of interactions based on if the right, left or both end overlap with a region
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param regions full path to regions file (bed format) or a GRanges object
#' @param chr column name of chromosome values
#' @param start column name of start positions
#' @param end column name of end positions
#' @param invert T/F if need those interactions that do NOT overlaps with any regions
#'
#' @return GenomicInteractions object filtered by regions, by default with additional columns regarding overlap on each node. If invert=T no additional columns.
#'
#' @importFrom magrittr `%>%`
#' @importFrom GenomicInteractions anchorOne anchorTwo
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps overlapsAny countOverlaps
#' @importFrom data.table fread
#'
#' @export
interactionsByRegions <- function(interactions,regions,chr=NULL,start=NULL,end=NULL, invert=F)
{
  ## Setting pipe operator from magrittr package
  `%>%` <- magrittr::`%>%`

  if (class(regions) == "GRanges")
  {
    regionsGR <- regions
  }
  if (class(regions) == "character")
  {
    ## Reading regions and transforming to Genomic Ranges
    if(!is.null(chr) & !is.null(start) & !is.null(end)) ## If the file has header
    {
      regions_df <- data.table::fread(regions, header=T)
      regionsGR <- GenomicRanges::makeGRangesFromDataFrame(regions_df, seqnames.field = chr, start.field = start, end.field = end, keep.extra.columns = T)
    }
    if(is.null(chr) & is.null(start) & is.null(end)) ## If the file has not header
    {
      chr="V1"
      start="V2"
      end="V3"
      regions_df <- data.table::fread(regions)
      regionsGR <- GenomicRanges::makeGRangesFromDataFrame(regions_df, seqnames.field = chr, start.field = start, end.field = end, keep.extra.columns = T)
    }
    else ## Message if the arguments are no correctly filled
    {
      stop("chr, start, end must be all filled if regions have header")
    }
  }
  if (invert == T)
  {
    interactions_regions <- unique(IRanges::subsetByOverlaps(interactions, regionsGR, invert = invert))
    message(paste("Integration with regions results in",length(interactions_regions),"interactions that do not overlap with any given region"))
    annotate_POEuce(interactions_regions)
  }
  else
  {
    ## Subseting ibed by overlap with regions
    interactions_regions <- unique(IRanges::subsetByOverlaps(interactions, regionsGR, invert = invert))
    message(paste("Integration with regions results in",length(interactions_regions),"interactions"))

    if (length(interactions_regions) != 0)
    {

      ## Annotating and calculating overlaps in anchor one
      interactions_regions$overlap_I <- IRanges::overlapsAny(anchorOne(interactions_regions), regionsGR)
      interactions_regions$n_overlap_I <- IRanges::countOverlaps(anchorOne(interactions_regions), regionsGR)

      ## Annotating and calculating overlaps in anchor two
      interactions_regions$overlap_II <- IRanges::overlapsAny(anchorTwo(interactions_regions), regionsGR)
      interactions_regions$n_overlap_II <- IRanges::countOverlaps(anchorTwo(interactions_regions), regionsGR)


      annotate_POEuce(interactions_regions)
      ## Annotating type of interactions depending of overlap with regions
      interactions_regions$int <- ifelse(interactions_regions$overlap_I & interactions_regions$overlap_II,
                                         paste(GenomicInteractions::anchorOne(interactions_regions)$node.class,GenomicInteractions::anchorTwo(interactions_regions)$node.class,"II", sep = "_"),
                                         ifelse(interactions_regions$overlap_I & !interactions_regions$overlap_II,
                                                paste(GenomicInteractions::anchorOne(interactions_regions)$node.class,GenomicInteractions::anchorTwo(interactions_regions)$node.class,"I", sep = "_"),
                                                ifelse(!interactions_regions$overlap_I & interactions_regions$overlap_II,
                                                       paste(GenomicInteractions::anchorTwo(interactions_regions)$node.class,GenomicInteractions::anchorOne(interactions_regions)$node.class,"I", sep = "_"),NA)))

    }
  }
    return(interactions_regions)

}
