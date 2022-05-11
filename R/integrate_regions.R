#' Filters GenomicInteractions object by integrating a regions file
#'
#' This function filters a GenomicInteractions object from load_interactions by integrating regions from a external regions file, and classifies the interactions by type of interactions based on if the right, left or both end overlap with a region
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param regions full path to regions file (bed format) or a GRanges object
#' @param chr column name of chromosome values
#' @param start column name of start positions
#' @param end column name of end positions
#' @param invert T/F if need those interactions that do NOT overlaps with any regions
#'
#' @return GenomicInteractions object filtered by regions with 'int' column containing type of interactions
#'
#' @importFrom magrittr `%>%`
#' @importFrom dplyr as_tibble group_by summarise n
#' @importFrom GenomicInteractions anchorOne anchorTwo
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom IRanges subsetByOverlaps
#' @importFrom data.table fread
#'
#' @export
integrate_regions <- function(interactions,regions,chr=NULL,start=NULL,end=NULL, invert=F)
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

  }
  else
  {
    ## Subseting ibed by overlap with regions
    interactions_regions <- unique(IRanges::subsetByOverlaps(interactions, regionsGR, invert = invert))
    message(paste("Integration with regions results in",length(interactions_regions),"interactions"))

    if (length(interactions_regions) != 0)
    {

      ## Annotating and calculating overlaps in anchor one
      a <- dplyr::as_tibble(GenomicRanges::findOverlaps(GenomicInteractions::anchorOne(interactions_regions), regionsGR)) %>% dplyr::group_by(queryHits) %>% dplyr::summarise(n_overlap=dplyr::n())
      interactions_regions$overlap_I <- F
      interactions_regions$overlap_I[a$queryHits] <- T
      interactions_regions$n_overlap_I <- 0
      interactions_regions$n_overlap_I[a$queryHits] <- a$n_overlap

      ## TODO use countoverlaps instead of findoverlaps

      ## Annotating and calculating overlaps in anchor two
      b <- dplyr::as_tibble(GenomicRanges::findOverlaps(GenomicInteractions::anchorTwo(interactions_regions), regionsGR)) %>% dplyr::group_by(queryHits) %>% dplyr::summarise(n_overlap=dplyr::n())
      interactions_regions$overlap_II <- F
      interactions_regions$overlap_II[b$queryHits] <- T
      interactions_regions$n_overlap_II <- 0
      interactions_regions$n_overlap_II[b$queryHits] <- b$n_overlap


      ## Annotating type of interactions depending of overlap with regions
      interactions_regions$int <- ifelse(interactions_regions$overlap_I & interactions_regions$overlap_II &
                                           (interactions_regions$gene_I == "." | interactions_regions$gene_II == ".") &
                                           (interactions_regions$gene_I != "uce" & interactions_regions$gene_II != "uce"), "P_OE_II",
                                         ifelse(interactions_regions$overlap_I & interactions_regions$overlap_II &
                                                  (interactions_regions$gene_I == "." | interactions_regions$gene_II == ".") &
                                                  (interactions_regions$gene_I == "uce" | interactions_regions$gene_II == "uce"), "uce_OE_II",
                                                ifelse(interactions_regions$overlap_I & interactions_regions$overlap_II &
                                                         !interactions_regions$gene_I %in% c(".","uce") & !interactions_regions$gene_II %in% c(".","uce"),"P_P_II",
                                                       ifelse(interactions_regions$overlap_I & interactions_regions$overlap_II & interactions_regions$gene_I != "." & interactions_regions$gene_II != "." &
                                                                (interactions_regions$gene_I == "uce" | interactions_regions$gene_II == "uce") &
                                                                !(interactions_regions$gene_I == "uce" & interactions_regions$gene_II == "uce"),"P_uce_II",
                                                              ifelse((interactions_regions$overlap_I | interactions_regions$overlap_II) & !(interactions_regions$overlap_I == interactions_regions$overlap_II) &
                                                                       !interactions_regions$gene_I %in% c(".","uce") & !interactions_regions$gene_II %in% c(".","uce"),"P_P_I",
                                                                     ifelse((interactions_regions$overlap_I & !interactions_regions$overlap_II &
                                                                               !interactions_regions$gene_I %in% c(".","uce") & interactions_regions$gene_II == ".") |
                                                                              (interactions_regions$overlap_II & !interactions_regions$overlap_I &
                                                                                 !interactions_regions$gene_II %in% c(".","uce") & interactions_regions$gene_I == "."),"P_OE_I",
                                                                            ifelse((interactions_regions$overlap_I & !interactions_regions$overlap_II &
                                                                                      !interactions_regions$gene_II %in% c(".","uce") & interactions_regions$gene_I == ".") |
                                                                                     (interactions_regions$overlap_II & !interactions_regions$overlap_I &
                                                                                        !interactions_regions$gene_I %in% c(".","uce") & interactions_regions$gene_II == "."),"OE_P_I",
                                                                                   ifelse((interactions_regions$overlap_I & !interactions_regions$overlap_II &
                                                                                             interactions_regions$gene_I == "uce" & !interactions_regions$gene_II %in% c(".","uce")) |
                                                                                            (interactions_regions$overlap_II & !interactions_regions$overlap_I &
                                                                                               interactions_regions$gene_II == "uce" & !interactions_regions$gene_I %in% c(".","uce")),"uce_P_I",
                                                                                          ifelse((interactions_regions$overlap_I & !interactions_regions$overlap_II &
                                                                                                    interactions_regions$gene_II == "uce" & !interactions_regions$gene_I %in% c(".","uce")) |
                                                                                                   (interactions_regions$overlap_II & !interactions_regions$overlap_I &
                                                                                                      interactions_regions$gene_I == "uce" & !interactions_regions$gene_II %in% c(".","uce")),"P_uce_I",
                                                                                                 ifelse((interactions_regions$overlap_I & !interactions_regions$overlap_II &
                                                                                                           interactions_regions$gene_I == "uce" & interactions_regions$gene_II == ".") |
                                                                                                          (interactions_regions$overlap_II & !interactions_regions$overlap_I &
                                                                                                             interactions_regions$gene_II == "uce" & interactions_regions$gene_I == "."),"uce_OE_I",
                                                                                                        ifelse((interactions_regions$overlap_I & !interactions_regions$overlap_II &
                                                                                                                  interactions_regions$gene_I == "." & interactions_regions$gene_II == "uce") |
                                                                                                                 (interactions_regions$overlap_II & !interactions_regions$overlap_I &
                                                                                                                    interactions_regions$gene_II == "." & interactions_regions$gene_I == "uce"),"OE_uce_I",
                                                                                                               ifelse((interactions_regions$overlap_I | interactions_regions$overlap_II) & !(interactions_regions$overlap_I == interactions_regions$overlap_II) &
                                                                                                                        interactions_regions$gene_I=="uce" & interactions_regions$gene_II=="uce","uce_uce_I",
                                                                                                                      ifelse(interactions_regions$overlap_I & interactions_regions$overlap_II &
                                                                                                                               interactions_regions$gene_I=="uce" & interactions_regions$gene_II=="uce","uce_uce_II",NA)))))))))))))



    }
  }
    return(interactions_regions)

}
