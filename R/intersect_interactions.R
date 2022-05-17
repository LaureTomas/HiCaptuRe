#' Computes the intersect between a list of interactions objects
#'
#' This function computes all the possible intersects for a given list of GenomicInteractions objects
#'
#' @param interactions_list list of GenomicInteractions objects from \code{\link{load_interactions}}
#' @param distance.boxplot logical, by default FALSE, plot a boxplot of log10 distance of each intersection in the upset plot
#' @param ... extra arguments for UpSetR
#'
#' @return list with 3 elements:
#' 1. 'upset_plot'  an upset plot of the intersections
#' 2. 'venn' a venn diagram of the intersections
#' 3. 'intersections' a list of GenomicInteractions objects for all the intersections
#'
#' @importFrom dplyr as_tibble
#' @importFrom stringr str_split_fixed
#' @importFrom GenomicInteractions GenomicInteractions calculateDistances
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom gplots venn
#' @importFrom UpSetR upset fromList
#'
#' @export
intersect_interactions <- function(interactions_list, distance.boxplot=F,...)
{
  la <- list()
  if (distance.boxplot)
  {
    steps <- length(interactions_list)+length(interactions_list)^2 + 1
  }
  if (!distance.boxplot)
  {
    steps <- length(interactions_list)+length(interactions_list)^2
  }
  p <- progressr::progressor(steps = steps)

  if(distance.boxplot)
  {

    for (i in 1:length(interactions_list))
    {
      p(sprintf("Preparing"))
      a <- interactions_list[[i]]
      a$dist <- GenomicInteractions::calculateDistances(a)
      aa <- paste(interactions_list[[i]]@regions[interactions_list[[i]]@anchor1],interactions_list[[i]]@regions[interactions_list[[i]]@anchor2],sep = "_")
      names(aa) <- a$dist
      la[[names(interactions_list)[i]]] <- aa
    }

    ## Based on the function FromList from UpSetR package
    ## Generate a dataframe with the intersection and the distance of each interaction
    p(sprintf("Boxplot"))

    elements1 <- unlist(la,use.names = T)
    names(elements1) <- gsub("\\.","",gsub(paste0(names(la),collapse = "|"),"",names(elements1)))
    elements <- unique(elements1)
    elements <- elements1[!duplicated(elements1)]
    data <- unlist(lapply(la, function(x) {
      x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)

    data <- data.frame(matrix(data, ncol = length(la), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(la)

    data_final <- suppressWarnings(cbind(data,as.numeric(names(elements))))
    colnames(data_final)[ncol(data_final)] <- "dist"
    data_final$log10dist <- log10(data_final$dist)

    p(sprintf("Upset"))
    uplot <- suppressWarnings(UpSetR::upset(data_final, nsets = length(interactions_list),boxplot.summary = "log10dist",...))
  }
  else
  {
    for (i in 1:length(interactions_list))
    {
      p(sprintf("Preparing"))
      aa <- paste(interactions_list[[i]]@regions[interactions_list[[i]]@anchor1],interactions_list[[i]]@regions[interactions_list[[i]]@anchor2],sep = "_")
      la[[names(interactions_list)[i]]] <- aa
    }

    p(sprintf("Upset"))
    uplot <- suppressWarnings(UpSetR::upset(UpSetR::fromList(la), nsets = length(interactions_list),...))
  }

  p(sprintf("Venn"))
  venn_plot <- gplots::venn(la,show.plot = F)
  p(sprintf("Intersections"))
  int <- attr(venn_plot,"intersections")

  for (i in 1:length(int))
  {
    p(sprintf("Intersections to GI"))
    int_name <- str_split(names(int[i]),":")[[1]][1]
    df <- dplyr::as_tibble(str_split_fixed(int[[i]],pattern = ":|-|_",n = 6))
    a1 <- GenomicRanges::makeGRangesFromDataFrame(df[,1:3], seqnames.field = "V1",start.field = "V2",end.field = "V3")
    a2 <- GenomicRanges::makeGRangesFromDataFrame(df[,4:6], seqnames.field = "V4",start.field = "V5",end.field = "V6")
    gi <- GenomicInteractions::GenomicInteractions(a1,a2)
    gi <- IRanges::subsetByOverlaps(interactions_list[[int_name]],gi)

    int[[i]] <- gi
  }


  message("\nTo see the venn diagram: plot(variable$venn)")
  return(list(intersections=int, upset_plot=uplot, venn = venn_plot))


}
