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
  if (length(interactions_list)< 6)
  {
    progress_bar = txtProgressBar(min=1, max=length(interactions_list)+length(interactions_list)^2, style = 3, char="=")
  }else
  {
    progress_bar = txtProgressBar(min=1, max=length(interactions_list), style = 3, char="=")
  }

  if(distance.boxplot)
  {

    for (i in 1:length(interactions_list))
    {
      a <- interactions_list[[i]]
      a$dist <- GenomicInteractions::calculateDistances(a)
      aa <- paste(interactions_list[[i]]@regions[interactions_list[[i]]@anchor1],interactions_list[[i]]@regions[interactions_list[[i]]@anchor2],sep = "_")
      names(aa) <- a$dist
      la[[names(interactions_list)[i]]] <- aa
      setTxtProgressBar(progress_bar, value = i)

    }

    ## Based on the function FromList from UpSetR package
    ## Generate a dataframe with the intersection and the distance of each interaction
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

    plot <- suppressWarnings(UpSetR::upset(data_final, nsets = length(interactions_list),boxplot.summary = "log10dist",...))
  }
  else
  {
    for (i in 1:length(interactions_list))
    {
      aa <- paste(interactions_list[[i]]@regions[interactions_list[[i]]@anchor1],interactions_list[[i]]@regions[interactions_list[[i]]@anchor2],sep = "_")
      la[[names(interactions_list)[i]]] <- aa
      setTxtProgressBar(progress_bar, value = i)
    }
    plot <- suppressWarnings(UpSetR::upset(UpSetR::fromList(la), nsets = length(interactions_list),...))
  }
  if(length(interactions_list) < 6)
  {
    plot2 <- gplots::venn(la,show.plot = F)
    int <- attr(gplots::venn(la,show.plot = F),"intersections")

    if (length(int) == 1)
    {
      is <- length(interactions_list)+length(interactions_list)^2
    }else
    {
      is <- seq(length(interactions_list),length(interactions_list)+length(interactions_list)^2,length.out=length(int)) ## for the progress bar
    }


    for (i in 1:length(int))
    {
      int_name <- str_split(names(int[i]),":")[[1]][1]
      df <- dplyr::as_tibble(str_split_fixed(int[[i]],pattern = ":|-|_",n = 6))
      a1 <- GenomicRanges::makeGRangesFromDataFrame(df[,1:3], seqnames.field = "V1",start.field = "V2",end.field = "V3")
      a2 <- GenomicRanges::makeGRangesFromDataFrame(df[,4:6], seqnames.field = "V4",start.field = "V5",end.field = "V6")
      gi <- GenomicInteractions::GenomicInteractions(a1,a2)
      gi <- IRanges::subsetByOverlaps(interactions_list[[int_name]],gi)

      int[[i]] <- gi
      setTxtProgressBar(progress_bar, value = is[i])
    }
    message("\nTo see the venn diagram do: plot(variable$venn)")
    return(list(upset_plot=plot,venn = plot2,intersections=int))
  }
  else
  {
    message("\nVenn diagram and intersections only can be generated with less than 6 sets")
    return(list(upset_plot = plot))
  }
}
