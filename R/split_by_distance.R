#' Computes the number of interactions by distance
#'
#' This function computes the number of interactions by distance for a given GenomicInteractions object
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param cell_type character variable with the name of the cell type or sample
#'
#' @return list with 2 tables: short_int_dist_table and long_int_dist_table, with the number of short interactions and long interactions respectivetly
#'
#' @importFrom magrittr `%>%`
#' @importFrom dplyr as_tibble summarise group_by n full_join
#' @importFrom GenomicInteractions calculateDistances
#' @importFrom utils read.table
#' @importFrom tibble add_column
#' @importFrom reshape2 melt
#'
#' @export
split_by_distance <- function(interactions, cell_type="cell_type")
  ## Introduce parameters to give P and OE annotations in dataframe (chr, start, end, annot)
  ## in list, see how dealt with name of element in list
  ## If short or long distances
{

  ## Setting pipe operator from magrittr package
  `%>%` <- magrittr::`%>%`

  ## Compute distance

  message("Computing Distance")

  interactions$distance <- GenomicInteractions::calculateDistances(interactions, method="midpoint")

  if(is.null(interactions$int))
  {
    stop("The interactions should have a 'int' column")
  }

  ## Computing short distances
  ## Compute total number of interactions by distance

  a1 <- dplyr::as_tibble(as.data.frame(interactions)[,c("distance","int")]) %>% dplyr::summarise("<100Kb"=sum(distance < 1*10^5, na.rm = T),
                                                                                       "100-200Kb"=sum(distance >= 1*10^5 & distance < 2*10^5, na.rm = T),
                                                                                       "200-300Kb"=sum(distance >= 2*10^5 & distance < 3*10^5, na.rm = T),
                                                                                       "300-400Kb"=sum(distance >= 3*10^5 & distance < 4*10^5, na.rm = T),
                                                                                       "400-500Kb"=sum(distance >= 4*10^5 & distance < 5*10^5, na.rm = T),
                                                                                       "500-600Kb"=sum(distance >= 5*10^5 & distance < 6*10^5, na.rm = T),
                                                                                       "600-700Kb"=sum(distance >= 6*10^5 & distance < 7*10^5, na.rm = T),
                                                                                       "700-800Kb"=sum(distance >= 7*10^5 & distance < 8*10^5, na.rm = T),
                                                                                       "800-900Kb"=sum(distance >= 8*10^5 & distance < 9*10^5, na.rm = T),
                                                                                       "900Kb-1Mb"=sum(distance >= 9*10^5 & distance < 1*10^6, na.rm = T),
                                                                                       ">1Mb"=sum(distance >= 1*10^6, na.rm = T),
                                                                                       # "Trans"=sum(is.na(distance)),
                                                                                       "int"="Total")

  ## Compute number of interactions by distance by interactions type

  a2 <- dplyr::as_tibble(as.data.frame(interactions)[,c("distance","int")]) %>% dplyr::group_by(int) %>% dplyr::summarise("<100Kb"=sum(distance < 1*10^5, na.rm = T),
                                                                                                                "100-200Kb"=sum(distance >= 1*10^5 & distance < 2*10^5, na.rm = T),
                                                                                                                "200-300Kb"=sum(distance >= 2*10^5 & distance < 3*10^5, na.rm = T),
                                                                                                                "300-400Kb"=sum(distance >= 3*10^5 & distance < 4*10^5, na.rm = T),
                                                                                                                "400-500Kb"=sum(distance >= 4*10^5 & distance < 5*10^5, na.rm = T),
                                                                                                                "500-600Kb"=sum(distance >= 5*10^5 & distance < 6*10^5, na.rm = T),
                                                                                                                "600-700Kb"=sum(distance >= 6*10^5 & distance < 7*10^5, na.rm = T),
                                                                                                                "700-800Kb"=sum(distance >= 7*10^5 & distance < 8*10^5, na.rm = T),
                                                                                                                "800-900Kb"=sum(distance >= 8*10^5 & distance < 9*10^5, na.rm = T),
                                                                                                                "900Kb-1Mb"=sum(distance >= 9*10^5 & distance < 1*10^6, na.rm = T),
                                                                                                                ">1Mb"=sum(distance >= 1*10^6, na.rm = T),
                                                                                                                # "Trans"=sum(is.na(distance)),
                                                                                                                "total_per_int"=dplyr::n())
  a <- dplyr::full_join(a1,a2, by=names(a1)) %>% tibble::add_column(cell_type=cell_type,
                                                                            ibed=length(interactions))

  short_int_dist_table <- dplyr::as_tibble(reshape2::melt(a, id.vars=c("int","total_per_int","cell_type","ibed")))


  ## Computing long distances

  ## Compute total number of interactions by distance

  a1 <- dplyr::as_tibble(as.data.frame(interactions)[,c("distance","int")]) %>% dplyr::summarise("<1Mb"=sum(distance < 1*10^6, na.rm = T),
                                                                                       "1-5Mb"=sum(distance >= 1*10^6 & distance < 5*10^6, na.rm = T),
                                                                                       "5-10Mb"= sum(distance >= 5*10^6 & distance < 10*10^6, na.rm = T),
                                                                                       "10-50Mb"=sum(distance >= 10*10^6 & distance < 50*10^6, na.rm = T),
                                                                                       "50-100Mb"=sum(distance >= 50*10^6 & distance < 100*10^6, na.rm = T),
                                                                                       ">100Mb"=sum(distance > 100*10^6, na.rm = T),
                                                                                       # "Trans"=sum(is.na(distance)),
                                                                                       "int"="Total")

  ## Compute number of interactions by distance by interactions type

  a2 <- dplyr::as_tibble(as.data.frame(interactions)[,c("distance","int")]) %>% dplyr::group_by(int) %>% dplyr::summarise("<1Mb"=sum(distance < 1*10^6, na.rm = T),
                                                                                                                "1-5Mb"=sum(distance >= 1*10^6 & distance < 5*10^6, na.rm = T),
                                                                                                                "5-10Mb"=sum(distance >= 5*10^6 & distance < 10*10^6, na.rm = T),
                                                                                                                "10-50Mb"=sum(distance >= 10*10^6 & distance < 50*10^6, na.rm = T),
                                                                                                                "50-100Mb"=sum(distance >= 50*10^6 & distance < 100*10^6, na.rm = T),
                                                                                                                ">100Mb"=sum(distance > 100*10^6, na.rm = T),
                                                                                                                # "Trans"=sum(is.na(distance)),
                                                                                                                "total_per_int"=n())
  a <- dplyr::full_join(a1,a2, by=names(a1)) %>% tibble::add_column(cell_type=cell_type,
                                                                            ibed=length(interactions))

  long_int_dist_table <- dplyr::as_tibble(reshape2::melt(a, id.vars=c("int","total_per_int","cell_type","ibed")))

  return(list(short_int_dist_table=short_int_dist_table,
              long_int_dist_table=long_int_dist_table))
}
