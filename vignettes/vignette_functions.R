## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----example_data-------------------------------------------------------------
ibed1 <- system.file("extdata", "ibed1_example.ibed", package="HiCaptuRe")
ibed2 <- system.file("extdata", "ibed2_example.ibed", package="HiCaptuRe")
annotation <- system.file("extdata", "annotation_example.txt", package="HiCaptuRe")
regions <- system.file("extdata", "regions_example.bed", package="HiCaptuRe")

## ----comment=F, message=F,library---------------------------------------------
library(HiCaptuRe)

## ----load_interactions--------------------------------------------------------
ibed1 <- load_interactions(file = ibed1)

## ----ibed1--------------------------------------------------------------------
ibed1

## ----slots--------------------------------------------------------------------
slotNames(ibed1)

## ----parameters---------------------------------------------------------------
ibed1@parameters

## ----annotate-----------------------------------------------------------------
ibed_annotated <- annotate_interactions(interactions = ibed1,
                                        annotation = annotation)
ibed_annotated

## ----annotate_parameters------------------------------------------------------
ibed_annotated@parameters

## ----interactionsByBaits------------------------------------------------------
baits_of_interest <- c("MTOR","APOA2","TP53")
ibed_byBaits <- interactionsByBaits(interactions = ibed_annotated,
                                    baits = baits_of_interest)

## ----interactionsByBaits_show-------------------------------------------------
ibed_byBaits

## ----ByBaits------------------------------------------------------------------
ibed_byBaits@ByBaits

## ----interactionsByRegions----------------------------------------------------
ibed_byRegions <- interactionsByRegions(interactions = ibed_annotated,
                                        regions = regions)

## ----interactionsByRegions_show-----------------------------------------------
ibed_byRegions

## ----byRegions----------------------------------------------------------------
ibed_byRegions@ByRegions

## ----ibed2--------------------------------------------------------------------
ibed2 <- load_interactions(file = ibed2)
ibed2_annotated <- annotate_interactions(interactions = ibed2, annotation = annotation)

## ----intersect_interactions---------------------------------------------------
interactions_list <- list(A = ibed_annotated, B = ibed2_annotated)
output <- intersect_interactions(interactions_list = interactions_list)

## ----intersections------------------------------------------------------------
lapply(output$intersections, function(x) x[1:2])

## ----upset, out.width="50%"---------------------------------------------------
output$upset

## ----venn, out.width="50%"----------------------------------------------------
output$venn

## ----distance_summary---------------------------------------------------------
dist_sum <- distance_summary(interactions = ibed_annotated,
                             breaks = seq(0, 10^6, 10^5),
                             sample = "ibed1")
dist_sum

## ----absolute, fig.show="hold", out.width="50%"-------------------------------
plots <- plot_distance_summary(distances = dist_sum, type_of_value = "absolute")
plots$int_dist
plots$total_dist

## ----by_int_type, out.width="50%"---------------------------------------------
plots <- plot_distance_summary(distances = dist_sum, type_of_value = "by_int_type")
plots$int_dist_norm_int


## ----by_total, fig.show="hold",  out.width="50%"------------------------------
plots <- plot_distance_summary(distances = dist_sum, type_of_value = "by_total")
plots$int_dist_norm_total
plots$total_dist_norm_total

## ----export, eval=F-----------------------------------------------------------
#  export_interactions(interactions = ibed_annotated,
#                      file = "/path/to/folder/ibed_annotated.ibed",
#                      type = "ibed")

