#' Interanl functions to access data held in a HiCaptuRe object.
#'
#' Use these functions to access data stored in each of the slots of a
#' HiCapture object.
#'
#' @param x A HiCaptuRe object
#' @name getters
#'
#' @return For 'getParameters', a named list of named vectors. For 'getByBaits', a list
#' of tibbles with bait-centric information. For 'getByRegions', a list of GRanges with
#' regions-centric information.
#'
#' @examples
#' ibed1 <- system.file("extdata", "ibed1_example.zip", package = "HiCaptuRe")
#' interactions <- load_interactions(ibed1, select_chr = "19")
#' getParameters(interactions)
#'
#' baits <- c("ENST00000332235", "ENST00000516525")
#' interactions_baits <- interactionsByBaits(interactions = interactions, baits = baits)
#' getByBaits(interactions_baits)
#'
#' regions <- GenomicRanges::GRanges(seqnames = 19, ranges = IRanges::IRanges(start = c(500000, 1000000), end = c(510000, 1100000)))
#' interactions_regions <- interactionsByRegions(interactions = interactions, regions = regions)
#' getByRegions(interactions_regions)
#'

## parameters

#' @rdname getters
#' @aliases getParameters
#' @export
setGeneric("getParameters", function(x) {
    standardGeneric("getParameters")
})
#' @rdname getters
#'
#' @export
setMethod("getParameters", c(x = "HiCaptuRe"), function(x) {
    return(x@parameters)
})

## byBaits

#' @rdname getters
#' @aliases getByBaits
#' @export
setGeneric("getByBaits", function(x) {
    standardGeneric("getByBaits")
})
#' @rdname getters
#'
#' @export
setMethod("getByBaits", c(x = "HiCaptuRe"), function(x) {
    return(x@ByBaits)
})

## byRegions

#' @rdname getters
#' @aliases getByRegions
#' @export
setGeneric("getByRegions", function(x) {
    standardGeneric("getByRegions")
})
#' @rdname getters
#'
#' @export
setMethod("getByRegions", c(x = "HiCaptuRe"), function(x) {
    return(x@ByRegions)
})
