#' Functions to set data held in a HiCaptuRe object.
#'
#' Use these functions to set data stored in each of the slots of a
#' HiCapture object.
#'
#' @param x A HiCaptuRe object
#' @param y Data to add (parameters, tibble with bait-centric information or GRanges with region-centric)
#' @name setters
#'
#' @return A HiCaptuRe object
#'

## parameters

#' @rdname setters
#' @aliases setParameters
#' @export
setGeneric("setParameters", function(x, y) {
    standardGeneric("setParameters")
})
#' @rdname setters
#' @export
setMethod("setParameters", c(x = "HiCaptuRe", y = "list"), function(x, y) {
    x@parameters <- y
    return(x)
})

## byBaits

#' @rdname setters
#' @aliases setByBaits
#' @export
setGeneric("setByBaits", function(x, y) {
    standardGeneric("setByBaits")
})
#' @rdname setters
#' @export
setMethod("setByBaits", c(x = "HiCaptuRe", y = "list"), function(x, y) {
    x@ByBaits <- y
    return(x)
})

## byRegions

#' @rdname setters
#' @aliases setByRegions
#' @export
setGeneric("setByRegions", function(x, y) {
    standardGeneric("setByRegions")
})
#' @rdname setters
#' @export
setMethod("setByRegions", c(x = "HiCaptuRe", y = "list"), function(x, y) {
    x@ByRegions <- y
    return(x)
})
