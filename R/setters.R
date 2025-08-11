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
#' @keywords internal


## parameters

#' @rdname setters
#' @aliases .setParameters
setGeneric(".setParameters", function(x, y) {
    standardGeneric(".setParameters")
})
#' @rdname setters
setMethod(".setParameters", c(x = "HiCaptuRe", y = "list"), function(x, y) .set_slot(x, "parameters", y))

## byBaits

#' @rdname setters
#' @aliases .setByBaits
setGeneric(".setByBaits", function(x, y) {
    standardGeneric(".setByBaits")
})
#' @rdname setters
setMethod(".setByBaits", c(x = "HiCaptuRe", y = "list"), function(x, y) .set_slot(x, "ByBaits", y))

## byRegions

#' @rdname setters
#' @aliases .setByRegions
setGeneric(".setByRegions", function(x, y) {
    standardGeneric(".setByRegions")
})
#' @rdname setters
setMethod(".setByRegions", c(x = "HiCaptuRe", y = "list"), function(x, y) .set_slot(x, "ByRegions", y))


#' @keywords internal
.set_slot <- function(x, slot, value) {
  methods::slot(x, slot) <- value
  methods::validObject(x)
  x
}
