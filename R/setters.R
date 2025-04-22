#' Functions to set data held in a HiCaptuRe object.
#'
#' Use these functions to set data stored in each of the slots of a
#' HiCapture object.
#'
#' @name setters
#' @param HiCaptuRe A HiCaptuRe object
#' @rdname setters
#'
#' @return A HiCaptuRe object
#'
#' @rdname setters
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
#' @rdname setters
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
#' @rdname setters
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
