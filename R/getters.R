#' Functions to access data held in a HiCaptuRe object.
#'
#' Use these functions to access data stored in each of the slots of a
#' HiCapture object.
#'
#' @name getters
#' @param HiCaptuRe A HiCaptuRe object
#' @rdname getters
#'
#' @return For 'getParameters', a named list of named vectors. For 'getByBaits', a list
#' of tibbles with bait-centric information. For 'getByRegions', a list of GRanges with
#' regions-centric information.
#'
#' @rdname getters
#' @export
setGeneric("getParameters", function(x) {
  standardGeneric("getParameters")
})
#' @rdname getters
#' @export
setMethod("getParameters", c(x = "HiCaptuRe"), function(x) {
  return(x@parameters)
})
#' @rdname getters
#' @export
setGeneric("getByBaits", function(x) {
  standardGeneric("getByBaits")
})
#' @rdname getters
#' @export
setMethod("getByBaits", c(x = "HiCaptuRe"), function(x) {
  return(x@ByBaits)
})
#' @rdname getters
#' @export
setGeneric("getByRegions", function(x) {
  standardGeneric("getByRegions")
})
#' @rdname getters
#' @export
setMethod("getByRegions", c(x = "HiCaptuRe"), function(x) {
  return(x@ByRegions)
})
#' @rdname getters
#' @export
setGeneric("getRegions", function(x) {
  standardGeneric("getRegions")
})
#' @rdname getters
#' @export
setMethod("getRegions", c(x = "HiCaptuRe"), function(x) {
  return(x@regions)
})
