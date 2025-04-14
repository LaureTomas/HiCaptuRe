#' HiCaptuRe Class
#'
#' @description
#'
#' A S4 class to represent interactions between genomic regions
#'
#' @slot parameters List of parameters used to create the object and subsequence analysis
#' @slot ByBaits List of Baits used by \link[HiCaptuRe]{interactionsByBaits}.
#' @slot ByRegions List of Regions used by \link[HiCaptuRe]{interactionsByRegions}
#'
#' @note This class contains a \link[GenomicInteractions]{GenomicInteractions} object inside therefore all methods available to it can be used. This type of object should be created through the function \link[HiCaptuRe]{load_interactions}
#'
#' @importFrom methods new
#'
#' @exportClass HiCaptuRe
setClass("HiCaptuRe",
  contains = c("GenomicInteractions"),
  slots = c(
    parameters = "list",
    ByBaits = "list",
    ByRegions = "list"
  )
)

setValidity("HiCaptuRe", function(object) {
  if (!(all(c("bait_1", "bait_2", "ID_1", "ID_2") %in% names(object@elementMetadata)))) {
    stop("Valid HiCaptuRe object must contain bait_1, bait_2, ID_1, ID_2 columns")
  }
  return(TRUE)
})

setGeneric("HiCaptuRe", function(genomicInteractions, parameters, ByBaits, ByRegions) {
  standardGeneric("HiCaptuRe")
})

setMethod("HiCaptuRe", c("GenomicInteractions", "list", "list", "list"), function(genomicInteractions, parameters, ByBaits, ByRegions) {
  new("HiCaptuRe", genomicInteractions, parameters = parameters, ByBaits = ByBaits, ByRegions = ByRegions)
})

setGeneric("getParameters", function(x) {
  standardGeneric("getParameters")
})

setMethod("getParameters", c(x = "HiCaptuRe"), function(x) {
  return(x@parameters)
})

setGeneric("setParameters", function(x, y) {
  standardGeneric("setParameters")
})

setMethod("setParameters", c(x = "HiCaptuRe", y = "list"), function(x, y) {
  x@parameters <- y
  return(x)
})

setGeneric("getByBaits", function(x) {
  standardGeneric("getByBaits")
})

setMethod("getByBaits", c(x = "HiCaptuRe"), function(x) {
  return(x@ByBaits)
})

setGeneric("setByBaits", function(x, y) {
  standardGeneric("setByBaits")
})

setMethod("setByBaits", c(x = "HiCaptuRe", y = "list"), function(x, y) {
  x@ByBaits <- y
  return(x)
})


setGeneric("getByRegions", function(x) {
  standardGeneric("getByRegions")
})

setMethod("getByRegions", c(x = "HiCaptuRe"), function(x) {
  return(x@ByRegions)
})

setGeneric("setByRegions", function(x, y) {
  standardGeneric("setByRegions")
})

setMethod("setByRegions", c(x = "HiCaptuRe", y = "list"), function(x, y) {
  x@ByRegions <- y
  return(x)
})
