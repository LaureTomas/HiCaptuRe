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
#' @importFrom S4Vectors elementMetadata
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
  if (!(all(c("bait_1", "bait_2", "ID_1", "ID_2") %in% names(S4Vectors::elementMetadata(object))))) {
    stop("Valid HiCaptuRe object must contain bait_1, bait_2, ID_1, ID_2 columns")
  }
  return(TRUE)
})

setGeneric("HiCaptuRe", function(genomicInteractions, parameters, ByBaits, ByRegions) {
  standardGeneric("HiCaptuRe")
})

setMethod("HiCaptuRe", c("GenomicInteractions", "list", "list", "list"), function(genomicInteractions, parameters, ByBaits, ByRegions) {
  if (!inherits(genomicInteractions, "GenomicInteractions")) {
    stop("First argument must be a GenomicInteractions object")
  }
  new("HiCaptuRe", genomicInteractions, parameters = parameters, ByBaits = ByBaits, ByRegions = ByRegions)
})

setMethod("show", "HiCaptuRe", function(object) {
  callNextMethod()  # call inherited show from GInteractions

  cat("\n📦 Slots in HiCaptuRe object:\n")
  slots <- slotNames(object)[1:3]
  for (s in slots) {
    slot_value <- slot(object, s)

    if (tolower(s) == "parameters" && is.list(slot_value)) {
      param_names <- names(slot_value)
      cat(sprintf("  • @%s(%d)       : %s\n", s, length(param_names), paste(param_names, collapse = ", ")))
    } else if (tolower(s) == "bybaits" && is.list(slot_value)) {
      if (length(slot_value)==0) {
        cat("  • @ByBaits(0)          : NULL\n")
      } else{
        # Gather bait counts for all tibbles
        details <- vapply(seq_along(slot_value), function(i) {
          df <- slot_value[[i]]
          if (tibble::is_tibble(df)) {
            sprintf("[[%d]] %d baits", i, nrow(df))
          } else {
            sprintf("[[%d]] not a tibble", i)
          }
        }, character(1))

        cat(sprintf("  • @ByBaits(%d)          : %s\n", length(slot_value), paste(details, collapse = "; ")))
      }

    } else if (tolower(s) == "byregions" && is.list(slot_value)) {
      if (length(slot_value)==0) {
        cat("  • @ByRegions(0)        : NULL\n")
      } else{
        # Gather range counts for all GRanges
        details <- vapply(seq_along(slot_value), function(i) {
          gr <- slot_value[[i]]
          if (inherits(gr, "GRanges")) {
            sprintf("[[%d]] %d regions", i, length(gr))
          } else {
            sprintf("[[%d]] not a GRanges", i)
          }
        }, character(1))

        cat(sprintf("  • @ByRegions(%d)        : %s\n", length(slot_value), paste(details, collapse = "; ")))
      }
    } else {
      # Default: try to show length or class
      len <- tryCatch(length(slot_value), error = function(e) NA)
      cat(sprintf("  • @%-15s : %s\n", s,
                  if (!is.na(len)) paste0("length = ", len) else paste("class =", class(slot_value)[1])))
    }
  }
})
