#' Annotates GenomicInteractions object from load_interactions
#'
#' This function annotate a GenomicInteractions object from load_interactions, based on a given annotations file
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param annotation full path to annotations file or a dataframe with 5 colums: chr, start, end, fragmentID, annotation
#' @param ... additional parameters for fread
#'
#' @return HiCaptuRe object annotated, with columns bait_1 and bait_2 substituted based on the given annotation
#'
#' @importFrom GenomicInteractions annotateInteractions anchorOne anchorTwo resetAnnotations annotateRegions
#' @importFrom GenomicRanges makeGRangesFromDataFrame split
#' @importFrom data.table fread
#' @importFrom S4Vectors elementMetadata
#' @importFrom methods is
#'
#' @export
annotate_interactions <- function(interactions, annotation, ...) {
  if (methods::is(annotation, "character")) {
    annotation_file <- normalizePath(annotation)
    ## Reading annotation file
    annotation <- data.table::fread(annotation, stringsAsFactors = F, ...)
  }
  if (any(class(annotation) == "data.frame")) {
    annotation_file <- deparse(substitute(annotation))

    ## Checking if file has 5 columns
    if (ncol(annotation) != 5) {
      stop("File has not 5 columns \nAnnotation must be chr start end fragmentID annotation")
    } else {
      ## Making a GRanges from the annotation and setting the NA as non-annotated
      cn <- colnames(annotation)
      annotation[[cn[5]]][is.na(annotation[[cn[5]]])] <- "non-annotated"
      annotationGR <- GenomicRanges::makeGRangesFromDataFrame(annotation, seqnames.field = cn[1], start.field = cn[2], end.field = cn[3], keep.extra.columns = T)

      ## Creating annotation list
      annot <- GenomicRanges::split(annotationGR[, -1], as.factor(annotation[[cn[5]]]))
      annotation.features <- list(annot = annot)

      suppressMessages(GenomicInteractions::annotateInteractions(interactions, annotation.features))

      interactions@elementMetadata[, "bait_1"] <- unlist(GenomicInteractions::anchorOne(interactions)@elementMetadata[, "annot.id"])
      interactions@elementMetadata[, "bait_2"] <- unlist(GenomicInteractions::anchorTwo(interactions)@elementMetadata[, "annot.id"])

      ## When using a annotation file with only baits the OE get a NA so we have to change that
      interactions@elementMetadata[is.na(interactions@elementMetadata[, "bait_2"]), "bait_2"] <- "."
      interactions@elementMetadata[is.na(interactions@elementMetadata[, "bait_1"]), "bait_1"] <- "."

      interactions <- annotate_BOE(interactions)

      cond <- ((interactions$ID_1 > interactions$ID_2) & interactions$int == "B_B") | ((interactions$ID_1 < interactions$ID_2) & interactions$int == "OE_B")

      a1 <- interactions@anchor1[cond]
      a2 <- interactions@anchor2[cond]

      interactions@anchor1[cond] <- a2
      interactions@anchor2[cond] <- a1

      cols <- sort(grep("_", colnames(S4Vectors::elementMetadata(interactions[cond]))[1:7], value = T))
      S4Vectors::elementMetadata(interactions[cond])[cols] <- S4Vectors::elementMetadata(interactions[cond])[cols[c(rbind(seq(2, length(cols), 2), seq(1, length(cols), 2)))]]

      interactions <- annotate_BOE(interactions)

      param <- getParameters(interactions)
      param$annotate <- c(annotation = annotation_file)
      interactions <- setParameters(interactions, param)
      interactions <- interactions[order(interactions$ID_1, interactions$ID_2)]

      return(interactions)
    }
  }
}


annotate_BOE <- function(interactions) {
  a1 <- GenomicInteractions::anchorOne(interactions)
  a1$bait_1 <- interactions@elementMetadata[, "bait_1"]
  a2 <- GenomicInteractions::anchorTwo(interactions)
  a2$bait_1 <- interactions@elementMetadata[, "bait_2"]

  a2$bait_1[is.na(a2$bait_1)] <- "."

  regions <- unique(c(a1, a2))

  regions$bait_1[is.na(regions$bait_1)] <- "non-annotated"
  bait <- regions[(regions$bait_1 != ".")]
  oe <- regions[regions$bait_1 == "."]

  baitl <- GenomicRanges::split(bait[, -1], as.factor(bait$bait_1))
  oel <- GenomicRanges::split(oe[, -1], as.factor(oe$bait_1))

  annotation.features <- list(B = baitl, OE = oel)
  GenomicInteractions::resetAnnotations(interactions)
  suppressMessages(GenomicInteractions::annotateInteractions(interactions, annotation.features))
  GenomicInteractions::annotateRegions(interactions, "fragmentID", unique(sort(c(interactions$ID_1, interactions$ID_2))))
  interactions$int <- paste(GenomicInteractions::anchorOne(interactions)$node.class, GenomicInteractions::anchorTwo(interactions)$node.class, sep = "_")

  return(interactions)
}
