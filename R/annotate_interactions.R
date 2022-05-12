#' Annotates GenomicInteractions object from load_interactions
#'
#' This function annotate a GenomicInteractions object from load_interactions, based on a given annotations file
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param annotation full path to annotations file
#' @param ... additional parameters for fread
#'
#' @return GenomicInteractions object annotated
#'
#' @importFrom magrittr `%>%`
#' @importFrom GenomicInteractions annotateInteractions anchorOne anchorTwo resetAnnotations
#' @importFrom GenomicRanges makeGRangesFromDataFrame split
#' @importFrom ddpcr quiet
#' @importFrom data.table fread
#'
#' @export
annotate_interactions <- function(interactions, annotation,...)
{

  ## Reading annotation file
  hind <- data.table::fread(annotation, stringsAsFactors = F,...)

  ## Checking if file has 5 columns
  if (ncol(hind)!=5)
  {
    stop("File has not 5 columns \nAnnotation must be chr start end fragID annotation")
  }
  else
  {
    message("Remember: Annotation must be chr start end fragID annotation")

    ## Making a GRanges from the annotation and setting the NA as non-annotated
    cn <- colnames(hind)
    hind[[cn[5]]][is.na(hind[[cn[5]]])] <- "non-annotated"
    hindGR <- GenomicRanges::makeGRangesFromDataFrame(hind, seqnames.field = cn[1], start.field = cn[2], end.field = cn[3], keep.extra.columns = T)

    ## Creating annotation list
    annot <- GenomicRanges::split(hindGR[,-1], as.factor(hind[[cn[5]]]))
    annotation.features = list(annot=annot)

    ddpcr::quiet(GenomicInteractions::annotateInteractions(interactions, annotation.features))

    interactions@elementMetadata[,"gene_I"] <- unlist(GenomicInteractions::anchorOne(interactions)@elementMetadata[,"annot.id"])
    interactions@elementMetadata[,"gene_II"] <- unlist(GenomicInteractions::anchorTwo(interactions)@elementMetadata[,"annot.id"])

    ## When using a annotation file with only baits the OE get a NA so we have to change that
    interactions@elementMetadata[is.na(interactions@elementMetadata[,"gene_II"]),"gene_II"] <- "."

    interactions <- annotate_POEuce(interactions)

    return(interactions)
  }

}


annotate_POEuce <- function(interactions)
{
  a1 <- GenomicInteractions::anchorOne(interactions)
  a1$gene_I <- interactions@elementMetadata[,"gene_I"]
  a2 <- GenomicInteractions::anchorTwo(interactions)
  a2$gene_I <- interactions@elementMetadata[,"gene_II"]

  a2$gene_I[is.na(a2$gene_I)] <- "."

  regions <- unique(c(a1,a2))

  regions$gene_I[is.na(regions$gene_I)] <- "non-annotated"
  prom <- regions[(regions$gene_I != ".") & (regions$gene_I != "uce")]
  oe <- regions[regions$gene_I == "."]
  uce <- regions[regions$gene_I == "uce"]

  proml <- GenomicRanges::split(prom[, -1], as.factor(prom$gene_I))
  oel <- GenomicRanges::split(oe[, -1], as.factor(oe$gene_I))
  ucel <- GenomicRanges::split(uce[, -1], as.factor(uce$gene_I))

  annotation.features = list(P = proml, OE = oel, uce = ucel)
  GenomicInteractions::resetAnnotations(interactions)
  suppressMessages(GenomicInteractions::annotateInteractions(interactions, annotation.features))
  interactions$int <- paste(GenomicInteractions::anchorOne(interactions)$node.class,GenomicInteractions::anchorTwo(interactions)$node.class, sep = "_")

  return(interactions)
}

