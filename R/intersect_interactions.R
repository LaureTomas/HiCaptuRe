#' Computes the intersect between a list of interactions objects
#'
#' This function computes all the possible intersects for a given list of HiCaptuRe objects
#'
#' @param interactions_list list of HiCaptuRe objects from \code{\link{load_interactions}}
#' @param distance.boxplot logical, by default FALSE, plot a boxplot of log10 distance of each intersection in the upset plot
#' @param ... extra arguments for UpSetR
#'
#' @return A list with:
#' \describe{
#'   \item{intersections}{A named list of intersected HiCaptuRe objects}
#'   \item{upset}{An UpSetR plot}
#'   \item{venn}{A ggVennDiagram plot (NULL if > 7 sets)}
#' }
#' @importFrom dplyr as_tibble coalesce reframe select group_by
#' @importFrom stringr str_split
#' @importFrom GenomicInteractions GenomicInteractions
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom gplots venn
#' @importFrom UpSetR upset fromList
#' @importFrom S4Vectors elementMetadata
#' @importFrom ggVennDiagram ggVennDiagram
#'
#' @examples
#' ibed1 <- system.file("extdata", "ibed1_example.zip", package="HiCaptuRe")
#' interactions1 <- load_interactions(ibed1, select_chr = "19")
#' ibed2 <- system.file("extdata", "ibed2_example.zip", package="HiCaptuRe")
#' interactions2 <- load_interactions(ibed2, select_chr = "19")
#' interactions_list <- list(ibed1 = ibed1, ibed2 = ibed2)
#' intersection <- intersect_interactions(interactions_list = interactions_list)
#'
#' @export
intersect_interactions <- function(interactions_list, distance.boxplot = FALSE, ...) {
  param <- lapply(interactions_list, function(x) getParameters(x)$digest)

  are_identical <- all(sapply(param[-1], function(x) identical(param[[1]], x)))

  if (!are_identical) {
    stop("HiCaptuRe objects in the interactions_list have different digest genome. Check load_interactions() and digest_genome() arguments.")
  }

  check_annot <- lapply(interactions_list, function(x) {
    as_tibble(elementMetadata(x@regions)) |>
      mutate(name = dplyr::coalesce(unlist(B.id), unlist(OE.id))) |>
      select(node.class, fragmentID, name) |>
      unique()
  })
  check_annot <- do.call(rbind, check_annot) |> unique()

  check_annot2 <- check_annot |>
    dplyr::select(-name) |>
    unique() |>
    dplyr::group_by(fragmentID) |>
    dplyr::reframe(n = n())
  if (any(check_annot2$n > 1)) {
    stop("HiCaptuRe objects in the interactions_list contain fragments with different B/OE classification. Check annotate_interactions().")
  }
  check_annot3 <- check_annot |>
    unique() |>
    dplyr::group_by(fragmentID) |>
    dplyr::reframe(n = n())
  if (any(check_annot3$n > 1)) {
    warning("HiCaptuRe objects in the interactions_list contain fragments with different annotation. Check annotate_interactions().")
  }

  original <- names(interactions_list)
  new_order <- order(unlist(lapply(interactions_list, length)))

  interactions_list <- interactions_list[new_order]

  la <- list()

  if (distance.boxplot) {
    for (i in 1:length(interactions_list))
    {
      aa <- paste(S4Vectors::elementMetadata(interactions_list[[i]])[["ID_1"]], S4Vectors::elementMetadata(interactions_list[[i]])[["ID_2"]], sep = "_")
      names(aa) <- S4Vectors::elementMetadata(interactions_list[[i]])[["distance"]]
      la[[names(interactions_list)[i]]] <- aa
    }

    ## Based on the function FromList from UpSetR package
    ## Generate a dataframe with the intersection and the distance of each interaction
    elements1 <- unlist(la, use.names = TRUE)
    names(elements1) <- gsub("\\.", "", gsub(paste0(names(la), collapse = "|"), "", names(elements1)))
    elements <- unique(elements1)
    elements <- elements1[!duplicated(elements1)]
    data <- unlist(lapply(la, function(x) {
      x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)

    data <- data.frame(matrix(data, ncol = length(la), byrow = FALSE))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(la)

    data_final <- suppressWarnings(cbind(data, as.numeric(names(elements))))
    colnames(data_final)[ncol(data_final)] <- "distance"
    data_final$log10dist <- log10(data_final$distance)

    uplot <- suppressWarnings(UpSetR::upset(data_final, nsets = length(interactions_list), boxplot.summary = "log10dist", mainbar.y.label = "# Interactions", sets.x.label = "# Interactions", ...))
  } else {
    for (i in 1:length(interactions_list))
    {
      aa <- paste(S4Vectors::elementMetadata(interactions_list[[i]])[["ID_1"]], S4Vectors::elementMetadata(interactions_list[[i]])[["ID_2"]], sep = "_")
      la[[names(interactions_list)[i]]] <- aa
    }

    uplot <- suppressWarnings(UpSetR::upset(UpSetR::fromList(la), nsets = length(interactions_list), mainbar.y.label = "# Interactions", sets.x.label = "# Interactions", ...))
  }

  venn_plot <- gplots::venn(la, show.plot = FALSE)

  if (length(la) < 8) {
    ggvenn <- ggVennDiagram::ggVennDiagram(la, label_percent_digit = 2)
  } else {
    ggvenn <- NULL
  }

  int <- attr(venn_plot, "intersections")

  for (i in 1:length(int))
  {
    int_name <- stringr::str_split(names(int[i]), ":")[[1]]

    for (j in 1:length(int_name))
    {
      sample <- int_name[j]

      if (j == 1) {
        m <- match(int[[i]], paste(S4Vectors::elementMetadata(interactions_list[[sample]])[["ID_1"]], S4Vectors::elementMetadata(interactions_list[[sample]])[["ID_2"]], sep = "_"))
        ints <- interactions_list[[sample]][m]

        cols <- names(GenomicRanges::mcols(ints))
        initial <- grep("CS.*", cols)[1] - 1
        final <- which(cols %in% c("counts", "int", "distance"))

        CS <- grep("CS.*", cols)
        names(GenomicRanges::mcols(ints))[CS] <- paste(cols[CS], sample, sep = "_")

        GenomicRanges::mcols(ints) <- GenomicRanges::mcols(ints)[, c(1:initial, CS, final)]
      } else {
        m <- match(int[[i]], paste(S4Vectors::elementMetadata(interactions_list[[sample]])[["ID_1"]], S4Vectors::elementMetadata(interactions_list[[sample]])[["ID_2"]], sep = "_"))
        ints2 <- interactions_list[[sample]][m]

        CS <- grep("CS.*", names(GenomicRanges::mcols(ints2)))
        names(GenomicRanges::mcols(ints2))[CS] <- paste(names(GenomicRanges::mcols(ints2))[CS], sample, sep = "_")

        count <- grep("counts", names(GenomicRanges::mcols(ints)))

        S4Vectors::elementMetadata(ints) <- cbind(S4Vectors::elementMetadata(ints)[1:(count - 1)], S4Vectors::elementMetadata(ints2)[CS], S4Vectors::elementMetadata(ints)[count:ncol(GenomicRanges::mcols(ints))])
      }
    }

    if (length(int_name) > 1) {
      sub_original <- original[sort(match(int_name, original))]

      cols <- names(GenomicRanges::mcols(ints))
      initial <- grep("CS_", cols)[1] - 1
      count <- grep("counts", cols)

      cs <- c()
      for (sample in sub_original)
      {
        cs <- c(cs, grep(paste0("CS_.*", sample, "$"), cols))
      }

      S4Vectors::elementMetadata(ints) <- S4Vectors::elementMetadata(ints)[, c(1:initial, cs, count:length(cols))]

      int_name <- paste(sub_original, collapse = ":")
    }

    int[[i]] <- ints
    names(int)[i] <- int_name
  }

  return(list(intersections = int, upset = uplot, venn = ggvenn))
}
