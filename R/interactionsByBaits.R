#' Filters HiCaptuRe object by baits
#'
#' This function filters a HiCaptuRe object from load_interactions by a set of given bait(s)
#'
#' @param interactions a HiCaptuRe object
#' @param baits character vector containing bait names
#' @param sep character separating baits names when several baits in same fragment
#'
#' @return The returned object includes a filtered set of interactions and updates the ByBaits slot with a tibble summarizing bait-wise interaction statistics (number of interactions, other ends, etc.). Baits that are not present in the interaction data will have empty statistics.
#'
#' @importFrom stringr str_split
#' @importFrom dplyr as_tibble group_by summarise filter select rename arrange tibble
#' @importFrom tidyr separate_rows
#'
#' @export
interactionsByBaits <- function(interactions, baits, sep = ",") {
  `%>%` <- magrittr::`%>%`

  one <- sapply(stringr::str_split(as.character(interactions$bait_1), sep), function(x) any(x %in% baits))
  two <- sapply(stringr::str_split(as.character(interactions$bait_2), sep), function(x) any(x %in% baits))
  interactions_baits <- interactions[one | two]

  if (length(interactions_baits) == 0) {
    warning("There is no interaction including the given bait(s). Could Interactions be annotated with a different bait nomenclature?")
    baits_final <- tibble(ID = NA, bait = baits, N_int = NA, NOE = NA, interactingID = NA, interactingAnnotation = NA, interactingDistance = NA)
  } else {
    baits_df <- dplyr::as_tibble(interactions_baits) %>%
      tidyr::separate_rows(bait_1, sep = sep) %>%
      tidyr::separate_rows(bait_2, sep = sep)

    baits1 <- baits_df %>%
      dplyr::filter(bait_1 %in% baits) %>%
      dplyr::select(ID_1, ID_2, bait_1, bait_2, distance)

    baits2 <- baits_df %>%
      dplyr::filter(bait_2 %in% baits) %>%
      dplyr::select(ID_2, ID_1, bait_2, bait_1, distance)

    colnames(baits2) <- colnames(baits1)

    baits_final <- dplyr::as_tibble(rbind(baits1, baits2)) %>%
      dplyr::group_by(ID_1, bait_1) %>%
      dplyr::summarise(
        N_int = length(unique(ID_2)),
        NOE = sum(bait_2 == "."),
        interactingID = paste(unique(ID_2), collapse = ","),
        interactingAnnotation = paste(unique(bait_2), collapse = ","),
        interactingDistance = paste(distance, collapse = ",")
      ) %>%
      dplyr::rename(
        "bait" = "bait_1",
        "ID" = "ID_1"
      )

    if (any(!baits %in% baits_final$bait)) {
      missing_baits <- baits[!baits %in% baits_final$bait]
      missing <- data.frame(ID = NA, bait = missing_baits, N_int = 0, NOE = 0, interactingID = NA, interactingAnnotation = NA, interactingDistance = NA)
      baits_final <- dplyr::as_tibble(rbind(baits_final, missing)) %>% dplyr::arrange(ID)
    }
  }

  bait_list <- getByBaits(interactions_baits)

  if (length(bait_list) == 0) {
    bait_list[[1]] <- baits_final
  } else {
    bait_list[[length(bait_list) + 1]] <- baits_final
  }

  interactions_baits <- setByBaits(interactions_baits, bait_list)

  param <- getParameters(interactions_baits)
  param[[paste0("ByBaits_", length(bait_list))]] <- c(interactions = deparse(substitute(interactions)), baits = paste(baits, collapse = ","), sep = sep)
  interactions_baits <- setParameters(interactions_baits, param)


  return(interactions_baits)
}
