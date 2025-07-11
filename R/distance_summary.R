#' Computes the number of interactions by distance
#'
#' This function computes the number of interactions by distance for a given GenomicInteractions object
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param breaks vector with breaks for split distances
#' @param sample character variable with the name of the sample
#'
#' @return list with 2 tables: short_int_dist_table and long_int_dist_table, with the number of short interactions and long interactions respectivetly
#'
#' @importFrom dplyr as_tibble reframe group_by n full_join mutate arrange
#' @importFrom GenomicInteractions calculateDistances is.cis
#' @importFrom tibble add_column
#'
#' @examples
#' ibed1 <- system.file("extdata", "ibed1_example.zip", package = "HiCaptuRe")
#' interactions1 <- load_interactions(ibed1, select_chr = "19")
#' df <- distance_summary(interactions = interactions1)
#'
#' @export
distance_summary <- function(interactions, breaks = seq(0, 10^6, 10^5), sample = "sample") {
    ## Set breaks
    breaks <- unique(c(0, breaks, Inf))

    interactions <- interactions[GenomicInteractions::is.cis(interactions)]

    if (is.null(interactions$int)) {
        stop("The interactions should have a 'int' column")
    }

    total <- dplyr::as_tibble(interactions) |>
        select(distance, int) |>
        dplyr::mutate(breaks = cut(distance, breaks = breaks)) |>
        dplyr::group_by(breaks) |>
        dplyr::reframe(value = n()) |>
        tibble::add_column(int = "Total", total_per_int = NA, sample = sample)

    per_int <- dplyr::as_tibble(interactions) |>
        select(distance, int) |>
        dplyr::mutate(breaks = cut(distance, breaks = breaks)) |>
        dplyr::group_by(int, breaks) |>
        dplyr::reframe(value = n())

    total_int <- dplyr::as_tibble(interactions) |>
        select(distance, int) |>
        dplyr::group_by(int) |>
        dplyr::reframe(total_per_int = n(), sample = sample)
    per_int_df <- dplyr::full_join(per_int, total_int, by = "int")

    df <- rbind(
        total[, sort(colnames(total))],
        per_int_df[, sort(colnames(per_int_df))]
    )
    df$HiCaptuRe <- length(interactions)

    df <- df[, c("int", "total_per_int", "sample", "HiCaptuRe", "breaks", "value")] |> dplyr::arrange(breaks)

    return(df)
}
