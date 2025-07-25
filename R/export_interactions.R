#' Export interactions in the desired output format
#'
#' This function exports interactions in different formats
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param file full path to desired output file (ibed, peakmatrix, washU, washUold, cytoscape, bedpe)
#' @param format type of output format (ibed, peakmatrix, washU, washUold, cytoscape, bedpe, seqmonk)
#' @param over.write TRUE/FALSE to over write the output file
#' @param cutoff Chicago score cutoff to export interactions
#' @param parameters TRUE/FALSE to also export the parameters of the given object

#' @return tibble object with the ibed table and save it in the desired output file
#'
#' @importFrom GenomicInteractions export.igraph
#' @importFrom dplyr as_tibble arrange bind_rows
#' @importFrom data.table fwrite
#' @importFrom igraph simplify as_edgelist
#' @importFrom stats setNames
#' @importFrom GenomeInfoDb seqlevelsStyle
#'
#'
#' @examples
#'
#' ibed1 <- system.file("extdata", "ibed1_example.zip", package = "HiCaptuRe")
#' interactions <- load_interactions(ibed1, select_chr = "19")
#' export_interactions(interactions = interactions, file = tempfile(), format = "ibed", over.write = TRUE)
#'
#' @export
export_interactions <- function(interactions, file, format = "ibed", over.write = FALSE, cutoff = 5, parameters = FALSE) {
    format <- match.arg(arg = format, choices = c("ibed", "peakmatrix", "washU", "washUold", "cytoscape", "bedpe", "seqmonk"), several.ok = FALSE)
    if (file.exists(file) & !over.write) {
        stop("File already exists. Use `over.write = TRUE` to overwrite.")
    }

    if (parameters) {
        .export_parameters(interactions, file)
    }

    params <- getParameters(interactions)
    type <- params$load[["format"]]
    pk <- params$peakmatrix2list

    if (format == "peakmatrix" & (type == "peakmatrix" & is.null(pk))) {
        m <- GenomicRanges::mcols(interactions)
        CS_m <- m[, grep("CS_", names(m))]
        interactions <- interactions[apply(CS_m, 1, function(x) any(x >= cutoff))]

        int_df <- dplyr::as_tibble(interactions)

        int_df <- int_df[, c(
            "seqnames1", "start1", "end1", "ID_1", "bait_1",
            "seqnames2", "start2", "end2", "ID_2", "bait_2",
            "distance", grep("CS_", colnames(int_df), value = TRUE)
        )]
        colnames(int_df) <- c(
            "baitChr", "baitStart", "baitEnd", "baitID", "baitName",
            "oeChr", "oeStart", "oeEnd", "oeID", "oeName",
            "dist", gsub("CS_", "", colnames(int_df)[12:ncol(int_df)])
        )
        data.table::fwrite(int_df, file = file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    } else if (format != "peakmatrix" & type == "peakmatrix") {
        if (!is.null(pk)) {
            if (is.list(interactions)) {
                int_list <- interactions
                files <- paste0(sub("\\.[^\\.]*$", "", file), "_", names(int_list), gsub(".*\\.", ".", basename(file)))
                export_function <- .export_dispatch(format)
                invisible(mapply(export_function, int_list, files))
            } else {
                interactions <- interactions[interactions$CS >= cutoff]
                export_function <- .export_dispatch(format)
                export_function(interactions, file)
            }
        } else {
            warning("Interactions input appears to be a peakmatrix. Automatically splitting by sample with `peakmatrix2list()`. And exporting as individual files")
            int_list <- peakmatrix2list(peakmatrix = interactions, cutoff = cutoff)
            files <- paste0(sub("\\.[^\\.]*$", "", file), "_", names(int_list), gsub(".*\\.", ".", basename(file)))
            export_function <- .export_dispatch(format)
            invisible(mapply(export_function, int_list, files))
        }
    } else {
        interactions <- interactions[interactions$CS >= cutoff]
        export_function <- .export_dispatch(format)
        export_function(interactions, file)
    }
}

.export_dispatch <- function(format) {
    switch(format,
        ibed = .export_ibed,
        washU = .export_washU,
        washUold = .export_washUold,
        cytoscape = .export_citoscape,
        bedpe = .export_bedpe,
        seqmonk = .export_seqmonk
    )
}


.export_parameters <- function(interactions, file) {
    file <- paste0(file, ".parameters")
    param <- getParameters(interactions)

    for (i in 1:length(param))
    {
        df <- as.data.frame(param[[i]])
        colnames(df) <- paste("#", toupper(names(param)[i]), "#")
        df <- rbind(df, "")
        rownames(df)[nrow(df)] <- ""
        if (i == 1) {
            data.table::fwrite(df, file = file, quote = FALSE, sep = "\t", col.names = TRUE)
        } else {
            suppressWarnings(data.table::fwrite(df, file = file, quote = FALSE, sep = "\t", col.names = TRUE, append = TRUE))
        }
    }
}


.export_ibed <- function(ints, file) {
    int_df <- dplyr::as_tibble(ints)[, c(
        "seqnames1", "start1", "end1", "bait_1",
        "seqnames2", "start2", "end2", "bait_2",
        "reads", "CS"
    )]
    colnames(int_df) <- c(
        "bait_chr", "bait_start", "bait_end", "bait_name",
        "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
        "N_reads", "score"
    )
    data.table::fwrite(int_df, file = file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}

.export_washU <- function(ints, file) {
    GenomeInfoDb::seqlevelsStyle(ints) <- "UCSC"
    int_df <- dplyr::as_tibble(ints)
    int_df <- dplyr::arrange(int_df, seqnames1, start1, end1, seqnames2, start2, end2)
    washU <- data.frame(
        paste(int_df$seqnames1, int_df$start1, int_df$end1, sep = "\t"),
        paste(int_df$seqnames2, ":", int_df$start2, "-", int_df$end2, ",", int_df$CS, sep = "")
    )
    colnames(washU) <- c("regionI", "regionIICS")
    data.table::fwrite(washU, file = file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

.export_washUold <- function(ints, file) {
    GenomeInfoDb::seqlevelsStyle(ints) <- "UCSC"
    int_df <- as.data.frame(ints)
    washU <- data.frame(
        paste(int_df$seqnames1, ":", int_df$start1, ",", int_df$end1, sep = ""),
        paste(int_df$seqnames2, ":", int_df$start2, ",", int_df$end2, sep = ""),
        int_df$CS
    )
    colnames(washU) <- c("regionI", "regionII", "CS")
    data.table::fwrite(washU, file = file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

.export_bedpe <- function(ints, file) {
    if (length(ints) == 0) {
        data.table::fwrite(data.frame(), file = file, col.names = FALSE, row.names = FALSE)
        return(invisible(NULL))
    }
    ints$name <- 1:length(ints)
    int_df <- dplyr::as_tibble(ints)[, c(
        "seqnames1", "start1", "end1",
        "seqnames2", "start2", "end2",
        "name", "CS", "strand1", "strand2"
    )]

    data.table::fwrite(int_df, file = file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

.export_citoscape <- function(ints, file) {
    net <- igraph::simplify(export.igraph(ints))
    nodes_edges <- igraph::as_edgelist(net)
    data.table::fwrite(nodes_edges, file = file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

.export_seqmonk <- function(ints, file) {
    if (length(ints) == 0) {
        data.table::fwrite(data.frame(), file = file, col.names = FALSE, row.names = FALSE)
        return(invisible(NULL))
    }
    ints$ID <- 1:length(ints)
    df1 <- dplyr::as_tibble(ints)[, c(
        "seqnames2", "start2", "end2", "bait_2",
        "reads", "CS", "ID"
    )]

    df2 <- dplyr::as_tibble(ints)[, c(
        "seqnames1", "start1", "end1", "bait_1",
        "reads", "CS", "ID"
    )]

    seqmonk <- dplyr::bind_rows(df1, stats::setNames(df2, names(df1))) |> dplyr::arrange(ID)
    data.table::fwrite(seqmonk[, -ncol(seqmonk)], file = file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}
