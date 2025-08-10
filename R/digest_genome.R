#' Digest genome for a specific restriction enzyme
#'
#' This function takes a genome installed and generates its digest for a given restriction enzyme
#'
#' @param genome character with the name of the genome version
#' @param RE_name name of the restriction enzyme
#' @param motif recognition motif of the restriction enzyme
#' @param cut_position cut position of the restriction enzyme inside the motif
#' @param select_chr a character vector containing the specific chromosomes to used from this genome, if NULL all chromosomes will be used
#' @param PAR_mask a logical value where the Y chromosome should exclude the Pseudoautosomical Regions (PAR) or not
#' @param PAR_file a full path to a file containing the coordinates of Y chromosome PAR with at least 3 colums with header: seqnames, start, end
#' @param ... extra arguments for read.table
#'
#' @return list object with 2 elements: a dataframe with the digested genome, and a the parameters used for the digestion.
#'
#' @note The package provides for a PAR coordinates file only for Homo sapiens for the genome version 38
#' @note The package provides the motives and cut positions for several restriction enzymes (HindII, MboI, DpnII, EcoRI, BamHI)
#'
#' @importFrom Biostrings replaceAt DNAStringSet matchPattern
#' @importFrom BSgenome getBSgenome
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevelsStyle seqnames seqlevels bsgenomeName
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom data.table rbindlist
#'
#'
#' @examples
#'
#' digest <- digest_genome(genome = "GRCh38", RE_name = "HindIII", select_chr = "19")
#'
#' @export
digest_genome <- function(genome = "GRCh38", RE_name = "HindIII", motif = NULL, cut_position = NULL, select_chr = c(seq_len(22), "X", "Y"), PAR_mask = TRUE, PAR_file = NULL, ...) {
    # Define enzyme database
    enzyme_db <- list(
        HindIII = list(motif = "AAGCTT", cut_position = 1),
        EcoRI   = list(motif = "GAATTC", cut_position = 1),
        BamHI   = list(motif = "GGATCC", cut_position = 1),
        MboI    = list(motif = "GATC", cut_position = 0),
        DpnII   = list(motif = "GATC", cut_position = 0)
    )

    # If user provides RE_name but no motif/cut_position, fill them in
    if (!is.null(RE_name)) {
        RE_name_upper <- toupper(RE_name)
        matched_enzyme <- names(enzyme_db)[toupper(names(enzyme_db)) == RE_name_upper]

        if (length(matched_enzyme) == 1) {
            enzyme <- enzyme_db[[matched_enzyme]]

            # If motif/cut_position not provided, fill them in
            if (is.null(motif)) {
                motif <- enzyme$motif
            } else if (motif != enzyme$motif) {
                warning(sprintf("Provided motif (%s) differs from known motif (%s) for enzyme %s.", motif, enzyme$motif, matched_enzyme))
            }

            if (is.null(cut_position)) {
                cut_position <- enzyme$cut_position
            } else if (cut_position != enzyme$cut_position) {
                warning(sprintf("Provided cut_position (%s) differs from known cut_position (%s) for enzyme %s.", cut_position, enzyme$cut_position, matched_enzyme))
            }
        } else {
            stop("Unknown restriction enzyme: ", RE_name)
        }
    }

    # Safety check
    if (is.null(motif) || is.null(cut_position)) {
        stop("Please provide both 'motif' and 'cut_position' or a valid 'RE_name'")
    }

    ## get genome
    genome_name <- genome
    genome <- suppressMessages(BSgenome::getBSgenome(genome))

    ## If Pseudoautosomal Regions (PAR) masked, read provided file for human
    if (PAR_mask) {
        if (is.null(PAR_file)) {
            PAR_file <- system.file("extdata", paste0("PAR_", gsub(" ", "_", metadata(genome)$organism), "_coordinates.txt"), package = "HiCaptuRe")
            if (file.exists(PAR_file)) {
                PAR <- read.delim(PAR_file, header = TRUE, ...)
            } else {
                stop_message <- paste("There is no PAR file provided for", metadata(genome)$organism, "\n It must be a headed file with seqnames,start,end columns")
                stop(stop_message)
            }
        } else {
            if (file.exists(PAR_file)) {
                PAR <- read.delim(PAR_file, header = TRUE, ...)
            } else {
                stop("The PAR file provided doesn't exist")
            }
        }
        # Harmonize styles
        PARGR <- GenomicRanges::makeGRangesFromDataFrame(PAR)

        if (!any(GenomeInfoDb::seqlevelsStyle(PARGR) %in% GenomeInfoDb::seqlevelsStyle(genome))) {
            warning("seqlevelsStyle of PAR changed to match seqlevelsStyle(genome)")
            GenomeInfoDb::seqlevelsStyle(PARGR) <- GenomeInfoDb::seqlevelsStyle(genome)
        }
    } ## PAR mask



    ## Select primary chromosomes
    chrs <- GenomeInfoDb::seqnames(genome)

    if (!is.null(select_chr)) {
        chrs <- chrs[chrs %in% select_chr]
    }
    if (length(chrs) == 0) {
        stop("Chromosomes selected in 'select_chr' are not present in this genome. Please try changing 'select_chr' or setting it to 'NULL'")
    }

    ## Digest genome by chromosomes
    digest <- lapply(chrs, function(chr) {
      chr_seq <- genome[[chr]]

      # Mask PAR regions (if requested and available for this chr)
      if (isTRUE(PAR_mask) && chr %in% GenomeInfoDb::seqlevels(PARGR)) {
        chr_PAR <- IRanges::subsetByOverlaps(
          PARGR,
          GenomicRanges::GRanges(chr, IRanges::IRanges(1, length(chr_seq)))
        )
        if (length(chr_PAR) > 0) {
          widths    <- IRanges::width(chr_PAR)
          mask_seqs <- Biostrings::DNAStringSet(strrep("N", widths))
          chr_seq   <- Biostrings::replaceAt(chr_seq, GenomicRanges::ranges(chr_PAR), mask_seqs)
        }
      }

      m <- Biostrings::matchPattern(motif, chr_seq)

      if (length(m) == 0) {
        starts <- 1
        ends   <- length(chr_seq)
      } else {
        correct <- IRanges::start(m) - 1 + cut_position
        starts <- c(1, correct + 1)
        ends   <- c(correct, length(chr_seq))

      }

      data.frame(
        seqnames = chr,
        start = starts,
        end   = ends,
        row.names = NULL,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }) |> data.table::rbindlist()

    if (!PAR_mask) {
        PAR_file <- "NULL"
    } else {
        PAR_file <- normalizePath(PAR_file)
    }

    digest$fragment_ID <- seq_len(nrow(digest))
    output <- list(
        digest = digest,
        parameters = c(
            "Genome" = genome_name,
            "Genome_Package" = GenomeInfoDb::bsgenomeName(genome),
            "Restriction_Enzyme" = RE_name,
            "Motif" = motif,
            "Cut_Position" = cut_position,
            "Selected_Chromosomes" = paste(select_chr, collapse = ","),
            "PAR_mask" = PAR_mask,
            "PAR_file" = PAR_file
        ),
        seqinfo = seqinfo(genome)[chrs]
    )
    return(output)
}
