#' Column name constants for format detection and export
#'
#' `.IBED_COLS` and `.PEAKMATRIX_COLS` store canonical column names for ibed and peakmatrix formats.
#' These are used in both `.detect_format()` and export functions to avoid repetition and ensure
#' consistency between import and export.
#'
#' @name constants
#'
#' @keywords internal
.IBED_COLS <- c(
    "bait_chr", "bait_start", "bait_end", "bait_name",
    "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
    "N_reads", "score"
)

.PEAKMATRIX_COLS <- c(
    "baitChr", "baitStart", "baitEnd", "baitID", "baitName",
    "oeChr", "oeStart", "oeEnd", "oeID", "oeName", "dist"
)
