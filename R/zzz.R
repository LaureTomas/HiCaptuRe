.onLoad <- function(libname, pkgname) {
  ## This line prevents data.table and tidyverse syntax being flagged up as a NOTE in R CMD check
  ## ("no visible binding for global variable")
  utils::globalVariables(c(
    ".SD", "B.id", "OE.id", "node.class", "ID", "ID_1", "ID_2", "annot", "bait_1", "bait_2", "breaks", "chr_1", "chr_2", "distance",
    "end1", "end2", "end_1", "end_2", "fragmentID", "int", "regionID", "seqnames1", "seqnames2",
    "start1", "start2", "start_1", "start_2", "total_per_int", "value", "regionID_1", "regionID_2", "name", "score_names"
  ))
  digest_genome <<- memoise::memoise(digest_genome)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(hicapture_attach_message())
}
