library(testthat)
library(HiCaptuRe)

test_that("load_interactions returns a valid HiCaptuRe object", {
    ibed1_file <- system.file("extdata", "ibed1_example.zip", package = "HiCaptuRe")
    ibed1 <- load_interactions(file = ibed1_file, genome = "BSgenome.Hsapiens.NCBI.GRCh38")
    expect_s4_class(ibed1, "HiCaptuRe")
    expect_true(length(ibed1) > 0)
    expect_true("GenomicInteractions" %in% class(as(ibed1, "GenomicInteractions")))
})

test_that("load_interactions errors on bad input", {
    bad <- tempfile(fileext = ".txt")
    writeLines("not\tvalid", bad)
    expect_error(load_interactions(bad), "Could not detect a valid format for the interaction file.")
})
