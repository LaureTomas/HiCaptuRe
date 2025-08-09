library(testthat)
library(HiCaptuRe)

test_that("load_interactions returns a valid HiCaptuRe object", {
    x <- .setup_chr19()
    expect_s4_class(x, "HiCaptuRe")
    expect_true(length(x) > 0)
    expect_true("GenomicInteractions" %in% class(as(x, "GenomicInteractions")))
})

test_that("load_interactions errors on bad input", {
    bad <- tempfile(fileext = ".txt")
    writeLines("not\tvalid", bad)
    expect_error(load_interactions(bad), "Could not detect a valid format for the interaction file.")
})
