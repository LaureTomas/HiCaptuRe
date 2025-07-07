library(testthat)
library(HiCaptuRe)

test_that("digest_genome returns GRanges fragments", {
  digest <- digest_genome(genome = "BSgenome.Hsapiens.NCBI.GRCh38", RE_name = "HindIII")
  expect_true(all(names(digest) == c("digest","parameters","seqinfo")))
  expect_equal(class(digest$digest), "data.frame")
  expect_true(all(c("seqnames", "start", "end") %in% colnames(digest$digest)))
})

test_that("digest_genome missing chromosomes", {
  expect_error(digest_genome(genome = "BSgenome.Hsapiens.NCBI.GRCh38", RE_name = "HindIII", select_chr = "chr1"),
               "Chromosomes selected in 'select_chr' are not present in this genome. Please try changing 'select_chr' or setting it to 'NULL'")
})
