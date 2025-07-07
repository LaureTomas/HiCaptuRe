library(testthat)
library(HiCaptuRe)
library(GenomicInteractions)
library(GenomicRanges)


test_that("HiCaptuRe constructor returns valid S4 object", {
  anchor1 <- GRanges(seqnames = "chr1", ranges = IRanges(1, 100), bait_1 = "A", ID_1 = 1)
  anchor2 <- GRanges(seqnames = "chr1", ranges = IRanges(200, 300), bait_2 = ".", ID_2 = 2, read = 1, CS = 5)

  gi <- GenomicInteractions(anchor1, anchor2)
  names(GenomicRanges::mcols(gi)) <- gsub(x = names(GenomicRanges::mcols(gi)), pattern = "anchor[1-2]\\.", "")

  hicapture <- HiCaptuRe:::HiCaptuRe(gi, parameters = list(digest = c(test="test"), load = c(file="test")), ByBaits = list(), ByRegions = list())

  expect_s4_class(hicapture, "HiCaptuRe")
  expect_s4_class(as(hicapture, "GenomicInteractions"), "GenomicInteractions")
  expect_equal(length(hicapture), 1)
  expect_true("counts" %in% colnames(mcols(hicapture)))
})

test_that("HiCaptuRe inherits required methods and accessors", {
  anchor1 <- GRanges(seqnames = "chr1", ranges = IRanges(1, 100), bait_1 = "A", ID_1 = 1)
  anchor2 <- GRanges(seqnames = "chr1", ranges = IRanges(200, 300), bait_2 = ".", ID_2 = 2, read = 1, CS = 5)

  gi <- GenomicInteractions(anchor1, anchor2)
  names(GenomicRanges::mcols(gi)) <- gsub(x = names(GenomicRanges::mcols(gi)), pattern = "anchor[1-2]\\.", "")

  hicapture <- HiCaptuRe:::HiCaptuRe(gi, parameters = list(digest = c(test="test"), load = c(file="test")), ByBaits = list(), ByRegions = list())

  expect_equal(as.character(seqnames(anchorOne(hicapture))[1]), "chr1")
  expect_equal(width(anchorOne(hicapture))[1], 100)
  expect_true(is.numeric(mcols(hicapture)$counts))
})

test_that("HiCaptuRe errors with invalid input", {
  anchor1 <- GRanges(seqnames = "chr1", ranges = IRanges(1, 100), )
  anchor2 <- GRanges(seqnames = "chr1", ranges = IRanges(200, 300), )

  badgi <- GenomicInteractions(anchor1, anchor2)
  names(GenomicRanges::mcols(badgi)) <- gsub(x = names(GenomicRanges::mcols(badgi)), pattern = "anchor[1-2]\\.", "")
  expect_error(HiCaptuRe:::HiCaptuRe(badgi, parameters = list(digest = c(test="test"), load = c(file="test")), ByBaits = list(), ByRegions = list()), "Valid HiCaptuRe object must contain bait_1, bait_2, ID_1, ID_2 columns")
})
