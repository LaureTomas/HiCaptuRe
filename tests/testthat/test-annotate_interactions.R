library(testthat)
library(HiCaptuRe)

test_that("annotate_interactions: data.frame input produces bait_* and preserves class/length", {
  x <- .setup_chr19()

  ann_df <- data.frame(
    chr = "19",
    start = 1L,
    end = 300000000L,
    fragmentID = "fragA",
    annotation = NA_character_,
    stringsAsFactors = FALSE
  )

  y <- annotate_interactions(x, ann_df)

  expect_s4_class(y, "HiCaptuRe")
  expect_identical(length(y), length(x))

  expect_true(all(c("bait_1", "bait_2") %in% colnames(S4Vectors::mcols(y))))
  expect_false(any(is.na(y$bait_1)))
  expect_false(any(is.na(y$bait_2)))

})

test_that("annotate_interactions: character file path input behaves like data.frame", {
  x <- .setup_chr19()

  ann_df <- data.table::data.table(
    chr = "19",
    start = 1L,
    end = 300000000L,
    fragmentID = "fragB",
    annotation = "bait"
  )

  tmp <- tempfile(fileext = ".tsv")
  data.table::fwrite(ann_df, tmp, sep = "\t")

  y <- annotate_interactions(x, tmp)

  expect_s4_class(y, "HiCaptuRe")
  expect_identical(length(y), length(x))
  expect_true(all(c("bait_1", "bait_2") %in% colnames(S4Vectors::mcols(y))))
  expect_false(any(is.na(y$bait_1)))
  expect_false(any(is.na(y$bait_2)))

})

test_that("annotate_interactions: input with wrong number of columns errors clearly", {
  x <- methods::new("HiCaptuRe") # minimal object just to hit the guard; will fail later anyway if it passes
  bad_annot <- data.frame(a = 1, b = 2, c = 3, d = 4)

  expect_error(
    annotate_interactions(x, bad_annot),
    "Annotation file must have exactly 5 columns: chr, start, end, fragmentID, annotation"
  )
})
