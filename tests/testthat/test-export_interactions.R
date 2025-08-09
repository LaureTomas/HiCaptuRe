library(testthat)
library(HiCaptuRe)

test_that("export_interactions: overwrite behavior", {
  x <- .setup_chr19()

  f <- tempfile(fileext = ".tsv")
  writeLines("stub", f)

  expect_error(
    export_interactions(x, f, format = "ibed", over.write = FALSE),
    "File already exists. Use `over.write = TRUE` to overwrite."
  )

  expect_invisible(export_interactions(x, f, format = "ibed", over.write = TRUE, cutoff = 0))
  expect_true(file.exists(f))

  tab <- data.table::fread(f)
  expect_true("bait_chr" %in% colnames(tab))
})

test_that("export_interactions: parameters file is written when parameters=TRUE", {
  x <- .setup_chr19()

  f <- tempfile(fileext = ".tsv")
  expect_invisible(export_interactions(x, f, format = "ibed", over.write = TRUE, parameters = TRUE, cutoff = 0))
  expect_true(file.exists(f))

  param <- paste0(f, ".parameters")
  expect_true(file.exists(param))
  expect_gt(file.info(param)$size, 0)
})

test_that("export_interactions: ibed respects cutoff", {
  x <- .setup_chr19()

  cutoff <- 0
  f <- tempfile(fileext = ".tsv")
  export_interactions(x, f, format = "ibed", over.write = TRUE, cutoff = cutoff)

  tab <- data.table::fread(f)
  expected_cols <- c(
    "bait_chr", "bait_start", "bait_end", "bait_name",
    "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name",
    "N_reads", "score"
  )
  expect_identical(colnames(tab), expected_cols)
  expect_gte(nrow(tab), 1)
})

test_that("export_interactions: bedpe exports empty file when all rows filtered out", {
  x <- setup_load_chr19()

  # set cutoff extremely high to drop everything
  f <- tempfile(fileext = ".bedpe")
  export_interactions(x, f, format = "bedpe", over.write = TRUE, cutoff = 1e9)

  expect_true(file.exists(f))
  # empty data.frame with col.names = FALSE -> empty file
  expect_equal(file.info(f)$size, 0)
})

test_that("export_interactions: seqmonk duplicates anchors and drops ID column in file", {
  x <- setup_load_chr19()

  # choose a cutoff that retains at least one row
  cutoff <- 0
  kept <- x[x$CS >= cutoff]
  expected_n <- 2L * length(kept)  # seqmonk writes both directions and drops ID column

  f <- tempfile(fileext = ".txt")
  export_interactions(x, f, format = "seqmonk", over.write = TRUE, cutoff = cutoff)

  tab <- data.table::fread(f)
  # 6 columns expected after dropping ID
  expect_equal(ncol(tab), 6L)
  expect_equal(nrow(tab), expected_n)
})

test_that("export_interactions: washU format writes two columns with UCSC style", {
  x <
