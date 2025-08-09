library(testthat)
library(HiCaptuRe)

test_that("peakmatrix2list: splits real peakmatrix example and preserves CS values", {
  x <- suppressWarnings(.setup_peakmatrix_chr19())

  lst <- peakmatrix2list(x)

  expect_true(length(lst) >= 1L)

  expect_true(all(vapply(lst, methods::is, logical(1), "HiCaptuRe")))

  m_all <- S4Vectors::mcols(x)
  cs_cols <- grep("^CS_", colnames(m_all), value = TRUE)
  expect_true(length(cs_cols) >= 1L)

  expect_identical(sort(names(lst)), sort(sub("^CS_", "", cs_cols)))

  # per-sample checks
  for (nm in names(lst)) {
    xi <- lst[[nm]]
    mi <- S4Vectors::mcols(xi)

    expect_true("CS" %in% colnames(mi))
    expect_false(any(grepl("^CS_", colnames(mi))))

    orig_col <- paste0("CS_", nm)
    keep_mask <- m_all[[orig_col]] >= 5

    p <- HiCaptuRe::getParameters(xi)
    expect_true("peakmatrix2list" %in% names(p))
    expect_match(p$peakmatrix2list[["split_by"]], paste0("^", orig_col, "$"))
  }
})

test_that("peakmatrix2list: error when already processed", {
  x <- suppressWarnings(.setup_peakmatrix_chr19())
  lst <- peakmatrix2list(x)
  xi  <- lst[[1]]

  expect_error(
    peakmatrix2list(xi),
    "already appears to have been processed by 'peakmatrix2list\\(\\)'"
  )
})

test_that("peakmatrix2list: errors when no CS_* columns present", {
  x <- suppressWarnings(.setup_peakmatrix_chr19())

  m <- S4Vectors::mcols(x)
  keep <- !grepl("^CS_", colnames(m))
  S4Vectors::mcols(x) <- m[, keep, drop = FALSE]

  expect_error(
    peakmatrix2list(x),
    "No CS_\\* columns found"
  )
})

test_that("peakmatrix2list: errors when format is not peakmatrix", {
  y <- .setup_chr19()

  params <- HiCaptuRe::getParameters(y)
  expect_false(identical(tolower(params$load[["format"]]), "peakmatrix"))

  expect_error(
    peakmatrix2list(y),
    "Input is not a peakmatrix, it is a "
  )
})
