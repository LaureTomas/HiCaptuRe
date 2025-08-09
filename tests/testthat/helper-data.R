.setup_chr19 <- function() {
  ibed1_file <- system.file("extdata", "ibed1_example.zip", package = "HiCaptuRe")
  x <- load_interactions(
    file = ibed1_file,
    genome = "BSgenome.Hsapiens.NCBI.GRCh38",
    select_chr = "19"
  )
  expect_s4_class(x, "HiCaptuRe")
  x
}


.setup_peakmatrix_chr19 <- function() {
  peakmatrix_file <- system.file("extdata", "peakmatrix_example.zip", package = "HiCaptuRe")
  x <- load_interactions(peakmatrix_file, select_chr = "19")
  expect_s4_class(x, "HiCaptuRe")
  params <- HiCaptuRe::getParameters(x)
  expect_identical(tolower(params$load[["format"]]), "peakmatrix")
  x
}
