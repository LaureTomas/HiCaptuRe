hicapture_attach_message <- function() {
  header <- cli::rule(
    left = cli::style_bold("Attaching HiCaptuRe package"),
    right = paste0("version ", utils::packageVersion("HiCaptuRe"))
  )
  info0 <- paste0(
    cli::col_br_yellow(cli::symbol$pointer), " ",
    cli::format_inline('If you want to see progress bar and messages in some functions, please set progressr::handlers(global = T) & progressr::handlers("progress")')
  )

  info <- paste0(
    cli::col_cyan(cli::symbol$info), " ",
    cli::format_inline("Please report any bug or ask for a feature in our {.href [GitHub](https://github.com/LaureTomas/HiCaptuRe)}")
  )

  paste(header, "", info0, info, sep = "\n")
}
