#' Generate a PDF report for the interaction object
#'
#' This function generates a PDF summary report of the GInteractions object
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param pdf_file file path to save the report
#' @param ... additional parameters for split_by_distance
#'
#' @return PDF report
#'
#' @import gridExtra
#'
#' @export
interactions_summary <- function(interactions, pdf_file, ...)
{

  dat <- as.data.frame(table(interactions$int))
  colnames(dat) <- c("category",  "count")

  dat$fraction <- dat$count/sum(dat$count)

  dat <- dat[order(dat$fraction), ]

  dat$ymax <- cumsum(dat$fraction)
  dat$ymin <- c(0, head(dat$ymax, n = -1))
  dat$label <- paste(dat$category, "\n", signif(100 * dat$fraction,
                                                3), "%")
  dat[dat$fraction < 0.03, "label"] <- NA
  dat$labelpos <- (dat$ymax + dat$ymin)/2
  dat$fraclabel <- paste(signif(100 * dat$fraction, 3), "%")
  plot1 <- ggplot(dat, aes_string(fill = "category", ymax = "ymax",
                              ymin = "ymin", xmax = 4, xmin = 2.3)) + geom_rect() +
    coord_polar(theta = "y") + xlim(c(0, 4)) + ylim(c(0,
                                                      1)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) +
    theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) +
    theme(axis.title = element_blank()) + theme(panel.border = element_blank()) +
    ggtitle("Interaction Classes") +
    geom_text(aes_string(x = (4 + 2.3)/2, y = "labelpos", label = "label")) + theme(legend.position = "none")


  df1 <- dat[order(dat$category),1:2]
  int_dist <- split_by_distance(interactions, ...)
  plots_dist <- plot_distance(int_dist)

  pdf(pdf_file,width = 13, height = 10)
  grid.arrange(grobs=list(plot1,tableGrob(df1,rows = NULL, theme = ttheme_default(base_size = 10)),
                          plots_dist[[1]],plots_dist[[2]],
                          plots_dist[[3]],plots_dist[[4]]), ncol=2)
  dev.off()
}
