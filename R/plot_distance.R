#' Generate the plot based on distance output
#'
#' This function plots the distances by different approaches for the output of distance_summary
#'
#' @param distances dataframes output from \code{\link{distance_summary}}
#' @param type_of_value an element from: absolute, by_int_type, by_total
#'
#' @details The type_of_value argument could be: Absolute: plots the Absolute number of interactions. By_int_type: plots the number of interactions divided by the total number of interactions of each type. By_total: plots the number of interactions divided by the total number of interactions
#'
#' @return list with the plots
#'
#' @importFrom ggplot2 ggplot
#'
#' @export
plot_distance_summary <- function(distances, type_of_value="absolute")
{

    if (type_of_value == "absolute")
    {
      plots <- list()

      distances1 <- distances[!(distances$int == "Total"),]
      distances1 <- distances1[!(distances1$breaks == "total_per_int"),]

      g <- ggplot2::ggplot(distances1, aes(x=int, y=value, fill=breaks)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="Number of interactions", fill="Distances",
             title=paste0("Number of interactions of ",unique(distances1$sample)," by type of interaction"),
             subtitle = "Absolute Values") + theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
        ggplot2::theme_minimal()

      plots[["int_dist"]] <- g

      distancesT <- distances[distances$int == "Total",]

      g <- ggplot2::ggplot(distancesT, aes(x=breaks, y=value, fill=breaks)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Distances",y="Number of interactions", fill="Distances",
             title = paste0("Total Number of interactions of ",unique(distances1$sample)),
             subtitle = "Absolute Values") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["total_dist"]] <- g

    }
    if (type_of_value == "by_int_type")
    {
      plots <- list()

      distances1 <- distances[!(distances$int == "Total"),]
      distances1 <- distances1[!(distances1$breaks == "total_per_int"),]

      g <- ggplot2::ggplot(distances1, aes(x=int, y=(value/total_per_int)*100, fill=breaks)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="% of interactions", fill="Distances",
             title=paste0("% of interactions of ",unique(distances1$sample)," by type of interaction"),
             subtitle = "Normalized by the total number of interaction of each type") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["int_dist_norm_int"]] <- g

    }

    if (type_of_value == "by_total")
    {
      plots <- list()

      distances1 <- distances[!(distances$int == "Total"),]
      distances1 <- distances1[!(distances1$breaks == "total_per_int"),]

      g <- ggplot2::ggplot(distances1, aes(x=int, y=(value/ibed)*100, fill=breaks)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="% of interactions", fill="Distances",
             title=paste0("% of interactions of ",unique(distances1$sample)," by type of interaction"),
             subtitle = "Normalized by the Total number of interactions") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["int_dist_norm_total"]] <- g

      distancesT <- distances[distances$int == "Total",]

      g <- ggplot2::ggplot(distancesT, aes(x=breaks, y=(value/ibed)*100, fill=breaks)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Distances",y="Number of interactions", fill="Distances",
                      title = paste0("Total Number of interactions of ",unique(distances1$sample)),
                      subtitle = "Absolute Values") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["total_dist_norm_total"]] <- g

    }
    return(plots)

}
