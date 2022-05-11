#' Generate the plot based on distance output
#'
#' This function plots the distances by different approaches for the output of split_by_distance
#'
#' @param distances list of dataframes output from \code{\link{split_by_distance}}
#' @param type_of_value an element from: absolute, by_int_type, by_total
#'
#' @details The type_of_value argument could be: Absolute: plots the Absolute number of interactions. By_int_type: plots the number of interactions divided by the total number of interactions of each type. By_total: plots the number of interactions divided by the total number of interactions
#'
#' @return list with the plots
#'
#' @importFrom ggplot2 ggplot
#'
#' @export
plot_distance <- function(distances, type_of_value="absolute")
{
  if (all(names(distances) %in% c("short_int_dist_table", "long_int_dist_table")))
  {

    if (type_of_value == "absolute")
    {
      plots <- list()

      short <- distances$short_int_dist_table
      short1 <- short[!(short$int == "Total"),]
      short1 <- short1[!(short1$variable == "total_per_int"),]

      g <- ggplot2::ggplot(short1, aes(x=int, y=value, fill=variable)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="Number of interactions", fill="Distances",
             title=paste0("Number of Short interactions of ",unique(short1$cell_type)," by type of interaction"),
             subtitle = "Absolute Values") + theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
        ggplot2::theme_minimal()

      plots[["short_int_dist"]] <- g

      shortT <- short[short$int == "Total",]

      g <- ggplot2::ggplot(shortT, aes(x=variable, y=value, fill=variable)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Distances",y="Number of interactions", fill="Distances",
             title = paste0("Total Number of Short interactions of ",unique(short1$cell_type)),
             subtitle = "Absolute Values") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["short_total_dist"]] <- g

      long <- distances$long_int_dist_table
      long1 <- long[!(long$int == "Total"),]
      long1 <- long1[!(long1$variable == "total_per_int"),]

      g <- ggplot2::ggplot(long1, aes(x=int, y=value, fill=variable)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="Number of interactions", fill="Distances",
             title = paste0("Number of Long interactions of ",unique(long1$cell_type)," by type of interaction"),
             subtitle = "Absolute Values") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["long_int_dist"]] <- g

      longT <- long[long$int == "Total",]

      g <- ggplot2::ggplot(longT, aes(x=variable, y=value, fill=variable)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Distances",y="Number of interactions", fill="Distances",
             title = paste0("Total Number of Long interactions of ",unique(long1$cell_type)),
             subtitle = "Absolute Values") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["long_total_dist"]] <- g

    }
    if (type_of_value == "by_int_type")
    {
      plots <- list()

      short <- distances$short_int_dist_table
      short1 <- short[!(short$int == "Total"),]
      short1 <- short1[!(short1$variable == "total_per_int"),]

      g <- ggplot2::ggplot(short1, aes(x=int, y=(value/total_per_int)*100, fill=variable)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="% of interactions", fill="Distances",
             title=paste0("% of Short interactions of ",unique(short1$cell_type)," by type of interaction"),
             subtitle = "Normalized by the total number of interaction of each type") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["short_int_dist_norm_int"]] <- g

      long <- distances$long_int_dist_table
      long1 <- long[!(long$int == "Total"),]
      long1 <- long1[!(long1$variable == "total_per_int"),]

      g <- ggplot2::ggplot(long1, aes(x=int, y=(value/total_per_int)*100, fill=variable)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="% of interactions", fill="Distances",
             title = paste0("Number of Long interactions of ",unique(long1$cell_type)," by type of interaction"),
             subtitle = "Normalized by the total number of interaction of each type") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["long_int_dist_norm_int"]] <- g

    }

    if (type_of_value == "by_total")
    {
      plots <- list()

      short <- distances$short_int_dist_table
      short1 <- short[!(short$int == "Total"),]
      short1 <- short1[!(short1$variable == "total_per_int"),]

      g <- ggplot2::ggplot(short1, aes(x=int, y=(value/ibed)*100, fill=variable)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="% of interactions", fill="Distances",
             title=paste0("% of Short interactions of ",unique(short1$cell_type)," by type of interaction"),
             subtitle = "Normalized by the Total number of interaction") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["short_int_dist_norm_int"]] <- g

      long <- distances$long_int_dist_table
      long1 <- long[!(long$int == "Total"),]
      long1 <- long1[!(long1$variable == "total_per_int"),]

      g <- ggplot2::ggplot(long1, aes(x=int, y=(value/ibed)*100, fill=variable)) +
        ggplot2::geom_bar(width = 1, size = 1, stat = "identity",position="dodge") +
        ggplot2::labs(x="Type of interactions",y="% of interactions", fill="Distances",
             title = paste0("Number of Long interactions of ",unique(long1$cell_type)," by type of interaction"),
             subtitle = "Normalized by the Total number of interaction") +
        ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = element_text(angle = 70, hjust = 1))

      plots[["long_int_dist_norm_int"]] <- g

    }
    return(plots)
  }
  else
  {
    stop("Input element is not ok")
  }
}
