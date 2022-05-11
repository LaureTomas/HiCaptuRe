#' Generate Cytoscape Network file
#'
#' This function creates the file needed to upload to Cytoscape
#'
#' @param interactions GenomicInteractions object from \code{\link{load_interactions}}
#' @param file full path to desired output file
#' @param type Type of network you want to obtained (Complete, P_P or P_OE)
#'
#' @return File needed to upload to Cytoscape
#'
#' @importFrom GenomicInteractions export.igraph
#' @importFrom data.table fwrite
#' @importFrom igraph simplify as_edgelist
#'
#' @export
interactions2cytoscape <- function(interactions, file, type="complete")
{
  if(type == "complete")
  {
    net <- igraph::simplify(export.igraph(interactions))
    nodes_edges <- igraph::as_edgelist(net)

    message(paste("Complete Network saved in",file))
  }
  if(type %in% c("P_P","P_OE"))
  {
    net <- igraph::simplify(GenomicInteractions::export.igraph(interactions[grepl(type,interactions$int)]))
    nodes_edges <- igraph::as_edgelist(net)

    message(paste(type, "Subnetwork saved in",file))
  }

  suppressMessages(data.table::fwrite(nodes_edges, file = file, row.names = F, col.names = F, quote = F, sep = "\t"))

  message("In Cytoscape do:")
  message("File > Import > Network from File ... (Ctrl + L)")
  message("Import Network From Table:")
  message("Advanced Options... > Delimiter: TAB > Uncheck 'Use first line as column names' > OK")
  message("Click first column > Source node (Green circle)")
  message("Click second column > Target node (Orange target)")
  message("Click OK and here you have your network")
}
