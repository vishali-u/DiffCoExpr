#' Create a co-expression network graph. Each gene will be a node, and if the
#' genes are co-expressed, there will be an edge connecting them.
#' 
#' @param edge_list A list of all the genes that should have an edge between
#'    them because they are co-expressed
#'    
#' @return A network graph where co-expressed genes are connected by an edge
#' 
#' @examples
#' # Using pbmc_srat.rds file in data folder
#' srat <- readRds(file = "data/pbmc_srat.rds")
#' expr_matrix <- getExpressionMatrix(srat = srat)
#' corr_matrix <- getCorrelationMatrix(expression_matrix = expr_matrix)
#' edge_list <- buildCoexpressionNetwork(correlation_matrix = corr_matrix)
#' network_graph <- plotCoexpressionNetwork(edge_list)
#' network_graph
#' 
#' @export
#' @importFrom igraph graph_from_data_frame layout_with_fr
plotCoexpressionNetwork <- function(edge_list) {
  network_graph <- igraph::graph_from_data_frame(edge_list, 
                                                 directed = FALSE)
  plot_graph <- plot(network_graph, 
                     vertex.label.cex = 0.5, 
                     layout = igraph::layout_with_fr(network_graph),
                     vertex.label.dist = 1.5,
                     vertex.size = 5)
  return(plot_graph)
}
