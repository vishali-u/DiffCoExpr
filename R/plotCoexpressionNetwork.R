#' Create a co-expression network graph. Each gene will be a node, and if the
#' genes are co-expressed, there will be an edge connecting them. If separate
#' communities/clusters of co-expressed genes are detected, each community
#' will be a different color. 
#' 
#' @param edgeList A list of all the genes that should have an edge between
#'    them because they are co-expressed
#'    
#' @return A network graph where co-expressed genes are connected by an edge
#' 
#' @examples
#' # Using an example coexpression network that was generated using the pbmc 
#' # data that is available in this package. This matrix only includes platelet
#' # cells.
#' coexprNetPath <- system.file("extdata", 
#'                              "coexpressionNetworkPlatelet.csv", 
#'                              package = "DiffCoExpr")
#' coexprNet <- read.csv(coexprNetPath)
#' # Delete the first column 
#' coexprNet <- coexprNet[, -1]
#' 
#' coexprNetworkGraph <- plotCoexpressionNetwork(edgeList = coexprNet)
#' coexprNetworkGraph
#' 
#' @export
#' @importFrom igraph graph_from_data_frame layout_with_fr
plotCoexpressionNetwork <- function(edgeList) {
  
  if (!is.data.frame(edgeList)) {
    stop("The input must be a dataframe.")
  }
  
  if (nrow(edgeList) == 0) {
    stop("The coexpression network you provided is empty. Provide a non-empty 
         network.")
  }
  
  networkGraph <- igraph::graph_from_data_frame(edgeList, 
                                                directed = FALSE)
  
  # Identify communities
  communities <- igraph::cluster_louvain(networkGraph)
  message(sprintf("There were %s communities detected.", length(communities)))
  
  # Assign a color to each community
  communityColors <- rainbow(length(communities))
  vertexColor <- communityColors[communities$membership]
  
  # Create the graph and a legend
  # Choose large graph layout
  layout <- igraph::layout_with_lgl(networkGraph) 
  plotGraph <- plot(networkGraph, 
                    vertex.label = NA, 
                    layout = layout,
                    vertex.color = adjustcolor(vertexColor, alpha.f = 0.7),
                    vertex.size = 4,
                    edge.color = adjustcolor("grey", alpha.f = 0.3))
  legend("topright", 
         legend = paste("Community", seq(length(communities))), 
         col = communityColors, 
         pch = 16)
  
  # Print a list of the genes in each community to the screen
  printCommunities(networkGraph = networkGraph,
                   communities = communities)
  
  return(plotGraph)
}

#' Print what genes are part of the same community
#' 
#' @param networkGraph an igraph object that contains some information about
#'    a coexpression network
#'    
#' @param communities a list of communities mapping community names to the
#'    genes that are in that community  
printCommunities <- function(networkGraph, communities) {
  for (i in 1:length(communities)) {
    cat("Community", i, ":\n")
    communityGenes <- communities[[i]]
    concatenatedGenes <- paste(communityGenes, collapse = ", ")
    cat(concatenatedGenes, "\n\n")
  }
}
