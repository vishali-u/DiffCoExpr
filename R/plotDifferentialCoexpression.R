#' Given two coexpression networks and expression matrices that were made
#' under two different 'conditions' (e.g. two different cell types) and two
#' co-expressed genes, plot the expression level of each gene across cells
#' under each condition as a scatter plot. The x-axis is expression of one gene 
#' and the y-axis is expression of the other gene. Each point represents the 
#' expression level of one gene in one cell under one condition vs. the 
#' expression level of a coexpressed gene in another cell under another 
#' condition. Genes that are highly coexpressed will have points close to the 
#' diagonal. Genes that are differentially coexpressed will have different 
#' patterns of scatter across the two different conditions. 
#' 
#' @param networkA An edge list for a coexpression network. Contains the
#'     correlation coefficient between pairs of genes. 

#' @param networkB An edge list for a different co-expression network. 
#'     Contains the correlation coefficient between pairs of genes. 
#' 
#' @param gene1 A gene name
#' 
#' @param gene2 Another gene name. Gene1 and Gene2 should be coexpressed.
#' 
#' @param expressionMatrixA The expression matrix that was used to create
#'     networkA
#'     
#' @param expressionMatrixB The expression matrix that was used to create
#'     networkB
#'     
#' @examples
#' # Using an example data that was generated using the pbmc data that is
#' # available in this package. 
#' 
#' # Condition A: Platelet cells
#' exprMatrixPlateletPath <- system.file("extdata", 
#'                                       "expressionMatrixPlatelet.csv", 
#'                                       package = "DiffCoExpr")
#' exprMatrixPlatelet <- read.csv(exprMatrixPlateletPath)
#' # Set the row names and remove the column containing gene names
#' rownames(exprMatrixPlatelet) <- exprMatrixPlatelet[, 1]
#' exprMatrixPlatelet <- exprMatrixPlatelet[, -1]
#' 
#' coexprNetPlateletPath <- system.file("extdata", 
#'                                      "coexpressionNetworkPlatelet.csv", 
#'                                      package = "DiffCoExpr")
#' coexprNetPlatelet <- read.csv(coexprNetPlateletPath)
#' # Delete the first column 
#' coexprNetPlatelet <- coexprNetPlatelet[, -1]
#' 
#' # Condition B: DC Cells
#' exprMatrixDCPath <- system.file("extdata", 
#'                                 "expressionMatrixDC.csv", 
#'                                 package = "DiffCoExpr")
#' exprMatrixDC <- read.csv(exprMatrixDCPath)
#' # Set the row names and remove the column containing gene names
#' rownames(exprMatrixDC) <- exprMatrixDC[, 1]
#' exprMatrixDC <- exprMatrixDC[, -1]
#' 
#' coexprNetDCPath <- system.file("extdata", 
#'                                "coexpressionNetworkDC.csv", 
#'                                package = "DiffCoExpr")
#' coexprNetDC <- read.csv(coexprNetDCPath)
#' # Delete the first column 
#' coexprNetDC <- coexprNetDC[, -1]
#' 
#' gene1 <- 'ADH5'
#' gene2 <- 'C4orf3'
#' 
#' plot <- plotDifferentialCoexpression(networkA = coexprNetPlatelet, 
#'                                      networkB = coexprNetDC, 
#'                                      gene1 = gene1, 
#'                                      gene2 = gene2,
#'                                      expressionMatrixA = exprMatrixPlatelet,
#'                                      expressionMatrixB = exprMatrixDC)
#' plot
#' 
#' @return a scatter plot of expression level of gene1 vs. expression level of
#'     gene 2 under two different conditions
#' 
#' @import ggplot2
#' @import dplyr
#' @export
plotDifferentialCoexpression <- function(networkA, 
                                         networkB, 
                                         gene1, 
                                         gene2,
                                         expressionMatrixA,
                                         expressionMatrixB) {
  # Some checks for invalid input
  checkForIncorrectInput(exprMatrix = expressionMatrixA, 
                         coexprNetwork = networkA, 
                         gene1 = gene1, 
                         gene2 = gene2)
  checkForIncorrectInput(exprMatrix = expressionMatrixB, 
                         coexprNetwork = networkB, 
                         gene1 = gene1, 
                         gene2 = gene2)
  
  networkA <- normalizeGenePairs(networkA)
  networkB <- normalizeGenePairs(networkB)
  
  # Merge the edge lists for both networks
  mergedList <- networkA %>%
    dplyr::inner_join(networkB, by = c("Gene1", "Gene2"), 
                      suffix = c("_A", "_B"))
  
  # Calculate the fold change
  mergedList$foldChange <- with(mergedList, 
                                 log2(Correlation_B / Correlation_A))
  
  coexpressedGenes <- checkGeneCoexpression(mergedList = mergedList,
                                            gene1 = gene1, 
                                            gene2 = gene2)
  message(sprintf("%s and %s are coexpressed in both networks. In networkA, the
                  correlation is %s. In networkB, the correlation is %s. 
                  The log2 fold change in correlation is %s.", 
                  gene1, gene2, 
                  coexpressedGenes$Correlation_A, 
                  coexpressedGenes$Correlation_B, 
                  coexpressedGenes$foldChange))
  
  message("Plotting expression levels...")
  
  # Extract expression levels for both genes and store each in a data frame
  exprA <- data.frame(
    t(expressionMatrixA[gene1, ]),
    t(expressionMatrixA[gene2, ])
  )
  colnames(exprA) <- c('ExpressionGene1', 'ExpressionGene2')
  
  exprB <- data.frame(
    t(expressionMatrixB[gene1, ]),
    t(expressionMatrixB[gene2, ])
  )
  colnames(exprB) <- c('ExpressionGene1', 'ExpressionGene2')

  exprA$Condition <- "Condition A"
  exprB$Condition <- "Condition B"
  
  # Create a ggplot
  plot <- 
    ggplot2::ggplot() +
    ggplot2::geom_point(data = exprA, ggplot2::aes(x = ExpressionGene1, 
                                                   y = ExpressionGene2, 
                                                   color = Condition),
                        size = 2) +
    ggplot2::geom_point(data = exprB, ggplot2::aes(x = ExpressionGene1, 
                                                   y = ExpressionGene2, 
                                                   color = Condition), 
                        size = 2) +
    ggplot2::labs(x = paste("Expression of", gene1),
                  y = paste("Expression of", gene2),
                  title = paste("Coexpression of", gene1, "and", gene2)) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c("Condition A" = "blue", 
                                           "Condition B" = "red")) +
    ggplot2::theme_gray()
  

  return(plot)
}

#' 'Normalize' the coexpression network so that the gene that comes first
#' alphabetically is always gene1. This makes sure that rows are removed 
#' incorrectly. For example, if network A has gene1 = A and gene2 = B while
#' network B has gene1 = B and gene2 = A, this pair should not be removed
#' from the differential coexpression analysis
#' 
#' @param network a dataframe containing where rows are gene pairs and the 
#'     corresponding correlation
#'     
#' @return the network but now Gene1 is always the gene that comes first
#'     alphabetically
normalizeGenePairs <- function(network) {
  network$Gene1 <- pmin(network$Gene1, network$Gene2)
  network$Gene2 <- pmax(network$Gene1, network$Gene2)
  return(network)
}

#' Check if the expression matrix and coexpression network are valid input
#' 
#' @param exprMatrix an expression matrix
#' @param coexprNetwork a coexpression network
#' @param gene1 a gene name
#' @param gene2 a different gene name
#' 
checkForIncorrectInput <- function(exprMatrix, coexprNetwork, gene1, gene2) {
  
  if (! (gene1 %in% rownames(exprMatrix))) {
    stop("Gene1 is not present in one of the expression matrices provided.")
  }
  
  if (! (gene2 %in% rownames(exprMatrix))) {
    stop("Gene2 is not present in one of the expression matrices provided.")
  }
  
  if (gene1 == gene2) {
    stop("Please input different genes for gene1 and gene2.")
  }
  
  if (nrow(coexprNetwork) == 0) {
    stop("One of the coexpression networks you provided is empty. Provide a 
          non-empty network.")
  }
  
  if (nrow(exprMatrix) == 0) {
    stop("One of the expression matrices you provided is empty. Provide a 
          non-empty expression matrix")
  }
  
  if (! all(sapply(exprMatrix, is.numeric))) {
    stop("Some of the values in the expression matrix you provided are not 
         numeric. Remove non-numeric values.")
  }
  
  if (! all(sapply(coexprNetwork$Correlation, is.numeric))) {
    stop("Some of the correlation values in the coexpression network you 
          provided are not numeric. Remove non-numeric values.")
  }
}

#' Given two genes, and a a dataframe that contains pairs of genes that are
#' coexpressed in two networks, check if the two genes are in the dataframe.
#' If the genes are not in the dataframe, this pair of genes is not 
#' coexpressed in both networks. 
#' 
#' @param mergedList a dataframe containing pairs of genes that are coexpressed
#'     in two networks, columns for the correlation in each network, and the
#'     log2 fold change of the correlation
#'     
#' @param gene1 a gene name
#' @param gene2 another gene name
#' 
#' @return the row in mergedList containing gene1 and gene2 if they are 
#'     coexpressed in both networks
#' 
#' @import dplyr
checkGeneCoexpression <- function(mergedList, gene1, gene2) {
  
  sortedGenes <- sort(c(gene1, gene2))
  normalizedGene1 <- sortedGenes[1]
  normalizedGene2 <- sortedGenes[2]
  
  # Check if the gene pair is in the merged list
  coexpressedPair <- mergedList %>%
    dplyr::filter(Gene1 == normalizedGene1 & Gene2 == normalizedGene2)
  
  if (nrow(coexpressedPair) == 0) {
    stop(sprintf("The genes %s and %s are not coexpressed in both networks.\n", 
                 gene1, gene2),
         "Pairs of genes that are coexpressed in both networks are:\n",
         paste(mergedList$Gene1, mergedList$Gene2, sep = "-", collapse = ", "))
  } else {
    return(coexpressedPair)
  }
}
