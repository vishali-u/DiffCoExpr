#' Given a table of differentially coexpressed genes, the expression matrices
#' that were used to generate the table of differentially coexpressed genes, and
#' two coexpressed genes, plot the expression level of each gene across cells
#' under each condition as a scatter plot. The x-axis is expression of one gene 
#' and the y-axis is expression of the other gene. Each point represents the 
#' expression level of one gene in one cell under one condition vs. the 
#' expression level of a coexpressed gene in another cell under another 
#' condition. Genes that are highly coexpressed will have points close to the 
#' diagonal. Genes that are differentially coexpressed will have different 
#' patterns of scatter across the two different conditions. 
#' 
#' @param diffCoexpTable A dataframe containing pairs of genes and their
#'     correlations in two networks: networkA and networkB. Also contains
#'     the log2FC
#' 
#' @param expressionMatrixA The expression matrix that was used to create
#'     networkA
#'     
#' @param expressionMatrixB The expression matrix that was used to create
#'     networkB
#'     
#' @param gene1 A gene name
#' 
#' @param gene2 Another gene name. Gene1 and Gene2 should be coexpressed. 
#'
#' @param conditionA an optional parameter to specify the condition/cell type
#'     that corresponds to expressionMatrixA
#'
#' @param conditionB an optional parameter to specify the condition/cell type
#'     that corresponds to expressionMatrixB   
#'
#' Note: this function assumes that the parameters are passed in the correct
#' order. That is, expressionMatrixA was used to generate the correlation
#' values stored in diffCoexpTable$CorrelationA, and expressionMatrixB was
#' used to generate the correlation values stored in 
#' diffCoexpTable$CorrelationB. Running the function in the incorrect order
#' most likely will not throw errors, but the results will not be meaningful. 
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
#' diffCoexpTable <- getDifferentialCoexpression(networkA = coexprNetPlatelet, 
#'                                               networkB = coexprNetDC)
#' 
#' gene1 <- 'ADH5'
#' gene2 <- 'C4orf3'
#'
#' plot <- plotDifferentialCoexpression(diffCoexpTable = diffCoexpTable,
#'                                      expressionMatrixA = exprMatrixPlatelet,
#'                                      expressionMatrixB = exprMatrixDC, 
#'                                      gene1 = gene1, 
#'                                      gene2 = gene2,
#'                                      conditionA = "platelets", 
#'                                      conditionB = "DC")
#' plot
#' 
#' @return a scatter plot of expression level of gene1 vs. expression level of
#'     gene 2 under two different conditions
#' 
#' @import ggplot2
#' @export
plotDifferentialCoexpression <- function(diffCoexpTable,
                                         expressionMatrixA,
                                         expressionMatrixB,
                                         gene1, 
                                         gene2,
                                         conditionA = NULL, 
                                         conditionB = NULL) {
  
  # --- Check for invalid input ----------------------------------------------

  checkForIncorrectInput(exprMatrix = expressionMatrixA, 
                         gene1 = gene1, 
                         gene2 = gene2)
  checkForIncorrectInput(exprMatrix = expressionMatrixB, 
                         gene1 = gene1, 
                         gene2 = gene2)
  
  if(nrow(diffCoexpTable) == 0) {
    stop("The differential coexpression data frame you provided is empty.
         Provide a non-empty data frame.")
  }
  
  requiredColumns <- c("Gene1", 
                       "Gene2",
                       "CorrelationA", 
                       "CorrelationB", 
                       "foldChange")
  
  if(! (all(requiredColumns %in% names(diffCoexpTable)))) {
    stop("The dataframe you provided for diffCoexpTable does not have all the
         needed columns. Please provide a dataframe with the columns Gene1, 
         Gene2, CorrelationA, CorrelationB, and foldChange")
  }
  
  # --- Find gene pair and expression levels ---------------------------------
  
  coexpressedGenes <- checkGeneCoexpression(mergedList = diffCoexpTable,
                                            gene1 = gene1, 
                                            gene2 = gene2)
  
  message(sprintf("%s and %s are coexpressed in both networks. In networkA, the
                  correlation is %s. In networkB, the correlation is %s. 
                  The log2 fold change in correlation is %s.", 
                  gene1, gene2, 
                  coexpressedGenes$CorrelationA, 
                  coexpressedGenes$CorrelationB, 
                  coexpressedGenes$foldChange))
  
  message("Plotting expression levels...")
  
  # --- Structure data for plotting -----------------------------------------
  
  # Convert expression matrices to dataframes
  expressionMatrixA <- as.data.frame(expressionMatrixA)
  expressionMatrixB <- as.data.frame(expressionMatrixB)
  
  # Extract expression levels for both genes and store each in a data frame
  exprAGene1 <- t(as.data.frame(expressionMatrixA[gene1, ]))
  exprAGene2 <- t(as.data.frame(expressionMatrixA[gene2, ]))
  exprBGene1 <- t(as.data.frame(expressionMatrixB[gene1, ]))
  exprBGene2 <- t(as.data.frame(expressionMatrixB[gene2, ]))
  
  exprA <- data.frame(
    ExpressionGene1 = exprAGene1[,1],
    ExpressionGene2 = exprAGene2[,1]
  )
  
  exprB <- data.frame(
    ExpressionGene1 = exprBGene1[,1],
    ExpressionGene2 = exprBGene2[,1]
  )
  
  # Define color mappings based on the conditions
  colorMappings <- c(
    "Network A" = "blue", 
    "Network B" = "red"
  )
  
  if (! is.null(conditionA)) {
    exprA$Network <- conditionA
    names(colorMappings)[1] <- conditionA
  } else {
    exprA$Network <- "Network A"
  }
  
  if (! is.null(conditionB)) {
    exprB$Network <- conditionB
    names(colorMappings)[2] <- conditionB
  } else {
    exprB$Network <- "Network B"
  } 
  
  # Create a ggplot
  plot <- 
    ggplot2::ggplot() +
    ggplot2::geom_point(data = exprA, ggplot2::aes(x = ExpressionGene1, 
                                                   y = ExpressionGene2, 
                                                   color = Network),
                        size = 2) +
    ggplot2::geom_point(data = exprB, ggplot2::aes(x = ExpressionGene1, 
                                                   y = ExpressionGene2, 
                                                   color = Network), 
                        size = 2) +
    ggplot2::labs(x = paste("Expression of", gene1),
                  y = paste("Expression of", gene2),
                  title = paste("Expression of", gene1, "vs.", gene2)) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = colorMappings) +
    ggplot2::theme_gray()
  

  return(plot)
}

#' Check if the expression matrix and coexpression network are valid input
#' 
#' @param exprMatrix an expression matrix
#' @param gene1 a gene name
#' @param gene2 a different gene name
#' 
#' @return no return value but throw an error if applicable
#' 
checkForIncorrectInput <- function(exprMatrix, gene1, gene2) {
  
  if (! (gene1 %in% rownames(exprMatrix))) {
    stop("Gene1 is not present in one of the expression matrices provided.")
  }
  
  if (! (gene2 %in% rownames(exprMatrix))) {
    stop("Gene2 is not present in one of the expression matrices provided.")
  }
  
  if (gene1 == gene2) {
    stop("Please input different genes for gene1 and gene2.")
  }
  
  if (nrow(exprMatrix) == 0) {
    stop("One of the expression matrices you provided is empty. Provide a 
          non-empty expression matrix")
  }
  
  if (! all(sapply(exprMatrix, is.numeric))) {
    stop("Some of the values in the expression matrix you provided are not 
         numeric. Remove non-numeric values.")
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
#' @importFrom dplyr filter
checkGeneCoexpression <- function(mergedList, gene1, gene2) {
  
  # Check if the gene pair is in the merged list
  coexpressedPair <- mergedList %>%
    dplyr::filter((Gene1 == gene1 & Gene2 == gene2) | 
                    (Gene1 == gene2 & Gene2 == gene1))
  
  if (nrow(coexpressedPair) == 0) {
    stop(sprintf("The genes %s and %s are not coexpressed in both networks.\n", 
                 gene1, gene2),
         "Pairs of genes that are coexpressed in both networks are:\n",
         paste(mergedList$Gene1, 
               mergedList$Gene2, 
               sep = " and ", 
               collapse = ", "))
  } else {
    # In case there was a duplicate, return just the first row
    return(coexpressedPair[1, ])
  }
}

# [END]
