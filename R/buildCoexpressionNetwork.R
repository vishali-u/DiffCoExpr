#' Identify the genes that have a positive correlation above a threshold 
#' correlation. Here the threshold correlation was chosen as the median 
#' correlation out of all the positive correlations.
#' 
#' @param correlation_matrix A matrix containing the correlation between all 
#'    pairs of genes
#' 
#' @return A data frame storing pairs of genes and the correlation between
#'     the genes
#'     
#' @examples
#' # Using pbmc_srat.rds file in data folder
#' srat <- readRds(file = "data/pbmc_srat.rds")
#' expr_matrix <- getExpressionMatrix(srat = srat)
#' corr_matrix <- getCorrelationMatrix(expression_matrix = expr_matrix)
#' edge_list <- buildCoexpressionNetwork(correlation_matrix = corr_matrix)
#'     
#' @export
#'
buildCoexpressionNetwork <- function(correlation_matrix) {
  # Only keep the gene pairs that have a positive correlation
  
  # Create a new matrix with NA where the correlation is less than or equal to 0
  positive_correlation_matrix <- correlation_matrix
  positive_correlation_matrix[positive_correlation_matrix <= 0] <- NA
  
  # Get the row and column indices of the positive correlations
  positive_pairs_indices <- which(correlation_matrix > 0, arr.ind = TRUE)
  
  # Create a data frame of gene pairs with positive correlations
  positive_correlation_pairs <- data.frame(
    Gene1 = rownames(correlation_matrix)[positive_pairs_indices[, 1]],
    Gene2 = colnames(correlation_matrix)[positive_pairs_indices[, 2]],
    Correlation = correlation_matrix[positive_pairs_indices]
  )
  
  # Convert the data frame back to a matrix
  # Get all unique gene names
  all_genes <- unique(c(positive_correlation_pairs$Gene1, 
                        positive_correlation_pairs$Gene2))
  
  # Create an empty matrix
  coexpression_matrix <- matrix(NA, 
                                nrow = length(all_genes), 
                                ncol = length(all_genes),
                                dimnames = list(all_genes, all_genes))
  
  # Add the correlations to the square matrix
  for (i in 1:nrow(positive_correlation_pairs)) {
    gene1 <- positive_correlation_pairs$Gene1[i]
    gene2 <- positive_correlation_pairs$Gene2[i]
    corr_value <- positive_correlation_pairs$Correlation[i]
    
    coexpression_matrix[gene1, gene2] <- corr_value
    coexpression_matrix[gene2, gene1] <- corr_value
  }
  
  # Get the median correlation from the coexpression matrix
  # Remember the diagnoal is all 1's, so extract the lower (or upper) 
  # triangular part of the matrix, excluding the diagonal
  lower_tri_values <- coexpression_matrix[lower.tri(coexpression_matrix)]
  
  # Calculate the median of these values
  median_correlation <- median(coexpression_matrix)
  
  coexpression_matrix[coexpression_matrix < median_correlation] <- NA
  
  # Construct a network graph from the co-expression matrix. Each node will be
  # a gene, and if the genes are co-expressed, there will be an edge 
  # connecting them
  
  # Get the gene names
  genes <- rownames(coexpression_matrix)
  
  # Create a dataframe for the edge list
  edge_list <- which(!is.na(coexpression_matrix), 
                     arr.ind = TRUE)
  edge_list <- data.frame(gene1 = genes[edge_list[, "row"]], 
                          gene2 = genes[edge_list[, "col"]], 
                          weight = coexpression_matrix[edge_list])
  # Remove self loops
  edge_list <- edge_list[edge_list$gene1 != edge_list$gene2, ]
  
  return(edge_list)
}




