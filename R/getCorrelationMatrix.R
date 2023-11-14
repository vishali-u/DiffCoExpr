#' Create (and visualize) the correlation between every pair of genes in the
#' expression matrix
#' 
#' @param expression_matrix A matrix containing the expression levels of genes
#'    per cell. The rows are genes and the columns are cells
#' @param visualize A boolean that will be used to determine if the function
#'    should also create a heatmap of the correlation matrix. Default is set
#'    to FALSE. 
#'    
#' @return A matrix where the rows and columns are genes, and the cell contains
#'    the correlation between the two genes
#' @examples
#' # Using saved Seurat object in the Data folder
#' srat <- readRds(file = "data/pbmc_srat.rds")
#' expr_matrix <- getExpressionMatrix(srat = srat)
#' corr_matrix <- getCorrelationMatrix(expression_matrix = expr_matrix)
#' corr_matrix
#' 
#' @references
#' Lemoine G., Scott-Boyer M., Ambroise B., Perin O., Droit A. (2021) GWENA: 
#' gene co-expression networks analysis and extended modules characterization 
#' in a single Bioconductor package. BMC Bioinformatics 22, 267.
#' \href{https://doi.org/10.1186/s12859-021-04179-4}{Link}.
#' 
#' @export
#' @import dplyr
#' @import ComplexHeatmap
getCorrelationMatrix <- function(expression_matrix, 
                                 visualize = FALSE) {
  
  # Filter out genes with zero variation
  # Calculate standard deviation for each gene
  gene_sd <- apply(expression_matrix, 1, sd)
  
  # Filter out genes with zero standard deviation
  expression_matrix <- expression_matrix[gene_sd > 0, ]
  
  # Filter out genes with low variation (out of the genes that have some
  # variation) by keeping the top 80% of genes based on variation
  # Convert matrix to data frame
  expression_matrix <- as.data.frame(expression_matrix)
  
  # Calculate the median expression for each gene
  variation <- lapply(expression_matrix, 
                      function(row) do.call(median, list(row)))
  
  top_pct_df <- data.frame(gene = names(variation), 
                           variation = unlist(variation), 
                           stringsAsFactors = FALSE)
  top_80 <- dplyr::top_frac(top_pct_df, 0.8, variation)
  expression_matrix <- expression_matrix[, top_80$gene]
  
  # Convert data frame back to matrix
  expression_matrix <- as.matrix(expression_matrix)
  
  # TODO: figure out how to build matrix without filtering outliers since
  # many rows could get filtered if data is noisy but if it does not work
  # maybe try filtering out genes that have a large proportion (maybe > 50 %)
  # of extreme expression levels
  # First filter out any outliers using z-scores
  z_scores <- scale(expression_matrix)
  # Check each row for extreme z scores and remove them
  extreme_rows <- apply(z_scores, 
                        1, 
                        function(row) any(row > 3 | row < -3))
  expression_matrix <- expression_matrix[!extreme_rows, ]

  # if too many rows got filtered out and there are less than 10 cells and 10
  # genes, stop running
  if (nrow(expression_matrix) < 10 || ncol(expression_matrix) < 10) {
    stop("After filtering, there is not enough data.")
  }
  
  # Compute the correlation matrix (for Pearson correlation)
  correlation_matrix <- cor(t(expression_matrix), method = "pearson")

  return(correlation_matrix)
}


