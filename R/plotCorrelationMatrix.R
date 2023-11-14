#' Plot a correlation matrix (columns and rows are genes and each cell contains
#' the correlation between the two genes)
#' 
#' @param correlation_matrix A matrix containing the correlation between all 
#'    pairs of genes
#'    
#' @return A heatmap where each row and column is a gene. If the genes are 
#'    highly correlated (i.e. co-expressed), the cell will be close to red. If
#'    the genes are not correlated (i.e. not co-expresseed), the cell will be
#'    close to blue.
#' 
#' @export
#' @import ComplexHeatmap
#' 
plotCorrelationMatrix <- function(correlation_matrix) {
  # TODO: use a different package or something because the gene names are way 
  # too small (maybe remove some outliers)
  # Plot a heatmap of the expression matrix to visualize differences in 
  # expression level across cells
  correlation_heatmap <- ComplexHeatmap::Heatmap(correlation_matrix, 
                                                 show_column_dend = FALSE,
                                                 show_row_dend = FALSE,
                                                 column_names_gp = grid::gpar(fontsize = 2),
                                                 row_names_gp = grid::gpar(fontsize = 2))
  return(correlation_heatmap)
} 