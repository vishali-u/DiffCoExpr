#' Plot an expression matrix (columns are cells and rows are genes)
#' 
#' @param expression_matrix A matrix containing gene counts per cell. Each
#'    row is a gene and each column is a cell
#'    
#' @return A heatmap where each row is a gene and each column is a cell. If the
#'    gene is highly expressed in the cell, the cell will be close to red. If
#'    the gene is lowly expressed in the cell, the cell will be close to blue.
#' 
#' @export
#' @import ComplexHeatmap
#' 
plotExpressionMatrix <- function(expression_matrix) {
  
  # TODO: use a different package maybe
  # Plot a heatmap of the expression matrix to visualize differences in 
  # expression level across cells
  expression_heatmap <- ComplexHeatmap::Heatmap(expression_matrix, 
                                                show_column_dend = FALSE,
                                                show_row_dend = FALSE,
                                                show_row_names = FALSE,
                                                show_column_names = FALSE, 
                                                row_title = "gene", 
                                                column_title = "cell",
                                                column_title_side = "bottom")
  return(expression_heatmap)
}