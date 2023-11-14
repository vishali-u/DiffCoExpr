#' Create (and visualize) the expression matrix for single-cell RNA-seq data
#' 
#' @param srat A Seurat object containing normalized counts from a scRNA_seq
#'    experiment.
#' @param cell_type A cell type that should be used to filter the expression
#'    matrix so all cells are of that cell type. Can also be a cluster if
#'    cell types were not used. Assumes cell type/cluster are in the Seurat
#'    object. Default is set to NULL, so all cells of all types are included
#'    in the expression matrix
#' 
#' @return A dense matrix where the rows are genes and the columns are cell
#'    barcodes. Each cell contains the log normalized counts of that gene in
#'    that cell
#' 
#' @examples
#' # Using saved Seurat object in the Data folder
#' TODO: removed the srat object from the data folder, must use the raw data and prepareData function
#' srat <- readRds(file = "data/pbmc_srat.rds")
#' expr_matrix <- get_expression_matrix(srat = srat)
#' expr_matrix
#' 
#' @references 
#' Butler, A. (2015). Seurat: Tools for Single Cell Genomics
#' R package version 4.4.0
#' \href{https://cran.r-project.org/web/packages/Seurat/}{Link}.
#' 
#' @export
#' @import Seurat
#' @import ComplexHeatmap
getExpressionMatrix <- function(srat, 
                                cell_type = NULL) {
  
  # Get a matrix of normalized values from the Seurat object
  data_matrix <- as.matrix(Seurat::GetAssayData(srat, 
                                                assay = "RNA", 
                                                slot = "data"))
  
  # Only use the variable features
  variable_genes <- intersect(row.names(data_matrix), 
                              Seurat::VariableFeatures(srat))
  
  # Extract the expression matrix for the selected genes
  expression_matrix <- data_matrix[variable_genes, ]
  
  # Filter by a cell type if applicable
  if (! is.null(cell_type)) {
    cells_to_keep <- Seurat::WhichCells(srat, 
                                        ident = cell_type)
    expression_matrix <- expression_matrix[, colnames(expression_matrix) %in% cells_to_keep]
    
  }
  
  return(expression_matrix)
  
}


