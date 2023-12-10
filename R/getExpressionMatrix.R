#' Create an expression matrix for single-cell RNA-seq data. Each row is a 
#' gene and each column is a cell.
#' 
#' @param srat A Seurat object containing normalized counts from a scRNA-seq
#'    experiment.
#' @param cellType A cell type that should be used to filter the expression
#'    matrix so all cells are of that cell type. Can also be a default cluster 
#'    name (e.g. c0, c1,...) if cell types were not used. Default is set to 
#'    NULL, so all cells of all types are included in the expression matrix.
#' 
#' @return A dense matrix where the rows are genes and the columns are cell
#'    barcodes. Each cell contains the log normalized counts of that gene in
#'    that cell.
#' 
#' @examples
#' # Using pbmc dataset available with package
#' geneMatrixPath <- system.file("extdata",
#'                               "filtered_gene_bc_matrices",
#'                                package =  "DiffCoExpr")
#' cellTypes <- system.file("extdata", 
#'                          "cell_types.csv", 
#'                           package = "DiffCoExpr")
#' 
#' pbmc <- prepareData(geneMatrixPath = geneMatrixPath,
#'                     cellTypesPath = cellTypes)
#'                     
#' exprMatrix <- getExpressionMatrix(srat = pbmc, 
#'                                   cellType = "Platelet")
#' exprMatrix
#' 
#' @references 
#' Butler, A. (2015). Seurat: Tools for Single Cell Genomics
#' R package version 4.4.0
#' \href{https://cran.r-project.org/web/packages/Seurat/}{Link}.
#' 
#' @export
#' @import Seurat
getExpressionMatrix <- function(srat, 
                                cellType = NULL) {
  
  if (! inherits(srat, "Seurat")) {
    stop("The input you provided for srat is not a Seurat object. You must pass 
         a seurat object containing normalized scRNA counts data.")
  }
  
  if (! is.null(cellType) && !(cellType %in% levels(Idents(srat)))) {
    stop("The cell type you provided is not a valid cell type in the Seurat
         object you provided. Re-run with a valid cell type or no cell type.")
  }
  
  # Get a matrix of normalized values from the Seurat object
  dataMatrix <- as.matrix(Seurat::GetAssayData(srat,
                                               assay = "RNA", 
                                               slot = "data"))
  
  # Only use the variable features
  variableGenes <- intersect(row.names(dataMatrix), 
                              Seurat::VariableFeatures(srat))
  
  # Extract the expression matrix for the selected genes
  expressionMatrix <- dataMatrix[variableGenes, ]
  
  # Filter by a cell type if applicable
  if (! is.null(cellType)) {
    cellsToKeep <- Seurat::WhichCells(srat,
                                      ident = cellType)
    expressionMatrix <- expressionMatrix[, colnames(expressionMatrix) %in% 
                                           cellsToKeep]
    
  }
  
  return(expressionMatrix)
  
}


