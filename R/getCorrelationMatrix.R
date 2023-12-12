#' Create a correlation matrix. Lowly expressed genes, and genes will too 
#' little variation across cells are filtered out. In the matrix, each
#' row and column are genes, and each cell contains the correlation between the
#' pair of genes. 
#' 
#' @param expressionMatrix A matrix containing the expression levels of genes
#'    per cell. The rows are genes and the columns are cells
#'    
#' @param minCellCount A minimum number of cells that should still be present
#'    in the data after filtering. If the number of cells after filtering is
#'    less than this value, stop running. The default is set to 5. This value
#'    must be greater than 2
#'
#' @param minGeneCount A minimum number of genes that should still be present
#'    in the data after filtering. If the number of genes after filtering is
#'    less than this value, stop running. The default is set to 5. This value
#'    must be greater than 2.
#'    
#' @param minPt A miminum percent to filter out genes. Any genes that have a
#'    variation that is less than the variation at minPt percentile will be
#'    removed. Value must be be greater than 0 and less than or equal to 1.
#'    Default is 0.20
#'    
#' @return A matrix where the rows and columns are genes, and the cell contains
#'    the correlation between the two genes
#'    
#' @examples
#' # Using an example expression matrix that was generated using the pbmc data
#' # that is available in this package. This matrix only includes platelet
#' # cells.
#' exprMatrixPath <- system.file("extdata", 
#'                               "expressionMatrixPlatelet.csv", 
#'                                package = "DiffCoExpr")
#' exprMatrix <- read.csv(exprMatrixPath)
#' # Set the row names and remove the column containing gene names
#' rownames(exprMatrix) <- exprMatrix[, 1]
#' exprMatrix <- exprMatrix[, -1]
#' 
#' corrMatrix <- getCorrelationMatrix(expressionMatrix = exprMatrix)
#' corrMatrix
#' 
#' @references
#' Lemoine G., Scott-Boyer M., Ambroise B., Perin O., Droit A. (2021) GWENA: 
#' gene co-expression networks analysis and extended modules characterization 
#' in a single Bioconductor package. BMC Bioinformatics 22, 267.
#' \href{https://doi.org/10.1186/s12859-021-04179-4}{Link}.
#' 
#' @import dplyr
#' @importFrom stats cor sd
#' @export
getCorrelationMatrix <- function(expressionMatrix, 
                                 minCellCount = 5,
                                 minGeneCount = 5,
                                 minPt = 0.20) {
  
  # --- Check some conditions to stop early if necessary --------------------
  
  if (! is.matrix(expressionMatrix) && ! is.data.frame(expressionMatrix)) {
    stop("Please provide a matrix or data.frame object for the 
         expressionMatrix.")
  }
  
  # Need at least 2 cells and 2 genes to get the correlation matrix
  if ((nrow(expressionMatrix) < 2) || (ncol(expressionMatrix) < 2)) {
    stop("There are fewer than 2 cells or fewer than 2 genes in this dataset.
         There is not enough data to get a correlation matrix")
  }
  
  # Check that minPt is in the correct range
  if (minPt < 0 || minPt > 1) {
    warning("The minPt value you provided is outside (0,1]. Using the default 
            minPt of 0.20 instead.")
    minPt <- 0.20
  }
  
  # Convert the expression matrix into a data frame if it is not a data frame
  # to do some calculations on the data
  if (is.matrix(expressionMatrix)) {
    expressionMatrix <- as.data.frame(expressionMatrix)
  }
  
  if (! all(sapply(expressionMatrix, is.numeric))) {
    stop("Some of the values in the matrix you provided are not numeric. Remove
         non-numeric values.")
  }
  
  # --- Filter out genes that show little to no variation --------------------
  
  # Filter out genes with zero variation
  # Calculate standard deviation for each gene and filter out genes with 0
  # standard variation
  geneSD <- apply(expressionMatrix, 1, stats::sd)
  geneSD <- geneSD[geneSD > 0]
  expressionMatrix <- expressionMatrix[names(geneSD), ]
  
  # Create a table mapping genes to their standard deviation (sorted by sd)
  geneVariationTable <- dplyr::tibble(gene = names(geneSD), 
                                      geneSD = geneSD) %>%
    dplyr::arrange(desc(geneSD))
  
  # Filter out genes with low variation (out of the genes that have some
  # variation) by filtering out the lower minPt % of genes based on variation
  cutoffValue <- quantile(geneVariationTable$geneSD, minPt)
  
  # Filter out genes below the cutoff
  topGenes <- geneVariationTable %>% 
    dplyr::filter(geneSD > cutoffValue) %>% 
    dplyr::pull(gene)
  expressionMatrix <- expressionMatrix[topGenes,]
  
  # Convert data frame back to matrix
  expressionMatrix <- as.matrix(expressionMatrix)
  
  # --- Make sure output is still large enough to continue analysis ----------
  
  # if too many rows got filtered out and there are less than 5 cells and 5
  # genes, stop running
  if (nrow(expressionMatrix) < minCellCount || 
      ncol(expressionMatrix) < minGeneCount) {
    stop("After filtering, there are fewer than 5 cells and/or 5 genes. This
         is not enough data for further analysis.")
  }
  
  # Compute the correlation matrix (for Pearson correlation)
  correlationMatrix <- stats::cor(t(expressionMatrix), method = "pearson")
  
  return(correlationMatrix)
}
