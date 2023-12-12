#' Identify the genes that have a positive correlation strictly higher than a  
#' threshold correlation.
#' 
#' @param correlationMatrix A matrix or data framer where each row and column 
#'    are genes, and each cell contains the correlation between the pairs of 
#'    genes.
#'    
#' @param thresholdCorrelation A percentile used to select only the genes with
#'    a high correlation. This value must be in the range [0, 1]. The default
#'    is 0.80, so only gene pairs that have a correlation higher than the 
#'    correlation at the 80th percentile will be used.
#' 
#' @return A data frame storing pairs of genes and the correlation between
#'     the genes that have a positive correlation above the threshold
#'     
#' @examples
#' # Using an example correlation matrix that was generated using the pbmc data
#' # that is available in this package. This matrix only includes platelet
#' # cells.
#' corrMatrixPath <- system.file("extdata", 
#'                               "correlationMatrixPlatelet.csv", 
#'                                package = "DiffCoExpr")
#' corrMatrix <- read.csv(corrMatrixPath)
#' # Set the row names and remove the column containing gene names
#' rownames(corrMatrix) <- corrMatrix[, 1]
#' corrMatrix <- corrMatrix[, -1]
#' 
#' coexprNetwork <- getCoexpressionNetwork(correlationMatrix = corrMatrix)
#' coexprNetwork
#' 
#' @importFrom stats quantile
#' @export
getCoexpressionNetwork <- function(correlationMatrix,
                                   thresholdCorrelation = 0.80) {
  
  # --- Check some conditions to stop early if necessary --------------------
  
  if (! is.matrix(correlationMatrix) && ! is.data.frame(correlationMatrix)) {
    stop("Please provide a matrix or data.frame object for correlationMatrix.")
  }
  
  if (! all(sapply(correlationMatrix, is.numeric))) {
    stop("Some of the values in the matrix you provided are not numeric. Remove
         non-numeric values.")
  }
  
  # Check that the input is a square matrix/data.frame
  if (nrow(correlationMatrix) != ncol(correlationMatrix)) {
    stop("The correlation matrix should be a square matrix.")
  }
  
  # Check that thresholdCorrelation is in the correct range
  if (thresholdCorrelation < 0 || thresholdCorrelation > 1) {
    warning("The thresholdCorrelation value you provided is outside (0,1]. 
            Using the default thresholdCorrelation of 0.80 instead.")
    thresholdCorrelation <- 0.80
  }
  
  # --- Filter out pairs of genes --------------------------------------------
  
  # Set correlation to NA where the correlation is less than or equal to 0
  correlationMatrix[correlationMatrix <= 0] <- NA
  
  # Get the row and column indices of the positive correlations
  positivePairIndices <- which(correlationMatrix > 0, arr.ind = TRUE)
  
  # Create a data frame of gene pairs with positive correlations
  positiveCorrelationPairs <- data.frame(
    Gene1 = rownames(correlationMatrix)[positivePairIndices[, 1]],
    Gene2 = colnames(correlationMatrix)[positivePairIndices[, 2]],
    Correlation = correlationMatrix[positivePairIndices]
  )
  
  # Each gene pair is listed twice (e.g. there is a row for Gene1 = A and
  # Gene2 = B and another row for Gene1 = B and Gene2 = A)
  # Take a subset of the gene pairs by only keeping the rows where Gene1
  # comes before Gene2 alphabetically
  positiveCorrelationPairs$Gene1 <- as.character(positiveCorrelationPairs$Gene1)
  positiveCorrelationPairs$Gene2 <- as.character(positiveCorrelationPairs$Gene2)
  positiveCorrelationPairs <- subset(positiveCorrelationPairs, 
                                     Gene1 < Gene2)
  
  # Only keep the pairs that have a correlation higher than the threshold
  cutoffValue <- stats::quantile(positiveCorrelationPairs$Correlation, 
                                 thresholdCorrelation)
  edgeList <- positiveCorrelationPairs[positiveCorrelationPairs$Correlation > 
                                         cutoffValue, ]
  
  return(edgeList)
}

# [END]