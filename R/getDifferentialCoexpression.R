#' Given two coexpression networks, identify the pairs that are coexpressed
#' in both networks. Return a dataframe of all the genes that are coexpressed
#' in both networks, and quantify differential coexpression using log fold 
#' change.
#' 
#' @param networkA An edge list for a coexpression network. Contains the
#'     correlation coefficient between pairs of genes. 
#'
#' @param networkB An edge list for a different co-expression network. 
#'     Contains the correlation coefficient between pairs of genes. 
#'     
#' @param thresholdLogFC A minimum log2 fold change to be used for filtering
#'     out gene pairs    
#'     
#' @examples
#' # Using an example data that was generated using the pbmc data that is
#' # available in this package. 
#' 
#' # Condition A: Platelet cells
#' coexprNetPlateletPath <- system.file("extdata", 
#'                                      "coexpressionNetworkPlatelet.csv", 
#'                                      package = "DiffCoExpr")
#' coexprNetPlatelet <- read.csv(coexprNetPlateletPath)
#' # Delete the first column 
#' coexprNetPlatelet <- coexprNetPlatelet[, -1]
#' 
#' # Condition B: DC Cells
#' coexprNetDCPath <- system.file("extdata", 
#'                                "coexpressionNetworkDC.csv", 
#'                                package = "DiffCoExpr")
#' coexprNetDC <- read.csv(coexprNetDCPath)
#' # Delete the first column 
#' coexprNetDC <- coexprNetDC[, -1]
#' 
#' diffCoexpressed <- getDifferentialCoexpression(networkA = coexprNetPlatelet, 
#'                                                networkB = coexprNetDC) 
#' diffCoexpressed
#' 
#' @return a dataframe where each row contains a pair of genes and the 
#'     correlation of the genes in each network
#' 
#' @importFrom dplyr inner_join mutate filter
#' @export
getDifferentialCoexpression <- function(networkA, 
                                        networkB, 
                                        thresholdLogFC = 0.20) {
  
  # --- Checks for invalid input -------------------------------------------
  
  if (nrow(networkA) == 0 || nrow(networkB) == 0) {
    stop("One of the coexpression networks you provided is empty. Provide a 
          non-empty network.")
  }
  
  if ((! all(sapply(networkA$Correlation, is.numeric))) || 
      (! all(sapply(networkB$Correlation, is.numeric)))) {
    stop("Some of the correlation values in the coexpression network you 
          provided are not numeric. Remove non-numeric values.")
  }
  
  # --- Compare networks ---------------------------------------------------
  
  # Merge the edge lists for both networks
  mergedList <- networkA %>%
    dplyr::inner_join(networkB, by = c("Gene1", "Gene2"), 
                      suffix = c("A", "B"))
  
  # Calculate the fold change
  mergedList <- mergedList %>%
    dplyr::mutate(foldChange = log2(CorrelationB / CorrelationA)) %>%
    dplyr::filter(abs(foldChange) > thresholdLogFC)
  
  return(mergedList)
}



