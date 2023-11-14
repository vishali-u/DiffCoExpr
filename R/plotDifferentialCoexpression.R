#' Visualize the distribution of differential co-expression that arise from
#' differences in co-expression levels between two conditions
#' 
#' @param network_A An edge list for a co-expression network. Contains the
#'     correlation coefficient between pairs of genes
#' @param network_B An edge list for a different co-expression network. 
#'     Contains the correlation coefficient between pairs of genes 
#'     
#' @return A histogram that illustrates the distribution of log fold changes.
#'     The fold change refers to the degree to which co-expression strength
#'     (e.g., correlation coefficient) changes between gene pairs across
#'     co-expression networks constructed at different conditions (e.g. cell
#'     types)
#'
#' @export
#' 
plotDifferentialCoexpression <- function(network_A, network_B) {
  # Merge the edge lists for both networks
  merged_list <- merge(network_A, network_B, 
                       by = c("gene1", "gene2"), 
                       all = TRUE,
                       suffixes = c("_A", "_B"))
  
  # Calculate the fold change
  merged_list$fold_change <- with(merged_list, 
                                  log2(weight_B / weight_A))

  # TODO: how to do a test for significance and multiple testing correction,
  # correlation coefficients are not normally distributed so cannot use most
  # parametric tests like t-test
  
  # TODO: if a test for significance is possible, change the histogram to a 
  # different graph (maybe scatterplot) that includes significance 
  # (such as graphing fold change vs p-value)
  fold_change_distribution <- hist(merged_list$fold_change, 
                                   main = "Distribution of Log Fold Changes", 
                                   xlab = "Log Fold Change", 
                                   breaks = 30)
  
  return(fold_change_distribution)
}