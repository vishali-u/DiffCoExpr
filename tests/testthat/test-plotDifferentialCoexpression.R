library(DiffCoExpr)

gene1 <- "Gene1"
gene2 <- "Gene2"

# Define data for testing
testMatrixA <- matrix(c(12, 22, 36, 47, 45, 66, 77, 84, 92, 10,
                        11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                        21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                        31, 32, 33, 34, 35, 36, 37, 38, 39, 45,
                        41, 42, 43, 44, 45, 46, 47, 48, 49, 52,
                        51, 52, 53, 54, 55, 56, 53, 58, 59, 60,
                        61, 62, 63, 64, 65, 66, 62, 68, 69, 70,
                        71, 72, 73, 74, 75, 76, 77, 71, 69, 80,
                        81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                        91, 92, 93, 94, 95, 96, 97, 98, 99, 100), 
                      nrow = 10, ncol = 10)
geneNames <- paste("Gene", 1:10, sep = "")
cellBarcodes <- paste("Cell", 1:10, sep = "")
rownames(testMatrixA) <- geneNames
colnames(testMatrixA) <- cellBarcodes

testMatrixB <- matrix(c(20, 12, 23, 44, 45, 66, 73, 82, 19, 10,
                        11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                        21, 22, 23, 24, 25, 26, 27, 28, 29, 32,
                        31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                        41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                        51, 52, 53, 54, 55, 56, 57, 58, 59, 63,
                        61, 62, 63, 64, 65, 66, 67, 68, 69, 50,
                        71, 72, 73, 74, 75, 76, 77, 78, 79, 70,
                        81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                        91, 92, 93, 94, 95, 96, 97, 98, 99, 100), 
                      nrow = 10, ncol = 10)
geneNames <- paste("Gene", 1:10, sep = "")
cellBarcodes <- paste("Cell", 1:10, sep = "")
rownames(testMatrixB) <- geneNames
colnames(testMatrixB) <- cellBarcodes

smallInvalidMatrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                               11, 12, 13, 14, 15, '16', 17, 18, 19, 20), 
                     nrow = 4, ncol = 5)
geneNames <- paste("Gene", 1:4, sep = "_")
cellBarcodes <- paste("Cell", 1:5, sep = "_")
rownames(smallInvalidMatrix) <- geneNames
colnames(smallInvalidMatrix) <- cellBarcodes

testDiffCoExpTable <- data.frame(Gene1 = c("Gene1", "Gene3"), 
                                 Gene2 = c("Gene2", "Gene4"),
                                 CorrelationA = c(0.8, 0.7), 
                                 CorrelationB = c(0.5, 0.6),
                                 foldChange = c(1.2, 3.4))

test_that("Valid Input", {

  plot <- plotDifferentialCoexpression(diffCoexpTable = testDiffCoExpTable,
                                       expressionMatrixA = testMatrixA, 
                                       expressionMatrixB = testMatrixB, 
                                       gene1 = gene1, 
                                       gene2 = gene2)
  expect_true(ggplot2::is.ggplot(plot))
})

test_that("Invalid Input", {
  
  # Check that function throws an error with a gene not present in the data
  expect_error(
    plotDifferentialCoexpression(diffCoexpTable = testDiffCoExpTable,
                                 expressionMatrixA = testMatrixA, 
                                 expressionMatrixB = testMatrixB, 
                                 gene1 = "invalid gene", 
                                 gene2 = gene2))
  
  # Check that function throws an error when trying to use the same gene
  expect_error(
    plotDifferentialCoexpression(diffCoexpTable = testDiffCoExpTable,
                                 expressionMatrixA = testMatrixA, 
                                 expressionMatrixB = testMatrixB, 
                                 gene1 = gene1, 
                                 gene2 = gene1))
  
  # Check that function throws an error when an empty dataframe is used for
  # either expression matrix
  emptyDf <- data.frame()
  expect_error(
    plotDifferentialCoexpression(diffCoexpTable = testDiffCoExpTable,
                                 expressionMatrixA = emptyDf, 
                                 expressionMatrixB = testMatrixB, 
                                 gene1 = gene1, 
                                 gene2 = gene2))
  expect_error(
    plotDifferentialCoexpression(diffCoexpTable = emptyDf,
                                 expressionMatrixA = testMatrixA, 
                                 expressionMatrixB = testMatrixB, 
                                 gene1 = gene1, 
                                 gene2 = gene2))
  
  # Check that function throws an error when the genes are found in the 
  # expression matrix but are not coexpressed
  expect_error(
    plotDifferentialCoexpression(diffCoexpTable = emptyDf,
                                 expressionMatrixA = testMatrixA, 
                                 expressionMatrixB = testMatrixB, 
                                 gene1 = "Gene7", 
                                 gene2 = gene2))
  
  # Check that the function throws an error when the expression matrix 
  # contains non-numeric values
  expect_error(
    plotDifferentialCoexpression(diffCoexpTable = testDiffCoExpTable,
                                 expressionMatrixA = smallInvalidMatrix, 
                                 expressionMatrixB = testMatrixB, 
                                 gene1 = gene1, 
                                 gene2 = gene2))
})