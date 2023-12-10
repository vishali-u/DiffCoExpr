library(DiffCoExpr)
set.seed(0)

test_that("Input that will not throw errors", {
  srat <- prepareData(
    geneMatrixPath = system.file("extdata",
                                 "filtered_gene_bc_matrices",
                                 package =  "DiffCoExpr"),
    cellTypesPath = system.file("extdata",
                                "cellTypes.csv",
                                package =  "DiffCoExpr")
  )
  
  filteredExprMatrix <- getExpressionMatrix(
    srat = srat,
    cellType = "B"
  )
  
  corrMatrix <- getCorrelationMatrix(
    expressionMatrix = filteredExprMatrix
  )
  
  # Check that the correlation matrix is a square matrix
  expect_identical(dim(corrMatrix)[1], dim(corrMatrix)[2])
  
  # Check that genes were filtered out
  expect_lt(dim(filteredExprMatrix)[1], dim(corrMatrix)[1])
  
  expect_warning(getCorrelationMatrix(expressionMatrix = filteredExprMatrix,
                                      minPt = 1.5))
})

test_that("Input that will throw errors", {
  
  # Test on a 10x10 matrix with some non-numeric values
  testMatrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                         11, 12, 13, 14, 15, '16', 17, 18, 19, 20,
                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                         51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                         61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
                         71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
                         81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
                         91, 92, 93, 94, 95, 96, 97, 98, 99, 100), 
                       nrow = 10, ncol = 10)
  geneNames <- paste("Gene", 1:10, sep = "_")
  cellBarcodes <- paste("Cell", 1:10, sep = "_")
  rownames(testMatrix) <- geneNames
  colnames(testMatrix) <- cellBarcodes
  expect_error(getCorrelationMatrix(expressionMatrix = testMatrix))
  
  # Test on an object that is not a data frame or a matrix
  expect_error(getCorrelationMatrix(expressionMatrix = "not a valid type"))
  
  # Test on a 1x2 matrix
  smallTestMatrix <- matrix(c(1, 2), nrow = 2, ncol = 1)
  geneNames <- paste("Gene", 1:1, sep = "_")
  cellBarcodes <- paste("Cell", 1:2, sep = "_")
  rownames(smallTestMatrix) <- geneNames
  colnames(smallTestMatrix) <- cellBarcodes
  expect_error(getCorrelationMatrix(expressionMatrix = smallTestMatrix))
  
})

