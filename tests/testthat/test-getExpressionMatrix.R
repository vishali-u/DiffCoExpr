library(DiffCoExpr)

# Create the Seurat object once to be used in both tests
srat <- prepareData(
  geneMatrixPath = system.file("extdata",
                               "filtered_gene_bc_matrices",
                               package =  "DiffCoExpr"),
  cellTypesPath = system.file("extdata",
                              "cell_types.csv",
                              package =  "DiffCoExpr")
)

test_that("Valid Input", {
  filteredExprMatrix <- getExpressionMatrix(
    srat = srat,
    cellType = "B"
  )
  
  exprMatrix <- getExpressionMatrix(
    srat = srat
  )
  
  # Check for the right number of B cells
  expect_equal(dim(filteredExprMatrix)[2], 344)
  
  # Check that the returned matrix has an approriate number of columns
  # Started with 2700 cells, so after filtering, there should be less than 2700
  expect_lt(dim(exprMatrix)[2], 2700)
})

test_that ("Invalid Input", {
  expect_error(getExpressionMatrix(srat = "", 
                                   cellType = "B"))
  expect_error(getExpressionMatrix(srat = srat, 
                                   cellType = "invalid cell type"))
})

