library(DiffCoExpr)

test_that("Expression matrix has right number of columns", {
  srat <- prepareData(
    gene_bc_matrix_path = "inst/extdata/pbmc/filtered_gene_bc_matrices",
    cell_types_path = "inst/extdata/pbmc/cell_types.csv"
  )
  
  filtered_expr_matrix <- getExpressionMatrix(
    srat = srat,
    cell_type = "B"
  )
  
  dimensions <- dim(filtered_expr_matrix)
  # Check for the right number of B cells
  expect_equal(dimensions[2], 344)
  
  expr_matrix <- getExpressionMatrix(
    srat = srat
  )
  
  # Started with 2700 cells, so after filtering, there should be less than 2700
  expect_lt(dim(expr_matrix)[2], 2700)
})