library(DiffCoExpr)

test_that("valid input", {
  srat <- prepareData(
    gene_bc_matrix_path = "inst/extdata/pbmc/filtered_gene_bc_matrices",
    cell_types_path = "inst/extdata/pbmc/cell_types.csv"
  )
  
  expr_matrix <- getExpressionMatrix(
    srat = srat,
    cell_type = "B"
  )
  
  corr_matrix <- getCorrelationMatrix(
    expression_matrix = expr_matrix
  )
  
  expect_identical(dim(corr_matrix)[1], dim(corr_matrix)[2])
  
})