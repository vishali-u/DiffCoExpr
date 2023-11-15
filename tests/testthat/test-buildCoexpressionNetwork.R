library(DiffCoExpr)

test_that("All edges have a positive correlation", {
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
  
  network <- buildCoexpressionNetwork(
    correlation_matrix = corr_matrix
  )
  
  expect_true(network$weight)
  
})