library(DiffCoExpr)

test_that("valid input", {
  srat <- prepareData(
    gene_bc_matrix_path = "inst/extdata/pbmc/filtered_gene_bc_matrices",
    cell_types_path = "inst/extdata/pbmc/cell_types.csv"
  )
  
  #expect_type(srat, "Seurat")
  expect_identical(levels(srat), 
                   c("Naive CD4+ T", "CD14+ Mono", "Memory CD4+", "B", 
                     "CD8+ T", "FCGR3A+ Mono", "NK", "DC", "Platelet"))
})
