library(DiffCoExpr)

test_that("valid input", {
  srat <- prepareData(
    geneMatrixPath = system.file("extdata",
                                 "filtered_gene_bc_matrices",
                                 package =  "DiffCoExpr"),
    cellTypesPath = system.file("extdata",
                                "cellTypes.csv",
                                package =  "DiffCoExpr")
  )
  
  expect_identical(levels(srat), 
                   c("Naive CD4+ T", "CD14+ Mono", "Memory CD4+", "B", 
                     "CD8+ T", "FCGR3A+ Mono", "NK", "DC", "Platelet"))
})

test_that("invalid input", {
  
  # Gene matrix path is invalid
  expect_error(prepareData(geneMatrixPath = "",
                           cellTypesPath = NULL))
  
  # Cell types path is invalid (does not exist and is not NULL)
  expect_error(prepareData(geneMatrixPath = 
                             system.file("extdata",
                                         "filtered_gene_bc_matrices",
                                          package =  "DiffCoExpr"), 
                           cellTypesPath = ""))
})


