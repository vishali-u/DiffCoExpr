library(DiffCoExpr)

test_that("Valid Input", {
  
  coexprNetPath <- system.file("extdata",
                               "coexpressionNetworkPlatelet.csv",
                               package = "DiffCoExpr")
  coexprNet <- read.csv(coexprNetPath)
  # Delete the first column
  coexprNet <- coexprNet[, -1]
  
  expect_success(plotCoexpressionNetwork(edgeList = coexprNet))
})

test_that("Invalid Input", {
  
  # Check that function fails on non-dataframe objects
  expect_error(plotCoexpressionNetwork(edgeList = "invalid type"))
  
  # Check that function fails on an empty dataframe
  emptyDf <- data.frame
  expect_error(plotCoexpressionNetwork(edgeList = emptyDf))
})