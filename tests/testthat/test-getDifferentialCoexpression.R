library(DiffCoExpr)

networkA <- data.frame(Gene1 = c("Gene1", "Gene3"), 
                       Gene2 = c("Gene2", "Gene4"),
                       Correlation = c(0.8, 0.7))
networkB <- data.frame(Gene1 = c("Gene1", "Gene5"),
                       Gene2 = c("Gene2", "Gene6"), 
                       Correlation = c(0.1, 0.8))
invalidNetwork <- data.frame(Gene1 = "Gene6", 
                             Gene2 = "Gene8", 
                             Correlation = "0.2")

test_that("Valid Input", {
  
  df <- getDifferentialCoexpression(networkA = networkA, 
                                    networkB = networkB)
  expect_equal(nrow(df), 1)
  
})

test_that("Invalid Input", {
  
  # Check that function throws an error when using an empty data frame
  expect_error(
    getDifferentialCoexpression(networkA = networkA, 
                                networkB = data.frame()))
  
  expect_error(
    getDifferentialCoexpression(networkA = data.frame(), 
                                networkB = networkB))
  
  # Check that function throws an error when using a network with non-numerical
  # values for the correlation
  expect_error(
    getDifferentialCoexpression(networkA = invalidNetwork, 
                                networkB = networkB))
  
})