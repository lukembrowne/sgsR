library(sgsR)
context("Spatial functions")

test_that("distance calculations are correct", {
  expect_equal(calcPairwiseDist(10, 10, 1), as.matrix(0))
})


