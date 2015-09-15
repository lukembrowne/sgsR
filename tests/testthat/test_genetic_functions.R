library(sgsR)
context("Genetic functions")

test_that("Allele frequencies", {
  expect_equal(calcAlleleFreqPop(c(0,1,2), c(1,2,3), Nallele = 4),
               c(1/6, 2/6, 2/6, 1/6 ))
  expect_equal(calcAlleleFreqCppInd(0,1, Nallele = 2),
               c(1/2, 1/2))
})

#
# ref_gen <- matrix(rep(1/3, 3), nrow = 1)
# alfreq1 <- matrix(c(1/2, 1/2, 0), nrow = 1)
# alfreq2 <- matrix(c(0, 0, 1), nrow = 1)
# Nloci = 1
# Nallele = 3
# Ngenecopies = 4
#
# test_that("Fij coefficient", {
#   expect_equal(calcFijPairwiseCpp(ref_gen, alfreq1, alfreq2,
#                                   Nloci, Nallele, Ngenecopies))
#
#
#
# })
#
#
