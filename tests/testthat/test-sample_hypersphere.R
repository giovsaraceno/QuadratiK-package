library(testthat)
test_that("Random Generation from PKBD works", {
   
   #------------------------------------------------------
   ## Generating data point on the Sphere
   ## and computing the densities
   size <- 100
   d = 3
   x_sp <- sample_hypersphere(d, size)
   
   expect_equal(dim(x_sp),c(size,d))
   expect_true(any(rowSums(x_sp^2) ==1))
   
   
})