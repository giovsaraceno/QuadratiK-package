library(testthat)
test_that("Poisson kernel-based tests work", {
   
   #------------------------------------------------------
   ## Uniformity test 
   ## 
   size <- 100
   d <- 3
   set.seed(012924)
   x_sp <- sample_hypersphere(d, n_points=size)
   
   unif_test <- pk.test(x_sp,rho=0.8)
   
   expect_true(class(unif_test)=="pk.test")
   expect_true(is.numeric(unif_test@Un))
   expect_true(is.numeric(unif_test@Vn))
   expect_false(unif_test@H0_Un)
   
   
})
