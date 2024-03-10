library(testthat)

## "Tests for pk.test function"

# Test 1: Error on Invalid x Input
test_that("Error on invalid x input", {
   
   expect_error(pk.test(x = "Not matrix", rho = 0.5),
                "x must be numeric", fixed=TRUE)
   
   expect_error(pk.test(x = rnorm(50), rho = 0.5),
                "x must be a matrix or a data.frame with dimension greater than 1.", fixed=TRUE)
})

# Test 2: Error on Invalid Quantile Input
test_that("Error on invalid Quantile input", {
   
   expect_error(pk.test(x = matrix(rnorm(100), ncol=2), rho = 0.5, 
                        Quantile = 1.1), "Quantile must be in (0,1]", 
                fixed=TRUE)
   expect_error(pk.test(x = matrix(rnorm(100), ncol=2), rho = 0.5, 
                        Quantile = -1), "Quantile must be in (0,1]", 
                fixed=TRUE)
})

# Test 3: Error on Invalid rho Input
test_that("Error on invalid rho input", {
   
   expect_error(pk.test(x = matrix(rnorm(100), ncol=2), rho = -0.1),
                "rho must be in (0,1)", fixed=TRUE)
   
   expect_error(pk.test(x = matrix(rnorm(100), ncol=2), rho = 2),
                "rho must be in (0,1)", fixed=TRUE)
})

# Test 4: Handling Vector x Input
test_that("Handle vector x input correctly", {
   
   expect_s4_class(pk.test(x = matrix(rnorm(100), ncol = 2), rho = 0.5), 
                   "pk.test")
   expect_s4_class(pk.test(x = data.frame(matrix(rnorm(100), ncol = 2)), 
                           rho = 0.5), "pk.test")
})

# Test 5: Main Functionality with Valid Inputs
test_that("Functionality with valid inputs", {
   
   size <- 100
   d <- 3
   set.seed(012924)
   x_sp <- sample_hypersphere(d, n_points=size)
   
   unif_test <- pk.test(x_sp,rho=0.8)
   
   expect_s4_class(unif_test, "pk.test")
   expect_true(is.numeric(unif_test@Un))
   expect_true(is.numeric(unif_test@Vn))
   expect_false(unif_test@H0_Un)
   
   # test show method
   output <- capture.output(show(unif_test))
   expect_true(any(grepl("H0 is rejected: ", output)))
   expect_true(any(grepl("Statistic Un: ", output)))
   expect_true(any(grepl("Statistic Vn: ", output)))
   
   # test summary method
   s <- summary(unif_test)
   expect_type(s$summary_tables, "list")
   expect_equal(nrow(s$test_results), 2)
   
})


