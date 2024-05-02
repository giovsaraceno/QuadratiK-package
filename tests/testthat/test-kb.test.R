#' 
#' Test for kb.test
#' 
#' 
#' 
#' @srrstats {G5.1, G5.5} data sets are generated using simple functions with 
#'                        fixed seed
#' @srrstats {G5.2,G5.2a,G5.2b} all the error and warning messages are tested
#' @srrstats {G5.4,G5.4a} correctness tested on simple cases
#' @srrstats {G5.8, G5.8a,G5.8b,G5.8c} edge conditions
#' 
#' @noRd
library(testthat)

## "Tests for kb.test function"

# Test 1: Verify Error on Invalid Method Input
test_that("Error on invalid method input", {
   set.seed(123)
   expect_error(kb.test(x = matrix(rnorm(100), ncol = 2), h=0.5, 
                        method = "invalid_method"), 
                "method must be one of 'bootstrap', 'permutation' or 
                     'subsampling'", fixed=TRUE)
})

# Test 2: Verify Error on Invalid b Input
test_that("Error on invalid b input", {
   set.seed(123)
   expect_error(kb.test(x = matrix(rnorm(100), ncol = 2), h=0.5, b = 10), 
                "b indicates the proportion used for the subsamples in the 
                     subsampling algoritm. It must be in (0,1].", fixed=TRUE)
})

# Test 3: Error on Invalid alternative Input
test_that("Error on invalid alternative input", {
   
   set.seed(123)
   expect_error(kb.test(x = matrix(rnorm(100), ncol= 2), 
                        y = matrix(rnorm(100), ncol= 2), h=0.5, 
                        alternative = "invalid"),
                "The algorithm for selecting the value of h can be performed
                    with respect to the following families of alternatives: 
                    'location', 'scale' or 'skewness'", fixed=TRUE)
})

# Test 4: Error on Invalid centeringType Input
test_that("Error on invalid centeringType input", {
   
   set.seed(123)
   expect_error(kb.test(x = matrix(rnorm(100), ncol= 2), 
                        y = matrix(rnorm(100), ncol= 2), h=0.5, 
                        centeringType = "invalid"),
                "centering must be chosen between 'Param' and 'Nonparam'", 
                fixed=TRUE)
})


# Test 5: Correct handling of vector x input
test_that("Handle vector x input correctly", {
   
   set.seed(123)
   # x is a vector
   result <- kb.test(x = rnorm(10), h=0.5)
   expect_s4_class(result, "kb.test")
   expect_equal(result@method, "Kernel-based quadratic distance Normality test")
   
   # x is a data.frame
   result <- kb.test(x = data.frame(matrix(rnorm(20),ncol=2)), h=0.5)
   expect_s4_class(result, "kb.test")
   
   # x is a matrix
   result <- kb.test(x = matrix(rnorm(20),ncol=2), h=0.5)
   expect_s4_class(result, "kb.test")
   
   # test summary method
   s <- summary(result)
   expect_true("matrix" %in% class(s$summary_tables))
   expect_equal(nrow(s$test_results), 2)
   
   # x is not numeric
   expect_error(kb.test(x = "invalid", h=0.5), "x must be numeric", fixed=TRUE)
   
   # x is not matrix or data.frame
   expect_error(kb.test(x = list(rnorm(10),rnorm(10)), h=0.5), 
                "x must be numeric", fixed=TRUE)
})


# Test 6: Testing main functionality: two-sample test
test_that("Functionality with valid inputs", {
   
   set.seed(123)
   dat <- generate_SN(d = 2, 100, 100, c(0,0),c(0,0), 1, 1, 0)
   x <- dat$X
   y <- dat$Y
   result <- kb.test(x=x, y=y, h=0.5, method = "subsampling", b = 0.5)
   expect_s4_class(result, "kb.test")
  expect_equal(result@method, "Kernel-based quadratic distance two-sample test")
   expect_true(is.numeric(result@Un))
   expect_false(result@H0_Un[1])
   expect_false(result@H0_Un[2])
   
   # test summary method
   s <- summary(result)
   expect_type(s$summary_tables, "list")
   expect_equal(nrow(s$test_results), 2)
   
   # Test parametric centering
   result <- kb.test(x, y, h=0.5, method = "bootstrap", centeringType = "Param")
   expect_s4_class(result, "kb.test")
   expect_equal(result@method,"Kernel-based quadratic distance two-sample test")
   
   ## Test all the methods for the CV computation
   result <- kb.test(x, y, h=0.5, method = "bootstrap")
   expect_s4_class(result, "kb.test")
   
   result <- kb.test(x, y, h=0.5, method = "permutation")
   expect_s4_class(result, "kb.test")
   
   # Test if y is a data.frame
   y <- data.frame(matrix(rnorm(200), ncol = 2))
   result <- kb.test(x, y, h=0.5, method = "subsampling", b = 0.5)
   expect_s4_class(result, "kb.test")
   
   # Test additional errors
   y <- matrix(rnorm(90), ncol = 3)
   expect_error(kb.test(x, y, h=0.5),
                "'x' and 'y' must have the same number of columns.", 
                fixed=TRUE)
   
})

# Test 7: Testing main functionality: k-sample test
test_that("Functionality with valid inputs", {
   set.seed(123)
   x <- matrix(rnorm(100), ncol = 2)
   y <- rep(c(1,2), each=25)
   result <- kb.test(x, y, h=0.5, method = "subsampling", b = 0.5)
   expect_s4_class(result, "kb.test")
   expect_equal(result@method, "Kernel-based quadratic distance k-sample test")
   expect_true(is.numeric(result@Un))
   expect_false(result@H0_Un[1])
   expect_false(result@H0_Un[2])
   
   # test show method
   output <- capture.output(show(result))
   expect_true(any(grepl("H0 is rejected: ", output)))
   expect_true(any(grepl("Test Statistic: ", output)))
   
   # test summary method
   s <- summary(result)
   expect_type(s$summary_tables, "list")
   expect_equal(nrow(s$test_results), 2)
   
   # Test all the methods for computing the CV
   result <- kb.test(x, y, h=0.5, method = "bootstrap")
   expect_s4_class(result, "kb.test")
   result <- kb.test(x, y, h=0.5, method = "permutation")
   expect_s4_class(result, "kb.test")
   
   # Test additional errors
   y <- rep(c(1,2), each=20)
   expect_error(kb.test(x, y, h=0.5),
                "'x' and 'y' must have the same number of rows.", fixed=TRUE)
   
})

# Test 8: Testing selection of h
test_that("Selection of h from kb.test", {
   
   set.seed(123)
   x <- matrix(rnorm(100), ncol = 2)
   y <- rep(c(1,2), each=25)
   
   result <- kb.test(x, method = "subsampling", mu_hat = c(0,0),
                     Sigma_hat = diag(2), b = 0.5)
   expect_s4_class(result, "kb.test")
   expect_equal(result@method, "Kernel-based quadratic distance Normality test")
   expect_equal(class(result@h$h_sel), "numeric")
   
   result <- kb.test(x, y, method = "subsampling", b = 0.5)
   expect_s4_class(result, "kb.test")
   expect_equal(class(result@h$h_sel), "numeric")
   
})


