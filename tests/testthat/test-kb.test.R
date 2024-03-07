library(testthat)

## "Tests for kb.test function"

# Test 1: Verify Error on Invalid Method Input
test_that("Error on invalid method input", {
   
   expect_error(kb.test(x = matrix(rnorm(100), ncol = 2), h=0.5, method = "invalid_method"), "method must be one of 'bootstrap', 'permutation' or 'subsampling'", fixed=TRUE)
})

# Test 2: Verify Error on Invalid b Input
test_that("Error on invalid b input", {
   expect_error(kb.test(x = matrix(rnorm(100), ncol = 2), h=0.5, b = 10), "b indicates the proportion used for the subsamples in the subsampling algoritm. It must be in (0,1]", fixed=TRUE)
})

# Test 3: Error on Invalid alternative Input
test_that("Error on invalid alternative input", {
   
   expect_error(kb.test(x = matrix(rnorm(100), ncol= 2), y = matrix(rnorm(100), ncol= 2), h=0.5, alternative = "invalid"),"The algorithm for selecting the value of h can be performed with respect to the following families of alternatives: 'location', 'scale' or 'skewness'", fixed=TRUE)
})

# Test 4: Error on Invalid centeringType Input
test_that("Error on invalid centeringType input", {
   
   expect_error(kb.test(x = matrix(rnorm(100), ncol= 2), y = matrix(rnorm(100), ncol= 2), h=0.5, centeringType = "invalid"),"centering must be chosen between 'Param' and 'Nonparam'", fixed=TRUE)
})


# Test 5: Correct handling of vector x input
test_that("Handle vector x input correctly", {
   result <- kb.test(x = rnorm(10), h=0.5)
   expect_s4_class(result, "kb.test")
   expect_equal(result@method, "Kernel-based quadratic distance Normality test")
})


# Test 6: Testing main functionality: two-sample test
test_that("Functionality with valid inputs", {
   x <- matrix(rnorm(100), ncol = 2)
   y <- matrix(rnorm(100), ncol = 2)
   result <- kb.test(x, y, h=0.5, method = "subsampling", b = 0.5)
   expect_s4_class(result, "kb.test")
   expect_equal(result@method, "Kernel-based quadratic distance two-sample test")
   expect_true(is.numeric(result@Dn))
   expect_false(result@H0)
   
})

# Test 7: Testing main functionality: k-sample test
test_that("Functionality with valid inputs", {
   x <- matrix(rnorm(100), ncol = 2)
   y <- rep(c(1,2), each=25)
   result <- kb.test(x, y, h=0.5, method = "subsampling", b = 0.5)
   expect_s4_class(result, "kb.test")
   expect_equal(result@method, "Kernel-based quadratic distance k-sample test")
   expect_true(is.numeric(result@Dn))
   expect_false(result@H0)
   
})


