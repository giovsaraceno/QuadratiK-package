library(testthat)

context("Tests for kb.test function")

# Test 1: Verify Error on Invalid Method Input
test_that("Error on invalid method input", {
   
   expect_error(kb.test(x = matrix(rnorm(100), ncol = 2), h=0.5, method = "invalid_method"),
                "method must be one of 'bootstrap', 'permutation' or 'subsampling'")
})

# Test 2: Verify Error on Invalid b Input
test_that("Error on invalid b input", {
   #expect_error(kb.test(x = matrix(rnorm(100), ncol = 2), h=0.5, b = -1), "b indicates the proportion used for the subsamples in the subsampling algoritm. It must be in (0,1]")
})

# Test 3: Error on Invalid alternative Input
test_that("Error on invalid alternative input", {
   #expect_error(kb.test(x = matrix(rnorm(100), ncol= 2), h=0.5, alternative = "invalid"),"alternatives must be one of 'skewness', 'location', 'scale'")
})


# Test 5: Correct handling of vector x input
test_that("Handle vector x input correctly", {
   expect_s4_class(kb.test(x = rnorm(10), h=0.5), "kb.test")
})


# Test 6: Testing main functionality with valid inputs
test_that("Functionality with valid inputs", {
   x <- matrix(rnorm(100), ncol = 2)
   result <- kb.test(x, h=0.5, method = "subsampling", b = 0.5, alternative = "skewness")
   expect_s4_class(result, "kb.test")
   expect_true(is.numeric(result@Dn))
   expect_false(result@H0)
   
})

# Add more tests for different scenarios, parameter combinations, edge cases etc.

test_that("Normal kernel-based tests work", {
   
   #------------------------------------------------------
   ## Normality test statistics
   ## 
   size <- 100
   d <- 2
   set.seed(081423)
   x <- matrix(rnorm(size*d),ncol=d,nrow=size)
   h <- 0.5
   
   norm_x <- kb.test(x=x,h=h)
   expect_true(class(norm_x)=="kb.test")
   expect_true(is.numeric(norm_x@Dn))
   expect_false(norm_x@H0)
   
   #------------------------------------------------------
   ## kernel-based quadratic distance Two-sample test
   
   y <- matrix(rnorm(size*d, mean=2),ncol=d,nrow=size)
   
   two_kb_sub <- kb.test(x,y,h=1.5, method = "subsampling", b=0.5)
   expect_true(class(two_kb_sub)=="kb.test")
   expect_true(two_kb_sub@H0)
   
   #------------------------------------------------------
   ## kernel-based quadratic distance k-sample test
   
   k <- 3
   x <- matrix(rnorm(size*d*k),ncol=d)
   y <- rep(1:k,each=size)
   
   k_kb_sub <- kb.test(x,y,h=1.5, method = "subsampling", b=0.5)
   expect_true(class(k_kb_sub)=="kb.test")
   expect_false(k_kb_sub@H0)
   
})
