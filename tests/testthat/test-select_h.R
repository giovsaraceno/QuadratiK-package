#' Test for select_h
#' 
#' 
#' 
#' @srrstats {G5.1, G5.5} data sets are generated using simple functions with 
#'                        fixed seed
#' @srrstats {G5.2,G5.2a,G5.2b} all the error and warning messages are tested
#' @srrstats {G5.8, G5.8a,G5.8b,G5.8c} edge conditions
#' 
#' @noRd
library(testthat)

## "Tests for select_h function"

# Test 1: Verify Error on Invalid Input
test_that("Error on invalid method input", {
   
   expect_error(select_h(x = matrix(rnorm(100), ncol = 2), 
                         alternative = "location", n_cores = "invalid"), 
                "n_cores must be a numeric value", fixed=TRUE)
   
   set.seed(123)
   expect_error(select_h(x = matrix(rnorm(100), ncol = 2), 
                        alternative = "invalid"), 
            "The alternative argument should be one of 'location', 'scale' or 
           'skewness'", fixed=TRUE)
   
   # x is not numeric
   expect_error(select_h(x = "invalid", alternative="skewness"), 
                "x must be numeric", fixed=TRUE)
   
   x <- matrix(rnorm(100), ncol = 2)
   x[1,] <- NA
   expect_error(select_h(x = x, alternative="skewness"), 
                "There are missing values in x!", fixed=TRUE)
   x[1,] <- Inf
   expect_error(select_h(x = x, alternative="skewness"), 
         "There are undefined values in x, that is Nan, Inf, -Inf", fixed=TRUE)
   
   x <- matrix(rnorm(100), ncol = 2)
   y <- matrix(rnorm(100), ncol = 2)
   y[1,] <- NA
   expect_error(select_h(x = x, y = y, alternative="skewness"), 
                "There are missing values in y!", fixed=TRUE)
   y[1,] <- Inf
   expect_error(select_h(x = x, y = y, alternative="skewness"), 
         "There are undefined values in y, that is Nan, Inf, -Inf", fixed=TRUE)
   
   x <- matrix(rnorm(100), ncol = 2)
   y <- matrix(rnorm(90), ncol = 3)
   expect_error(select_h(x, y, alternative="skewness"),
                "'x' and 'y' must have the same number of columns.", 
                fixed=TRUE)
   
   x <- matrix(rnorm(100), ncol = 2)
   y <- rep(c(1,2), each=20)
   expect_error(select_h(x, y, alternative="skewness"),
                "'x' and 'y' must have the same number of rows.", fixed=TRUE)

   y <- rep(c(1,2), each=25)
   expect_error(select_h(x, y, alternative="skewness", delta_dim=c(1,2,1)),
   "delta_dim must be 1 or a numeric vector of length equal to the 
              number of columns of pooled.", fixed=TRUE)
   
   x <- rnorm(50)
   y <- rnorm(50, 10)
   expect_error(select_h(x, y, alternative="skewness", 
                         h_values = c("a","b","c")), 
                "h_values must be a numeric vector", fixed=TRUE)
   expect_error(select_h(x, y, alternative="skewness", 
                         delta = c("a","b","c")), 
                "delta must be a numeric vector", fixed=TRUE)
   
})


## # Test 1: test for select_h
test_that("Select h", {
   set.seed(123)
   # normality
   result <- select_h(x = as.data.frame(matrix(rnorm(20),ncol=2)), 
                      alternative="location")
   expect_equal(class(result$h_sel), "numeric")
   expect_equal(class(result$power), "data.frame")
   
   result <- select_h(x = as.data.frame(matrix(rnorm(20),ncol=2)), 
                      alternative="scale", n_cores = 2)
   expect_equal(class(result$h_sel), "numeric")
   expect_equal(class(result$power), "data.frame")

   # two-sample
   result <- select_h(x = matrix(rnorm(20),ncol=2),
                      y = as.data.frame(matrix(rnorm(20),ncol=2)), 
                      alternative="location")
   expect_equal(class(result$h_sel), "numeric")
   expect_equal(class(result$power), "data.frame")
   
   result <- select_h(x = matrix(rnorm(20),ncol=2),
   y = matrix(rnorm(20),ncol=2), alternative="skewness")
   expect_equal(class(result$h_sel), "numeric")
   expect_equal(class(result$power), "data.frame")
   
   result <- select_h(x = matrix(rnorm(20),ncol=2),
                      y = matrix(rnorm(20),ncol=2), alternative="scale")
   expect_equal(class(result$h_sel), "numeric")
   expect_equal(class(result$power), "data.frame")

   # k-sample
   result <- select_h(x = matrix(rnorm(30),ncol=2), y = rep(c(1,2,3),each=5),
   alternative="scale")
   expect_equal(class(result$h_sel), "numeric")
   expect_equal(class(result$power), "data.frame")
   
   result <- select_h(x = matrix(rnorm(30),ncol=2), y = rep(c(1,2,3),each=5),
                      alternative="location")
   expect_equal(class(result$h_sel), "numeric")
   expect_equal(class(result$power), "data.frame")
})


