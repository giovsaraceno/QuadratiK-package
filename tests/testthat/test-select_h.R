library(testthat)
test_that("Select h alhorithm works", {
   
   #------------------------------------------------------
   ## Selection of h for two-sample test
   ## 
   size <- 100
   d <- 2
   set.seed(012924)
   x <- matrix(rnorm(size),ncol=d)
   y <- matrix(rnorm(size),ncol=d)
   h_sel <- select_h(x,y,"skewness")
   
   expect_true(is.numeric(h_sel$h_sel))
   expect_gt(max(h_sel$power), 0.5)
   
   
   #------------------------------------------------------
   ## Selection of h for k-sample test
   ## 
   size <- 150
   d <- 2
   k <- 3
   set.seed(012924)
   x <- matrix(rnorm(size*d*k),ncol=d)
   y <- rep(1:k,each=size)
   h_sel <- select_h(x,y,"location")
   
   expect_true(is.numeric(h_sel$h_sel))
   expect_gt(max(h_sel$power), 0.5)
   

})
