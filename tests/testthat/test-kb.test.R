library(testthat)
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
