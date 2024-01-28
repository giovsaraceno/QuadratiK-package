test_that("Normal kernel-based tests work", {
   
   size <- 100
   d <- 2
   set.seed(081423)
   x <- matrix(rnorm(size*d),ncol=d,nrow=size)
   y <- matrix(rnorm(size*d),ncol=d,nrow=size)
   z <- rbind(x,y)
   h <- 0.5
   H <- h^2*diag(d)
   #------------------------------------------------------
   ## Two-sample statistics
   
   # expect_true(is.numeric(stat2sample(x,y,h,matrix(0,nrow=1),diag(1),"Nonparam")))
   # expect_true(is.numeric(stat2sample(x,y,h,matrix(colMeans(z),nrow=1),cov(z),"Param")))
   # ## Critical value
   
   ## kernel-based quadratic distance Two-sample test
   #two_kb_sub <- kb.test(x,y,h=0.5, method = "subsampling", b=0.9)
   #expect_true(is.numeric(two_kb_sub@Dn))
   #expect_false(two_kb_sub@H0)
   
   #-------------------------------------------------------
   ## Test of Normality
   #norm_x <- kb.test(x, h=0.5)
   #expect_false(norm_x@H0)
   
})
