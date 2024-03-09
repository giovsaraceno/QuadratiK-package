library(testthat)

#------------------------------------------------------
## Clustering on the Sphere
## 

# Test 1: Test for summary_stat
test_that("summary_stat test", {
   
   dat<-rbind(matrix(rnorm(50),ncol=2), matrix(rnorm(50,4),ncol=2), matrix(rnorm(50,2),ncol=2))
   y <- rep(c(1,2,3),each=25)
   pkbd_res<- pkbc(dat, c(2,3))
   
   expect_error(summary_stat(pkbd_res, 4), "The provided pkbc object does not contain results for the requested 
           number of clusters")
   
   
})
