library(testthat)
#------------------------------------------------------
## Clustering on the Sphere
## 

# Test 1: Test for summary_stat
test_that("extract_stat test", {
   
   # dimension = 2
   dat<-rbind(matrix(rnorm(50),ncol=2), 
              matrix(rnorm(50,4),ncol=2), 
              matrix(rnorm(50,2),ncol=2))
   y <- rep(c(1,2,3),each=25)
   pkbd_res<- pkbc(dat, c(2,3))
   
   expect_error(extract_stats(pkbd_res, 4), 
   "The provided pkbc object does not contain results for the requested
           number of clusters")
   res1 <- extract_stats(pkbd_res, 3)
   expect_equal(length(res1), 2)
   
   # dimension = 3
   dat<-rbind(matrix(rnorm(60),ncol=3), 
              matrix(rnorm(60,4),ncol=3), 
              matrix(rnorm(60,2),ncol=3))
   y <- rep(c(1,2,3),each=20)
   pkbd_res<- pkbc(dat, c(2,3))
   
   res1 <- extract_stats(pkbd_res, 2)
   res2 <- extract_stats(pkbd_res, 3)
   expect_equal(length(res1), 3)
   expect_equal(length(res2), 3)
   
   # dimension = 4
   dat<-rbind(matrix(rnorm(60),ncol=4), 
              matrix(rnorm(60,4),ncol=4), 
              matrix(rnorm(60,2),ncol=4))
   y <- rep(c(1,2,3),each=15)
   pkbd_res<- pkbc(dat, c(2,3))
   
   res1 <- extract_stats(pkbd_res, 2)
   res2 <- extract_stats(pkbd_res, 3)
   expect_equal(length(res1), 4)
   expect_equal(length(res2), 4)
   
})
