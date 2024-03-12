library(testthat)
#------------------------------------------------------
## Clustering on the Sphere
## 

# Test 1: Test for true.labels
test_that("Clustering algorithm works", {
   
   dat<-rbind(matrix(rnorm(50),ncol=2), matrix(rnorm(50,4),ncol=2))
   y <- rep(c(1,2),each=25)
   
   val1 <- validation(pkbc(dat, 2), true_label = NULL, elbow.plot = FALSE)
   expect_equal(nrow(val1$metrics),5)
   expect_null(val1$elbow)
   expect_equal(length(val1$IGP), 2)
   
   val2 <- validation(pkbc(dat, c(2,3,4)), true_label = y, elbow.plot = TRUE)
   expect_equal(nrow(val2$metrics),8)
   expect_equal(length(val2$IGP), 4)
   
})
