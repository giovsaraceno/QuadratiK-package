library(testthat)

context("Tests for my custom clustering function")

#------------------------------------------------------
## Clustering on the Sphere
## 
size <- 100
groups<-c(rep(1, size), rep(2, size),rep(3,size))
rho=0.8
data1<-rpkb(size, c(1,0,0),rho,method='rejvmf')
data2<-rpkb(size, c(0,1,0),rho,method='rejacg')
data3<-rpkb(size, c(-1,0,0),rho,method='rejpsaw')
dat<-rbind(data1$x,data2$x, data3$x)

# Test 1: Verify Error on Invalid nClust
test_that("Error is thrown for invalid nClust", {
   expect_error(pkbc(dat, nClust = 0), "nClust must be greater than 0")
})

# Test 2: Test for valid input
test_that("Function works for valid input", {
   
   # Assuming your function is named 'my_clustering_function'
   # and it returns a list with a specific structure.
   result <- pkbc(dat, nClust = 3)
   expect_is(result, "list")
   expect_true(all(result$postProbs >= 0 & result$postProbs <= 1))
   expect_length(result$LogLik, 1)
})

# Test 3: Test for stopping rule
test_that("Function respects the stopping rule", {
   
   result_loglik <- pkbc(dat, nClust = 3, stoppingRule = 'loglik')
   result_max <- pkbc(dat, nClust = 3, stoppingRule = 'max')
   expect_not_equal(result_loglik, result_max)
   
})


test_that("Clustering algorithm works", {
   
   
   pkbd_res<- pkbc(dat, 3)
   
   expect_true(class(pkbd_res)== "pkbc")
   expect_type(pkbd_res@res_k, "list")

})
