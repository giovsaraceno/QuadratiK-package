library(testthat)

## "Tests for my custom clustering function"

#------------------------------------------------------
## Clustering on the Sphere
## 
size <- 100
groups <- c(rep(1, size), rep(2, size),rep(3,size))
rho <- 0.8
data1<-rpkb(size, c(1,0,0),rho,method='rejvmf')
data2<-rpkb(size, c(0,1,0),rho,method='rejacg')
data3<-rpkb(size, c(-1,0,0),rho,method='rejpsaw')
dat<-rbind(data1$x,data2$x, data3$x)

# Test 1: Verify Error on Invalid nClust
test_that("Error is thrown for invalid nClust", {
   dat <- matrix(rnorm(100),ncol=2)
   expect_error(pkbc(dat, nClust = 0), "Values in the input parameter nClust must be 
                        greater than 0")
})

# Test 2: Test for valid input
test_that("Function works for valid input", {
   
   result <- pkbc(dat, nClust = 3)
   expect_s4_class(result, "pkbc")
   expect_true(all(result@res_k$postProbs >= 0 & result@res_k$postProbs <= 1))
   expect_type(result@res_k[[3]]$LogLik, "double")
})

# Test 3: Test for stopping rule
test_that("Function respects the stopping rule", {
   
   result_loglik <- pkbc(dat, nClust = 3, stoppingRule = 'loglik')
   result_max <- pkbc(dat, nClust = 3, stoppingRule = 'max')
   expect_true(class(result_loglik)== "pkbc")
   expect_true(class(result_max)== "pkbc")
   
   expect_error(pkbc(dat, nClust = 3, stoppingRule = 'prova'))
   
})


test_that("Clustering algorithm works", {
   
   
   pkbd_res<- pkbc(dat, 3)
   
   expect_true(class(pkbd_res)== "pkbc")
   expect_type(pkbd_res@res_k, "list")

})
