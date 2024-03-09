library(testthat)

#------------------------------------------------------
## Clustering on the Sphere
## 

# Test 1: Verify Error on Invalid inputs
test_that("Error is thrown for invalid inputs", {
   dat <- matrix(rnorm(100),ncol=2)
   
   #Invalid nClust
   expect_error(pkbc(dat, nClust = 0), "Values in the input parameter nClust must be 
                        greater than 0")
   expect_error(pkbc(dat), "Input parameter nClust is required. Provide one specific 
                     value or a set of possible values.")
   
   #Invalid maxIter
   expect_error(pkbc(dat, nClust=2,maxIter=0), "Input parameter maxIter must be greater than 0")
   
   #Invalid initMethod
   expect_error(pkbc(dat, nClust=2,initMethod="Invalid"), 'Unrecognized value "Invalid" in input 
                parameter initMethod.')
   
   #Invalid numInit
   expect_error(pkbc(dat, nClust=2,numInit=0), "Input parameter numInit must be greater than 0")
})

# Test 2: Verify Error on Invalid data
test_that("Error is thrown for invalid data", {
   dat <- matrix(rnorm(6),ncol=2)
   
   #Invalid nClust
   expect_error(pkbc(dat, nClust = 4), 'Only 3 unique observations. When initMethod = "sampleData", must have more
                                 than numClust unique observations.')
   
})

# Test 3: Test for valid input
test_that("Function works for valid input", {
   
   dat<-rbind(matrix(rnorm(50),ncol=2), matrix(rnorm(50,4),ncol=2))
   result <- pkbc(dat, nClust = 2)
   expect_s4_class(result, "pkbc")
   expect_true(all(result@res_k$postProbs >= 0 & result@res_k$postProbs <= 1))
   expect_type(result@res_k[[2]]$LogLik, "double")
})

# Test 4: Test for stopping rule
test_that("Function respects the stopping rule", {
   
   dat<-rbind(matrix(rnorm(50),ncol=2), matrix(rnorm(50,4),ncol=2))
   result_loglik <- pkbc(dat, nClust = 3, stoppingRule = 'loglik')
   result_max <- pkbc(dat, nClust = 3, stoppingRule = 'max')
   expect_true(class(result_loglik)== "pkbc")
   expect_true(class(result_max)== "pkbc")
   
   expect_error(pkbc(dat, nClust = 3, stoppingRule = 'prova'))
   
})

# Test 5: Test for clustering algorithm
test_that("Clustering algorithm works", {
   
   dat<-rbind(matrix(rnorm(50),ncol=2), matrix(rnorm(50,4),ncol=2))
   pkbd_res<- pkbc(dat, 3)
   
   expect_true(class(pkbd_res)== "pkbc")
   expect_type(pkbd_res@res_k, "list")

})
