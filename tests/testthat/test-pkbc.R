#' Test for pkbc
#' 
#' Clustering on the Sphere
#' 
#' @srrstats {G5.1, G5.5} data sets are generated using simple functions with 
#'                        fixed seed
#' @srrstats {G5.2,G5.2a,G5.2b} all the error and warning messages are tested
#' @srrstats {G5.4,G5.4a} correctness tested on simple cases
#' @srrstats {G5.8, G5.8a,G5.8b,G5.8c} edge conditions
#' @srrstats {UL7.1} test demonstrate appropriate inputs for clustering
#' @srrstats {UL7.4} test for prediction function for pkbc object
#' 
#' @noRd
library(testthat)
# Test 1: Verify Error on Invalid inputs
test_that("Error is thrown for invalid inputs", {
   set.seed(123)
   dat <- matrix(rnorm(100),ncol=2)
   
   #Invalid nClust
   expect_error(pkbc(dat, nClust = 0), 
   "Values in the input parameter nClust must be 
                  greater than 0", fixed=TRUE)
   expect_error(pkbc(dat), 
   "Input parameter nClust is required. Provide one specific 
               value or a set of possible values.", fixed=TRUE)
   
   #Invalid maxIter
   expect_error(pkbc(dat, nClust=2,maxIter=0), 
                "Input parameter maxIter must be greater than 0")
   
   #Invalid initMethod
   expect_error(pkbc(dat, nClust=2,initMethod="Invalid"), 
   'Unrecognized value "Invalid" in input 
          parameter initMethod.')
   
   #Invalid numInit
   expect_error(pkbc(dat, nClust=2,numInit=0), 
                "Input parameter numInit must be greater than 0")
})

# Test 2: Verify Error on Invalid data
test_that("Error is thrown for invalid data", {
   set.seed(123)
   dat <- matrix(rnorm(6),ncol=2)
   
   #Invalid nClust
   expect_error(pkbc(dat, nClust = 4), 
   'Only 3 unique observations. When initMethod = "sampleData", must have more
                           than numClust unique observations.', fixed=TRUE)
   
})

# Test 3: Test for valid input
test_that("Function works for valid input", {
   
   set.seed(123)
   dat<-rbind(matrix(rnorm(50),ncol=2), matrix(rnorm(50,4),ncol=2))
   result <- pkbc(dat, nClust = 2)
   expect_s4_class(result, "pkbc")
   expect_true(all(result@res_k$postProbs >= 0 & result@res_k$postProbs <= 1))
   expect_type(result@res_k[[2]]$LogLik, "double")
})

# Test 4: Test for stopping rule
test_that("Function respects the stopping rule", {
   
   set.seed(123)
   dat<-rbind(matrix(rnorm(50),ncol=2), matrix(rnorm(50,4),ncol=2))
   result_loglik <- pkbc(dat, nClust = 3, stoppingRule = 'loglik')
   result_max <- pkbc(dat, nClust = 3, stoppingRule = 'max')
   result_memb <- pkbc(dat, nClust = 3, stoppingRule = 'membership')
   expect_true(class(result_loglik)== "pkbc")
   expect_true(class(result_max)== "pkbc")
   expect_true(class(result_memb)== "pkbc")
   
   expect_error(pkbc(dat, nClust = 3, stoppingRule = 'prova'))
   
})

# Test 5: Test for clustering algorithm
test_that("Clustering algorithm works", {
   
   set.seed(123)
   dat<-rbind(matrix(rnorm(50),ncol=2), matrix(rnorm(50,4),ncol=2))
   pkbd_res<- pkbc(dat, 3)
   
   expect_true(class(pkbd_res)== "pkbc")
   expect_type(pkbd_res@res_k, "list")

})
