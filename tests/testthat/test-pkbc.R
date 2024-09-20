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
   expect_error(pkbc(dat, nClust="invalid"), 
                "nClust must be a signle value or a numeric vector of possible
                  values", fixed=TRUE)
   
   #Invalid maxIter
   expect_error(pkbc(dat, nClust=2,maxIter=0), 
                "Input parameter maxIter must be greater than 0")
   
   #Invalid initMethod
   expect_error(pkbc(dat, nClust=2,initMethod="Invalid"), 
   'Unrecognized value Invalid in input 
          parameter initMethod.')
   
   #Invalid numInit
   expect_error(pkbc(dat, nClust=2,numInit=0), 
                "Input parameter numInit must be greater than 0")
})

# Test 2: Verify Error on Invalid data
test_that("Error is thrown for invalid data", {
   set.seed(123)
   dat <- matrix(rnorm(6),ncol=2)
   
   #Not enough observations
   expect_error(pkbc(dat, nClust = 4), 
   'Only 3 unique observations. When initMethod = "sampleData", must have more
                           than nClust unique observations.', fixed=TRUE)
   
   # data is a vector
   expect_error(pkbc(rep(c(1,2),each=50), nClust = 2), 
                'dat must be a matrix or a data.frame', fixed=TRUE)
   # data is a character matrix 
   expect_error(pkbc(matrix(rep(c("A","B"),each=100),ncol=2), nClust = 2), 
                'dat must be a numeric matrix or data.frame', fixed=TRUE)
   
   # NA in the data
   dat <- matrix(rnorm(100),ncol=2)
   dat[1,] <- NA 
   expect_error(pkbc(dat = dat, nClust = 2), 
                'There are missing values in the data set!', fixed=TRUE)
   
   # Inf or Nan in the data
   dat <- matrix(rnorm(100),ncol=2)
   dat[1,] <- Inf 
   expect_error(pkbc(dat = dat, nClust = 2), 
            'There are undefined values, that is Nan, Inf, -Inf', fixed=TRUE)
   dat <- matrix(rnorm(100),ncol=2)
   dat[1,] <- NaN 
   expect_error(pkbc(dat = dat, nClust = 2), 
                'There are missing values in the data set!', fixed=TRUE)
   
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
   dat <- as.data.frame(dat)
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
   
   # Generated well separated clusters and show that the clustering algorithm
   # works properly in this simple case
   set.seed(123)
   x <- rpkb(50, c(1,0,0),0.99, method = "rejacg")$x
   y <- rpkb(50, c(0,0,1),0.99, method = "rejacg")$x
   z <- rpkb(50, c(-1,0,0),0.99, method = "rejacg")$x
   dat<-rbind(x,y,z)
   label <- rep(c(3,1,2),each=50)
   pkbd_res<- pkbc(dat, 3)
   
   expect_true(class(pkbd_res)== "pkbc")
   expect_type(pkbd_res@res_k, "list")
   require(mclust)
 expect_equal(mclust::adjustedRandIndex(pkbd_res@res_k[[3]]$finalMemb,label),1)
   
   # Tests for stats_clusters, correct input
   expect_error(stats_clusters(pkbd_res, "invalid"), 
                'k must be an integer', fixed=TRUE)
   # Tests for stats_clusters
   stats_res <- stats_clusters(pkbd_res, 3)
   expect_equal(stats_res[[1]][1,4], mean(dat[,1]))
   expect_equal(stats_res[[2]][1,4], mean(dat[,2]))
   expect_equal(stats_res[[3]][1,4], mean(dat[,3]))
   
   # Tests for pkbc_validation: valid input
   expect_error(pkbc_validation(pkbd_res,label[-1]),
                "true_label must have the same length of finalMemb.")
   expect_error(pkbc_validation(pkbd_res,matrix(label)),
                "true_label must be a factor or a vector")
   
   # Tests for pkbc_validation
   cluster_valid <- pkbc_validation(pkbd_res,label)
   expect_equal(cluster_valid$metrics[2,1],1)
   expect_gt(cluster_valid$metrics[1,1],0.9)
   expect_equal(cluster_valid$IGP[[3]],c(1,1,1))
   
   # Tests for plot method
   pkbd_res<- pkbc(dat, c(2,3,4))
   expect_silent(plot(pkbd_res,3))
   expect_silent(plot(pkbd_res,k=3,true_label=label))
})

# Test 6: Test for plot method
test_that("plot method for the clustering algorithm", {
   
   # d=2
   set.seed(123)
   x <- rpkb(50, c(1,0),0.95, method = "rejacg")$x
   y <- rpkb(50, c(0,1),0.95, method = "rejacg")$x
   dat<-rbind(x,y)
   label <- rep(c(1,2),each=50)
   pkbd_res<- pkbc(dat, c(2,3))
   # Tests for plot method
   expect_silent(plot(pkbd_res,k = 2))
   
   # d>3
   set.seed(123)
   x <- rpkb(50, c(1,0,0,0),0.95, method = "rejacg")$x
   y <- rpkb(50, c(0,0,1,0),0.95, method = "rejacg")$x
   dat<-rbind(x,y)
   label <- rep(c(1,2),each=50)
   pkbd_res<- pkbc(dat, c(2,3))
   # Tests for plot method
   expect_silent(plot(pkbd_res,k = 2, pca_res = TRUE))
   expect_silent(plot(pkbd_res,k = 2,true_label=label))
})

# Test 7: Test for predict method
test_that("plot method for the clustering algorithm", {
   
   set.seed(123)
   x <- rpkb(50, c(1,0),0.99, method = "rejacg")$x
   y <- rpkb(50, c(-1,0),0.99, method = "rejacg")$x
   dat<-rbind(x,y)
   label <- rep(c(1,2),each=50)
   pkbd_res<- pkbc(dat, c(2,3))
   
   # Correct newdata input
   newdat <- "invalid"
   expect_error(predict(pkbd_res, k = 2, newdat), 
                'newdata must be a matrix or data.frame.', fixed=TRUE)
   
   newdat <- matrix("invalid", ncol=2, nrow=50)
   expect_error(predict(pkbd_res, k = 2, newdat), 
                'newdata must be numeric', fixed=TRUE)
   
   newdat <- matrix(rnorm(150),ncol=3)
   expect_error(predict(pkbd_res, k = 2, newdat), 
         'newdata must have the same number of variables as the training data.', 
                fixed=TRUE)
   
   
   # Tests for predict method
   expect_equal(predict(pkbd_res,k=2),pkbd_res@res_k[[2]]$finalMemb)
   
   newdat <- as.data.frame(rbind(rpkb(50, c(1,0),0.99, method = "rejacg")$x,
                   rpkb(50, c(-1,0),0.99, method = "rejacg")$x))
   expect_equal(predict(pkbd_res, k=2, newdat)$Memb,rep(c(2,1),each=50))
   

})

# Test 8: Test for show and summary methods
test_that("show and summary methods for the clustering algorithm", {
   
   set.seed(123)
   x <- rpkb(50, c(1,0),0.99, method = "rejacg")$x
   y <- rpkb(50, c(-1,0),0.99, method = "rejacg")$x
   dat<-rbind(x,y)
   label <- rep(c(1,2),each=50)
   pkbd_res<- pkbc(dat, c(2,3))
   output <- capture.output(show(pkbd_res))
   
   # Verify that the output contains expected components
   expect_true(any(grepl("Poisson Kernel-Based Clustering on the Sphere", 
                         output)))
   expect_true(any(grepl("Available components", output)))
   expect_true(any(grepl("Input Parameters", output)))
   expect_true(any(grepl("Considered possible number of clusters:  2 3", 
                         output)))
   expect_true(any({
   grepl("Available components for each value of number of clusters:", output)
      }))
   
   out_sum <- capture.output(summary(pkbd_res))
   
   # Verify that the output contains the expected components
   expect_true(any(grepl("Poisson Kernel-Based Clustering on the Sphere", 
                         out_sum)))
   expect_true(any(grepl("Summary:", out_sum)))
   expect_true(any(grepl("LogLik", out_sum)))
   expect_true(any(grepl("WCSS", out_sum)))
   expect_true(any(grepl("Results for 2 clusters:", out_sum)))
   expect_true(any(grepl("Estimated Mixing Proportions", out_sum)))
   expect_true(any(grepl("Clustering table", out_sum)))
   expect_true(any(grepl("Results for 3 clusters:", out_sum)))
   
   
})
