#' Test for rpkb
#' 
#' 
#' 
#' @srrstats {G5.1, G5.5} data sets are generated using simple functions with 
#'                        fixed seed
#' @srrstats {G5.2,G5.2a,G5.2b} all the error and warning messages are tested
#' @srrstats {PD4.0} the numeric output are tested
#' @srrstats {PD4.2} The function is tested for all the methods (distributions)
#'                   available
#' @srrstats {PD4.3} different values of the \code{uniroot} function are tested
#' @srrstats {PD4.4} The function uniroot satisfies this standard.
#' @srrstats {PD4.1} Tests that generated points are in expected range
#' 
#' 
#' @noRd
library(testthat)
test_that("Random Generation from PKBD works", {

   
   # Check the errors
   set.seed(123)
   expect_error(rpkb(10, mu = c(1,0,0), rho = 2), 
                'Input argument rho must be within [0,1)',fixed=TRUE)
   expect_error(rpkb(10, mu = c(0,0,0), rho = 0.8), 
                'Input argument mu cannot be a vector of zeros',fixed=TRUE)
   expect_error(rpkb(10, mu = 2, rho = 0.8), 
                'mu must have length >= 2',fixed=TRUE)
   
   expect_error(rpkb(-10, mu = c(1,0,0), rho = 0.8), 
                'n must be a positive integer',fixed=TRUE)
   expect_error(rpkb(10, mu = c(1,0,0), rho = 0.8, method = "invalid"), 
                'Unknown method',fixed=TRUE)
   
   #------------------------------------------------------
   ## Generating data point on the Sphere
   ## "rejvmf"
   size <- 100
   rho <- 0.8
   mu <- c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejvmf")
   expect_equal(dim(pkbd_dat$x),c(size,length(mu)))
   expect_true(is.matrix(pkbd_dat$x))
   expect_equal(rowSums(pkbd_dat$x^2), rep(1,size))
   
   ## "rejacg"
   size <- 100
   rho <- 0.8
   mu <- c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejacg")
   expect_equal(dim(pkbd_dat$x),c(size,length(mu)))
   expect_true(is.matrix(pkbd_dat$x))
   expect_equal(rowSums(pkbd_dat$x^2), rep(1,size))
   
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejacg", max.iter=20)
   expect_lt(pkbd_dat$beta$iter, 20)
   
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejacg", tol.eps=1e-7)
   expect_lt(pkbd_dat$beta$estim.prec, 1e-7)
   
   ## "rejpsaw"
   size <- 100
   rho <- 0.8
   mu <- c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejpsaw")
   expect_equal(dim(pkbd_dat$x),c(size,length(mu)))
   expect_true(is.matrix(pkbd_dat$x))
   expect_equal(rowSums(pkbd_dat$x^2), rep(1,size))
   
   # dimension = 2
   mu <- c(1,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejpsaw")
   expect_equal(dim(pkbd_dat$x),c(size,length(mu)))
   expect_true(is.matrix(pkbd_dat$x))
   expect_equal(rowSums(pkbd_dat$x^2), rep(1,size))
   
   # dimension > 3 
   mu <- c(1,0,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejpsaw")
   expect_equal(dim(pkbd_dat$x),c(size,length(mu)))
   expect_true(is.matrix(pkbd_dat$x))
   expect_equal(rowSums(pkbd_dat$x^2), rep(1,size))
   
   #------------------------------------------------------
   ## Test numerical values
   rho <- 0.9
   mu <- rnorm(3)
   mu <- mu/sqrt(sum(mu^2))
   size <- 100
   
   x <- rpkb(size, mu = mu, rho = rho, method = "rejpsaw")
   # Test if the mean of generated points is close to the true mean
   expect_lt(max(abs(colMeans(x$x)-mu)),1e-1)
   
   x <- rpkb(size, mu = mu, rho = 0.99, method = "rejpsaw")
   # Test if the mean of generated points is close to the true mean
   # higher rho, less variability
   expect_lt(max(abs(colMeans(x$x)-mu)),1e-2)
   
   
})

test_that("Random Generation from PKBD compared to wrapped Cauchy", {
   
   require(circular)
   
   # Parameters 
   n <- 1000  
   location <- circular(0)  
   rho <- 0.5  
   # Generate data from wrapped cauchy
   wc <- circular::rvonmises(n, mu = location, kappa = -log(rho))
   # Generate data from pkbd 
   pkbd <- rpkb(n, mu = c(1, 0), rho = 0.5)$x
   
   # Convert Cartesian coordinates to angles for comparison
   pkbd_angles <- atan2(pkbd[,2], pkbd[,1]) %% (2*pi)
   
   # Perform Kolmogorov-smirnov test 
   # (common for comparing two circular distributions)
   ks_test <- ks.test(pkbd_angles, wc)
   expect_lt(ks_test$p.value,0.05)

   
})
