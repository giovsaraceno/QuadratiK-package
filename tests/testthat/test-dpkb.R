#' 
#' Test for dpkb
#' 
#' 
#' 
#' @srrstats {G5.1, G5.5} data sets are generated using simple functions with 
#'                        fixed seed
#' @srrstats {G5.2,G5.2a,G5.2b} all the error and warning messages are tested
#' @srrstats {G5.4,G5.4a} correctness tested on simple cases
#' @srrstats {G5.8, G5.8a,G5.8b,G5.8c} edge conditions
#' 
#' @noRd
library(testthat)
test_that("Density of PKBD", {
   
   set.seed(123)
   size <- 100
   rho <- 0.8
   mu <- c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejvmf")
   
   # Check the errors
   # mu
   expect_error(dpkb(pkbd_dat$x, mu = c(1), rho = 0.8), 
                'mu must have length >= 2',fixed=TRUE)
   expect_error(dpkb(pkbd_dat$x, mu = c(0,0,0), rho = 0.8), 
                'Input argument mu cannot be a vector of zeros',fixed=TRUE)
   
   expect_error(dpkb(rnorm(10), mu = c(0,0,0), rho = 0.8), 
                "x must be a matrix or a data.frame",fixed=TRUE)
   expect_error(dpkb(matrix(rnorm(10),ncol=1), mu = c(0,0,0), rho = 0.8), 
                'x must have dimension >= 2',fixed=TRUE)
   expect_error(dpkb(pkbd_dat$x, mu = c(1,0), rho = 0.8), 
                'number of rows of x must be equal to the length of mu',
                fixed=TRUE)
   
   expect_error(dpkb(pkbd_dat$x, mu = c(1,0,0), rho = 2), 
                'Input argument rho must be within [0,1)',fixed=TRUE)
   
   #------------------------------------------------------
   ## Generating data point on the Sphere
   ## and computing the densities
   den <- dpkb(pkbd_dat$x, mu, rho)
   
   expect_equal(dim(den),c(size,1))
   expect_false(any(den<0))
   
   den <- dpkb(data.frame(pkbd_dat$x), mu, rho)
   
   expect_equal(dim(den),c(size,1))
   expect_false(any(den<0))
   
   log_den <- dpkb(pkbd_dat$x, mu, rho, logdens = TRUE)
   
   expect_equal(dim(log_den),c(size,1))
   
   
})
