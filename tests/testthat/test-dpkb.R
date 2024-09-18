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
#' @srrstats {PD4.0} numeric outputs are tested
#' @srrstats {PD4.1} Tests that density values are as expected for selected 
#'                   inputs
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

test_that("Test numerical values", {
   #------------------------------------------------------
   ## Test numerical values
   rho <- 0.8
   mu <- c(1,0,0)
   d <- 3
   x <- matrix(rnorm(3),nrow=1)
   # For rho=0 the density corresponds to the normalizing constant for d 
   # independently of x
   expect_equal(as.numeric(dpkb(x, mu, rho=0)), 1/(2*pi^(d/2)*gamma(d/2)^(-1)))
   
   # For x = mu the product of x*mu =1 in the poisson density. The remaining 
   # depends only on rho and d
   expect_equal(as.numeric(dpkb(matrix(mu,nrow=1), mu, rho=0.5)), 0.47746483)
   expect_equal(as.numeric(dpkb(matrix(mu,nrow=1), mu, rho=0.5, 
                                logdens = TRUE)), log(0.47746483))
   
})

test_that("tests - Density of PKBD", {
   
   set.seed(123)
   size <- 100
   rho <- 0.8
   mu <- c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejvmf")$x
   
   # Check the errors
   # x
   pkbd_dat[1,] <- NA
   expect_error(dpkb(pkbd_dat, mu = mu, rho = rho), 
                'There are missing values in x!',fixed=TRUE)
   pkbd_dat[1,] <- Inf
   expect_error(dpkb(pkbd_dat, mu = mu, rho = rho), 
         'There are undefined values in x, that is Nan, Inf, -Inf',fixed=TRUE)
   
   expect_error(dpkb(rnorm(10), mu = c(0,0,0), rho = 0.8), 
                "x must be a matrix or a data.frame",fixed=TRUE)
   expect_error(dpkb(matrix(rnorm(10),ncol=1), mu = c(0,0,0), rho = 0.8), 
                'x must have dimension >= 2',fixed=TRUE)
   
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejvmf")
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

