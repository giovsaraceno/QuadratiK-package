library(testthat)
test_that("Random Generation from PKBD works", {

   
   # Check the errors

   expect_error(rpkb(10, mu = c(1,0,0), rho = 2), 'Input argument rho must be within [0,1)',fixed=TRUE)
   expect_error(rpkb(10, mu = c(0,0,0), rho = 0.8), 'Input argument mu cannot be a vector of zeros',fixed=TRUE)
   expect_error(rpkb(-10, mu = c(1,0,0), rho = 0.8), 'n must be a positive integer',fixed=TRUE)
   expect_error(rpkb(10, mu = c(1,0,0), rho = 0.8, method = "invalid"), 'Unknown method',fixed=TRUE)
   
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
   
})
