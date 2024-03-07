library(testthat)
test_that("Random Generation from PKBD works", {
   
   #------------------------------------------------------
   ## Generating data point on the Sphere
   ## "rejvmf"
   size <- 100
   rho=0.8
   mu = c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejvmf")
   expect_equal(dim(pkbd_dat$x),c(size,length(mu)))
   expect_true(is.matrix(pkbd_dat$x))
   expect_equal(rowSums(pkbd_dat$x^2), rep(1,size))
   
   ## "rejacg"
   size <- 100
   rho=0.8
   mu = c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejacg")
   expect_equal(dim(pkbd_dat$x),c(size,length(mu)))
   expect_true(is.matrix(pkbd_dat$x))
   expect_equal(rowSums(pkbd_dat$x^2), rep(1,size))
   
   ## "rejpsaw"
   size <- 100
   rho=0.8
   mu = c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejpsaw")
   expect_equal(dim(pkbd_dat$x),c(size,length(mu)))
   expect_true(is.matrix(pkbd_dat$x))
   expect_equal(rowSums(pkbd_dat$x^2), rep(1,size))
   
})
