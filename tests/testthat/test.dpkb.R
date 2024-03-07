library(testthat)
test_that("Random Generation from PKBD works", {
   
   #------------------------------------------------------
   ## Generating data point on the Sphere
   ## and computing the densities
   size <- 100
   rho=0.8
   mu = c(1,0,0)
   pkbd_dat <- rpkb(size, mu = mu, rho = rho, method = "rejvmf")
   den <- dpkb(pkbd_dat$x, mu, rho)
   
   expect_equal(dim(den),c(size,1))
   expect_false(any(den<0))
   
   
})
