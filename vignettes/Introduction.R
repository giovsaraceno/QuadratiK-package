## ----eval = FALSE-------------------------------------------------------------
#  install.packages("QuadratiK")

## ----eval = FALSE-------------------------------------------------------------
#  library(devtools)
#  devtools::install_github('giovsaraceno/QuadratiK-package')

## -----------------------------------------------------------------------------
library(QuadratiK)

## ----eval=FALSE---------------------------------------------------------------
#  ?kb.test
#  ?select_h

## -----------------------------------------------------------------------------
x <- matrix(rnorm(100), ncol = 2)
# Does x come from a multivariate standard normal distribution?
kb.test(x, h=0.4)

## -----------------------------------------------------------------------------
x <- matrix(rnorm(100,4), ncol = 2)
# Does x come from a multivariate standard normal distribution?
kb.test(x, mu_hat = c(4,4), Sigma_hat = diag(2), h = 0.4)

## -----------------------------------------------------------------------------
x <- matrix(rnorm(100), ncol = 2)
y <- matrix(rnorm(100,mean = 5), ncol = 2)
# Do x and y come from the same distribution?
kb.test(x, y, h = 0.4)

## -----------------------------------------------------------------------------
x1 <- matrix(rnorm(100), ncol = 2)
x2 <- matrix(rnorm(100), ncol = 2)
x3 <- matrix(rnorm(100, mean = 5), ncol = 2)
y <- rep(c(1, 2, 3), each = 50)
# Do x1, x2 and x3 come from the same distribution?
x <- rbind(x1, x2, x3)
kb.test(x, y, h = 0.4)

## ----eval=FALSE---------------------------------------------------------------
#  ?pk.test

## -----------------------------------------------------------------------------
# Generate points on the sphere from the uniform ditribution 
x <- sample_hypersphere(d = 3, n_points = 100)
# Does x come from the uniform distribution on the sphere?
pk.test(x, rho = 0.9)

## ----eval=FALSE---------------------------------------------------------------
#  ?dpkb
#  ?rpkb

## -----------------------------------------------------------------------------
mu <- c(1,0,0)
rho <- 0.9
x <- rpkb(n = 100, mu = mu, rho = rho)
head(x$x)
dens_x <- dpkb(x$x, mu = mu, rho = rho)
head(dens_x)

## ----eval=FALSE---------------------------------------------------------------
#  ?pkbc

## -----------------------------------------------------------------------------
# Generate 3 samples from the PKBD with different location directions
x1 <- rpkb(n = 100, mu = c(1,0,0), rho = rho)
x2 <- rpkb(n = 100, mu = c(-1,0,0), rho = rho)
x3 <- rpkb(n = 100, mu = c(0,0,1), rho = rho)
x <- rbind(x1$x, x2$x, x3$x)
# Perform the clustering algorithm
# Serch for 2, 3 or 4 clusters
cluster_res <- pkbc(dat = x, nClust = c(2, 3, 4))
summary(cluster_res)

## -----------------------------------------------------------------------------
# Predict the membership of new data with respect to the clustering results
x_new <- rpkb(n = 10, mu = c(1,0,0), rho = rho)
memb_mew <- predict(cluster_res, k = 3, newdata = x_new$x)
memb_mew$Memb


## -----------------------------------------------------------------------------
# Compute measures for evaluating the clustering results
val_res <- pkbc_validation(cluster_res)
val_res

## ----fig.show = "hold", fig.width=6, fig.height=6-----------------------------
# Plot method for the pkbc object:
# - scatter plot of data points on the sphere
# - elbow plot for helping the choice of the number of clusters
plot(cluster_res)

