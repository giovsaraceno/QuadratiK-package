## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
#  #' srr tags
#  #'
#  #'
#  #' @srrstats {G1.5} two-sample test example in the associated paper

## ----message=FALSE------------------------------------------------------------
library(sn)
library(mvtnorm)
library(QuadratiK)
n <- 200
d <- 4
skewness_y <- 0.5
set.seed(2468)
x_2 <- rmvnorm(n, mean = rep(0,d))
y_2 <- rmsn(n=n, xi=0, Omega = diag(d), alpha=rep(skewness_y,d))

