## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
#  #' srr tags
#  #'
#  #'
#  #' @srrstats {G1.5} k-sample test example in the associated paper

## -----------------------------------------------------------------------------
library(mvtnorm)
library(QuadratiK)
sizes <- rep(200,3)
eps <- 1
set.seed(2468)
x1 <- rmvnorm(sizes[1], mean = c(0,sqrt(3)*eps/3))
x2 <- rmvnorm(sizes[2], mean = c(-eps/2,-sqrt(3)*eps/6))
x3 <- rmvnorm(sizes[3], mean = c(eps/2,-sqrt(3)*eps/6))
x <- rbind(x1, x2, x3)
y <- as.factor(rep(c(1,2,3), times=sizes))

