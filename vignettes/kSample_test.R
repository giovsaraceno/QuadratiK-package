## -----------------------------------------------------------------------------
library(mvtnorm)
library(QuadratiK)
sizes <- rep(200,3)
eps = 1
set.seed(2468)
x1 <- rmvnorm(sizes[1], mean = c(0,sqrt(3)*eps/3))
x2 <- rmvnorm(sizes[2], mean = c(-eps/2,-sqrt(3)*eps/6))
x3 <- rmvnorm(sizes[3], mean = c(eps/2,-sqrt(3)*eps/6))
x <- rbind(x1, x2, x3)
y <- as.factor(rep(c(1,2,3), times=sizes))

## -----------------------------------------------------------------------------
h=1.5
set.seed(2468)
k_test <- kb.test(x=x, y=y, h=h)
k_test

## -----------------------------------------------------------------------------
summary_ktest <- summary(k_test)
summary_ktest$summary_tables

## -----------------------------------------------------------------------------
#k_test_h <- kb.test(x=x, y=y)

## ----eval=FALSE, fig.width=6, fig.height=4------------------------------------
#  set.seed(2468)
#  h_k <- select_h(x=x, y=y, alternative="skewness")

## ----eval=FALSE---------------------------------------------------------------
#  h_k$h_sel

