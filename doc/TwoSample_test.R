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

## -----------------------------------------------------------------------------
h <- 2
set.seed(2468)
two_test <- kb.test(x=x_2, y=y_2, h=h)

## ----fig.width=6, fig.height=8------------------------------------------------
summary_two <- summary(two_test)

## -----------------------------------------------------------------------------
summary_two$summary_tables

## -----------------------------------------------------------------------------
# two_test_h <- kb.test(x=x_2, y=y_2)

## ----eval=FALSE, fig.width=6, fig.height=4------------------------------------
#  set.seed(2468)
#  h_test2 <- select_h(x=x_2, y=y_2, alternative="skewness")

## ----eval=FALSE---------------------------------------------------------------
#  h_test2$h_sel

