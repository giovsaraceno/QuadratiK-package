## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
#  #' srr tags
#  #'
#  #'
#  #' @srrstats {G1.5} two-sample test example in the associated paper

## ----eval=FALSE---------------------------------------------------------------
#  help(kb.test)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 

## ----message=FALSE------------------------------------------------------------
library(sn)
library(mvtnorm)
library(QuadratiK)
n <- 100
d <- 4
skewness_y <- 0.5
set.seed(2468)
x_2 <- rmvnorm(n, mean = rep(0, d))
y_2 <- rmsn(n = n, xi = 0, Omega = diag(d), alpha = rep(skewness_y, d))

## -----------------------------------------------------------------------------
set.seed(2468)
two_test <- kb.test(x = x_2, y = y_2)
two_test

## ----fig.width=6, fig.height=4------------------------------------------------
two_test@h$h_sel
two_test@h$power.plot

## ----eval=FALSE---------------------------------------------------------------
#  help(select_h)

## ----fig.width=6, fig.height=8------------------------------------------------
summary_two <- summary(two_test)

## -----------------------------------------------------------------------------
summary_two$summary_tables

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(2468)
#  two_test_h <- select_h(x = x_2, y = y_2, alternative = "skewness")

## -----------------------------------------------------------------------------
x_pool <- rbind(x_2, y_2)
y_memb <- rep(c(1, 2), each = n)
h <- two_test@h$h_sel
set.seed(2468)
kb.test(x = x_pool, y = y_memb, h = h)

