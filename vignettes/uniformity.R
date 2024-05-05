## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
#  #' srr tags
#  #'
#  #'
#  #' @srrstats {G1.5} Uniformity test example in the associated paper

## -----------------------------------------------------------------------------
library(QuadratiK)
n <- 200
d <- 3
set.seed(2468)
z <- matrix(rnorm(n * d), n, d)
dat_sphere <- z/sqrt(rowSums(z^2))

