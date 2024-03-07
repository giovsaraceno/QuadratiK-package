## -----------------------------------------------------------------------------
library(QuadratiK)
n <- 200
d <- 3
set.seed(2468)
z <- matrix(rnorm(n * d), n, d)
dat_sphere <- z/sqrt(rowSums(z^2))

## -----------------------------------------------------------------------------
rho = 0.7
set.seed(2468)
res_unif <- pk.test(x=dat_sphere, rho=rho)

res_unif

## ----fig.width=6, fig.height=8------------------------------------------------
summary_unif <- summary(res_unif)

