## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
# #' srr tags
# #'
# #'
# #' @srrstats {G1.5} Uniformity test example in the associated paper

## ----message=FALSE------------------------------------------------------------
library(QuadratiK)
n <- 200
d <- 3
set.seed(2468)
z <- matrix(rnorm(n * d), n, d)
dat_sphere <- z/sqrt(rowSums(z^2))

## -----------------------------------------------------------------------------
rho <- 0.7
set.seed(2468)
res_unif <- pk.test(x = dat_sphere, rho = rho)

show(res_unif)

## ----fig.width=6, fig.height=8------------------------------------------------
summary_unif <- summary(res_unif)

## ----warning=FALSE, message=FALSE---------------------------------------------
# Load necessary libraries
library(movMF)        
library(sphunif)      

## -----------------------------------------------------------------------------
set.seed(2468)
# Define the mean directions of the 4 von Mises-Fisher distributions
means <- rbind(
  c(1, 0),    
  c(0, 1),    
  c(-1, 0),   
  c(0, -1)    
)
# Define the concentration parameter (kappa)
kappa <- 5
# Generate 100 samples from a mixture of 4 von Mises-Fisher distributions
samples <- matrix(rmovMF(100, theta = kappa * means), ncol = 2)

## -----------------------------------------------------------------------------
# Run the pk.test from the QuadratiK package to test the data
pk_test_result <- pk.test(samples, rho = 0.8)

# Run the Bingham and Ajne tests from the sphunif package
other_test_result <- unif_test(samples, type = c("Bingham", "Ajne"))

pk_test_result
other_test_result


