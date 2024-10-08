## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
#  #' srr tags
#  #'
#  #'
#  #' @srrstats {PD4.2} The results of the random sampling using the different
#  #'                   distributions are compared graphically, by plotting on
#  #'                   the sphere the generated observations.
#  #'

## ----message=FALSE------------------------------------------------------------
library(QuadratiK)

## -----------------------------------------------------------------------------
mu <- c(0,0,1)
d <- 3
n <- 1000
rho <- 0.8

## -----------------------------------------------------------------------------
n <- 1000
set.seed(2468)
# Generate observations using the rejection algorithm with von-Mises 
# distribution envelopes
dat1 <- rpkb(n = n, rho=rho, mu=mu, method="rejvmf")
# Generate observations using the rejection algorithm with angular central 
# Gaussian distribution envelopes
dat2 <- rpkb(n = n, rho=rho, mu=mu, method="rejacg")
# Generate observations using the projected Saw distribution
dat3 <- rpkb(n = n, rho=rho, mu=mu, method="rejpsaw")

## -----------------------------------------------------------------------------
summary(dat1)

## -----------------------------------------------------------------------------
x <- rbind(dat1$x, dat2$x, dat3$x)

## ----fig.width=6, fig.height=8------------------------------------------------
library(rgl)
# Legend information
classes <- c("rejvmf", "rejacg", "rejpsaw")
# Fix a color for each method
colors <- c("red", "blue", "green")
labels <- factor(rep(colors, each = 1000))
# Element needed for the Legend position in the following plot
offset <- 0.25
# Coordinates for legend placement
legend_x <- max(x[,1]) + offset  
legend_y <- max(x[,2]) + offset
legend_z <- seq(min(x[,3]), length.out = length(classes), by = offset)

open3d()
# Create the legend
for (i in seq_along(classes)) {
   text3d(legend_x, legend_y, legend_z[i], texts = classes[i], adj = c(0, 0.5))
   points3d(legend_x-0.1, legend_y, legend_z[i], col = colors[i], size = 5)
}
title3d("", line = 3, cex = 1.5, font=2, add=TRUE)
# Plot the sampled observations colored with respect to the used method
plot3d(x[,1], x[,2], x[,3], col = labels, size = 5, add=TRUE)
title3d("", line = 3, cex = 1.5, font=2, add=TRUE)
# Add a Sphere as background
rgl.spheres(0 , col = "transparent", alpha = 0.2)
# Rotate and zoom the generated 3 dimensional plot to facilitate visualization
view3d(theta = 10, phi = -25, zoom = 0.5)
# rglwidget is added in order to display the generated figure into the html 
# replication file.
rglwidget()
close3d()


