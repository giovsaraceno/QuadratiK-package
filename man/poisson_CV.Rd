% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/critical_value.R
\name{poisson_CV}
\alias{poisson_CV}
\title{Compute the critical value for the Poisson KBQD tests for Uniformity}
\usage{
poisson_CV(d, size, rho, B, Quantile)
}
\arguments{
\item{d}{the dimension of generated samples.}

\item{size}{the number of observations to be generated.}

\item{rho}{the concentration parameter for the Poisson kernel.}

\item{B}{the number of replications.}

\item{Quantile}{the quantile of the distribution use to select the critical
value.}
}
\value{
the critical value for the specified dimension, size and level.
}
\description{
This function computes the empirical critical value for the U-statistics for
testing uniformity on the sphere based on the centered poisson kernel.
}
\details{
For each replication, a sample of d-dimensional observations from the uniform
distribution on the Sphere are generated and the Poisson kernel-based
U-statistic is computed. After B iterations, the critical value is selected
as the \code{Quantile} of the empirical distribution of the computed test
statistics.
}
\references{
Ding Yuxin, Markatou Marianthi, Saraceno Giovanni (2023). “Poisson
Kernel-Based Tests for Uniformity on the d-Dimensional Sphere.”
Statistica Sinica. doi: doi:10.5705/ss.202022.0347
}
\keyword{internal}
