% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{pk.test-class}
\alias{pk.test-class}
\title{An S4 class for Poisson kernel-based quadratic distance tests.}
\description{
A class to represent the results of Poisson kernel-based
quadratic distance tests for Uniformity on the sphere.
}
\section{Slots}{

\describe{
\item{\code{method}}{Description of the test.}

\item{\code{x}}{Matrix of data}

\item{\code{Un}}{The value of the U-statistic.}

\item{\code{CV_Un}}{The critical value for Un computed through replications.}

\item{\code{H0_Un}}{A logical value indicating whether or not the null hypothesis is
rejected according to Un.}

\item{\code{Vn}}{The value of the V-statistic.}

\item{\code{CV_Vn}}{The critical value for Vn computed following the asymptotic
distribution.}

\item{\code{H0_Vn}}{A logical value indicating whether or not the null hypothesis is
rejected according to Vn.}

\item{\code{rho}}{The concentration parameter of the Poisson kernel.}

\item{\code{B}}{Number of replications.}

\item{\code{var_Un}}{exact variance of the kernel-based U-statistic.}
}}

\examples{
# create a pk.test object
d=3
size=100
x_sp <- sample_hypersphere(d, n_points=size)
pk.test(x_sp,rho=0.8)

}
\seealso{
\code{\link[=pk.test]{pk.test()}}
}
