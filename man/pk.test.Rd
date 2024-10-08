% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pk.test.R
\name{pk.test}
\alias{pk.test}
\alias{pk.test,ANY-method}
\alias{show,pk.test-method}
\title{Poisson kernel-based quadratic distance test of Uniformity on the sphere}
\usage{
pk.test(x, rho, B = 300, Quantile = 0.95)

\S4method{pk.test}{ANY}(x, rho, B = 300, Quantile = 0.95)

\S4method{show}{pk.test}(object)
}
\arguments{
\item{x}{A numeric d-dim matrix of data points on the Sphere S^(d-1).}

\item{rho}{Concentration parameter of the Poisson kernel function.}

\item{B}{Number of Monte Carlo iterations for critical value estimation of Un
(default: 300).}

\item{Quantile}{The quantile to use for critical value estimation,
0.95 is the default value.}

\item{object}{Object of class \code{pk.test}}
}
\value{
An S4 object of class \code{pk.test} containing the results of the
Poisson kernel-based tests. The object contains the following slots:
\itemize{
\item \code{method}: Description of the test performed.
\item \code{x} Data matrix.
\item \code{Un} The value of the U-statistic.
\item \code{CV_Un} The empirical critical value for Un.
\item \code{H0_Vn} A logical value indicating whether or not the null
hypothesis is rejected according to Un.
\item \code{Vn} The value of the V-statistic Vn.
\item \code{CV_Vn} The critical value for Vn computed following the
asymptotic distribution.
\item \code{H0_Vn} A logical value indicating whether or not the null
hypothesis is rejected according to Vn.
\item \code{rho} The value of concentration parameter used for the Poisson
kernel function.
\item \code{B} Number of replications for the critical value of the
U-statistic Un.
}
}
\description{
This function performs the kernel-based quadratic distance goodness-of-fit
tests for Uniformity for spherical data \code{x} using the Poisson kernel
with concentration parameter \code{rho}. \cr
The Poisson kernel-based test for uniformity exhibits excellent results
especially in the case of multimodal distributions, as shown in the example
of the \href{../doc/uniformity.html}{Uniformity test on the Sphere vignette}.
}
\details{
Let \eqn{x_1, x_2, ..., x_n} be a random sample with empirical distribution
function \eqn{\hat F}.
We test the null hypothesis of uniformity on the
\eqn{d}-dimensional sphere, i.e. \eqn{H_0:F=G}, where \eqn{G} is the uniform
distribution on the \eqn{d}-dimensional sphere \eqn{\mathcal{S}^{d-1}}.
We compute the U-statistic estimate of the sample KBQD (Kernel-Based
Quadratic Distance)
\deqn{U_{n}=\frac{1}{n(n-1)}\sum_{i=2}^{n}\sum_{j=1}^{i-1}K_{cen}
(\mathbf{x}_{i}, \mathbf{x}_{j}),}
then the first test statistic is given as
\deqn{T_{n}=\frac{U_{n}}{\sqrt{Var(U_{n})}},}
with
\deqn{Var(U_{n})= \frac{2}{n(n-1)}
\left[\frac{1+\rho^{2}}{(1-\rho^{2})^{d-1}}-1\right],}
and the V-statistic estimate of the KBQD
\deqn{V_{n} = \frac{1}{n}\sum_{i=1}^{n}\sum_{j=1}^{n}K_{cen}
(\mathbf{x}_{i}, \mathbf{x}_{j}),}
where \eqn{K_{cen}} denotes the Poisson kernel \eqn{K_\rho} centered with
respect to the uniform distribution on the \eqn{d}-dimensional sphere, that
is
\deqn{K_{cen}(\mathbf{u}, \mathbf{v}) = K_\rho(\mathbf{u}, \mathbf{v}) -1}
and
\deqn{K_\rho(\mathbf{u}, \mathbf{v}) = \frac{1-\rho^{2}}{\left(1+\rho^{2}-
2\rho (\mathbf{u}\cdot \mathbf{v})\right)^{d/2}},}
for every \eqn{\mathbf{u}, \mathbf{v} \in \mathcal{S}^{d-1} 
\times \mathcal{S}^{d-1}}.

The asymptotic distribution of the V-statistic is an infinite combination
of weighted independent chi-squared random variables with one degree of
freedom. The cutoff value is obtained using the Satterthwaite approximation
\eqn{c \cdot \chi_{DOF}^2}, where \deqn{c=\frac{(1+\rho^{2})-
(1-\rho^{2})^{d-1}}{(1+\rho)^{d}-(1-\rho^{2})^{d-1}}} and \deqn{DOF(K_{cen}
)=\left(\frac{1+\rho}{1-\rho} \right)^{d-1}\left\{ 
\frac{\left(1+\rho-(1-\rho)^{d-1} \right )^{2}}
{1+\rho^{2}-(1-\rho^{2})^{d-1}}\right \}.}.
For the \eqn{U}-statistic the cutoff is determined empirically:
\itemize{
\item Generate data from a Uniform distribution on the d-dimensional sphere;
\item Compute the test statistics for \code{B} Monte Carlo(MC) replications;
\item Compute the 95th quantile of the empirical distribution of the test
statistic.
}
}
\note{
A U-statistic is a type of statistic that is used to estimate a population
parameter. It is based on the idea of averaging over all possible \emph{distinct}
combinations of a fixed size from a sample.
A V-statistic considers all possible tuples of a certain size, not just
distinct combinations and can be used in contexts where unbiasedness is not
required.
}
\examples{
# create a pk.test object
x_sp <- sample_hypersphere(3, n_points=100)
unif_test <- pk.test(x_sp,rho=0.8)
unif_test

}
\references{
Ding, Y., Markatou, M. and Saraceno, G. (2023). “Poisson Kernel-Based Tests
for Uniformity on the d-Dimensional Sphere.” Statistica Sinica.
doi:10.5705/ss.202022.0347
}
\seealso{
\linkS4class{pk.test}
}
