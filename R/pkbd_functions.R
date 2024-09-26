#'
#' The Poisson kernel-based Distribution (PKBD)
#' 
#' @description 
#' The Poisson kernel-based densities are based on the normalized Poisson kernel
#' and are defined on the \eqn{d}-dimensional unit sphere. Given a vector 
#' \eqn{\mathbf{\mu} \in \mathcal{S}^{d-1}}, where \eqn{\mathcal{S}^{d-1}= 
#' \{x \in \mathbb{R}^d : ||x|| = 1\}}, and a parameter \eqn{\rho} such that 
#' \eqn{0 < \rho < 1}, the probability density function of a \eqn{d}-variate 
#' Poisson kernel-based density is defined by:
#' \deqn{f(\mathbf{x}|\rho, \mathbf{\mu}) = \frac{1-\rho^2}{\omega_d 
#' ||\mathbf{x} - \rho \mathbf{\mu}||^d},}
#' where \eqn{\mu} is a vector orienting the center of the distribution, 
#' \eqn{\rho} is a parameter to control the concentration of the distribution 
#' around the vector \eqn{\mu} and it is related to the variance of the 
#' distribution. Recall that, for \eqn{x = (x_1, \ldots, x_d) \in \mathbb{R}^d},
#' \eqn{||x|| = \sqrt{x_1^2 + \ldots + x_d^2}}. Furthermore, \eqn{\omega_d =
#' 2\pi^{d/2} [\Gamma(d/2)]^{-1}} is the surface area of the unit sphere in
#' \eqn{\mathbb{R}^d} (see Golzy and Markatou, 2020). When \eqn{\rho \to 0}, 
#' the Poisson kernel-based density tends to the uniform density on the sphere.
#' Connections of the PKBDs to other distributions are discussed in detail in 
#' Golzy and Markatou (2020). Here we note that when \eqn{d=2}, PKBDs reduce to 
#' the wrapped Cauchy distribution. Additionally, with precise choice of the 
#' parameters \eqn{\rho} and \eqn{\mu} the two-dimensional PKBD becomes a 
#' two-dimensional projected normal distribution. However, the connection with 
#' the \eqn{d}-dimensional projected normal distributions does not carry beyond 
#' \eqn{d=2}.  
#' Golzy and Markatou (2020) proposed an acceptance-rejection method for 
#' simulating data from a PKBD using von Mises-Fisher envelops (\code{rejvmf} 
#' method). Furthermore Sablica, Hornik and Leydold (2023) proposed new ways for
#' simulating from the PKBD, using angular central Gaussian envelops 
#' (\code{rejacg}) or using the projected Saw distributions (\code{rejpsaw}).
#'  
#' @param x Matrix (or data.frame) of data point on the sphere 
#'          \eqn{\mathcal{S}^{d-1}}, with \eqn{d \ge 2}.
#' @param mu Location parameter with same length as the rows of \code{x}.
#' @param rho Concentration parameter, with \eqn{0 \le} \code{rho} \eqn{< 1}.
#' @param logdens Logical; if 'TRUE', densities are returned in logarithmic
#'                scale.
#'
#' @return 
#' \code{dpkb} gives the density value;
#' \code{rpkb} generates random observations from the PKBD. 
#' 
#' @details
#' This function \code{dpkb()} computes the density value for a given point 
#' \code{x} from the Poisson kernel-based distribution with mean direction
#' vector \code{mu} and concentration parameter \code{rho}.
#'
#' @examples
#' # Generate some data from pkbd density
#' pkbd_dat <- rpkb(10, c(0.5,0), 0.5)
#' 
#' # Calculate the PKBD density values
#' dens_val <- dpkb(pkbd_dat$x, c(0.5,0.5),0.5)
#'
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G2.0,G2.0a} check input mu 
#' @srrstats {G2.7,G2.8} different input x
#' @srrstats {G2.13,G2.14,G2.14a,G2.15,G2.16} error for NA, Nan, Inf, -Inf
#' 
#' @export
dpkb <- function(x, mu, rho, logdens = FALSE) {
   # validate input
   if (length(mu) < 2) {
      stop('mu must have length >= 2')
   }
   # If x is data.frame, convert into a matrix
   if(is.data.frame(x)) {
      x <- as.matrix(x)
   } else if(!is.matrix(x)){
      stop("x must be a matrix or a data.frame")
   }
   if(any(is.na(x))){
      stop("There are missing values in x!")
   } else if(any(is.infinite(x) |is.nan(x))){
      stop("There are undefined values in x, that is Nan, Inf, -Inf")
   }
   p <- ncol(x)
   if (p < 2) {
      stop('x must have dimension >= 2')
   }
   if (length(mu) != p) {
      stop('number of rows of x must be equal to the length of mu')
   }
   if (rho >= 1 | rho < 0) {
      stop('Input argument rho must be within [0,1)')
   }
   # norm input
   if (sum(abs(mu)) == 0) {
      stop('Input argument mu cannot be a vector of zeros')
   }
   mu <- mu / sqrt(sum(mu^2))
   x <- x / sqrt(rowSums(x^2))
   
   # calculate log(density)
   logretval <- (log(1 - rho^2) - log(2) - p/2*log(pi) + lgamma(p/2) -
                    p/2*log(1 + rho^2 - 2*rho * (x %*% mu)))
   if (logdens) {
      return(logretval)
   } else {
      return(exp(logretval))
   }
}
#'
#' @rdname dpkb
#'
#' @param n number of observations.
#' @param mu location vector parameter with length indicating the dimension of 
#'          generated points.
#' @param rho Concentration parameter, with \eqn{0 \le} \code{rho} \eqn{< 1}.
#' @param method string that indicates the method used for sampling 
#'          observations. The available methods are 
#' \itemize{
#'    \item \code{'rejvmf'} acceptance-rejection algorithm using 
#'                         von Mises-Fisher envelops (Algorithm in Table 2 of 
#'                         Golzy and Markatou 2020);
#'    \item \code{'rejacg'} using angular central Gaussian envelops 
#'                         (Algorithm in Table 1 of Sablica et al. 2023);
#'    \item \code{'rejpsaw'} using projected Saw distributions 
#'                         (Algorithm in Table 2 of Sablica et al. 2023).
#' }
#' @param max.iter the maximum number of iterations (for 'rejacg' method).
#' @param tol.eps the desired accuracy of convergence tolerance 
#'               (for 'rejacg' method).
#' 
#' @details
#' The number of observations generated is determined by \code{n} for 
#' \code{rpkb()}. This function returns a list with the matrix of generated 
#' observations \code{x}, the number of tries \code{numTries} and the number of
#' acceptances \code{numAccepted}.
#' 
#' A limitation of the \code{rejvmf} is that the method does not ensure the
#' computational feasibility of the sampler for \eqn{\rho} approaching 1.
#' 
#' If the chosen method is 'rejacg', the function \code{uniroot}, from the 
#' \code{stat} package, is used to estimate the beta parameter. In this case, 
#' the complete results are provided as output. 
#' 
#' @note
#' If the required packages (`movMF` for \code{rejvmf} method, and
#' `Tinflex` for \code{rejpsaw}) are not installed, the function will display a
#' message asking the user to install the missing package(s).
#' 
#' @references
#' Golzy, M. and Markatou, M. (2020) Poisson Kernel-Based Clustering on the 
#' Sphere: Convergence Properties, Identifiability, and a Method of Sampling, 
#' Journal of Computational and Graphical Statistics, 29:4, 758-770, 
#' DOI: 10.1080/10618600.2020.1740713.
#' 
#' Sablica L., Hornik K. and Leydold J. (2023) "Efficient sampling from the PKBD
#' distribution", Electronic Journal of Statistics, 17(2), 2180-2209.
#'
#' @srrstats {G1.0} Reference section reports the related literature.
#' @srrstats {G1.4} roxigen2 is used.
#' @srrstats {G2.0,G2.0a} check input mu. 
#' @srrstats {PD1.0} references on choice and usage of the probability 
#'                   distributions are provided.
#' @srrstats {PD3.2} all parameters of uniroot function are specified
#'                   (for method 'rejacg').
#' @srrstats {PD3.3} the complete results of uniroot function are provided as 
#'                   output (for method 'rejacg').
#' @srrstats {PD3.1} the three algorithms for random generation are implemented
#'                   as separate functions.
#' 
#' @export
rpkb <- function(n, mu, rho, method = 'rejacg', 
                 tol.eps = .Machine$double.eps^0.25, 
                 max.iter = 1000) {
   

   # nocov start
   if (!requireNamespace("movMF", quietly = TRUE)) {
      install <- readline(prompt = "'movMF' is required for 'rejvmf' method.
                           Would you like to install it now? (yes/no): ")
      if (tolower(install) == "yes") {
         install.packages("movMF")
      } else {
         message("'rejvmf' cannot be performed without 'movMF'.")
         return(NULL)
      }
   }
   
   if (!requireNamespace("Tinflex", quietly = TRUE)) {
      install <- readline(prompt = "'Tinflex' is required for 'rejpsaw' method. 
                           Would you like to install it now? (yes/no): ")
      if (tolower(install) == "yes") {
         install.packages("Tinflex")
      } else {
         message("'rejpsaw' cannot be performed without 'Tinflex'.")
         return(NULL)
      }
   }
   # nocov end
   if (rho >= 1 | rho < 0) {
      stop('Input argument rho must be within [0,1)')
   }
   # norm input
   if (length(mu) < 2) {
      stop('mu must have length >= 2')
   }
   if (sum(abs(mu)) == 0) {
      stop('Input argument mu cannot be a vector of zeros')
   }
   if (!is.numeric(n) | n < 0) {
      stop('n must be a positive integer')
   }
   # norm input
   if (!method%in%c('rejvmf','rejacg','rejpsaw')) {
      stop('Unknown method')
   }
   
   mu <- mu / sqrt(sum(mu^2))
   p <- length(mu)
   
   if (method == 'rejvmf') {
      
      res <- rejvmf(n, rho, mu, p)
      
   } else if (method == 'rejacg') {
      
      res <- rejacg(n, rho, mu, p, tol.eps, max.iter)
      
   } else if (method == 'rejpsaw') {
      
      res <- rejpsaw(n, rho, mu, p)
      
   }
   
   return(res)
}
#'
#' Random generation with vonMises distribution evelopes.
#'
#' @param n number of observations.
#' @param rho is the concentration parameter, with 0 <= rho < 1.
#' @param mu location vector parameter with length indicating the dimension of 
#'          generated points.
#' @param p dimension.
#' 
#' @references
#' Golzy, M. and Markatou, M. (2020) Poisson Kernel-Based Clustering on the 
#' Sphere: Convergence Properties, Identifiability, and a Method of Sampling, 
#' Journal of Computational and Graphical Statistics, 29:4, 758-770, 
#' DOI: 10.1080/10618600.2020.1740713.
#'
#' @importFrom stats runif
#' 
#' @noRd
#' @keywords internal
rejvmf <- function(n, rho, mu, p) {
   
   kappa<-p*rho/(1+rho^2)
   M1 <- (1 + rho) / ((1 - rho)^(p-1)*exp(kappa))
   retvals <- matrix(nrow = n, ncol = p)
   numAccepted <- 0
   numTries <- 0
   while (numAccepted < n) {
      numTries <- numTries + 1
      theta<-kappa*mu
      yx<- movMF::rmovMF(1,theta,1)
      v<- yx%*%mu
      f1<- (1-rho**2)/ (( 1 + rho**2 - 2*rho*v)^(p/2))
      u <- runif(1)
      if (u <= f1/(M1 * exp(kappa*v))) {
         numAccepted <- numAccepted + 1
         retvals[numAccepted,] <- yx
      }
   }
   
   res <- list(x = retvals, numAccepted = numAccepted, numTries = numTries)
   
   return(res)
}
#'
#' Random generation with angular central gaussian evelopes.
#'
#' @param n number of observations.
#' @param rho is the concentration parameter, with 0 <= rho < 1.
#' @param mu location vector parameter with length indicating the dimension of 
#'          generated points.
#' @param p dimension.
#' 
#' @references
#' Sablica L., Hornik K. and Leydold J. (2023) "Efficient sampling from the PKBD
#' distribution", Electronic Journal of Statistics, 17(2), 2180-2209.
#' 
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats uniroot
#' 
#' @noRd
#' @keywords internal
rejacg <- function(n, rho, mu, p, tol.eps, max.iter){
      
   lambda <- (2*rho)/(1+rho^2)
   beta_lambda <- lambda/(2 - lambda)
   
   C_d_lambda <- function(beta){
      return(-4*(p - 1)*beta^3 + (4*p - lambda^2*(p - 2)^2)*beta^2 + 
                2*p*(p - 2)*lambda^2*beta - p^2*lambda^2)
   }
   
   beta_star <- uniroot(C_d_lambda, interval = c(beta_lambda, 1), 
                        tol = tol.eps, maxiter = max.iter)
   beta_res <- beta_star
   beta_star <- beta_star$root
   
   b1 <- beta_star/(1 - beta_star)
   b2 <- -1 + 1/sqrt(1 - beta_star)
   
   retvals <- matrix(nrow = n, ncol = p)
   numAccepted <- 0
   numTries <- 0
   while (numAccepted < n) {
      
      numTries <- numTries + 1
      # Step 5
      u <- runif(1)
      # Step 6-7
      z <- rnorm(p)
      # Step 8
      q <- (t(mu)%*%z + b2*t(mu)%*%z)/sqrt(t(z)%*%z + b1*(t(mu)%*%z)^2)
      
      # Step 9
      M1 <- p/2 * (-log(1 - lambda*q) +log(1 - beta_star*q^2) - 
                      log(2/(1 + sqrt(1 - lambda^2/beta_star))))
      
      yx <- (z + b2*t(mu)%*%z%*%mu)/as.numeric(sqrt(t(z)%*%z + 
                                                       b1*(t(mu)%*%z)^2))
      
      if (log(u) <= M1) {
         numAccepted <- numAccepted + 1
         retvals[numAccepted,] <- yx
      }
   }
   
   res <- list(x = retvals, numAccepted = numAccepted, numTries = numTries, 
               beta = beta_res)
   
   return(res)
}
#'
#' Random generation with projected Saw distribution.
#'
#' @param n number of observations.
#' @param rho is the concentration parameter, with 0 <= rho < 1.
#' @param mu location vector parameter with length indicating the dimension of 
#'          generated points.
#' @param p dimension.
#' 
#' @references
#' Sablica L., Hornik K. and Leydold J. (2023) "Efficient sampling from the PKBD
#' distribution", Electronic Journal of Statistics, 17(2), 2180-2209.
#' 
#' 
#' @noRd
#' @keywords internal
rejpsaw <- function(n, rho, mu, p){

   retvals <- vapply(1:n, function(x){
      lambda <- (2*rho)/(1+rho^2)
      # pSaw: log-density and derivatives
      lpdf <- function(t) {-p*log(1-lambda*t) + (p-3)*(log(1+t) + log(1-t))}
      dlpdf <- function(t) {p*lambda/(1-lambda*t) - (p-3)*(2*t/(1-t^2))}
      d2lpdf <- function(t) {(p*lambda^2)/((1-lambda*t)^2) - 
            2*(p-3)*(1+t^2)/((1-t^2)^2)}
      
      t1 <- p*lambda/(p-3 - sqrt((p-3)^2 - p*(p-6)*lambda^2))
      t2 <- p*lambda/(p-3 + sqrt((p-3)^2 - p*(p-6)*lambda^2))
      if(p>3) ib <- c(-1,t2)
      if(p<3) ib <- c(t1,0.9999999)
      if(p==3) ib <- c(t1+0.000001,t2-0.000001)
      
      ## Create generator object for pSaw distribution
      pSaw <- Tinflex::Tinflex.setup.C(lpdf, dlpdf, d2lpdf, 
                                       ib=ib,cT=1, rho=1.05)
      ## Print data about generator object.
      #print(gen)
      
      # Step 1
      W <- Tinflex::Tinflex.sample(pSaw, n=1)
      # Step 2
      y <- sample_hypersphere(p,1)
      # Step 3
      y <- y - mu%*%t(y)%*%t(mu)
      # Step 4
      y <- y/as.numeric(sqrt(y%*%t(y)))
      # Step 5
      
      xy <- W*mu + sqrt(1 - W^2)*y
      return(xy)
      
   }, numeric(p))
   
   retvals <- t(retvals)
   numAccepted <- NA
   numTries <- NA
   
   res <- list(x = retvals, numAccepted = numAccepted, numTries = numTries)
   return(res)   
}