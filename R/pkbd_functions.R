#'
#' The Poisson kernel-based Distribution (PKBD)
#' 
#' Density function and random number generation from the Poisson kernel-based Distribution
#' with mean direction vector \code{mu} and concentration parameter \code{rho}.
#'  
#' @param x Matrix (or data.frame) with number of columns >=2.
#' @param mu Location parameter with same length as the rows of x. Normalized to length one.
#' @param rho Concentration parameter, with 0 <= rho < 1.
#' @param logdens Logical; if 'TRUE', densities d are given as log(d).
#'
#' @return 
#' \code{dpkb} gives the density value.
#' \code{rpkb} generates random observations from the PKBD. 
#' 
#' The number of observations generated is determined by \code{n} for \code{rpkb}.
#' This function returns a list with the matrix of generated observations \code{x},
#' the number of tries \code{numTries} and the number of acceptances \code{numAccepted}.
#'
#' 
#' @examples
#' # Generate some data from pkbd density
#' pkbd_dat <- rpkb(10, c(0.5,0), 0.5)
#' 
#' # Calculate the PKBD density values
#' dens_val <- dpkb(pkbd_dat$x, c(0.5,0.5),0.5)
#'
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
   p <- ncol(x)
   if (p < 2) {
      stop('vectors must have length >= 2')
   }
   if (length(mu) != p) {
      stop('vectors and mu must have the same length')
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
#' @param mu location vector parameter with length indicating the dimension of generated points.
#' @param rho is the concentration parameter, with 0 <= rho < 1.
#' @param method string that indicates the method used for sampling observations. The available methods are 
#' \itemize{
#'    \item \code{'rejvmf'} acceptance-rejection algorithm using von Mises-Fisher envelops (Algorithm in Table 2 of Golzy and Markatou 2020);
#'    \item \code{'rejacg'} using angular central Gaussian envelops (Algorithm in Table 1 of Sablica et al. 2023);
#'    \item \code{'rejpsaw'} using projected Saw distributions (Algorithm in Table 2 of Sablica et al. 2023).
#' }
#' 
#' @references
#' Golzy, M., Markatou, M. (2020) Poisson Kernel-Based Clustering on the Sphere:
#' Convergence Properties, Identifiability, and a Method of Sampling, Journal of
#' Computational and Graphical Statistics, 29:4, 758-770, 
#' DOI: 10.1080/10618600.2020.1740713.
#' 
#' Sablica L., Hornik K., Leydold J. (2023) "Efficient sampling from the PKBD 
#' distribution", Electronic Journal of Statistics, 17(2), 2180-2209.
#'
#' @import movMF
#' @import Tinflex
#'
#' @export
rpkb <- function(n, mu, rho, method = 'rejvmf') {
   
   if (rho >= 1 | rho < 0) {
      stop('Input argument rho must be within [0,1)')
   }
   # norm input
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
      kappa<-p*rho/(1+rho^2)
      M1 <- (1 + rho) / ((1 - rho)^(p-1)*exp(kappa))
      retvals <- matrix(nrow = n, ncol = p)
      numAccepted <- 0
      numTries <- 0
      while (numAccepted < n) {
         numTries <- numTries + 1
         theta<-kappa*mu
         yx<- rmovMF(1,theta,1)
         v<- yx%*%mu
         f1<- (1-rho**2)/ (( 1 + rho**2 - 2*rho*v)^(p/2))
         u <- runif(1)
         if (u <= f1/(M1 * exp(kappa*v))) {
            numAccepted <- numAccepted + 1
            retvals[numAccepted,] <- yx
         }
      }
      
   } else if (method == 'rejacg') {
      
      lambda <- (2*rho)/(1+rho^2)
      beta_lambda <- lambda/(2 - lambda)
      
      C_d_lambda <- function(beta){
         return(-4*(p - 1)*beta^3 + (4*p - lambda^2*(p - 2)^2)*beta^2 + 2*p*(p - 2)*lambda^2*beta - p^2*lambda^2)
      }
      
      beta_star <- uniroot(C_d_lambda, interval = c(beta_lambda, 1))
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
         M1 <- p/2 * (-log(1 - lambda*q) +log(1 - beta_star*q^2) - log(2/(1 + sqrt(1 - lambda^2/beta_star))))
         
         yx <- (z + b2*t(mu)%*%z%*%mu)/as.numeric(sqrt(t(z)%*%z + b1*(t(mu)%*%z)^2))
         
         if (log(u) <= M1) {
            numAccepted <- numAccepted + 1
            retvals[numAccepted,] <- yx
         }
      }
      
   } else if (method == 'rejpsaw') {
      
      retvals <- sapply(1:n, function(x){
         lambda <- (2*rho)/(1+rho^2)
         # pSaw: log-density and derivatives
         lpdf <- function(t) {-p*log(1-lambda*t) + (p-3)*(log(1+t) + log(1-t))}
         dlpdf <- function(t) {p*lambda/(1-lambda*t) - (p-3)*(2*t/(1-t^2))}
         d2lpdf <- function(t) {(p*lambda^2)/((1-lambda*t)^2) - 2*(p-3)*(1+t^2)/((1-t^2)^2)}
         
         t1 <- p*lambda/(p-3 - sqrt((p-3)^2 - p*(p-6)*lambda^2))
         t2 <- p*lambda/(p-3 + sqrt((p-3)^2 - p*(p-6)*lambda^2))
         if(p>3) ib=c(-1,t2)
         if(p<3) ib=c(t1,0.9999999)
         if(p==3) ib=c(t1+0.000001,t2-0.000001)
         
         ## Create generator object for pSaw distribution
         pSaw <- Tinflex.setup.C(lpdf, dlpdf, d2lpdf, ib=ib,cT=1, rho=1.05)
         ## Print data about generator object.
         #print(gen)
         
         # Step 1
         W <- Tinflex.sample(pSaw, n=1)
         # Step 2
         y <- sample_hypersphere(p,1)
         # Step 3
         y <- y - mu%*%t(y)%*%t(mu)
         # Step 4
         y <- y/as.numeric(sqrt(y%*%t(y)))
         # Step 5
         
         xy <- W*mu + sqrt(1 - W^2)*y
         return(xy)
         
      })
      
      retvals <- t(retvals)
      numAccepted <- NA
      numTries <- NA
      
   }
   return(list(x = retvals, numAccepted = numAccepted, numTries = numTries))
}
