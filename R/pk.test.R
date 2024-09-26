#' 
#' Poisson kernel-based quadratic distance test of Uniformity on the sphere
#' 
#' @description
#'  This function performs the kernel-based quadratic distance goodness-of-fit 
#' tests for Uniformity for spherical data \code{x} using the Poisson kernel 
#' with concentration parameter \code{rho}. \cr
#' The Poisson kernel-based test for uniformity exhibits excellent results 
#' especially in the case of multimodal distributions, as shown in the example
#' of the \href{../doc/uniformity.html}{Uniformity test on the Sphere vignette}.
#' 
#' @details
#' Let \eqn{x_1, x_2, ..., x_n} be a random sample with empirical distribution
#' function \eqn{\hat F}. 
#' We test the null hypothesis of uniformity on the 
#' \eqn{d}-dimensional sphere, i.e. \eqn{H_0:F=G}, where \eqn{G} is the uniform
#' distribution on the \eqn{d}-dimensional sphere \eqn{\mathcal{S}^{d-1}}.
#' We compute the U-statistic estimate of the sample KBQD (Kernel-Based 
#' Quadratic Distance)
#' \deqn{U_{n}=\frac{1}{n(n-1)}\sum_{i=2}^{n}\sum_{j=1}^{i-1}K_{cen}
#' (\mathbf{x}_{i}, \mathbf{x}_{j}),}
#' then the first test statistic is given as
#' \deqn{T_{n}=\frac{U_{n}}{\sqrt{Var(U_{n})}},}
#' with
#' \deqn{Var(U_{n})= \frac{2}{n(n-1)}
#' \left[\frac{1+\rho^{2}}{(1-\rho^{2})^{d-1}}-1\right],}
#' and the V-statistic estimate of the KBQD  
#' \deqn{V_{n} = \frac{1}{n}\sum_{i=1}^{n}\sum_{j=1}^{n}K_{cen}
#' (\mathbf{x}_{i}, \mathbf{x}_{j}),}
#' where \eqn{K_{cen}} denotes the Poisson kernel \eqn{K_\rho} centered with
#' respect to the uniform distribution on the \eqn{d}-dimensional sphere, that 
#' is
#' \deqn{K_{cen}(\mathbf{u}, \mathbf{v}) = K_\rho(\mathbf{u}, \mathbf{v}) -1} 
#' and
#' \deqn{K_\rho(\mathbf{u}, \mathbf{v}) = \frac{1-\rho^{2}}{\left(1+\rho^{2}-
#' 2\rho (\mathbf{u}\cdot \mathbf{v})\right)^{d/2}},}
#' for every \eqn{\mathbf{u}, \mathbf{v} \in \mathcal{S}^{d-1} 
#' \times \mathcal{S}^{d-1}}.
#' 
#' The asymptotic distribution of the V-statistic is an infinite combination
#' of weighted independent chi-squared random variables with one degree of 
#' freedom. The cutoff value is obtained using the Satterthwaite approximation 
#' \eqn{c \cdot \chi_{DOF}^2}, where \deqn{c=\frac{(1+\rho^{2})-
#' (1-\rho^{2})^{d-1}}{(1+\rho)^{d}-(1-\rho^{2})^{d-1}}} and \deqn{DOF(K_{cen}
#' )=\left(\frac{1+\rho}{1-\rho} \right)^{d-1}\left\{ 
#' \frac{\left(1+\rho-(1-\rho)^{d-1} \right )^{2}}
#' {1+\rho^{2}-(1-\rho^{2})^{d-1}}\right \}.}.
#' For the \eqn{U}-statistic the cutoff is determined empirically:
#' -  Generate data from a Uniform distribution on the d-dimensional sphere;
#' - Compute the test statistics for \code{B} Monte Carlo(MC) replications;
#' - Compute the 95th quantile of the empirical distribution of the test
#'   statistic.
#' 
#' @seealso \linkS4class{pk.test}
#' 
#' @note
#' A U-statistic is a type of statistic that is used to estimate a population
#' parameter. It is based on the idea of averaging over all possible *distinct*
#' combinations of a fixed size from a sample. 
#' A V-statistic considers all possible tuples of a certain size, not just
#' distinct combinations and can be used in contexts where unbiasedness is not
#' required.
#'
#' @param x A numeric d-dim matrix of data points on the Sphere S^(d-1).
#' @param rho Concentration parameter of the Poisson kernel function.
#' @param B Number of Monte Carlo iterations for critical value estimation of Un
#'          (default: 300).
#' @param Quantile The quantile to use for critical value estimation, 
#'                 0.95 is the default value.
#'
#' @return An S4 object of class \code{pk.test} containing the results of the 
#' Poisson kernel-based tests. The object contains the following slots:
#'\itemize{
#'   \item \code{method}: Description of the test performed.
#'   \item \code{x} Data matrix.
#'   \item \code{Un} The value of the U-statistic.
#'   \item \code{CV_Un} The empirical critical value for Un.
#'   \item \code{H0_Vn} A logical value indicating whether or not the null 
#'                      hypothesis is rejected according to Un.
#'   \item \code{Vn} The value of the V-statistic Vn.
#'   \item \code{CV_Vn} The critical value for Vn computed following the 
#'                      asymptotic distribution.
#'   \item \code{H0_Vn} A logical value indicating whether or not the null 
#'                      hypothesis is rejected according to Vn.
#'   \item \code{rho} The value of concentration parameter used for the Poisson
#'                    kernel function.
#'   \item \code{B} Number of replications for the critical value of the 
#'                  U-statistic Un.    
#'}
#'
#'
#' @references
#' Ding, Y., Markatou, M. and Saraceno, G. (2023). “Poisson Kernel-Based Tests 
#' for Uniformity on the d-Dimensional Sphere.” Statistica Sinica. 
#' doi:10.5705/ss.202022.0347
#' 
#' @examples
#' # create a pk.test object
#' x_sp <- sample_hypersphere(3, n_points=100)
#' unif_test <- pk.test(x_sp,rho=0.8)
#' unif_test
#'
#' @importFrom stats qchisq
#'
#' @useDynLib QuadratiK
#' @importFrom Rcpp sourceCpp
#'
#' @srrstats {G1.0} Reference section reports the related literature
#' @srrstats {G1.3} description of parameter
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @export
setGeneric("pk.test",function(x, rho, B = 300, Quantile = 0.95){
   standardGeneric("pk.test")
})
#' @rdname pk.test
#' 
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G2.7,G2.8} different types of input are considered
#' @srrstats {G2.13,G2.14,G2.14a,G2.15,G2.16} error for NA, Nan, Inf, -Inf
#' 
#' @export
setMethod("pk.test", signature(x = "ANY"),
          function(x,  rho, B = 300, Quantile = 0.95)
          {
             
             if(!is.numeric(x) & !is.data.frame(x)){
                stop("x must be numeric")
             } else if(is.data.frame(x)) {
                x <- as.matrix(x)
             } else if(!is.matrix(x)){
                stop("x must be a matrix or a data.frame with dimension greater
                     than 1.")
             }
             if(any(is.na(x))){
                stop("There are missing values in x!")
             } else if(any(is.infinite(x) |is.nan(x))){
                stop("There are undefined values in x, that is Nan, Inf, -Inf")
             }
             
             if(Quantile<=0 | Quantile>1){
                stop("Quantile must be in (0,1].")
             }
             if(!is.numeric(rho) | (rho<=0 | rho>1)){
                stop("rho must be in (0,1).")
             }
             
             
             METHOD <- "Poisson Kernel-based quadratic distance test of 
                        Uniformity on the Sphere"
             
             d<- ncol(x)
             n <- nrow(x)
             
             pk <- statPoissonUnif(x, rho)
             
             #x.prod <- x%*%t(x)
             #pk <- ((1-rho^2)/ (( 1 + rho^2 - 2*rho*x.prod)^(d/2))) - 1
             
             dof <- DOF(d, rho)
             qu_q <- qchisq(Quantile,df=dof$DOF)
             CV_Vn <- dof$Coefficient*qu_q
             
             var_Un <- (2/(n*(n-1)))*((1+rho^2)/((1-rho^2)^(d-1)) -1)
             
             CV_Un <- poisson_CV(d=d, size=n, rho=rho, B=B, Quantile=Quantile )
             CV_Un <- CV_Un/sqrt(var_Un)
             
             res <- new("pk.test", Un = pk[1]/sqrt(var_Un), CV_Un = CV_Un, 
                        Vn = pk[2], CV_Vn = CV_Vn, method = METHOD, x = x, B= B,
                        rho= rho, H0_Un = pk[1]/sqrt(var_Un) > CV_Un, 
                        H0_Vn = pk[2] > CV_Vn, var_Un=var_Un)
             return(res)
          })
#' @rdname pk.test
#'
#' @param object Object of class \code{pk.test}
#' 
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @export
setMethod("show", "pk.test",
          function(object) {
             cat( "\n", object@method, "\n")
             cat("Selected consentration parameter rho: ", object@rho, "\n")
             
             cat("\n")
             cat("U-statistic:\n")
             cat("\n")
             cat("H0 is rejected: ", object@H0_Un , "\n")
             cat("Statistic Un: ", object@Un, "\n")
             cat("Critical value: ", object@CV_Un,"\n")
             
             cat("\n")
             cat("V-statistic:\n")
             cat("\n")
             cat("H0 is rejected: ", object@H0_Vn, "\n")
             cat("Statistic Vn: ", object@Vn, "\n")
             cat("Critical value: ", object@CV_Vn,"\n")
             
             cat("\n")
          })
#' Summarizing kernel-based quadratic distance results
#'
#' \code{summary} method for the class \code{pk.test}
#'
#' @param object Object of class \code{pk.test}
#' 
#' @return List with the following components:
#' \itemize{
#'    \item \code{summary_tables} Table of computed descriptive statistics per 
#'                                variable.
#'    \item \code{test_results} Data frame with the results of the performed 
#'                              Poisson kernel-based test.
#'    \item \code{qqplots} Figure with qq-plots for each variable against the 
#'                         uniform distribution.
#' }
#' 
#' @seealso [pk.test()] and \linkS4class{pk.test} for additional details.
#'
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' 
#' @examples
#' # create a pk.test object
#' x_sp <- sample_hypersphere(3, n_points=100)
#' unif_test <- pk.test(x_sp,rho=0.8)
#' summary(unif_test)
#' 
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @importFrom stats IQR
#' @importFrom stats median
#' @importFrom stats sd
#' @importFrom stats runif
#' 
#' @name summary.pk.test
#' @rdname summary.pk.test
#' @aliases summary,pk.test-method
#' @export
#'
setMethod("summary", "pk.test", function(object) {
   
   dat_x <- as.data.frame(object@x)
   
   figure <- NA
   
   plot_list <- list()
   stats <- list()
   for(i in seq_len(ncol(dat_x))) {
      
      unif_data <- runif(nrow(dat_x),-1,1)
      probs <- seq(0, 1, length.out = nrow(dat_x))
      # qq_df <- data.frame(
      #    x = quantile(unif_data, probs = seq(0, 1, length.out = nrow(dat_x))),
      #    sample_quantiles = quantile(dat_x[,i], probs = probs))
      x <- quantile(unif_data, probs = seq(0, 1, length.out = nrow(dat_x))) 
      sample_quantiles <- quantile(dat_x[,i], probs = probs)
      
      pl <- ggplot(mapping=aes(x = x, y = sample_quantiles)) +
         geom_line(col="blue") +
         theme_minimal()+
         geom_abline(slope = 1, intercept = 0,col="red") +
         ggtitle(paste("QQ Plot against Uniform - ",names(dat_x)[i])) +
         xlab("Theoretical Quantiles") +
         ylab("Sample Quantiles")
      
      
      stats_step <- data.frame(matrix(c(mean(dat_x[,i]),sd(dat_x[,i]),
                                        median(dat_x[,i]),IQR(dat_x[,i]),
                                        min(dat_x[,i]),max(dat_x[,i])),
                                      nrow=6,ncol=1,byrow=TRUE))
      colnames(stats_step) <- c(paste(names(dat_x)[i]))
      rownames(stats_step) <- c("mean", "sd", "median", "IQR", "min", "max")
      
      stats[[i]] <- stats_step
      
      plot_list[[length(plot_list) + 1]] <- list(pl)
      
      
   }
   stats <- do.call(cbind, stats)
   colnames(stats) <- c(names(dat_x))
   rownames(stats) <- c("mean", "sd", "median", "IQR", "min", "max")
   
   plot_list <- do.call(c, plot_list)
   figure <- ggarrange(plotlist = plot_list, ncol = 1)
   
   # Print main results of the test
   cat( "\n", object@method, "\n")
   test_results <- data.frame(
      Test_Statistics = c(object@Un, object@Vn),
      Critical_Value = c(object@CV_Un,object@CV_Vn),
      Reject_H0 = c(object@H0_Un, object@H0_Vn)
   )
   print(test_results)
   print(figure)
   return(list(summary_tables = stats, 
               test_results = test_results, 
               qqplots = figure))
})
