#' Select the value of the kernel tuning parameter
#'
#' This function computes the kernel bandwidth of the Gaussian kernel for the two-sample and k-sample kernel-based quadratic distance (KBQD) tests.
#'
#' @param x Data set of observations from X.
#' @param y Data set of observations from Y.
#' @param alternative Family of alternative chosen for selecting h, between "location", "scale" and "skewness".
#' @param method The method used for critical value estimation ("subsampling", "bootstrap", or "permutation").
#' @param b The size of the subsamples used in the subsampling algorithm .
#' @param B The number of iterations to use for critical value estimation, B = 150 as default.
#' @param delta_dim Vector of coefficient of alternative with respect to each dimension
#' @param delta Vector of parameter values indicating chosen alternatives
#' @param h_values Values of the tuning parameter used for the selection
#' @param Nrep Number of bootstrap/permutation/subsampling replications.
#' @param n_cores Number of cores used to parallel the h selection algorithm (default:2).
#' @param Quantile The quantile to use for critical value estimation, 0.95 is the default value.
#' @param power.plot Logical. If TRUE, it is displayed the plot of power for values in h_values and delta.
#' 
#' @return A list with the following attributes:
#' \itemize{
#'    \item \code{h_sel} the selected value of tuning parameter h;
#'    \item \code{power} matrix of power values computed for the considered values of \code{delta} and \code{h_values};
#'    \item \code{power.plot} power plots (if \code{power.plot} is \code{TRUE}).
#' }
#' @details
#' The function performs the selection of the optimal value for the tuning parameter \eqn{h} of the normal kernel function, for the two-sample and k-sample KBQD tests. It performs a small simulation study, generating sample according to the family of \code{alternative} specified, for the chosen values of \code{h_values} and \code{delta}.
#' 
#' @references
#' Markatou, M., Saraceno, G., Chen, Y. (2023). “Two- and k-Sample Tests Based on
#' Quadratic Distances.” Manuscript, (Department of Biostatistics, University at
#' Buffalo)
#' 
#' @examples
#' # Select the value of h using the mid-power algorithm
#' \donttest{
#' x <- matrix(rnorm(100),ncol=2)
#' y <- matrix(rnorm(100),ncol=2)
#' h_sel <- select_h(x,y,"skewness")
#' h_sel
#' }
#'
#' @importFrom moments skewness
#' @importFrom sn rmsn
#' @import doParallel
#' @import foreach
#' @import stats
#' @import rlecuyer
#' @import ggplot2
#' @import RcppEigen
#'
#' @useDynLib QuadratiK
#'
#' @export
select_h <- function(x, y=NULL, alternative=NULL, method="subsampling", b=0.8, B=100, delta_dim=1, delta=NULL, h_values=NULL,Nrep=50, n_cores=2, Quantile=0.95, power.plot=TRUE) {
   
   
   # Convert vectors to a single column matrix
   if(is.vector(x)) {
      x <- matrix(x, ncol = 1)
   } else if(is.data.frame(x)) {
      x <- as.matrix(x)
   } else if(!is.matrix(x)){
      stop("x must be a matrix or a data.frame")
   }
   
   if(!is.null(y)){
      if(is.vector(y) | is.factor(y)) {
         y <- matrix(as.numeric(y), ncol = 1)
      } else if(is.data.frame(y)) {
         y <- as.matrix(y)
      } else if(!is.matrix(y)){
         stop("y must be a vector or a matrix")
      }
   } #else {
      #stop("y must be provided.")
   #}
   
   if(!(alternative %in% c("location", "scale", "skewness"))){
      
      stop("The alternative argument should be one of 'location', 'scale' or 'skewness'")
   }
   
   n <- nrow(x)
   d <- ncol(x)

   if (is.null(h_values)) h_values <- seq(0.4, 3.3, 0.4)
   
   if(is.null(delta)){
      if(alternative=='location'){
         delta <- c(0.2,0.3,0.4)
      } else if(alternative=='scale'){
         delta <- 1 + c(0.1,0.3,0.5)
      } else if(alternative=='skewness'){
         delta <- c(0.2, 0.3, 0.6)
      }
   }
   
   if(!is.null(y)){
      
      d_y <- ncol(y)
      if(d_y > 1){
         
         # Check that they have the same number of columns (features):
         if(ncol(x) != ncol(y)) {
            stop("'x' and 'y' must have the same number of columns.")
         }
         
         m <- nrow(y)
         pooled <- rbind(x, y)
         
      } else if(d_y==1){
         
         # if(nrow(x) != length(y)) {
         #    stop("Number of rows of 'x' and the length of 'y' must coincide.")
         # }
         if(nrow(x) != nrow(y)) {
          stop("'x' and 'y' must have the same number of rows.")
         }
         
         K <- length(unique(y))
         nk <- round(n/K)
         pooled <- x
      }
      
   } else {
      pooled <- x
   }
   
   if(length(delta_dim)==1){
      delta_dim <- rep(1, d)
   }else{
      if (!is.numeric(delta_dim) || length(delta_dim) != d) {
         stop("delta_dim must be 1 or a numeric vector of length equal to the number of columns of pooled.")
      }
   }
   
   # Compute estimates of mean, covariance and skewness from the pooled sample
   mean_dat <- colMeans(pooled)
   S_dat <- cov(pooled)
   S_dat <- diag(diag(S_dat), nrow = d, ncol = d)
   skew_data <- skewness(pooled)
   
   # Define the objective function for the alternative of the two-sample test
   objective_2 <- function(h,k) {
      
      dk <- delta_dim*delta[k]
      
      if(alternative=='location'){
         xnew <- sn::rmsn(n, xi = mean_dat, Omega = S_dat, alpha = skew_data)
         ynew <- sn::rmsn(m, xi = mean_dat + dk, Omega = S_dat, alpha = skew_data)
      } else if(alternative=='scale'){
         xnew <- sn::rmsn(n, xi = mean_dat, Omega = S_dat, alpha = skew_data)
         ynew <- sn::rmsn(m, xi = mean_dat , Omega = dk*S_dat, alpha = skew_data)
      } else if(alternative=='skewness'){
         xnew <- sn::rmsn(n, xi = mean_dat, Omega = S_dat, alpha = skew_data)
         ynew <- sn::rmsn(m, xi = mean_dat, Omega = S_dat, alpha = skew_data+dk)
      }
      
      #x <- as.matrix(xnew[sample(n, replace = FALSE),])
      #y <- as.matrix(ynew[sample(m, replace = FALSE),])
      
      STATISTIC <- stat2sample(xnew, ynew, h, matrix(0,nrow=1),diag(1),"Nonparam")
      CV <- compute_CV(B, Quantile, pooled, n, m, h, method, b)
      
      return(STATISTIC < CV)
   }
   # Define the objective function for the alternative k-sample test
   objective_k <- function(h,k) {
      
      dk <- delta_dim*delta[k]
      
      if(alternative=='location'){
         mean_tilde <- mean_dat + dk
         S_tilde <- S_dat
         skew_tilde <- skew_data
      } else if(alternative=='scale'){
         mean_tilde <- mean_dat
         S_tilde <- S_dat*dk
         skew_tilde <- skew_data
      } else if(alternative=='skewness'){
         mean_tilde <- mean_dat
         skew_tilde <- skew_data + dk
         S_tilde <- S_dat
      }
      
      xnew <- rbind(
         sn::rmsn(nk*(K-1), xi = mean_dat, Omega = S_dat, alpha = skew_data),
         sn::rmsn(nk, xi = mean_tilde, Omega = S_tilde, alpha = skew_tilde))
      ynew <- matrix(rep(1:K, each=nk),ncol=1)
      
      
      sizes_new <- as.vector(table(ynew))
      cum_size_new <- c(0,cumsum(sizes_new))
      STATISTIC <- stat_ksample_cpp(xnew, ynew, h, sizes_new, cum_size_new)
      CV <- cv_ksample(xnew, ynew, h, B, b, Quantile, method)
      
      return(c(STATISTIC < CV))
   }
   # Define the objective function for the alternative normality test
   objective_norm <- function(h,k) {
      
      dk <- delta_dim*delta[k]
      
      if(alternative=='location'){
         mean_tilde <- mean_dat + dk
         S_tilde <- S_dat
         skew_tilde <- skew_data
      } else if(alternative=='scale'){
         mean_tilde <- mean_dat
         S_tilde <- S_dat*dk
         skew_tilde <- skew_data
      } else if(alternative=='skewness'){
         mean_tilde <- mean_dat
         skew_tilde <- skew_data + dk
         S_tilde <- S_dat
      }
      
      xnew <- sn::rmsn(n, xi = mean_tilde, Omega = S_tilde, alpha = skew_tilde)
      
      STATISTIC <- kbNormTest(xnew, h, matrix(mean_dat,nrow = 1), S_dat, centeringType = "Param")
      CV <- normal_CV(d,n,h,matrix(mean_dat,nrow = 1),S_dat,B,Quantile)
      
      return(c(STATISTIC < CV))
   }
   
   num_cores <- as.numeric(n_cores)
   registerDoParallel(cores=n_cores)
   
   D <- length(delta)
   
   k_values <- 1:D
   rep_values <- 1:Nrep
   
   params <- expand.grid(Rep=rep_values, h = h_values)
   params <- split(params, seq(nrow(params)))
   
   res <- data.frame(Rep=numeric(), delta=numeric(), h=numeric(), power=numeric())
   pars = NULL
   
   for(k in k_values){
      results <- foreach(pars = params, .combine = rbind, .packages=c("sn", "moments", "stats", "rlecuyer","QuadratiK")) %dopar% {
         h <- as.numeric(pars$h)
         if(is.null(y)){
            
            objective_result <- objective_norm(h,k)
         
         } else if(d_y > 1){
            
            objective_result <- objective_2(h, k)
            
         } else if(d_y == 1){
            
            objective_result <- objective_k(h, k)
         }
         data.frame(Rep=pars$Rep, delta=delta[k], h=h, score=objective_result)
      }
      
      results$score <- 1 - results$score
      results_mean <- aggregate(score ~ h + delta , results, mean)
      names(results_mean)[3] <- "power"
      
      res <- rbind(res,results_mean)
      
      # Select the minimum h where power > 0.5
      # If no such h exists, select the minimum h with the maximum power
      min_h_power_gt_05 <- subset(results_mean, power >= 0.5)
      if (nrow(min_h_power_gt_05) > 0) {
         min_h <- min_h_power_gt_05$h[1]
         break
      } else {
         min_h <- NULL
      }
   }
   
   if(is.null(min_h)){
      results_mean <- results_mean[order(-results_mean$power, results_mean$h), ]
      min_h <- results_mean$h[1]
   }
   
   stopImplicitCluster()
   results <- list(h_sel = min_h, power = res)
   
   pl <- ggplot(res, aes(x = h, y = power)) +
      geom_line(aes(col = as.factor(delta)), linewidth = 0.9, alpha=.9) +
      labs(y="Power")+
      theme_minimal()+
      theme_light()+
      labs(color=expression(delta))+
      theme(legend.title = element_text(size=16),
            legend.text = element_text(size = 18),
            plot.title = element_text(size = 16),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.text.x = element_text(size = 11),
            axis.text.y = element_text(size = 11),
            strip.text = element_text(size = 14)) +
      scale_color_brewer(palette='Set1')

   results$power.plot <- pl
   
   if(power.plot){
      print(pl)
   }
   return(results)
}
