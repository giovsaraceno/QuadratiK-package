#' Select the value of the kernel tuning parameter
#'
#' This function computes the kernel bandwidth of the Gaussian kernel for the 
#' normality, two-sample and k-sample kernel-based quadratic distance (KBQD) 
#' tests.
#'
#' @param x Data set of observations from X.
#' @param y Numeric matrix or vector of data values. Depending on the input 
#'          \code{y}, the selection of h is performed for the corresponding 
#'          test.
#' \itemize{
#'    \item if \code{y} = NULL, the function performs the tests for normality on
#'          \code{x}.
#'    \item if \code{y} is a data matrix, with same dimensions of \code{x}, the 
#'          function performs the two-sample test between \code{x} and \code{y}.
#'    \item if \code{y} is a numeric or factor vector, indicating the group 
#'          memberships for each observation, the function performs the k-sample
#'          test.
#' }
#' @param alternative Family of alternative chosen for selecting h, between 
#'                    "location", "scale" and "skewness".
#' @param method The method used for critical value estimation 
#'               ("subsampling", "bootstrap", or "permutation").
#' @param b The size of the subsamples used in the subsampling algorithm .
#' @param B The number of iterations to use for critical value estimation, 
#'          B = 150 as default.
#' @param delta_dim Vector of coefficient of alternative with respect to each 
#'                  dimension
#' @param delta Vector of parameter values indicating chosen alternatives
#' @param h_values Values of the tuning parameter used for the selection
#' @param Nrep Number of bootstrap/permutation/subsampling replications.
#' @param n_cores Number of cores used to parallel the h selection algorithm. 
#'                If this is not provided, the function will detect the 
#'                available cores.
#' @param Quantile The quantile to use for critical value estimation, 0.95 is 
#'                 the default value.
#' @param power.plot Logical. If TRUE, it is displayed the plot of power for 
#'                   values in h_values and delta.
#' 
#' @return A list with the following attributes:
#' \itemize{
#'    \item \code{h_sel} the selected value of tuning parameter h;
#'    \item \code{power} matrix of power values computed for the considered 
#'                       values of \code{delta} and \code{h_values};
#'    \item \code{power.plot} power plots (if \code{power.plot} is \code{TRUE}).
#' }
#' 
#' @details
#' The function performs the selection of the optimal value for the tuning 
#' parameter \eqn{h} of the normal kernel function, for normality test, the 
#' two-sample and k-sample KBQD tests. It performs a small simulation study,
#' generating samples according to the family of \code{alternative} specified, 
#' for the chosen values of \code{h_values} and \code{delta}.
#' 
#' We consider target alternatives \eqn{F_\delta(\hat{\mathbf{\mu}},
#' \hat{\mathbf{\Sigma}}, \hat{\mathbf{\lambda}})}, where 
#' \eqn{\hat{\mathbf{\mu}}, \hat{\mathbf{\Sigma}}} and 
#' \eqn{\hat{\mathbf{\lambda}}} indicate the location,
#' covariance and skewness parameter estimates from the pooled sample. 
#' - Compute the estimates of the mean \eqn{\hat{\mu}}, covariance matrix
#'  \eqn{\hat{\Sigma}} and skewness \eqn{\hat{\lambda}} from the pooled sample.
#' - Choose the family of alternatives \eqn{F_\delta = F_\delta(\hat{\mu}
#' ,\hat{\Sigma}, \hat{\lambda})}. \cr \cr
#' *For each value of \eqn{\delta} and \eqn{h}:*
#' - Generate \eqn{\mathbf{X}_1,\ldots,\mathbf{X}_{k-1}  \sim F_0}, for 
#' \eqn{\delta=0};
#' - Generate \eqn{\mathbf{X}_k \sim F_\delta};
#' - Compute the \eqn{k}-sample test statistic between \eqn{\mathbf{X}_1, 
#' \mathbf{X}_2, \ldots, \mathbf{X}_k} with kernel parameter \eqn{h};
#' - Compute the power of the test. If it is greater than 0.5, 
#' select \eqn{h} as optimal value. 
#' - If an optimal value has not been selected, choose the \eqn{h} which
#'  corresponds to maximum power. 
#'
#' The available \code{alternative} are \cr
#' *location* alternatives, \eqn{F_\delta = 
#' SN_d(\hat{\mu} + \delta,\hat{\Sigma}, \hat{\lambda})},with 
#' \eqn{\delta = 0.2, 0.3, 0.4}; \cr
#' *scale* alternatives, 
#' \eqn{F_\delta = SN_d(\hat{\mu} ,\hat{\Sigma}*\delta, \hat{\lambda})}, 
#' \eqn{\delta = 0.1, 0.3, 0.5}; \cr
#' *skewness* alternatives, 
#' \eqn{F_\delta = SN_d(\hat{\mu} ,\hat{\Sigma}, \hat{\lambda} + \delta)}, 
#' with \eqn{\delta = 0.2, 0.3, 0.6}. \cr
#' The values of \eqn{h = 0.6, 1, 1.4, 1.8, 2.2} and \eqn{N=50} are set as 
#' default values. \cr
#' The function \code{select_h()} allows the user to 
#' set the values of \eqn{\delta} and \eqn{h} for a more extensive grid search. 
#' We suggest to set a more extensive grid search when computational resources 
#' permit.
#' 
#' @note Please be aware that the `select_h()` function may take a significant 
#' amount of time to run, especially with larger datasets or when using an 
#' larger number of parameters in \code{h_values} and \code{delta}. Consider 
#' this when applying the function to large or complex data.
#' 
#' @references
#' Markatou, M. and Saraceno, G. (2024). “A Unified Framework for
#' Multivariate Two- and k-Sample Kernel-based Quadratic Distance 
#' Goodness-of-Fit Tests.” \cr 
#' https://doi.org/10.48550/arXiv.2407.16374
#' 
#' Saraceno, G., Markatou, M., Mukhopadhyay, R. and Golzy, M. (2024). 
#' Goodness-of-Fit and Clustering of Spherical Data: the QuadratiK package 
#' in R and Python. \cr
#' https://arxiv.org/abs/2402.02290.
#' 
#' @seealso The function \code{select_h} is used in the [kb.test()] function.
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
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel detectCores
#' @importFrom parallel stopCluster
#' @import foreach 
#' @importFrom stats cov
#' @importFrom stats aggregate
#' @importFrom stats power
#' @import ggplot2
#' @import RcppEigen
#' @import rlecuyer
#' @importFrom Rcpp sourceCpp
#'
#' @useDynLib QuadratiK, .registration = TRUE
#'
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G2.0, G2.0a} input y, delta_dim, B, b
#' @srrstats {G2.1, G2.1a} input y, delta_dim, delta, h_values
#' @srrstats {G2.2, G2.4b} input y
#' @srrstats {G2.3, G2.3a, G2.3b} input alternative, method
#' @srrstats {G2.7,G2.8} different types of input are considered
#' @srrstats {G2.13,G2.14,G2.14a,G2.15,G2.16} error for NA, Nan, Inf, -Inf
#' 
#' @export
select_h <- function(x, y=NULL, alternative=NULL, method="subsampling", b=0.8, 
                     B=100, delta_dim=1, delta=NULL, h_values=NULL,Nrep=50, 
                     n_cores=2, Quantile=0.95, power.plot=TRUE) {
   
   # Convert vectors to a single column matrix
   if(!is.numeric(x) & !is.data.frame(x)){
      stop("x must be numeric")
   }
   # Convert vectors to a single column matrix
   if(is.vector(x)) {
      x <- matrix(x, ncol = 1)
   } else if(is.data.frame(x)) {
      x <- as.matrix(x)
   }
   
   if(any(is.na(x))){
      stop("There are missing values in x!")
   } else if(any(is.infinite(x) |is.nan(x))){
      stop("There are undefined values in x, that is Nan, Inf, -Inf")
   }
   
   if(!is.null(y)){
      if(is.vector(y) | is.factor(y)) {
         y <- matrix(as.numeric(y), ncol = 1)
      } else if(is.data.frame(y)) {
         y <- as.matrix(y)
      }
      if(any(is.na(y))){
         stop("There are missing values in y!")
      } else if(any(is.infinite(y) |is.nan(y))){
         stop("There are undefined values in y, that is Nan, Inf, -Inf")
      }
   }
   
   if(!(alternative %in% c("location", "scale", "skewness"))){
      
      stop("The alternative argument should be one of 'location', 'scale' or 
           'skewness'")
   }
   
   n <- nrow(x)
   d <- ncol(x)

   if (is.null(h_values)) {
      h_values <- seq(0.4, 3.3, 0.4)
   } else if(!(is.numeric(h_values) & is.vector(h_values))){
      stop("h_values must be a numeric vector")
   }
   
   if(is.null(delta)){
      if(alternative=='location'){
         delta <- c(0.2,0.3,0.4)
      } else if(alternative=='scale'){
         delta <- 1 + c(0.1,0.3,0.5)
      } else if(alternative=='skewness'){
         delta <- c(0.2, 0.3, 0.6)
      }
   } else if(!(is.numeric(delta) & is.vector(delta))){
      stop("delta must be a numeric vector.")
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
         stop("delta_dim must be 1 or a numeric vector of length equal to the 
              number of columns of pooled.")
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
         xnew <- rmsn(n, xi = mean_dat, Omega = S_dat, alpha = skew_data)
         ynew <- rmsn(m, xi = mean_dat + dk, Omega = S_dat, 
                          alpha = skew_data)
      } else if(alternative=='scale'){
         xnew <- rmsn(n, xi = mean_dat, Omega = S_dat, alpha = skew_data)
         ynew <- rmsn(m, xi = mean_dat , Omega = dk*S_dat, 
                          alpha = skew_data)
      } else if(alternative=='skewness'){
         xnew <- rmsn(n, xi = mean_dat, Omega = S_dat, alpha = skew_data)
         ynew <- rmsn(m, xi = mean_dat, Omega = S_dat, alpha = skew_data+dk)
      }
      
      #x <- as.matrix(xnew[sample(n, replace = FALSE),])
      #y <- as.matrix(ynew[sample(m, replace = FALSE),])
      
      STATISTIC <- stat2sample(xnew, ynew, h, rep(0,d),
                               diag(d),"Nonparam",compute_variance=FALSE)
      CV <- compute_CV(B, Quantile, pooled, n, m, h, method, b,
                       compute_variance=FALSE)
      cv <- CV$cv
      return(c(STATISTIC[1] < cv[1]))
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
         rmsn(nk*(K-1), xi = mean_dat, Omega = S_dat, alpha = skew_data),
         rmsn(nk, xi = mean_tilde, Omega = S_tilde, alpha = skew_tilde))
      ynew <- matrix(rep(1:K, each=nk),ncol=1)
      
      
      sizes_new <- as.vector(table(ynew))
      cum_size_new <- c(0,cumsum(sizes_new))
      STATISTIC <- stat_ksample_cpp(xnew, ynew, h, sizes_new,
                                    cum_size_new,compute_variance=FALSE)
      CV <- cv_ksample(xnew, ynew, h, B, b, Quantile, method,
                       compute_variance=FALSE)
      cv <- CV$cv
      return(c(STATISTIC[1] < cv[1]))
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
      
      STATISTIC <- kbNormTest(xnew, h, mean_dat, S_dat)
      CV <- normal_CV(d,n,h,mean_dat,S_dat,B,Quantile)
      
      return(c(STATISTIC[1] < CV))
   }
   
   # chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
   # if (nzchar(chk) && chk == "TRUE") {
   #    # nocov start
   #    num_cores <- 2
   #    # nocov end
   # } elseif(is.null(n_cores)) {
   #    num_cores <- detectCores()
   # } else 
   if (is.numeric(n_cores)) {
      num_cores <- as.numeric(n_cores)
   } else {
      stop("n_cores must be a numeric value")
   }
   
   cl <- parallel::makeCluster(num_cores)
   doParallel::registerDoParallel(cl)
   # parallel::clusterExport(cl, varlist = c("kbNormTest", "normal_CV",
   #                                       "cv_ksample", "stat_ksample_cpp",
   #                                       "compute_CV","stat2sample"))
   
   D <- length(delta)
   
   k_values <- 1:D
   rep_values <- 1:Nrep
   
   params <- expand.grid(Rep=rep_values, h = h_values)
   params <- split(params, seq(nrow(params)))
   
   res <- data.frame(delta=numeric(),
                     h=numeric(), power=numeric())
   
   i <- NULL
   for(k in k_values){
   
      results <- foreach::foreach(i = seq_along(params), 
                         .combine = rbind,
                         .packages=c("sn", "moments", "stats",
                                     "rlecuyer")) %dopar% {

         h <- as.numeric(params[[i]]$h)

         if(is.null(y)){
            
            objective_result <- objective_norm(h,k)
         
         } else if(d_y > 1){
            
            objective_result <- objective_2(h, k)
            
         } else if(d_y == 1){
            
            objective_result <- objective_k(h, k)
         }
         
         data.frame(delta=delta[k], 
                    h=h, power=objective_result)
      }
      
      results$power <- 1 - results$power
      results_mean <- aggregate(power ~ h + delta , data = results, FUN = mean) 
      
      res <- rbind(res,results_mean)
      
      # Select the minimum h where power > 0.5
      # If no such h exists, select the minimum h with the maximum power
      min_h_power_gt_05 <- subset(res, power >= 0.5)
      if (nrow(min_h_power_gt_05) > 0) {
         min_h <- min_h_power_gt_05$h[1]
         break
      } else {
         min_h <- NULL
      }
   }
   
   if(is.null(min_h)){
      res <- res[order(-res$power, res$h), ]
      min_h <- res$h[1]
   }

   parallel::stopCluster(cl)
   
   
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
   
   results <- list(h_sel = min_h, power = res, power.plot = pl)
   
   if(power.plot){
      print(pl)
   }
   return(results)
   
}
