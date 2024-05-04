#' Compute the critical value for two-sample KBQD tests
#'
#' This function computes the critical value for two-sample kernel tests with
#' centered Gaussian kernel
#' using one of three methods: bootstrap, permutation, or subsampling.
#'
#' @param B the number of bootstrap/permutation/subsampling samples to generate.
#' @param Quantile the quantile of the bootstrap/permutation/subsampling 
#'                 distribution to use as the critical value.
#' @param data_pool a matrix containing the data to be used in the test.
#' @param size_x the number of rows in the \code{data_pool} matrix 
#'               corresponding to group X.
#' @param size_y the number of rows in the \code{data_pool} matrix 
#'               corresponding to group Y.
#' @param h the tuning parameter for the kernel test.
#' @param method the method to use for computing the critical value 
#'               (one of "bootstrap", "permutation", or "subsampling").
#' @param b the subsampling block size (only used if \code{method} is 
#'          "subsampling").
#'
#' @return the critical value for the specified method and significance level.
#' 
#' @references
#' Markatou Marianthi, Saraceno Giovanni, Chen Yang (2023). “Two- and k-Sample 
#' Tests Based on Quadratic Distances.” Manuscript, (Department of 
#' Biostatistics, University at Buffalo)
#'
#'
#' @useDynLib QuadratiK
#' @importFrom stats quantile
#' @import RcppEigen
#'
#' @srrstats {G1.4a} roxigen2 is used
#' 
#' @keywords internal
compute_CV<-function(B, Quantile, data_pool, size_x, size_y, h, method, b=1){
   
   Results <- matrix(rep(0,2*B),ncol=2)
   for(i in 1:B) {
      
      if(method=="bootstrap"){
         Ind_x<-sample(1:(size_x+size_y), size_x, replace = T)
         Ind_y<-sample(1:(size_x+size_y), size_y, replace = T)
         data_x_star<-as.matrix(data_pool[Ind_x, ])
         data_y_star<-as.matrix(data_pool[Ind_y, ])
      } else if(method=="permutation"){
         Ind_x<-sample(1:(size_x+size_y), size_x, replace = F)
         Ind_y<-setdiff(1:(size_x+size_y), Ind_x)
         data_x_star<-as.matrix(data_pool[Ind_x, ])
         data_y_star<-as.matrix(data_pool[Ind_y, ])
      } else if(method=="subsampling"){
         Ind_x<-sample(1:size_x, round(size_x*b), replace = F)
         Ind_y<-sample((size_x+1):(size_x+size_y), round(size_y*b), replace = F)
         newind <- c(Ind_x,Ind_y)
         newsample <- sample(newind,length(newind),replace = F)
         m <- length(newsample)
         data_x_star<-as.matrix(data_pool[newsample[1:round(size_x*b)], ])
         data_y_star<-as.matrix(data_pool[newsample[(round(size_x*b)+1):m], ])
      }
      
      Results[i,] <-   stat2sample(data_x_star,data_y_star,h,
                                  rep(0,ncol(data_pool)),diag(1),"Nonparam")[1:2]
   }
   
   cv_res <- apply(Results,2,function(x) as.numeric(quantile(x,Quantile,na.rm=T)))
   return(list(cv=cv_res))
}
#' 
#' Compute the critical value for the Poisson KBQD tests for Uniformity
#' 
#' This function computes the empirical critical value for the U-statistics for 
#' testing uniformity on the sphere based on the centered poisson kernel.  
#'
#' @param d the dimension of generated samples.
#' @param size the number of observations to be generated.
#' @param rho the concentration parameter for the Poisson kernel.
#' @param B the number of replications.
#' @param Quantile the quantile of the distribution use to select the critical 
#'                 value.
#'
#' @return the critical value for the specified dimension, size and level.
#'
#' @details 
#' For each replication, a sample of d-dimensional observations from the uniform
#' #' distribution on the Sphere are generated and the Poisson kernel-based 
#' U-statistic is computed. After B iterations, the critical value is selected 
#' as the \code{Quantile} of the empirical distribution of the computed test 
#' statistics.
#' 
#' @references
#' Ding Yuxin, Markatou Marianthi, Saraceno Giovanni (2023). “Poisson 
#' Kernel-Based Tests for Uniformity on the d-Dimensional Sphere.” 
#' Statistica Sinica. doi: doi:10.5705/ss.202022.0347
#' 
#' @useDynLib QuadratiK
#' @importFrom stats quantile
#' @import RcppEigen
#'
#' @srrstats {G1.4a} roxigen2 is used
#'  
#' @keywords internal
poisson_CV<-function(d, size, rho, B, Quantile){
   
   Results <- rep(0,B)
   for(i in 1:B) {
      
      dat <- sample_hypersphere(d, size)
      pk <- statPoissonUnif(dat,rho)
      
      Results[i] <- pk[1]
      
   }
   return(as.numeric(quantile(Results, Quantile)))
   
}
#' 
#' Compute the critical value for the KBQD tests for multivariate Normality 
#' 
#' This function computes the empirical critical value for the Normality test
#' based on the KBQD tests using the centered Gaussian kernel.
#'
#' @param d the dimension of generated samples.
#' @param size the number of observations to be generated.
#' @param h the concentration parameter for the Gaussian kernel.
#' @param mu_hat Mean vector for the reference distribution.
#' @param Sigma_hat Covariance matrix of the reference distribution.
#' @param B the number of replications.
#' @param Quantile the quantile of the distribution use to select the critical 
#'                 value
#'
#' @return the critical value for the specified dimension, size and level.
#' 
#' @details 
#' For each replication, a sample from the d-dimensional Normal distribution 
#' with mean vector \code{mu_hat} and covariance matrix \code{Sigma_hat} is 
#' generated and the KBQD test U-statistic for Normality is computed. 
#' After B iterations, the critical value is selected as the \code{Quantile} 
#' of the empirical distribution of the computed test statistics.
#' 
#' @useDynLib QuadratiK
#' @import mvtnorm
#' @import RcppEigen
#'
#' @srrstats {G1.4a} roxigen2 is used
#' 
#' @keywords internal
normal_CV<-function(d, size, h, mu_hat, Sigma_hat, B = 150, Quantile=0.95){
   
   Results <- rep(0,B)
   for(i in 1:B) {
      
      dat <- rmvnorm(n = size, mean=mu_hat, sigma=Sigma_hat)
      dat <- matrix(dat, ncol = d, byrow = TRUE)
      kbnorm <- kbNormTest(dat, h, mu_hat, Sigma_hat)
      Results[i] <- kbnorm[1]
      
   }
   return(as.numeric(quantile(Results, Quantile)))
   
}
#' Compute the critical value for the KBQD k-sample tests
#'
#' This function computes the empirical critical value for the k-sample KBQD 
#' tests using the centered Gaussian kernel, with bootstrap, permutation, or 
#' subsampling.
#'
#' @param x matrix containing the observations to be used in the k-sample test
#' @param y vector indicating the sample for each observation
#' @param h the tuning parameter for the test using the Gaussian kernel
#' @param B the number of bootstrap/permutation/subsampling samples to generate
#' @param b the subsampling block size (only used if \code{method} is 
#'          "subsampling")
#' @param Quantile the quantile of the bootstrap/permutation/subsampling 
#'                 distribution to use as the critical value
#' @param method the method to use for computing the critical value 
#'               (one of "bootstrap", "permutation")
#'
#' @return a vector of two critical values corresponding to different 
#' formulation of the k-sample test statistics.
#'
#' @useDynLib QuadratiK
#' @importFrom stats quantile
#' @import RcppEigen
#'
#' @srrstats {G1.4a} roxigen2 is used
#' 
#' @keywords internal
cv_ksample <- function(x, y, h, B=150, b=0.9, Quantile =0.95, 
                       method="subsampling"){
   
   sizes <- as.vector(table(y))
   cum_size <- c(0,cumsum(sizes))
   K <- length(unique(y))
   n <- nrow(x)
   d <- ncol(x)
   
   Results <- matrix(rep(0,2*B),ncol=2)
   
   for(i in 1:B) {
      
      if(method=="bootstrap"){
         
         ind_k <- unlist(lapply(1:K, function(k) sample(1:sizes[k], sizes[k], 
                                                      replace=T) + cum_size[k]))
         ind_k <- sample(ind_k,length(ind_k),replace = F)
         
      } else if(method=="permutation"){
         
         ind_k <- sample(1:n, n, replace=F)
         
      } else if(method=="subsampling"){
         
         ind_k <- unlist(lapply(1:K, function(j) {
            sample(1:sizes[j], round(sizes[j]*b), replace=F) + cum_size[j]}))
         ind_k <- sample(ind_k,length(ind_k),replace = F)
      }
      
      data_k <- as.matrix(x[ind_k,])
      y_ind <- y[ind_k,]
      
      sizes_sub <- as.vector(table(y_ind))
      cum_size_sub <- c(0,cumsum(sizes_sub))
      
      Results[i,] <-   stat_ksample_cpp(as.matrix(data_k), c(y_ind), h, 
                                        sizes_sub, cum_size_sub)[1:2]
   }
   
   cv_k <- apply(Results,2,function(x) as.numeric(quantile(x,Quantile,na.rm=T)))
   return(list(cv=cv_k))
}
