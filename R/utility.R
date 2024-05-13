#'
#' Degrees of freedom (DOF) for the Poisson kernel
#' 
#' Compute the Degrees of Freedom (DOF) of the Poisson Kernel given the
#' dimension d and concentration parameter rho
#'
#' @param d the number of dimensions
#' @param rho concentration parameter
#'
#' @return a list containing the DOF and the coefficient c of the asymptotic
#' distribution
#' 
#' @srrstats {G1.4a} roxigen2 is used
#' @keywords internal
DOF <- function(d, rho){
   num_c <- (1+rho^2)/((1-rho^2)^(d-1)) -1
   den_c <- (1+rho)/((1-rho)^(d-1)) -1
   DOF <- (den_c^2)/num_c
   c <- num_c/den_c
   result <- list("DOF"=DOF,"Coefficient"=c)
   return(result)
}
#'
#' Degrees of freedom (DOF) for the Normal kernel
#' 
#' Compute the Degrees of Freedom (DOF) of the normal Kernel centered with 
#' respect to the standard normal distribution, given the dimension d and the
#' bandwidth parameter h.
#'
#' @param Sigma_h covariance matrix of the gaussian kernel 
#' @param V Covariance matrix of the tested distribution G
#'
#' @return a list containing the DOF and the coefficient c of the asymptotic
#' distribution
#' 
#' @srrstats {G1.4a} roxigen2 is used
#' @keywords internal
DOF_norm <- function(Sigma_h, V){
   
   num_dof <- det(Sigma_h)^(-1/2) - det(Sigma_h + 2*V)^(-1/2)
   den_dof <- det(Sigma_h)^(-1/2) * det(Sigma_h + 4*V)^(-1/2) 
               -2 * det(Sigma_h + V)^(-1/2) * det(Sigma_h + 3*V)^(-1/2)
               + det(Sigma_h + 2*V)^(-1/2)
   
   dof <- num_dof^2/den_dof
   
   const <- den_dof/num_dof
   
   # h2 <- h^2
   # dof_num <- 1 - (h2/(h2+2))^(d/2)
   # dof_den <- ((h2+2)/(h2+4))^(d/2) - 2* (h2/(h2+1))^(d/2)*(h2/(h2+3))^(d/2) +
   #    (h2/(h2+2))^(d)
   # dof <- dof_num^2/dof_den
   # 
   # c_num <- ((h2+2)/(h2*(h2+4)))^(d/2) - 2* (h2/(h2+1))^(d/2)*(1/(h2+3))^(d/2) +
   #    (h/(h2+2))^(d)
   # const <- c_num/dof_num
   # 
   result <- list("DOF"=dof,"Coefficient"=const)
   return(result)
}
#'
#' Exact variance of normality test 
#' 
#' Compute the exact variance of kernel test for normality under the null 
#' hypothesis that G=N(0,I).
#'
#' @param Sigma_h covariance matrix of the gaussian kernel 
#' @param V Covariance matrix of the tested distribution G
#' @param n sample size
#'
#' @return the value of computed variance
#' 
#' @srrstats {G1.4a} roxigen2 is used
#' @keywords internal
var_norm <- function(Sigma_h, V, n){
   
   d <- nrow(Sigma_h)
   
   res <- det(Sigma_h)^(-1/2) * det(Sigma_h + 4*V)^(-1/2) -2 * det(Sigma_h + V)^(-1/2) * det(Sigma_h + 3*V)^(-1/2) + det(Sigma_h + 2*V)^(-1)
   
   res <- 2/(n*(n-1)) * 1/(2*pi)^(d) * res
   
   # h2 <- h^2
   # h_const <- ((h2+2)/(h2^2*(h2+4)))^(d/2) -2/((h2+1)*(h2+3))^(d/2) + (h2+2)^(-d)
   # res <- 2/(n*(n-1)) * 1/(2*pi)^(d)  * h_const
   # 
   return(res)
}
#'
#' Generate random sample from the hypersphere
#'
#' Generate random sample from the uniform distribution on the hypersphere
#' 
#' @param d Number of dimensions.
#' @param n_points Number of sampled observations.
#'
#' @return Data matrix with the sampled observations.
#'
#' @examples
#' x_sp <- sample_hypersphere(3,100)
#'
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @export
sample_hypersphere <- function(d, n_points=1) {
   z <- matrix(rnorm(n_points * d), n_points, d)
   norm <- sqrt(rowSums(z^2))
   return(sweep(z, 1, norm, "/"))
}
#' Generate two samples data from skew-normal distributions
#'
#' This function generates data from skew-normal distributions with the 
#' specified parameters of means and covariance matrices.
#'
#' @param d number of dimensions.
#' @param size_x the number of observations for sample X
#' @param size_y the number of observations for sample Y
#' @param mu_x the mean of X
#' @param mu_y the mean of Y
#' @param sigma_x the standard deviation of X
#' @param sigma_y the standard deviation of Y
#' @param skewness_y the skewness of Y (the skewness of X is set to zero).
#'
#' @return a list containing the generated X and Y data sets.
#'
#' @importFrom sn rmsn
#'
#' @useDynLib QuadratiK
#' 
#' @srrstats {G1.4a} roxigen2 is used 
#' @keywords internal
generate_SN<-function(d, size_x, size_y, mu_x, mu_y, 
                      sigma_x, sigma_y, skewness_y){
   
   muX <- rep(mu_x, d)
   SigmaX <- diag(rep(sigma_x,d))
   data_x <- rmsn(n=size_x, xi=muX, Omega = SigmaX, alpha=rep(0,d))
   
   muY <- rep(mu_y, d)
   SigmaY <- diag(rep(sigma_y,d))
   alphaY <- rep(skewness_y, d)
   data_y <- rmsn(n=size_y, xi=muY, Omega = SigmaY, alpha=alphaY)
   
   colnames(data_x)<-NULL; colnames(data_y)<-NULL
   return(list("X"=data_x, "Y"=data_y))
}
#'
#' QQ-plot of given two samples using ggplot2
#'
#' @param sample1 matrix of observations in sample 1
#' @param sample2 matrix of observations in sample 2
#' @param main_title title of the generated plot
#'
#' @return QQ-plot of given samples
#'
#' @import ggplot2
#' 
#' @srrstats {G1.4a} roxigen2 is used
#' 
#' @keywords internal
compare_qq <- function(sample1, sample2, main_title) {
   # Compute quantiles
   quantiles1 <- quantile(sample1, 
                          probs = seq(0, 1, length.out = length(sample1)))
   quantiles2 <- quantile(sample2, 
                          probs = seq(0, 1, length.out = length(sample2)))
   df <- data.frame(q1 = quantiles1, q2 = quantiles2)
   
   with(df, {pl <- ggplot(df, aes(x=q1,y=q2))+
      geom_abline(slope=1, col="red") +
      geom_line(col="blue", linewidth=0.9)+
      ggtitle(main_title) +
      theme_minimal()+
      theme(legend.title = element_text(size=14),
            legend.text = element_text(size = 16),
            plot.title = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text.x = element_text(size = 11),
            axis.text.y = element_text(size = 11),
            strip.text = element_text(size = 14)) +
      scale_color_brewer(palette='Set1')
   
      return(pl)})
}
#'
#' Compute and display some descriptive statistics for the two sample tests
#'
#' @param var1 vector of observations of a given variable from sample 1
#' @param var2 vector of observations of a given variable from sample 2
#' @param var_name Name of the variable displayed
#' @param eps precision of displayed statistics
#'
#' @import ggplot2
#'
#' @return Computed statistics with a plot
#' 
#' @srrstats {G1.4a} roxigen2 is used
#' 
#' @keywords internal
compute_stats <- function(var1, var2, var_name,eps=3) {
   
   overall <- c(var1, var2)
   stats <- data.frame(matrix(c(mean(var1),mean(var2),mean(overall),sd(var1),
                                sd(var2),sd(overall),median(var1),median(var2),
                                median(overall),IQR(var1),IQR(var2),
                                IQR(overall),min(var1),min(var2),min(overall),
                                max(var1),max(var2),max(overall)),
                                 nrow=6,ncol=3,byrow=TRUE))
   colnames(stats) <- c("Group 1", "Group 2", "Overall")
   rownames(stats) <- c("mean", "sd", "median", "IQR", "min", "max")
   
   pl <- ggplot() +
      ggpp::annotate('table', x = 0.5, y = 0.5, 
                     label = data.frame(Stat = rownames(stats),stats), 
                     hjust = 0.5, vjust = 0.5) +
      theme_void() +
      ggtitle(paste(var_name))+
      scale_color_brewer(palette='Set1')
   
   
   return(list(plots=pl,stats=stats))
}
