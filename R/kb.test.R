#'
#' @title 
#' Kernel-based quadratic distance (KBQD) Goodness-of-Fit tests
#'
#' @description
#' This function performs the kernel-based quadratic distance goodness-of-fit 
#' tests. It includes tests for multivariate normality, two-sample tests and
#' \eqn{k}-sample tests.
#' 
#'
#' @param x Numeric matrix or vector of data values.
#' @param y Numeric matrix or vector of data values. Depending on the input 
#'          \code{y}, the corresponding test is performed.
#' \itemize{
#'    \item if \code{y} = NULL, the function performs the tests for normality on
#'          \code{x}
#'    \item if \code{y} is a data matrix, with same dimensions of \code{x}, the 
#'          function performs the two-sample test between \code{x} and \code{y}.
#'    \item if \code{y} is a numeric or factor vector, indicating the group 
#'          memberships for each observation, the function performs the k-sample
#'          test.
#' }
#' @param h Bandwidth for the kernel function. If a value is not provided, the 
#'          algorithm for the selection of an optimal h is performed 
#'          automatically. See the function \code{\link{select_h}} for more 
#'          details.
#' @param method The method used for critical value estimation ("subsampling", 
#'               "bootstrap", or "permutation")(default: "subsampling").
#' @param B The number of iterations to use for critical value estimation 
#'          (default: 150).
#' @param b The size of the subsamples used in the subsampling algorithm  
#'          (default: 0.8).
#' @param Quantile The quantile to use for critical value estimation, 0.95 is 
#'                 the default value.
#' @param mu_hat Mean vector for the reference distribution.
#'                
#' @param Sigma_hat Covariance matrix of the reference distribution.
#' @param centeringType String indicating the method used for centering the 
#'                      normal kernel ('Param' or 'Nonparam').
#' @param K_threshold maximum number of groups allowed. Default is 10. It is a 
#'                    control parameter. Change in case of more than 10 samples.
#' @param alternative Family of alternative chosen for selecting h, between 
#'                    "location", "scale" and "skewness" (only if \code{h} 
#'                    is not provided).
#'
#' @details 
#' The function \code{kb.test} performs the kernel-based quadratic
#' distance tests using the Gaussian kernel with bandwidth parameter \code{h}.
#' Depending on the shape of the input \code{y} the function performs the tests
#' of multivariate normality, the non-parametric two-sample tests or the 
#' k-sample tests.
#' 
#' The quadratic distance between two probability distributions \eqn{F} and
#' \eqn{G} is
#' defined as \deqn{d_{K}(F,G)=\iint K(x,y)d(F-G)(x)d(F-G)(y),}
#' where \eqn{G} is a distribution whose goodness of fit we wish to assess and 
#' \eqn{K} denotes the Normal kernel defined as 
#' \deqn{ K_{{h}}(\mathbf{s}, \mathbf{t}) = (2 \pi)^{-d/2} 
#' \left(\det{\mathbf{\Sigma}_h}\right)^{-\frac{1}{2}}  
#' \exp\left\{-\frac{1}{2}(\mathbf{s} - \mathbf{t})^\top 
#' \mathbf{\Sigma}_h^{-1}(\mathbf{s} - \mathbf{t})\right\},}
#' for every \eqn{\mathbf{s}, \mathbf{t} \in \mathbb{R}^d \times 
#' \mathbb{R}^d}, with covariance matrix \eqn{\mathbf{\Sigma}_h=h^2 I} and
#' tuning parameter \eqn{h}. \cr
#' \itemize{
#'    \item **Test for Normality**: \cr
#'    Let \eqn{x_1, x_2, ..., x_n} be a random sample with empirical 
#'    distribution function \eqn{\hat F}. We test the null hypothesis of 
#'    normality, i.e. \eqn{H_0:F=G=\mathcal{N}_d(\mu, \Sigma)}. 
#'    
#'    We consider the U-statistic estimate of the sample KBQD
#'    \deqn{U_{n}=\frac{1}{n(n-1)}\sum_{i=2}^{n}\sum_{j=1}^{i-1}
#'    K_{cen}(\mathbf{x}_{i}, \mathbf{x}_{j}),}
#'    then the first test statistics is 
#'    \deqn{T_{n}=\frac{U_{n}}{\sqrt{Var(U_{n})}},}
#'    with \eqn{Var(U_n)} computed exactly following Lindsay et al.(2014),
#'    and the V-statistic estimate \deqn{V_{n} = \frac{1}{n}\sum_{i=1}^{n}
#'    \sum_{j=1}^{n}K_{cen}(\mathbf{x}_{i}, \mathbf{x}_{j}),}
#'    where \eqn{K_{cen}} denotes the Normal kernel \eqn{K_h} with parametric 
#'    centering with respect to the considered normal distribution 
#'    \eqn{G = \mathcal{N}_d(\mu, \Sigma)}.
#'    
#'    The asymptotic distribution of the V-statistic is an infinite combination
#'    of weighted independent chi-squared random variables with one degree of 
#'    freedom. The cutoff value is obtained using the Satterthwaite 
#'    approximation \eqn{c \cdot \chi_{DOF}^2}, where \eqn{c} and \eqn{DOF} 
#'    are computed exactly following the formulas in Lindsay et al.(2014).
#' 
#'    For the \eqn{U}-statistic the cutoff is determined empirically:
#'    -  Generate data from the considered normal distribution ;
#'    - Compute the test statistics for \code{B} Monte Carlo(MC) replications;
#'    - Compute the 95th quantile of the empirical distribution of the test
#'      statistic.
#'     
#'    \item **k-sample test**: \cr
#'    Consider \eqn{k} random samples of i.i.d. observations 
#'    \eqn{\mathbf{x}^{(i)}_1,
#'    \mathbf{x}^{(i)}_{2},\ldots, \mathbf{x}^{(i)}_{n_i} \sim F_i}, 
#'    \eqn{i = 1, \ldots, k}.  
#'    We test if the samples are generated from the same *unknown* distribution,
#'    that is \eqn{H_0: F_1 = F_2 = \ldots = F_k} versus 
#'    \eqn{H_1: F_i \not = F_j}, for some \eqn{1 \le i \not = j \le k}. \cr
#'    We construct a matrix distance \eqn{\hat{\mathbf{D}}}, with 
#'    off-diagonal elements 
#'    \deqn{\hat{D}_{ij} = \frac{1}{n_i n_j} \sum_{\ell=1}^{n_i}
#'    \sum_{r=1}^{n_j}K_{\bar{F}}(\mathbf{x}^{(i)}_\ell,\mathbf{x}^{(j)}_r), 
#'    \qquad \mbox{ for }i \not= j}
#'    and in the diagonal
#'    \deqn{\hat{D}_{ii} = \frac{1}{n_i (n_i -1)} \sum_{\ell=1}^{n_i}
#'    \sum_{r\not= \ell}^{n_i} K_{\bar{F}}(\mathbf{x}^{(i)}_\ell,
#'    \mathbf{x}^{(i)}_r), \qquad \mbox{ for }i = j,}
#'    where \eqn{K_{\bar{F}}} denotes the Normal kernel \eqn{K_h}
#'    centered non-parametrically with respect to 
#'    \deqn{\bar{F} = \frac{n_1 \hat{F}_1 + \ldots + n_k \hat{F}_k}{n}, 
#'    \quad \mbox{ with } n=\sum_{i=1}^k n_i.}
#'    
#'    We compute the trace statistic
#'    \deqn{\mathrm{trace}(\hat{\mathbf{D}}_n) =  \sum_{i=1}^{k}\hat{D}_{ii}}
#'    and \eqn{D_n}, derived considering all the possible pairwise comparisons 
#'    in the *k*-sample null hypothesis, given as
#'    \deqn{D_n = (k-1) \mathrm{trace}(\hat{\mathbf{D}}_n) 
#'    - 2 \sum_{i=1}^{k}\sum_{j> i}^{k}\hat{D}_{ij}.}
#'    
#'    We compute the empirical critical value by employing numerical techniques 
#'    such as the bootstrap, permutation and subsampling algorithms:
#'    - Generate k-tuples, of total size \eqn{n_B}, from the pooled sample 
#'    following one of the sampling methods;
#'    - Compute the k-sample test statistic;
#'    - Repeat \code{B} times;
#'    - Select the \eqn{95^{th}} quantile of the obtained values. 
#'    
#'    \item **Two-sample test**: \cr
#'    Let \eqn{x_1, x_2, ..., x_{n_1} \sim F} and 
#'    \eqn{y_1, y_2, ..., y_{n_2} \sim G} be
#'    random samples from the distributions \eqn{F} and \eqn{G}, respectively.
#'    We test the null hypothesis that the two samples are generated from 
#'    the same *unknown* distribution, that is \eqn{H_0: F=G} vs 
#'    \eqn{H_1:F\not=G}. The test statistics coincide with the \eqn{k}-sample 
#'    test statistics when \eqn{k=2}.
#' }
#' 
#' ### Kernel centering
#' The arguments \code{mu_hat} and \code{Sigma_hat} indicate the normal model 
#' considered for the normality test, that is \eqn{H_0: F = N(}\code{mu_hat},
#' \code{Sigma_hat}). 
#' For the two-sample and \eqn{k}-sample tests, \code{mu_hat} and 
#' \code{Sigma_hat} can 
#' be used for the parametric centering of the kernel, in the case we want to 
#' specify the reference distribution, with \code{centeringType = "Param"}. 
#' This is the default method when the test for normality is performed.
#' The normal kernel centered with respect to 
#'  \eqn{G \sim N_d(\mathbf{\mu}, \mathbf{V})} can be computed as
#'  \deqn{K_{cen(G)}(\mathbf{s}, \mathbf{t}) = 
#'  K_{\mathbf{\Sigma_h}}(\mathbf{s}, \mathbf{t}) - 
#'  K_{\mathbf{\Sigma_h} + \mathbf{V}}(\mathbf{\mu}, \mathbf{t}) 
#' -  K_{\mathbf{\Sigma_h} + \mathbf{V}}(\mathbf{s}, \mathbf{\mu}) +
#'  K_{\mathbf{\Sigma_h} + 2\mathbf{V}}(\mathbf{\mu}, \mathbf{\mu}).} 
#' We consider the non-parametric centering of the kernel with respect to 
#' \eqn{\bar{F}=(n_1 F_1 + \ldots n_k F_k)/n} where \eqn{n=\sum_{i=1}^k n_i}, 
#' with \code{centeringType = "Nonparam"}, for the two- and \eqn{k}-sample 
#' tests. 
#' Let \eqn{\mathbf{z}_1,\ldots, \mathbf{z}_n} denote the pooled sample. For any
#'  \eqn{s,t \in \{\mathbf{z}_1,\ldots, \mathbf{z}_n\}}, it is given by 
#' \deqn{K_{cen(\bar{F})}(\mathbf{s},\mathbf{t}) =    K(\mathbf{s},\mathbf{t}) -
#'  \frac{1}{n}\sum_{i=1}^{n} K(\mathbf{s},\mathbf{z}_i) - 
#'  \frac{1}{n}\sum_{i=1}^{n} K(\mathbf{z}_i,\mathbf{t}) + 
#'  \frac{1}{n(n-1)}\sum_{i=1}^{n} \sum_{j \not=i}^{n} 
#'  K(\mathbf{z}_i,\mathbf{z}_j).}
#'  
#'
#' @seealso \linkS4class{kb.test} for the class definition.
#' 
#' @return An S4 object of class \code{kb.test} containing the results of the 
#' kernel-based quadratic distance tests, based on the normal kernel. The object
#' contains the following slots:
#' \itemize{
#'   \item \code{method}: Description of the kernel-based quadratic 
#'   distance test performed.
#'   \item \code{x} Data list of samples X (and Y).
#'   \item \code{Un} The value of the U-statistic.
#'   \item \code{H0_Un} A logical value indicating whether or not the null 
#'   hypothesis is rejected according to Un.
#'   \item \code{CV_Un} The critical value computed for the test Un.
#'   \item \code{Vn} The value of the V-statistic (if available).
#'   \item \code{H0_Vn} A logical value indicating whether or not the null 
#'   hypothesis is rejected according to Vn (if available).
#'   \item \code{CV_Vn} The critical value computed for the test Vn 
#'   (if available).
#'   \item \code{h} List with the value of bandwidth parameter used for the 
#'   normal kernel function. If \code{select_h} is used, the matrix of computed 
#'   power values and the corresponding power plot are also provided. 
#'   \item \code{B} Number of bootstrap/permutation/subsampling replications.
#'   \item \code{var_Un} exact variance of the kernel-based U-statistic.
#'   \item \code{cv_method} The method used to estimate the critical value 
#'   (one of "subsampling", "permutation" or "bootstrap").
#'   
#' }
#' 
#' @note
#' For the two- and \eqn{k}-sample tests, the slots \code{Vn}, \code{H0_Vn} and 
#' \code{CV_Vn} are empty, while the computed statistics are both reported in
#' slots \code{Un}, \code{H0_Un} and \code{CV_Un}.
#' 
#' A U-statistic is a type of statistic that is used to estimate a population
#' parameter. It is based on the idea of averaging over all possible *distinct*
#' combinations of a fixed size from a sample. 
#' A V-statistic considers all possible tuples of a certain size, not just
#' distinct combinations and can be used in contexts where unbiasedness is not
#' required.
#'
#' @references
#' Markatou, M. and Saraceno, G. (2024). “A Unified Framework for 
#' Multivariate Two- and k-Sample Kernel-based Quadratic Distance 
#' Goodness-of-Fit Tests.” \cr
#' https://doi.org/10.48550/arXiv.2407.16374
#'
#' Lindsay, B.G., Markatou, M. and Ray, S. (2014) "Kernels, Degrees of Freedom,
#' and Power Properties of Quadratic Distance Goodness-of-Fit Tests", Journal 
#' of the American Statistical Association, 109:505, 395-410, 
#' DOI: 10.1080/01621459.2013.836972
#'
#'
#' @examples
#' # create a kb.test object
#' x <- matrix(rnorm(100),ncol=2)
#' y <- matrix(rnorm(100),ncol=2)
#' 
#' # Normality test
#' my_test <- kb.test(x, h=0.5)
#' my_test
#' 
#' # Two-sample test
#' my_test <- kb.test(x,y,h=0.5, method="subsampling",b=0.9,
#'                      centeringType = "Nonparam")
#' my_test
#' 
#' # k-sample test
#' z <- matrix(rnorm(100,2),ncol=2)
#' dat <- rbind(x,y,z)
#' group <- rep(c(1,2,3),each=50)
#' my_test <- kb.test(x=dat,y=group,h=0.5, method="subsampling",b=0.9)
#' my_test
#'
#' @import Rcpp
#' @import RcppEigen
#'
#' @useDynLib QuadratiK
#' @importFrom Rcpp sourceCpp
#'
#' @srrstats {G1.0} Reference section reports the related literature
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G2.0a,G2.1,G2.3b} The code considers the different inputs
#' 
#' @export
setGeneric("kb.test",function(x, y=NULL, h = NULL, method = "subsampling", 
                              B = 150, b = NULL, Quantile = 0.95, 
                              mu_hat = NULL, Sigma_hat = NULL, 
                              centeringType="Nonparam", 
                              K_threshold=10, alternative="skewness"){
   
   standardGeneric("kb.test")
})
#' @rdname kb.test
#' 
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G2.0,G2.3a, G2.4b} The code considers the different inputs
#' @srrstats {G2.7,G2.8} different types of input are considered
#' @srrstats {G2.13,G2.14,G2.14a,G2.15,G2.16} error for NA, Nan, Inf, -Inf
#' 
#' @export
setMethod("kb.test", signature(x = "ANY"),
          function(x, y=NULL, h=NULL, method="subsampling", B = 150,
                   b = 0.9, Quantile = 0.95, mu_hat = NULL,
                   Sigma_hat = NULL, centeringType="Nonparam",
                   K_threshold=10, alternative="skewness"){
             
             
             if(!(method%in%c("bootstrap", "permutation", "subsampling"))){
                stop("method must be one of 'bootstrap', 'permutation' or 
                     'subsampling'")
             }
             if(b<=0 | b>1){
                stop("b indicates the proportion used for the subsamples in the
                     subsampling algoritm. It must be in (0,1].")
             }
             
             if(!(alternative%in%c("skewness", "location", "scale"))){
               stop("The algorithm for selecting the value of h can be performed
                    with respect to the following families of alternatives: 
                    'location', 'scale' or 'skewness'")
             }
             
             if(!(centeringType%in%c("Param", "Nonparam"))){
                stop("centering must be chosen between 'Param' and 'Nonparam'")
             }
             
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
             
             size_x <- nrow(x)
             k <- ncol(x)
             
             if (is.null(h) & is.null(y)){
                
                #stop("A value of the tuning parameter h must be provided to 
                #perform the kernel-based quadratic distance Normality tests")
                h_best <- select_h(x=x, alternative=alternative, method=method,
                                   b=b, B=B, power.plot=FALSE)
                h <- h_best$h_sel
             
             } else if (is.null(h)& !(is.null(y))){
                
                h_best <- select_h(x=x, y=y, alternative=alternative, 
                                   method=method, b=b, B=B, power.plot=FALSE)
                h <- h_best$h_sel

             } else {
                h_best <- list(h_sel=h)
             }
             
             if(is.null(y)){
                
                
                METHOD <- "Kernel-based quadratic distance Normality test"
                
                # Compute the estimates of mean and covariance from the data
                if(is.null(mu_hat)){
                   mu_hat <- rep(0,k)
                } else {
                   x <- x - mu_hat
                   mu_hat <- rep(0,k)
                }
                if(is.null(Sigma_hat)){
                   Sigma_hat <- diag(k)
                }
                
                STATISTIC <- kbNormTest(x, h, mu_hat, Sigma_hat)
                CV_Un <- normal_CV(k, size_x, h, mu_hat, Sigma_hat, B, Quantile)
                Sigma_h <- h^2*diag(k)
                dof <- DOF_norm(Sigma_h, Sigma_hat)
                qu_q <- qchisq(Quantile,df=dof$DOF)
                CV_Vn <- dof$Coefficient*qu_q
                
                var_Un <- var_norm(Sigma_h, Sigma_hat, size_x)
                CV_Un <- CV_Un/sqrt(var_Un)
                
                H0_Un <- (STATISTIC[1]/sqrt(var_Un) > CV_Un)
                
                res <- new("kb.test", Un = STATISTIC[1]/sqrt(var_Un), 
                           Vn = STATISTIC[2], CV_Un = CV_Un, CV_Vn = CV_Vn, 
                           H0_Un = H0_Un, H0_Vn = STATISTIC[2] > CV_Vn, 
                           method = METHOD, data = list(x = x), 
                           B= B, h= h_best, var_Un = var_Un)
                
             } else {
                
                K <- length(unique(y))
                
                if(K > K_threshold){ 
                   #Here we consider a maximum number of groups. However this 
                   #threshold can be set to a higher value by the user if needed
                   
                   # Check that they have the same number of columns (features):
                   if(!is.null(y) && ncol(x) != ncol(y)) {
                      stop("'x' and 'y' must have the same number of columns.")
                   }
                   
                   size_y <- nrow(y)
                   data_pool <- rbind(x, y)
                   
                   
                   METHOD <- "Kernel-based quadratic distance two-sample test"
                   
                   if(centeringType == "Param"){
                      
                      # Compute the estimates of mean and covariance from the 
                      # data
                      if(is.null(mu_hat)){
                         mu_hat <- colMeans(data_pool)
                      }
                      if(is.null(Sigma_hat)){
                         Sigma_hat <- cov(data_pool)
                      }
                      
                      STATISTIC <- stat2sample(x, y, h, mu_hat, 
                                               Sigma_hat, "Param",
                                               compute_variance = TRUE)
                      
                   } else if(centeringType == "Nonparam"){
                      
                      STATISTIC <- stat2sample(x, y, h, rep(0,k),
                                               diag(k),"Nonparam",
                                               compute_variance = TRUE)
                   }
                   
                   CV <- compute_CV(B, Quantile, data_pool, size_x, size_y, h, 
                                    method, b, compute_variance = TRUE)
                   STATISTIC[1] <- STATISTIC[1]/sqrt(STATISTIC[3])
                   STATISTIC[2] <- STATISTIC[2]/sqrt(STATISTIC[4])
                   CV$cv[1] <- CV$cv[1]/sqrt(STATISTIC[3])
                   CV$cv[2] <- CV$cv[2]/sqrt(STATISTIC[4])
                   
                   H0 <- (STATISTIC[1:2] > CV$cv)
                   
                   res <- new("kb.test", Un = STATISTIC[1:2], CV_Un = CV$cv, 
                              H0_Un = H0, method = METHOD, 
                              data = list(x = x, y = y), cv_method = method, 
                              B= B, h= h_best, var_Un= STATISTIC[3:4])
                } else {
                   
                   
                   # Check that they have the same number of rows (observations)
                   if(!is.null(y) && nrow(x) != nrow(y)) {
                      stop("'x' and 'y' must have the same number of rows.")
                   }

                   METHOD <- "Kernel-based quadratic distance k-sample test"
                   
                   sizes <- as.vector(table(y))
                   cum_size <- c(0,cumsum(sizes))
                   STATISTIC <- stat_ksample_cpp(x, c(y), h, sizes, cum_size,
                                                 compute_variance = TRUE)
                   
                   CV <- cv_ksample(x, y, h, B, b, Quantile, method,
                                    compute_variance = TRUE)
                   
                   STATISTIC[1] <- STATISTIC[1]/sqrt(STATISTIC[3])
                   STATISTIC[2] <- STATISTIC[2]/sqrt(STATISTIC[4])
                   CV$cv[1] <- CV$cv[1]/sqrt(STATISTIC[3])
                   CV$cv[2] <- CV$cv[2]/sqrt(STATISTIC[4])
                   
                   H0 <- (STATISTIC[1:2] > CV$cv)
                   
                   res <- new("kb.test", Un = STATISTIC[1:2], CV_Un = CV$cv, 
                              H0_Un = H0, method = METHOD,
                              data = list(x = x, y = y), cv_method = method, 
                              B= B, h= h_best, var_Un = STATISTIC[3:4])
                }
             }
             
             return(res)
          })
#' @rdname kb.test
#'
#' @param object Object of class \code{kb.test}
#'
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @export
setMethod("show", "kb.test",
 function(object) {
    cat( "\n", object@method, "\n")
    
    if(length(object@Vn)==0){
       
       cat("U-statistic\t Dn \t\t Trace \n")
       cat("------------------------------------------------\n")
       cat("Test Statistic:\t", object@Un[1], "\t", object@Un[2], "\n")
       cat("Critical Value:\t", object@CV_Un[1], "\t", object@CV_Un[2], "\n")
       cat("H0 is rejected:\t", object@H0_Un[1], "\t\t", object@H0_Un[2], "\n")
       
       cat("CV method: ", object@cv_method, "\n")
    } else {
       
       cat("\t\tU-statistic\tV-statistic\n")
       cat("------------------------------------------------\n")
       cat("Test Statistic:\t", object@Un, "\t", object@Vn, "\n")
       cat("Critical Value:\t", object@CV_Un, "\t", object@CV_Vn, "\n")
       cat("H0 is rejected:\t", object@H0_Un, "\t\t", object@H0_Vn, "\n")
       
    }
    
    cat("Selected tuning parameter h: ", object@h$h_sel, "\n")
    
    cat("\n")
 })
#'
#' Summarizing kernel-based quadratic distance results
#'
#' \code{summary} method for the class \code{kb.test}
#'
#' @param object Object of class \code{kb.test}
#'
#' @return List with the following components:
#' \itemize{
#'    \item \code{summary_tables} Table of computed descriptive statistics per 
#'    variable (and per group if available).
#'    \item \code{test_results} Data frame with the results of the performed 
#'    kernel-based quadratic distance test.
#'    \item \code{qqplots} Figure with qq-plots for each variable.
#' }
#' @seealso [kb.test()] and \linkS4class{kb.test} for more details.
#'
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#'
#'@examples
#' # create a kb.test object
#' x <- matrix(rnorm(100),ncol=2)
#' # Normality test
#' my_test <- kb.test(x, h=0.5)
#' summary(my_test)
#' 
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @importFrom stats IQR
#' @importFrom stats median
#' @importFrom stats sd
#' @importFrom stats qqnorm
#' 
#' @name summary.kb.test
#' @rdname summary.kb.test
#' @aliases summary,kb.test-method
#' @export
setMethod("summary", "kb.test", function(object) {
   
   if(object@method=="Kernel-based quadratic distance k-sample test"){
      
      x <- object@data$x
      y <- object@data$y
      k <- length(unique(y))
      
      stats <- list()
      for(i in seq_len(ncol(x))){
         res <- rbind(as.numeric(by(x[,i],y,mean)),
                      as.numeric(by(x[,i],y,sd)),
                      as.numeric(by(x[,i],y,median)),
                      as.numeric(by(x[,i],y,IQR)),
                      as.numeric(by(x[,i],y,min)),
                      as.numeric(by(x[,i],y,max))
         )
         res <- as.data.frame(res)
         res <- cbind(res, c(mean(x[,i]),sd(x[,i]),median(x[,i]),
                             IQR(x[,i]),min(x[,i]),max(x[,i])))
         rownames(res) <- c("mean","sd", "median", "IQR", "min", "max")
         colnames(res) <- c(paste("Group ",seq(1,k),sep=""),"Overall")
         stats[[i]] <- res
      }
      
      figure <- NA
      
      
   }
   
   if(object@method=="Kernel-based quadratic distance two-sample test"){
      
      sample1 <- as.data.frame(object@data$x)
      sample2 <- as.data.frame(object@data$y)
      colnames(sample2) <- colnames(sample1)
      plot_list <- lapply(names(sample1), function(name) {
         list(
            compare_qq(sample1[[name]], sample2[[name]], name)
            #compute_stats(sample1[[name]], sample2[[name]], name)$plots
         )
      })
      plot_list <- do.call(c, plot_list)
      figure <- ggarrange(plotlist = plot_list, ncol = 1)
      # figure <- ggarrange(plotlist = plot_list, ncol = 2, 
      #                     nrow = length(plot_list) / 2, widths=c(1,1.3))
       
      stats <- lapply(names(sample1), function(name) {
         
         compute_stats(sample1[[name]], sample2[[name]], name)$stats
         
      })
      print(figure)
      
   }
   
   if(object@method=="Kernel-based quadratic distance Normality test"){
      
      dat_x <- as.data.frame(object@data$x)
      
      plot_list <- list()
      stats <- list()
      for(i in seq_len(ncol(dat_x))) {
         
         # qq_df <- data.frame(x = sort(qqnorm(dat_x[,i], plot = FALSE)$x), 
         #             sample_quantiles = quantile(dat_x[,i], 
         #                         probs = seq(0, 1, length.out = nrow(dat_x))))
         x <- sort(qqnorm(dat_x[,i], plot = FALSE)$x)
         sample_quantiles <- quantile(dat_x[,i], 
                                 probs = seq(0, 1, length.out = nrow(dat_x)))
         
         pl <- ggplot(mapping=aes(x = x, y = sample_quantiles)) +
            geom_line(col="blue") +
            theme_minimal()+
            geom_abline(slope = 1, intercept = 0,col="red") +
            ggtitle(paste("QQ Plot against Normal - ",names(dat_x)[i])) +
            xlab("Theoretical Quantiles") +
            ylab("Sample Quantiles")
         
         stats_step <- data.frame(matrix(c(mean(dat_x[,i]),sd(dat_x[,i]),
                                           median(dat_x[,i]),IQR(dat_x[,i]),
                                           min(dat_x[,i]),max(dat_x[,i])),
                                             nrow=6,ncol=1,byrow=TRUE))
         colnames(stats_step) <- c(paste(names(dat_x)[i]))
         rownames(stats_step) <- c("mean", "sd", "median", "IQR", "min", "max")
         
         stats[length(stats) +1] <- stats_step
         
   #    pl_stat <- ggplot() +
   # geom_table_npc(data = data.frame(Stat = rownames(stats_step), stats_step),
   #                      aes(npcx = 0.5, npcy = 0.5, 
   #        label = list(data.frame(Stat = rownames(stats_step), stats_step))),
   #                      hjust = 0.5, vjust = 0.5) +
   #    # annotate('table', x = 0.5, y = 0.5, 
   #    #          label = data.frame(Stat = rownames(stats_step),stats_step),
   #    #          hjust = 0.5, vjust = 0.5) +
   #       theme_void() +
   #       ggtitle("")+
   #       scale_color_brewer(palette='Set1')
   #    
      #plot_list[[length(plot_list) + 1]] <- list(pl,pl_stat)
         plot_list[[length(plot_list) + 1]] <- list(pl)
         
      }
      plot_list <- do.call(c, plot_list)
      figure <- ggarrange(plotlist = plot_list,ncol = 1)
      print(figure)
      
      stats <- do.call(cbind, stats)
      colnames(stats) <- c(names(dat_x))
      rownames(stats) <- c("mean", "sd", "median", "IQR", "min", "max")
      
      
   }
   # Print main results of the test
   cat( "\n", object@method, "\n")
   if(length(object@Vn)==0){
      test_results <- data.frame(
         Statistic = c("Dn", "Trace"),
         Value = object@Un,
         Critical_Value = object@CV_Un,
         Reject_H0 = object@H0_Un
      )
   } else {
      test_results <- data.frame(
         Statistic = c("Un", "Vn"),
         Value = c(object@Un,object@Vn),
         Critical_Value = c(object@CV_Un,object@CV_Vn),
         Reject_H0 = c(object@H0_Un,object@H0_Vn)
      )
   }
   print(test_results)
   
   return(list(summary_tables = stats, 
               test_results = test_results, 
               qqplots = figure))
})
