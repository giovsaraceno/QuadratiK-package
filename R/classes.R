#' @rdname kb.test-class
#'
#' @title An S4 class for kernel-based distance tests with normal kernel
#'
#' @description A class to represent the results of Gaussian kernel-based 
#' quadratic distance tests. This includes the normality test, the two-sample 
#' test statistics and the k-sample tests.
#'
#' @slot method String indicating the normal kernel-based quadratic distance 
#'              test performed.
#' @slot Un The value of the test U-statistics.
#' @slot Vn The value of the test V-statistic.
#' @slot H0_Un A logical value indicating whether or not the null hypothesis is 
#'          rejected according to U-statistics.
#' @slot H0_Vn A logical value indicating whether or not the null hypothesis is 
#'          rejected according to Vn.
#' @slot data List of samples X (and Y).
#' @slot CV_Un The critical value computed for the test Un.
#' @slot CV_Vn The critical value computed for the test Vn.
#' @slot cv_method The method used to estimate the critical value (one of 
#' "subsampling", "permutation" or "bootstrap").
#' @slot h A list with the value of bandwidth parameter used for the Gaussian 
#'         kernel. If the function \code{select_h} is used, then also the matrix
#'         of computed power values and the resulting power plot are provided. 
#' @slot B Number of bootstrap/permutation/subsampling replications.
#' @slot var_Un exact variance of the kernel-based U-statistic.
#'
#' @examples
#' # create a kb.test object
#' x <- matrix(rnorm(100),ncol=2)
#' y <- matrix(rnorm(100),ncol=2)
#' # Normality test
#' kb.test(x, h=0.5)
#' 
#' # Two-sample test
#' kb.test(x,y,h=0.5, method="subsampling",b=0.9)
#' 
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @export
setClass("kb.test",
         slots = list(
            method = "character",
            Un = "numeric",
            Vn = "numeric",
            CV_Un = "numeric",
            CV_Vn = "numeric",
            H0_Un = "logical",
            H0_Vn = "logical",
            data = "list",
            cv_method = "character",
            B = "numeric",
            h = "list",
            var_Un = "numeric"
         )
)
#' @rdname pk.test-class
#'
#' @title An S4 class for Poisson kernel-based quadratic distance tests.
#'
#' @description A class to represent the results of Poisson kernel-based 
#'              quadratic distance tests for Uniformity on the sphere.
#'
#' @slot method The method used for the test ("Poisson Kernel-based quadratic 
#'              distance test of Uniformity on the Sphere").
#' @slot x Matrix of data
#' @slot Un The value of the U-statistic.
#' @slot CV_Un The critical value for Un computed through replications.
#' @slot H0_Un A logical value indicating whether or not the null hypothesis is 
#'             rejected according to Un.
#' @slot Vn The value of the V-statistic.
#' @slot CV_Vn The critical value for Vn computed following the asymptotic 
#'             distribution.
#' @slot H0_Vn A logical value indicating whether or not the null hypothesis is 
#'             rejected according to Vn.
#' @slot rho The concentration parameter of the Poisson kernel.
#' @slot B Number of replications.
#' @slot var_Un exact variance of the kernel-based U-statistic.
#'
#' @examples
#' # create a pk.test object
#' d=3
#' size=200
#' x_sp <- sample_hypersphere(d, n_points=size)
#' pk.test(x_sp,rho=0.8)
#'
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @export
setClass("pk.test",
         slots = list(
            method = "character",
            Un = "numeric",
            CV_Un = "numeric",
            H0_Un = "logical",
            Vn = "numeric",
            CV_Vn = "numeric",
            H0_Vn = "logical",
            x = "matrix",
            B = "numeric",
            rho = "numeric",
            var_Un = "numeric"
         )
)
#' @rdname pkbc-class
#'
#' @title A S4 class for the clustering algorithm on the sphere based on
#' Poisson kernel-based distributions.
#'
#' @description A class to represent the results of Poisson kernel-based 
#'              clustering procedure for spherical observations.
#'
#' @slot res_k List of objects with the results of the clustering algorithm for 
#'             each value of possible number of clusters considered.
#' @slot input List of input data
#'
#' @details See the function \code{pkbc} for more details.
#' 
#' @examples 
#' data("wireless")
#' res <- pkbc(as.matrix(wireless[,-8]),4)
#' 
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {UL4.0, UL4.1, UL4.2} class object created for the clustering 
#'            results. It allows to create an object without running the 
#'            algorithm. All the input control parameters can be extracted from
#'            the result object.
#' 
#' @export
setClass("pkbc",
         slots = list(
            res_k = "list",
            input = "list"
         )
)
