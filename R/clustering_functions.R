#'
#' Poisson kernel-based clustering on the sphere
#'
#' @description
#' The function \code{pkbc()} performs the Poisson kernel-based clustering 
#' algorithm on the sphere proposed by Golzy and Markatou (2020). 
#' The proposed algorithm is based on a mixture, with \eqn{M} components, of 
#' Poisson kernel-based densities on the hypersphere \eqn{\mathcal{S}^{d-1}}
#' given by
#' \deqn{f(x|\Theta) = \sum_{j=1}^M \alpha_j f_j(x|\rho_j, \mu_j)}
#' where \eqn{\alpha_j}'s are the mixing proportions and \eqn{f_j(x|\rho_j, 
#' \mu_j)}'s denote the probability density function of a \eqn{d}-variate 
#' Poisson kernel-based density given as
#' \deqn{f(\mathbf{x}|\rho, \mathbf{\mu}) = \frac{1-\rho^2}{\omega_d 
#' ||\mathbf{x} - \rho \mathbf{\mu}||^d}.}
#' The parameters \eqn{\alpha_j, \mu_j, \rho_j} are estimated through a 
#' iterative reweighted EM algorithm. \cr
#' The proposed clustering algorithm exhibits excellent results when 
#' (1) the clusters are not well separated; (2) the data points
#' are fairly well concentrated around the vectors \eqn{\mu_j} of each cluster;
#' (3) the percentage of noise in the data increases.
#' 
#' @param dat Data matrix or data.frame of data points on the sphere to be 
#'            clustered. The observations in \code{dat} are normalized by 
#'            dividing with the length of the vector to ensure
#'            that they lie on the \eqn{d}-dimensional sphere. Note that 
#'            \eqn{d > 1}.
#' @param nClust Number of clusters. It can be a single value or a numeric 
#'               vector.
#' @param maxIter The maximum number of iterations before a run is terminated.
#' @param stoppingRule String describing the stopping rule to be used within 
#'                     each run. Currently must be either \code{'max'}, 
#'                     \code{'membership'}, or \code{'loglik'}.
#' @param initMethod String describing the initialization method to be used.
#'                   Currently must be \code{'sampleData'}.
#' @param numInit Number of initialization.
#'
#' @details 
#' We set all concentration parameters equal to 0.5 and all mixing proportions
#' to be equal. \cr
#' The initialization method \code{'sampleData'} indicates that observation 
#' points are randomly chosen as initializers of the centroids \eqn{\mu_j}.
#' This random starts strategy has a chance of not obtaining initial 
#' representatives from the underlying clusters, then the clustering is 
#' performed \code{numInit} times and the random start with the highest
#' likelihood is chosen as the final estimate of the parameters.
#' 
#' The possible \code{stoppingRule} for each iteration are: \cr
#' - \code{'loglik'} run the algorithm until the change in log-likelihood from 
#' one iteration to the next is less than a given threshold (1e-7) \cr
#' - \code{'membership'} run the algorithm until the membership is unchanged for
#' all points from one iteration to the next \cr
#' - \code{'max'} reach a maximum number of iterations \code{maxIter}
#' 
#' The obtained estimates are used for assigning final memberships, identifying
#' the \code{nClust} clusters, according to the following rule
#' \deqn{P(x_i, \Theta) = \arg\max_{j \in \{1, \ldots, k\}} \{ \frac{\alpha_j 
#' f_j(x_i|\mu_j, \rho_j)}{f(x_i, \Theta)}\}.}
#' The number of clusters \code{nClust} must be provided as input to the
#' clustering algorithm.
#' 
#' @seealso [dpkb()] and [rpkb()] for more information on the Poisson 
#'          kernel-based distribution. \cr
#'          \linkS4class{pkbc} for the class definition.
#'          
#' @note
#' The clustering algorithm is tailored for data points on the sphere
#' \eqn{\mathcal{S}^{d-1}}, but it can also be performed on spherically
#' transformed observations, i.e. data points on the Euclidean space 
#' \eqn{\mathbb{R}^d} that are normalized such that they lie on the 
#' corresponding \eqn{d}-dimensional sphere \eqn{\mathcal{S}^{d-1}}.
#'
#' @return An S4 object of class \code{pkbc} containing the results of the 
#' clustering procedure based on Poisson kernel-based distributions. The object 
#' contains the following slots:
#'
#'   \code{res_k}: List of results of the Poisson kernel-based clustering 
#'   algorithm for each value of number of clusters specified in \code{nClust}. 
#'   Each object in the list contains:
#'   \itemize{
#'      \item \code{postProbs} Posterior probabilities of each observation for 
#'      the indicated clusters.
#'      \item \code{LogLik} Maximum value of log-likelihood function
#'      \item \code{wcss} Values of within-cluster sum of squares computed with 
#'      Euclidean distance and cosine similarity, respectively.
#'      \item \code{params} List of estimated parameters of the mixture model
#'           \itemize{
#'              \item \code{mu} estimated centroids
#'              \item \code{rho} estimated concentration parameters rho
#'              \item \code{alpha} estimated mixing proportions
#'           }
#'      \item \code{finalMemb} Vector of final memberships
#'      \item \code{runInfo} List of information of the EM algorithm iterations
#'           \itemize{
#'              \item \code{lokLikVec} vector of log-likelihood values
#'              \item \code{numIterPerRun} number of E-M iterations per run
#'           }
#'   }
#'   \code{input}: List of input information.
#'
#' @references
#' Golzy, M. and Markatou, M. (2020) Poisson Kernel-Based Clustering on the 
#' Sphere: Convergence Properties, Identifiability, and a Method of Sampling, 
#' Journal of Computational and Graphical Statistics, 29:4, 758-770, 
#' DOI: 10.1080/10618600.2020.1740713.
#'
#' @examples
#' #We generate three samples of 100 observations from 3-dimensional
#' #Poisson kernel-based densities with rho=0.8 and different mean directions
#' size<-100
#' groups<-c(rep(1, size), rep(2, size),rep(3,size))
#' rho<-0.8
#' set.seed(081423)
#' data1<-rpkb(size, c(1,0,0),rho)
#' data2<-rpkb(size, c(0,1,0),rho)
#' data3<-rpkb(size, c(0,0,1),rho)
#' dat<-rbind(data1$x,data2$x, data3$x)
#'
#' #Perform the clustering algorithm with number of clusters k=3.
#' pkbd<- pkbc(dat=dat, nClust=3)
#' show(pkbd)
#'
#' @srrstats {G2.0a} Documentation of input nClust
#' @srrstats {G1.3} description of parameter
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G2.0,G2.2,G2.3a} The code considers the different types of input
#' 
#' @srrstats {UL1.0,UL1.1} format of input data is specified, with checks 
#' @srrstats {UL1.4} assumptions are taken into consideration
#' 
#' @export
setGeneric("pkbc",function(dat, 
                           nClust = NULL,
                           maxIter = 300,
                           stoppingRule = "loglik",
                           initMethod = "sampleData",
                           numInit = 10){
   
   standardGeneric("pkbc")
})
#' @rdname pkbc
#' 
#' @srrstats {G2.0} input nClust
#' @srrstats {G2.7,G2.8} different types of tabular input are considered
#' @srrstats {G2.13,G2.14,G2.14a,G2.15,G2.16} error for NA, Nan, Inf, -Inf
#' 
#' @srrstats {UL2.1} normalization of input dat
#' @srrstats {UL2.3} inputs are all checked
#' @srrstats {UL7.0} error messages are inserted for inappropriate inputs
#' 
#' @importFrom methods new
#' @importFrom stats uniroot
#' 
#' 
#' @export
setMethod("pkbc", signature(dat = "ANY"),
    function(dat,
             nClust = NULL,
             maxIter = 300,
             stoppingRule = "loglik",
             initMethod = "sampleData",
             numInit = 10){
       # Constant defining threshold by which log likelihood must change 
       # to continue iterations, if applicable.
       LL_TOL <- 1e-7
       
       # validate input
       if(is.null(nClust)) {
          stop("Input parameter nClust is required. Provide one specific 
               value or a set of possible values.")
       } else if(is.vector(nClust) & is.numeric(nClust)) {
          if(any(nClust < 1)){
             stop('Values in the input parameter nClust must be 
                  greater than 0')
          }
       } else {
          stop("nClust must be a signle value or a numeric vector of possible
                  values")
       }
       if (maxIter < 1) {
          stop("Input parameter maxIter must be greater than 0")
       }
       if ( !(stoppingRule %in% c("max", "membership", "loglik")) ) {
          stop(paste("Unrecognized value ", stoppingRule, " in input 
          parameter stoppingRule.", sep=''))
       }
       if (!(initMethod %in% c("sampleData"))) {
          stop(paste("Unrecognized value ", initMethod, " in input 
          parameter initMethod.", sep=''))
       }
       if (numInit < 1) {
          stop("Input parameter numInit must be greater than 0")
       }
       
       # set options for stopping rule
       checkMembership <- stoppingRule == "membership"
       checkLoglik <- stoppingRule == "loglik"
       
       if(is.data.frame(dat)){
          dat <- as.matrix(dat)
       } else if(!is.matrix(dat)){
          stop("dat must be a matrix or a data.frame")
       }
       if(!is.numeric(dat)) {
          stop("dat must be a numeric matrix or data.frame")
       }
       
       if(any(is.na(dat))) {
          stop("There are missing values in the data set!")
       } else if(any(is.infinite(dat) |is.nan(dat))){
          stop("There are undefined values, that is Nan, Inf, -Inf")
       }
       # Normalize the data
       dat <- dat/sqrt(rowSums(dat^2))
       
       # Initialize variables of interest
       numVar <- ncol(dat)
       numData <- nrow(dat)
       
       res_k <- list()
       for(numClust in nClust){
          
          logLikVec <- rep(-Inf, numInit)
          numIterPerRun <- rep(-99, numInit)
          alpha_best <- rep(-99, numClust)
          rho_best <- rep(-99, numClust)
          mu_best <- matrix(nrow = numClust, ncol = numVar)
          normprobMat_best <- matrix(-99, nrow=numData, ncol=numClust)
          if (initMethod == "sampleData") {
             uniqueData <- unique(dat)
             numUniqueObs <- nrow(uniqueData)
             if (numUniqueObs < numClust) {
                stop(paste("Only", numUniqueObs, "unique observations.",
                           'When initMethod = "sampleData", must have more
                           than nClust unique observations.'
                ))
             }
          }
          log_w_d <- (numVar / 2) * (log(2) + log(pi)) - lgamma(numVar / 2)
          
          # Begin initializations
          for (init in 1:numInit) {
             
             # initialize parameters
             alpha_current <- rep(1/numClust, numClust)
             rho_current <- rep(0.5,numClust)
             
             # Initializing mu
             if (initMethod == 'sampleData') {
                # sampling w/o replacement from unique data points
                mu_current <- matrix(t(uniqueData[sample(numUniqueObs, 
                                 size=numClust, replace=FALSE) ,]), 
                                 ncol=numVar, byrow=TRUE)
             }
             
             currentIter <- 1
             membCurrent <- rep.int(0, times = numData)
             loglikCurrent <- -Inf
             
             # Begin EM iterations
             while (currentIter <= maxIter) {
                v_mat <- dat %*% t(mu_current)
                alpha_mat_current <- matrix(
                   alpha_current,
                   nrow = numData,
                   ncol = numClust,
                   byrow = TRUE
                )
                rho_mat_current <- matrix(
                   rho_current,
                   nrow = numData,
                   ncol = numClust,
                   byrow = TRUE
                )
                log_probMat_denom <- log(1 + rho_mat_current^2 - 
                                            2*rho_mat_current*v_mat)
                # eq (11) in ICMLD16 paper
                log_probMat <- log(1 - (rho_mat_current^2)) - log_w_d - 
                                    (numVar / 2) * log_probMat_denom
                ######### E step done#####################################
                ########## M step now#####################################
                # Denominator of eq (18) in ICMLD16 paper
                probSum <- matrix(
                   exp(log_probMat) %*% alpha_current,
                   nrow = numData,
                   ncol = numClust
                )
                # eq (18) in ICMLD16 paper
                log_normProbMat_current <- log(alpha_mat_current) + 
                                             log_probMat - log(probSum)
                # denominator of eq (20) in ICMLD16 paper
                log_weightMat <- log_normProbMat_current - log_probMat_denom
                
                # eq (19) in ICMLD16 paper
                alpha_current <- colSums(exp(log_normProbMat_current)) / numData
                # numerator of fraction in eq (21) of ICMLD16 paper
                mu_num_sum_MAT <- t(exp(log_weightMat)) %*% dat
                norm_vec <- function(x) sqrt(sum(x^2))
                mu_denom <- apply(mu_num_sum_MAT, 1, norm_vec)
                # eq (21) of ICMLD16 paper without sign function
                mu_current <- mu_num_sum_MAT / mu_denom
                for(h in 1:numClust){
                   # update rho params
                   sum_h_weightMat <- sum(exp(log_weightMat[,h]))
                   alpha_current_h <- alpha_current[h]
                   mu_denom_h <- mu_denom[h]
                   root_func <- function(x) {
                      (-2*numData*x*alpha_current_h)/(1 - x^2) + 
                         numVar*mu_denom_h -numVar*x*sum_h_weightMat
                   }
                   rho_current[h] <- (uniroot(
                      f = root_func,
                      interval = c(0,1),
                      f.lower = root_func(0),
                      f.upper = root_func(1),
                      tol = 0.001
                   ))$root
                } # for h
                
                # Terminate iterations if maxIter is reached. Redundant 
                # with while condition, but this ensures that the stored 
                # numIter below is accurate.
                if (currentIter >= maxIter) {
                   break
                }
                
                # Terminate iterations if membership is unchanged, if applicable
                if (checkMembership) {
                   # calculate and store group membership
                   membPrevious <- membCurrent
                   membCurrent <-apply(exp(log_normProbMat_current),1,which.max)
                   if (all(membPrevious == membCurrent)) {
                      break
                   }
                }
                # Terminate iterations if log-likelihood is unchanged, if 
                # applicable.
                if (checkLoglik) {
                   loglikPrevious <- loglikCurrent
                   # Calculate log likelihood for the current iteration
                   loglikCurrent <- sum(log( exp(log_probMat)%*%alpha_current))
                   if (abs(loglikPrevious - loglikCurrent) < LL_TOL) {
                      break
                   }
                }
                # Update counter to NEXT iteration number.
                currentIter <- currentIter + 1
             } # for currentIter
             # Compare the current log likelihood to the best stored log 
             # likelihood; if it exceeds the stored, then store it and all 
             # associated estimates.
             # alpha, rho, mu, normprobMat
             loglikCurrent <- sum(log( exp(log_probMat) %*% alpha_current ))
             if (loglikCurrent > max(logLikVec)) {
                alpha_best <- alpha_current
                rho_best <- rho_current
                mu_best <- mu_current
                normprobMat_best <- exp(log_normProbMat_current)
             }
             logLikVec[init] <- loglikCurrent
             numIterPerRun[init] <- currentIter - 1
             # Store run info
          } # for init
          # Calculate final membership
          memb_best <- apply(normprobMat_best, 1, which.max)
          
          # Calculate the within-cluster sum of squares
          wcss <- vapply(1:numClust, function(cl) {
             idx <- which(memb_best==cl)
             
             if(length(idx)>0){
                
                dat_memb <- matrix(t(dat[idx,]),ncol=numVar, byrow=TRUE)
                mu_memb <- matrix(rep(as.numeric(mu_best[cl,]), 
                                      times = nrow(dat_memb)), 
                                       ncol=numVar, 
                                       byrow=TRUE)
                return(sum(apply(dat_memb - mu_memb, 1, norm, "2")**2))
                
             } else {
                return(0)
             }
          }, numeric(1))
          wcss <- sum(wcss)
          
          # Calculate the within-cluster sum of squares using cosine 
          # similarity
          wcss_cos <- vapply(1:numClust, function(cl) {
             idx <- which(memb_best==cl)
             if(length(idx)>0){
                dat_memb <- matrix(t(dat[idx,]),ncol=numVar, byrow=TRUE)
                mu_memb <- matrix(rep(mu_best[cl,], nrow(dat_memb)), 
                                  ncol=numVar,byrow=TRUE)
                
                return(sum(dat_memb * mu_memb))
             } else {
                return(0)
             }
          }, numeric(1))
          
          wcss_cos <- sum(wcss_cos)
          
          res_k[[numClust]] <- list(postProbs = normprobMat_best,
                                    LogLik = max(logLikVec),
                                    wcss = c(wcss,wcss_cos),
                                    params = list(mu = mu_best, 
                                                  rho = rho_best, 
                                                  alpha = alpha_best),
                                    finalMemb = memb_best,
                                    runInfo = list(
                                       logLikVec = logLikVec,
                                       numIterPerRun = numIterPerRun
                                    ))
          
       }
       
       
       # Store and return results object.
       results <- new("pkbc",
                      res_k = res_k,
                      input = list(
                         dat = dat,
                         nClust = nClust,
                         maxIter = maxIter,
                         stoppingRule = stoppingRule,
                         initMethod = initMethod,
                         numInit = numInit
                      ))
       return(results)
    })
#' @rdname pkbc
#'
#' @param object Object of class \code{pkbc}
#'
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {UL4.3,UL4.3a} this function reports the print method
#' 
#' @export
setMethod("show", "pkbc", function(object) {
   cat("Poisson Kernel-Based Clustering on the Sphere (pkbc) Object\n")
   cat("------------------------------------------------------------\n\n")
   
   cat("Available components:\n")
   # Print input parameters
   cat("Input Parameters:\n")
   print(names(object@input))
   cat("\n\n")
   
   cat("Considered possible number of clusters: ", object@input$nClust, "\n\n")
   
   cat("Available components for each value of number of clusters: \n")
   
   i <- object@input$nClust[1]
   print(names(object@res_k[[i]]))
   
   invisible(object)
})
#'
#' Summarizing PKBD mixture Fits
#' 
#' Summary method for class "pkbc"
#' 
#' @param object Object of class \code{pkbc}
#' 
#' @return Display the logLikelihood values and within cluster sum of squares 
#'         (wcss) for all the values of number of clusters provided. For each of
#'         these values the estimated mixing proportions are showed together 
#'         with a table with the assigned memberships.
#' 
#' @examples
#' dat <- rbind(matrix(rnorm(100),2),matrix(rnorm(100,5),2))
#' res <- pkbc(dat,2:4)
#' summary(res)
#' 
#' @seealso [pkbc()] for the clustering algorithm \cr
#'          \linkS4class{pkbc} for the class object definition. 
#'
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {UL4.4} summary method for pkbc object
#' 
#' @name summary.pkbc
#' @rdname summary.pkbc
#' @aliases summary,pkbc-method
#' @export
setMethod("summary", "pkbc", function(object) {
   cat("Poisson Kernel-Based Clustering on the Sphere (pkbc) Results\n")
   cat("------------------------------------------------------------\n\n")
   
   summaryMatrix <- do.call(rbind, lapply(object@res_k[object@input$nClust], 
            function(res) {
               c(LogLik = res$LogLik, wcss = sum(res$wcss))
            }))
   rownames(summaryMatrix) <- names(object@res_k[object@input$nClust])
   colnames(summaryMatrix) <- c("LogLik", "WCSS")
   cat("Summary:\n")
   print(summaryMatrix)
   cat("\n")
   
   # Detailing each cluster's results
   for (numClust in object@input$nClust) {
      res <- object@res_k[[numClust]]
      
      cat(sprintf("Results for %s clusters:\n", numClust))
      
      cat("Estimated Mixing Proportions (alpha):\n")
      print(res$params$alpha)
      cat("\n")
      
      cat("Clustering table:\n")
      print(table(res$finalMemb))
      cat("\n\n")
   }
   
   invisible(object)
})
#' 
#' @title Descriptive statistics for the clusters identified by the Poisson 
#' kernel-based clustering.
#'
#' @description
#' Method for objects of class \code{pkbc} which computes some 
#' descriptive for each variable with respect to the detected groups.
#' 
#' @param object Object 
#' @param ... possible additional inputs
#' 
#' @seealso [pkbc()] for the clustering algorithm \cr
#'          \linkS4class{pkbc} for the class object definition.
#'          
#' @export
setGeneric("stats_clusters",function(object,...){
   
   standardGeneric("stats_clusters")
})
#' 
#' Descriptive statistics for the identified clusters
#'
#' Method for objects of class \code{pkbc} which computes descriptive 
#' statistics for each variable with respect to the detected groups. 
#'
#' @param object Object of class \code{pkbc}.
#' @param k Number of clusters to be used.
#'
#' @details 
#' The function computes mean, standard deviation, median, 
#' inter-quantile range, minimum and maximum for each variable in the data set 
#' given the final membership assigned by the clustering algorithm. 
#' 
#' @examples
#' #We generate three samples of 100 observations from 3-dimensional
#' #Poisson kernel-based densities with rho=0.8 and different mean directions
#' dat<-matrix(rnorm(300),ncol=3)
#'
#' #Perform the clustering algorithm
#' pkbc_res<- pkbc(dat, 3)
#' stats_clusters(pkbc_res, 3)
#' 
#'
#' @return List with computed descriptive statistics for each dimension. 
#'
#'
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {UL3.2} true label can be provided as a separate input
#' @srrstats {UL3.4} the function computes summary statistics with respect to 
#'                   the identified clusters.
#'                   
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats IQR 
#' 
#' @name stats_clusters
#' @rdname stats_clusters 
#' @aliases stats_clusters,pkbc-method                 
#' @export
setMethod("stats_clusters", "pkbc", function(object, k){
   
   if(!(is.numeric(k) & length(k)==1)){
      stop("k must be an integer")
   } else if(!(k %in% object@input$nClust)){
      stop("The provided pkbc object does not contain results for the requested
           number of clusters")
   } 
   
   x <- object@input$dat
   y <- as.factor(object@res_k[[k]]$finalMemb)
   
   metrics <- list()
   
   for(i in seq_len(ncol(x))){
      res <- rbind(as.numeric(by(x[,i],y,mean)),
                   as.numeric(by(x[,i],y,sd)),
                   as.numeric(by(x[,i],y,median)),
                   as.numeric(by(x[,i],y,IQR)),
                   as.numeric(by(x[,i],y,min)),
                   as.numeric(by(x[,i],y,max))
      )
      res <- cbind(res, c(mean(x[,i]), sd(x[,i]), median(x[,i]), 
                          IQR(x[,i]), min(x[,i]), max(x[,i])))
      res <- as.data.frame(res)
      rownames(res) <- c("mean","sd", "median", "IQR", "min", "max")
      colnames(res) <- c(paste("Group ",seq(1,k),sep=""), "Overall")
      metrics[[length(metrics)+1]] <- res
   }
   
   return(metrics)
})
#' 
#' 
#' @title Plotting method for Poisson kernel-based clustering
#' 
#' @description
#' Plots for a pkbc object.
#'  
#' @param x Object of class \code{pkbc}
#' @param ... Additional arguments that can be passed to the plot function
#' @param k number of considered clusters. If it is not provided the scatter 
#'          plot is displayed for each value of number of clusters present in 
#'          the \code{x} object 
#' @param true_label factor or vector of true membership to clusters (if 
#'                   available). It must have the same length of final 
#'                   memberships.
#' @param pca_res Logical. If TRUE the results from PCALocantore are also 
#'                reported (when dimension is greater than 3). 
#'
#' @details 
#' - scatterplot: If dimension is equal to 2 or 3, points are displayed on the 
#' circle and sphere, respectively. If dimension if greater than 3, the 
#' spherical Principal Component procedure proposed by Locantore et al. (1999),
#' is applied for dimensionality reduction and the first three principal
#' components are normalized and displayed on the sphere. For d > 3, the
#' complete results from the \code{PcaLocantore} function (package \code{rrcov})
#' are returned if \code{pca_res=TRUE}.
#' - elbow plot: the within cluster sum of squares (wcss) is computed using the 
#' Euclidean distance (left) and the cosine similarity (right). 
#' 
#' @note
#' The elbow plot is commonly used as a graphical method for choosing the
#' *appropriate* number of clusters. Specifically, plotting the wcss versus the
#' number of clusters, the suggested number of clusters correspond to the point
#' in which the plotted line has the greatest change in slope, showing 
#' an elbow.  
#' 
#' @seealso [pkbc()] for the clustering algorithm \cr
#'          \linkS4class{pkbc} for the class object definition.
#' 
#' @examples
#' dat<-matrix(rnorm(300),ncol=3)
#' pkbc_res<- pkbc(dat, 3)
#' plot(pkbc_res, 3)
#' 
#' @references
#' Locantore, N., Marron, J.S., Simpson, D.G. et al. (1999) "Robust principal 
#' component analysis for functional data." Test 8, 1–73. 
#' https://doi.org/10.1007/BF02595862
#' 
#' 
#' @srrstats {UL6.1} this function includes a plot method
#' @srrstats {UL6.0,UL6.2} plot method for pkbc object
#' 
#' @name plot.pkbc
#' @rdname plot.pkbc
#' @aliases plot.pkbc
#' @aliases plot,pkbc,ANY-method
#' @return The scatter-plot(s) and the elbow plot.
#' 
#' @export
setMethod("plot", c(x = "pkbc"), 
          function(x, k = NULL, true_label=NULL, pca_res=FALSE, ...) {
             
             if(is.null(k)){
                for (k in x@input$nClust) {
                   
                   scatterplotMethod(x, k, true_label, pca_res)
                }
             } else {
                
                scatterplotMethod(x, k, true_label, pca_res)
             }
             
             # Elbow plot
             elbowMethod(x)

})
#' Scatter-plot of data points colored by the final membership.
#' 
#' @param object Object of class \code{pkbc}
#' @param k Number of clusters to be used.
#' @param true_label factor or vector of true membership to clusters (if 
#'                   available). It must have the same length of final 
#'                   memberships.
#' @param pca_res Logical. If TRUE the results from PCALocantore when dimension
#'                os greater than 3 are also reported. 
#'
#' 
#' 
#' @return 
#' If dimension is equal to 2 or 3, points are displayed on the circle and 
#' sphere, respectively. If dimension if greater than 3, the spherical Principal
#' Component procedure proposed by Locantore et al., (1999) is applied for 
#' dimensionality reduction and the first three principal components are 
#' normalized and displayed on the sphere.For d > 3, the complete results from 
#' the \code{PcaLocantore} function (package \code{rrcov}) are returned if 
#' pca_res=TRUE.
#' 
#' @importFrom scatterplot3d scatterplot3d
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
#' @importFrom graphics par
#' @import ggplot2
#' @importFrom rrcov PcaLocantore
#' 
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @keywords internal
#' @noRd
scatterplotMethod <- function(object, k, true_label = NULL, pca_res = FALSE) {
   x <- object@input$dat
   y <- as.factor(object@res_k[[k]]$finalMemb)
   
   if (ncol(x) == 2) {
      # 2D Scatter plot using ggplot2
      df <- data.frame(V1 = x[,1], V2 = x[,2], clusters = as.factor(y))
      with(df, {
         pl <- ggplot(df, aes(x = V1, y = V2, color = clusters)) +
            geom_point(size = 2) +  # Smaller points
            theme_minimal() +
            labs(color = "Cluster") +
            theme(legend.position = "right")  # Add legend
         print(pl)
      })
      
   } else if (ncol(x) == 3) {
      # Define layout for two plots side by side if true_label is provided
      if (!is.null(true_label)) {
         par(mfrow = c(1, 2), mar = c(3, 3, 2, 1) + 0.1, oma = c(0, 0, 0, 0))  
      }
      
      # 3D Scatter plot for clustering result
      colors <- rainbow(length(unique(y)))
      color_labels <- colors[as.numeric(y)]
      
      s3d <- scatterplot3d(x[,1], x[,2], x[,3], color = color_labels, pch = 16, 
                           cex.symbols = 0.7, main = "Final Membership", 
                           xlab = "X", ylab = "Y", zlab = "Z",
                           xlim = range(x[,1]), ylim = range(x[,2]), 
                           zlim = range(x[,3]))
      
      legend("topright", legend = levels(y), col = colors, pch = 16, 
             pt.cex = 1.2, cex = 0.5, bty = "n", title = "Clusters")
      
      # If true labels are provided, plot them in a separate panel
      if (!is.null(true_label)) {
         colors_true <- rainbow(length(unique(true_label)))
         color_labels_true <- colors_true[as.numeric(true_label)]
         
         s3d_true <- scatterplot3d(x[,1], x[,2], x[,3], 
                                   color = color_labels_true, pch = 16, 
                                   cex.symbols = 0.7,
                                   main = "True Labels", xlab = "X", 
                                   ylab = "Y", zlab = "Z",
                                   xlim = range(x[,1]), ylim = range(x[,2]), 
                                   zlim = range(x[,3]))
         
         legend("topright", legend = levels(as.factor(true_label)), 
                col = colors_true, pch = 16, pt.cex = 1.2, cex = 0.5, 
                bty = "n", title = "True Labels")
      }
      
      # Reset layout if it was changed
      if (!is.null(true_label)) {
         par(mfrow = c(1, 1))  # Reset layout
      }
      
   } else {
      # For data with more than 3 dimensions, perform PCA and plot the 
      # first 3 components
      pca_result <- PcaLocantore(x)
      pca_data <- data.frame(PC1 = pca_result@scores[,1], 
                             PC2 = pca_result@scores[,2], 
                             PC3 = pca_result@scores[,3])
      pca_data <- pca_data / sqrt(rowSums(pca_data^2))
      pca_data <- data.frame(pca_data, Cluster = y)
      
      if (!is.null(true_label)) {
         par(mfrow = c(1, 2), mar = c(3, 3, 2, 1) + 0.1, oma = c(0, 0, 0, 0)) 
      }
      
      colors <- rainbow(length(unique(y)))
      color_labels <- colors[as.numeric(pca_data$Cluster)]
      
      s3d <- scatterplot3d(pca_data$PC1, pca_data$PC2, pca_data$PC3, 
                           color = color_labels, pch = 16, cex.symbols = 0.7,
                           main = "Final Membership", xlab = "PC1", 
                           ylab = "PC2", zlab = "PC3")
      
      legend("topright", legend = levels(y), col = colors, pch = 16, 
             pt.cex = 1.2, cex = 0.5, bty = "n", title = "Clusters")
      
      if (!is.null(true_label)) {
         colors_true <- rainbow(length(unique(true_label)))
         color_labels_true <- colors_true[as.numeric(true_label)]
         
         s3d_true <- scatterplot3d(pca_data$PC1, pca_data$PC2, pca_data$PC3, 
                                   color = color_labels_true, pch = 16, 
                                   cex.symbols = 0.7,
                                   main = "True Labels", xlab = "PC1", 
                                   ylab = "PC2", zlab = "PC3")
         
         legend("topright", legend = levels(as.factor(true_label)), 
                col = colors_true, pch = 16, pt.cex = 1.2, cex = 0.5, 
                bty = "n", title = "True Labels")
      }
      
      if (pca_res) {
         return(pca_result)
      }
      
      # Reset layout if it was changed
      if (!is.null(true_label)) {
         par(mfrow = c(1, 1))  # Reset layout
      }
   }
}
#'
#' Elbow plot
#'
#'
#' @param object Object of class \code{pkbc}
#'
#' @return Returns a figure with two panels representing the elbow plot for the
#'         within sum of squares computed with the Euclidean distance and the 
#'         cosine similarity.
#' 
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' 
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @keywords internal
#' @noRd
elbowMethod <- function(object){
   
   wcss <- matrix(ncol=2)
   
   for(k in object@input$nClust){
      wcss <- rbind(wcss, object@res_k[[k]]$wcss)
   }
   wcss <- wcss[-1,]
   wcss <- matrix(wcss,ncol=2)
   colnames(wcss) <- c("Euclidean","cosine_similarity")
   rownames(wcss) <- paste0("k=",object@input$nClust)
   
   wcss_values <- data.frame(k = object@input$nClust, wcss = wcss[,1])
   pl <- ggplot(wcss_values, aes(x = k, y = wcss)) +
      geom_line() +
      geom_point() +
      labs(title = "Euclidean", x = "Number of clusters", 
           y = "Within-cluster sum of squares (WCSS)") +
      theme_minimal()
   
   wcss_values <- data.frame(k = object@input$nClust, wcss = wcss[,2])
   pl_cos <- ggplot(wcss_values, aes(x = k, y = wcss)) +
      geom_line() +
      geom_point() +
      labs(title = "Cosine Similarity", x = "Number of clusters", 
           y = "Within-cluster sum of squares (WCSS)") +
      theme_minimal()
   
   fig <- ggarrange(plotlist = list(pl, pl_cos), nrow=1)
   print(fig)

}
#' @title
#' Cluster spherical observations using a mixture of Poisson kernel-based 
#' densities
#' 
#' @description
#' Obtain predictions of membership for spherical observations based on a 
#' mixture of Poisson kernel-based densities estimated by \code{pkbc}
#' 
#' @param object Object of class \code{pkbc}
#' @param k Number of clusters to be used.
#' @param newdata a data.frame or a matrix of the data. If missing the 
#'                clustering data obtained from the \code{pkbc} object are 
#'                classified.
#' 
#' @return Returns a list with the following components
#' \itemize{
#'    \item Memb: vector of predicted memberships of \code{newdata}
#'    \item Probs: matrix where entry (i,j) denotes the probability that
#'                 observation i belongs to the k-th cluster. 
#' }
#' 
#' @examples
#' # generate data
#' dat <- rbind(matrix(rnorm(100),ncol=2),matrix(rnorm(100,5),ncol=2))
#' res <- pkbc(dat,2)
#' 
#' # extract membership of dat
#' predict(res,k=2)
#' # predict membership of new data
#' newdat <- rbind(matrix(rnorm(10),ncol=2),matrix(rnorm(10,5),ncol=2))
#' predict(res, k=2, newdat)
#'  
#' @seealso [pkbc()] for the clustering algorithm \cr
#'          \linkS4class{pkbc} for the class object definition.
#' 
#' @srrstats {G1.a} roxygen2 is used
#' @srrstats {UL3.3} prediction function for pkbc object
#' 
#' @name predict.pkbc
#' @rdname predict.pkbc
#' @aliases predict,pkbc-method
#' @export
setMethod("predict", signature(object="pkbc"), 
          function(object, k, newdata=NULL){
  
   if(is.null(newdata)){
      return(object@res_k[[k]]$finalMemb)
   }
   if (!is.matrix(newdata) && !is.data.frame(newdata)) {
      stop("newdata must be a matrix or data.frame.")
   }
   
   if (is.data.frame(newdata)) {
      newdata <- as.matrix(newdata)
   }
   
   if(!is.numeric(newdata)){
      stop("newdata must be numeric")
   }
             
   # Ensure that x has the same number of columns as the training data
   if (ncol(newdata) != ncol(object@input$dat)) {
   stop("newdata must have the same number of variables as the training data.")
   }
   
   numVar <- ncol(newdata)
   numData <- nrow(newdata)
   
   # Normalize newdata as in the training process
   x <- newdata / sqrt(rowSums(newdata^2))
   
   # Extract model parameters from the pkbc object
   mu <- object@res_k[[k]]$params$mu
   rho <- object@res_k[[k]]$params$rho
   alpha <- object@res_k[[k]]$params$alpha
   
   log_w_d <- (numVar / 2) * (log(2) + log(pi)) - lgamma(numVar / 2)
   
   v_mat <- x %*% t(mu)
   alpha_mat <- matrix(
      alpha,
      nrow = numData,
      ncol = k,
      byrow = TRUE
   )
   rho_mat <- matrix(
      rho,
      nrow = numData,
      ncol = k,
      byrow = TRUE
   )
   log_probMat_denom <- log(1 + rho_mat^2 - 2*rho_mat*v_mat)
   log_probMat <- log(1 - (rho_mat^2)) - log_w_d - 
                  (numVar / 2) * log_probMat_denom
   probSum <- matrix(
      exp(log_probMat) %*% alpha,
      nrow = numData,
      ncol = k
   )
   log_normProbMat <- log(alpha_mat) + log_probMat - log(probSum)
   memb <-apply(exp(log_normProbMat),1,which.max)
   
   res <- list()
   res$Memb <- memb
   res$Probs <- exp(log_normProbMat)
   return(res)
   
}) 
#'
#' Validation of Poisson kernel-based clustering results
#'
#' Method for objects of class \code{pkbc} which computes evaluation measures 
#' for clustering results.
#' The following evaluation measures are computed: 
#' In-Group Proportion (Kapp and Tibshirani (2007)). If true label are 
#' provided, ARI, Average Silhouette Width (Rousseeuw (1987)), Macro-Precision 
#' and Macro-Recall are computed.
#'
#' @param object Object of class \code{pkbc}
#' @param true_label factor or vector of true membership to clusters (if 
#'                   available). It must have the same length of final 
#'                   memberships.
#'
#' @details   
#' The IGP is a statistical measure that quantifies the proportion of 
#' observations within a group that belong to the same predefined category or 
#' class. It is often used to assess the homogeneity of a group by evaluating 
#' how many of its members share the same label. A higher IGP indicates that the
#' group is more cohesive, while a lower proportion suggests greater diversity 
#' or misclassification within the group (Kapp and Tibshirani 2007).
#' 
#' The Adjusted Rand Index (ARI) is a statistical measure used in data 
#' clustering analysis. It quantifies the similarity between two partitions of 
#' a dataset by comparing the assignments of data points to clusters. The ARI 
#' value ranges from 0 to 1, where a value of 1 indicates a perfect match 
#' between the partitions and a value close to 0 indicates a random assignment 
#' of data points to clusters.
#' 
#' Each cluster can represented by a so-called silhouette which is based on the
#' comparison of its tightness and separation. The average silhouette width 
#' provides an evaluation of clustering validity, and might be used to select 
#' an *appropriate* number of clusters (Rousseeuw 1987). 
#' 
#' Macro Precision is a metric used in multi-class classification that 
#' calculates the precision for each class independently and then takes the 
#' average of these values. Precision for a class is defined as the proportion 
#' of true positive predictions out of all predictions made for that class. 
#' 
#' Macro Recall is similar to Macro Precision but focuses on recall. Recall for 
#' a class is the proportion of true positive predictions out of all actual 
#' instances of that class. Macro Recall is the average of the recall values 
#' computed for each class.
#'
#' @note
#' Note that Macro Precision and Macro Recall depend on the assigned labels, 
#' while the ARI measures the similarity between partition up to label 
#' switching.   
#' 
#' If the required packages (`mclust` for ARI, `clusterRepro` for IGP, and
#' `cluster` for ASW) are not installed, the function will display a message
#' asking the user to install the missing package(s).
#' 
#' @return List with the following components:
#' \itemize{
#'    \item \code{metrics} Table of computed evaluation measures for each value
#'                         of number of clusters in the \code{pkbc} object. The
#'                         number of cluster is indicated as column name. 
#'    \item \code{IGP} List of in-group proportions for each value of number of 
#'                     clusters specified.
#' }
#' 
#' @seealso [pkbc()] for the clustering algorithm \cr
#'          \linkS4class{pkbc} for the class object definition.
#'
#' @references
#' Kapp, A.V. and Tibshirani, R. (2007) "Are clusters found in one dataset present 
#' in another dataset?", Biostatistics, 8(1), 9–31, 
#' https://doi.org/10.1093/biostatistics/kxj029
#' 
#' Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. Journal of Computational and Applied 
#' Mathematics, 20, 53–65.
#'
#' @importFrom stats dist
#' @importFrom utils install.packages
#'
#' @examples
#' #We generate three samples of 100 observations from 3-dimensional
#' #Poisson kernel-based densities with rho=0.8 and different mean directions
#' 
#' size<-20
#' groups<-c(rep(1, size), rep(2, size),rep(3,size))
#' rho<-0.8
#' set.seed(081423)
#' data1<-rpkb(size, c(1,0,0),rho,method='rejvmf')
#' data2<-rpkb(size, c(0,1,0),rho,method='rejvmf')
#' data3<-rpkb(size, c(1,0,0),rho,method='rejvmf')
#' data<-rbind(data1$x,data2$x, data3$x)
#'
#' #Perform the clustering algorithm
#' pkbc_res<- pkbc(data, 3)
#' pkbc_validation(pkbc_res)
#' 
#' 
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G2.0, G2.0a, G2.1, G2.1a,G2.2} input true_label
#' @srrstats {UL3.2} true label can be provided as a separate input 
#' 
#' @export
pkbc_validation <- function(object, true_label=NULL){
   
   # nocov start
   if (!requireNamespace("mclust", quietly = TRUE)) {
      install <- readline(prompt = "'mclust' is required for ARI. Would you 
                          like to install it now? (yes/no): ")
      if (tolower(install) == "yes") {
         install.packages("mclust")
      } else {
         message("ARI will not be computed without 'mclust'.")
         return(NULL)
      }
   }
   
   if (!requireNamespace("cluster", quietly = TRUE)) {
      install <- readline(prompt = "'cluster' is required for ASW. Would you 
                          like to install it now? (yes/no): ")
      if (tolower(install) == "yes") {
         install.packages("cluster")
      } else {
         message("ASW will not be computed without 'cluster'.")
         return(NULL)
      }
   }
   
   if (!requireNamespace("clusterRepro", quietly = TRUE)) {
      install <- readline(prompt = "'clusterRepro' is required for IGP. 
                           Would you like to install it now? (yes/no): ")
      if (tolower(install) == "yes") {
         install.packages("clusterRepro")
      } else {
         message("IGP will not be computed without 'clusterRepro'.")
         return(NULL)
      }
   }
   # nocov end
   
   x <- object@input$dat
   x <- x/sqrt(rowSums(x^2))
   
   if(!is.null(true_label)){
      if((is.factor(true_label) | is.vector(true_label)) &!is.list(true_label)){
         if(length(true_label)!=nrow(x)){
            stop("true_label must have the same length of finalMemb.")
         }
      } else {
         stop("true_label must be a factor or a vector")
      }
   }
   
   #test_res <- matrix(nrow = 4)
   if(is.null(true_label)){
      metrics <- matrix(nrow=1)
   } else {
      metrics <- matrix(nrow=4)
   }
   igp_k <- list()
   
   # For each number of cluster considered
   for(k in object@input$nClust){
      
      # Compute the In-Group Proportion
      if(k>1){
         igp_k[[k]] <- clusterRepro::IGP.clusterRepro(as.data.frame(t(x)), 
                              as.data.frame(t(object@res_k[[k]]$params$mu)))$IGP
      }
      
      # Compute the Average Silhouette Width
      sil <- mean(cluster::silhouette(x =object@res_k[[k]]$finalMemb, 
                                      dist =dist(x))[,3])
      
      
      # If true labels are provided
      if(!is.null(true_label)){
         
         # Compute the Adjusted Rand Index
         ari <- mclust::adjustedRandIndex(object@res_k[[k]]$finalMemb, 
                                          true_label)
         
         m <- max(length(unique(true_label)),
                  length(unique(object@res_k[[k]]$finalMemb)))
         conf_mat <- table(factor(true_label,levels=c(1:m)), 
                           factor(object@res_k[[k]]$finalMemb,levels=c(1:m)))
         
         # Compute Macro-Precision and Macro-Recall
         precision <- recall <- numeric(nrow(conf_mat))
         for (i in seq_len(ncol(conf_mat))) {
            TP <- conf_mat[i, i]
            FP <- sum(conf_mat[, i]) - TP
            FN <- sum(conf_mat[i, ]) - TP
            
            precision[i] <- TP / (TP + FP)
            recall[i] <- TP / (TP + FN)
         }
         
         macroPrecision <- mean(precision,na.rm=T)
         macroRecall <- mean(recall,na.rm=T)
         
         metrics <- cbind(metrics, c(sil, ari, macroPrecision, 
                                     macroRecall))
      } else {
         ari <- NA
         metrics <- cbind(metrics, c(sil))
      }
      
   }
   metrics <- metrics[,-1]
   if(is.null(true_label)){
      metrics <- as.data.frame(matrix(metrics,nrow=1),colnames=NULL)
   } else {
      metrics <- as.data.frame(metrics,colnames=NULL)
   }
   
   if(!is.null(true_label)){
      rownames(metrics) <- c("ASW", "ARI","Macro_Precision", 
                             "Macro_Recall")
   } else {
      rownames(metrics) <- c("ASW")
   }
   colnames(metrics) <- object@input$nClust

   results <- list(metrics = metrics, IGP = igp_k)
   
   return(results)
}
