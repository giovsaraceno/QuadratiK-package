#'
#' Poisson kernel-based clustering on the sphere
#'
#' The function \code{pkbc} performs the Poisson kernel-based clustering 
#' algorithm on the sphere based on the Poisson kernel-based densities.
#'
#' @param dat Data matrix or data.frame of data points on the sphere to be 
#'            clustered. The observations in dat are normalized to ensure that
#'            they lie on the d-simensional sphere. Note that d > 1.
#' @param nClust Number of clusters. It can be a single value or a numeric 
#'               vector.
#' @param maxIter The maximum number of iterations before a run is terminated.
#' @param stoppingRule String describing the stopping rule to be used within 
#'                     each run. Currently must be either: 
#'                \code{'max'} (until the change in the log-likelihood is less 
#'                than a given threshold (1e-7)), 
#'                \code{'membership'} (until the membership is unchanged), or 
#'                \code{'loglik'} (based on a maximum number of iterations).
#' @param initMethod String describing the initialization method to be used.
#'                   Currently must be \code{'sampleData'}.
#' @param numInit Number of initializations.
#'
#' @details The function estimates the parameter of a mixture of Poisson
#' kernel-based densities. The obtained estimates are used for assigning final 
#' memberships, identifying the \code{nClust} clusters.
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
#' Golzy, M., Markatou, M. (2020) Poisson Kernel-Based Clustering on the Sphere:
#' Convergence Properties, Identifiability, and a Method of Sampling, Journal of
#' Computational and Graphical Statistics, 29:4, 758-770, 
#' DOI: 10.1080/10618600.2020.1740713.
#'
#' @examples
#' #We generate three samples of 100 observations from 3-dimensional
#' #Poisson kernel-based densities with rho=0.8 and different mean directions
#' size<-100
#' groups<-c(rep(1, size), rep(2, size),rep(3,size))
#' rho<-0.8
#' set.seed(081423)
#' data1<-rpkb(size, c(1,0,0),rho,method="rejvmf")
#' data2<-rpkb(size, c(0,1,0),rho,method="rejvmf")
#' data3<-rpkb(size, c(0,0,1),rho,method="rejvmf")
#' dat<-rbind(data1$x,data2$x, data3$x)
#'
#' #Perform the clustering algorithm with number of clusters k=3.
#' pkbd<- pkbc(dat, 3)
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
                           nClust=NULL,
                           maxIter = 300,
                           stoppingRule = 'loglik',
                           initMethod = 'sampleData',
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
             nClust=NULL,
             maxIter = 300,
             stoppingRule = 'loglik',
             initMethod = 'sampleData',
             numInit = 10){
       # Constant defining threshold by which log likelihood must change 
       # to continue iterations, if applicable.
       LL_TOL <- 1e-7
       
       # validate input
       if(is.null(nClust)){
          stop("Input parameter nClust is required. Provide one specific 
               value or a set of possible values.")
       } else if(is.vector(nClust) & is.numeric(nClust)){
          if(any(nClust < 1)){
             stop('Values in the input parameter nClust must be 
                  greater than 0')
          }
       } else {
          stop("nClust must be a signle value or a numeric vector of possible
                  values")
       }
       if (maxIter < 1) {
          stop('Input parameter maxIter must be greater than 0')
       }
       if ( !(stoppingRule %in% c('max', 'membership', 'loglik')) ) {
          stop(paste('Unrecognized value "', stoppingRule, '" in input 
          parameter stoppingRule.', sep=''))
       }
       if ( !(initMethod %in% c('sampleData')) ) {
          stop(paste('Unrecognized value "', initMethod, '" in input 
          parameter initMethod.', sep=''))
       }
       if (numInit < 1) {
          stop('Input parameter numInit must be greater than 0')
       }
       
       # set options for stopping rule
       checkMembership <- stoppingRule == 'membership'
       checkLoglik <- stoppingRule == 'loglik'
       
       if(is.data.frame(dat)){
          dat <- as.matrix(dat)
       } else if(!is.matrix(dat)){
          stop("dat must be a matrix or a data.frame")
       }
       if(!is.numeric(dat)){
          stop("dat must be a numeric matrix or data.frame")
       }
       
       if(any(is.na(dat))){
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
          normprobMat_best <- matrix(-99,nrow=numData, ncol=numClust)
          if (initMethod == 'sampleData') {
             uniqueData <- unique(dat)
             numUniqueObs <- nrow(uniqueData)
             if (numUniqueObs < numClust) {
                stop(paste('Only', numUniqueObs, 'unique observations.',
                           'When initMethod = "sampleData", must have more
                           than numClust unique observations.'
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
#' @seealso [pkbc()]
#'
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {UL4.4} summary method for pkbc object
#' 
#' @rdname summary.pkbc
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
#' @details The function computes mean, standard deviation, median, 
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
#' @return List with computed descriptive statistics for each variable. 
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
#' @rdname stats_clusters                  
#' @export
setMethod("stats_clusters", "pkbc", function(object, k){
   
   if(!(k %in% object@input$nClust)){
      stop("The provided pkbc object does not contain results for the requested
           number of clusters")
   } else if(!(is.numeric(k) & length(k)==1)){
      stop("k must be an integer")
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
#' @rdname plot.pkbc
#' 
#' @title Plotting method for Poisson kernel-based clustering
#' 
#' @description
#' Plots for a pkbc object.
#'  
#' @param x Object of class \code{pkbc}
#' @param true_label factor or vector of true membership to clusters (if 
#'                   available). It must have the same length of final 
#'                   memberships.
#' @param pca_res Logical. If TRUE the results from PCALocantore when dimension
#'                is greater than 3 are also reported. 
#'
#' @details 
#' - scatterplot: If dimension is equal to 2 or 3, points are displayed on the 
#' circle and sphere, respectively. If dimension if greater than 3, the 
#' spherical Principal Component procedure proposed by Locantore et al., (1999)
#' is applied for dimensionality reduction and the first three principal
#' components are normalized and displayed on the sphere. For d > 3, the
#' complete results from the \code{PcaLocantore} function (package \code{rrcov})
#' are returned if pca_res=TRUE.
#' - elbow plot: the within cluster sum of squares (wcss) is computed using the 
#' Euclidean distance and the cosine similarity. 
#' 
#' @examples
#' \donttest{
#' dat<-matrix(rnorm(300),ncol=3)
#' pkbc_res<- pkbc(dat, 3)
#' stats_clusters(pkbc_res, 3)
#' }
#' 
#' @references
#' Locantore, N., Marron, J.S., Simpson, D.G. et al. (1999) "Robust principal 
#' component analysis for functional data." Test 8, 1–73. 
#' https://doi.org/10.1007/BF02595862
#' 
#' @return One of the following plot:
#' - scatterplot of data points colored by final membership
#' - elbow plot
#' 
#' @srrstats {UL6.1} this function includes a plot method
#' @srrstats {UL6.0,UL6.2} plot method for pkbc object
#' 
#' @export
setMethod("plot", signature(x="pkbc"), 
          function(x, true_label=NULL, pca_res=FALSE) {
repeat {
   # Display plot options
   cat("Select a plot option:\n")
   cat("0: Exit\n")
   cat("1: Scatterplot\n")
   cat("2: Elbow plot\n")
   
   # Read user choice
   choice <- as.integer(readline(prompt="Enter your choice: "))
   
   # Check for exit condition
   if (choice == 0) {
      cat("Exiting plot menu.\n")
      break
   }
   
   # Validate input
   if (choice %in% c(1,2)) {
      if (choice == 1) {
         
   k <- as.integer(readline(prompt="Enter the number of clusters to display: "))

         if(k %in% x@input$nClust){
            
            # Call the scatterplot function
            scatterplotMethod(x, k, true_label,pca_res)
         } else {
            cat("Invalid number of clusters. Please try again.\n")
         }
      } else if(choice == 2){
         
         elbowMethod(x)
         
      }
   } else {
      cat("Invalid choice. Please try again.\n")
   }
}
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
#' @importFrom rgl plot3d
#' @importFrom rgl layout3d
#' @importFrom rgl title3d
#' @importFrom rgl next3d
#' @importFrom rgl rgl.spheres
#' @importFrom rgl open3d
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom rrcov PcaLocantore
#' 
#' @srrstats {G1.4} roxigen2 is used
#' 
#' @keywords internal
#' @noRd
scatterplotMethod <- function(object, k, true_label=NULL, pca_res=FALSE) {
   x <- object@input$dat
   y <- as.factor(object@res_k[[k]]$finalMemb)
   if (ncol(x) == 2) {
      
      df <- data.frame(V1 = x[,1], V2 = x[,2], clusters = as.factor(y))
      with(df, {pl <- ggplot(df, aes(x = V1, y = V2, color = clusters)) +
         geom_point() +
         theme_minimal() +
         labs(color = "Cluster") 
      print(pl)})
      
   } else if (ncol(x) == 3) {
      
      if(!is.null(true_label)){
         
         my_palette <- colorRampPalette(c("orange","blue", "green", "red")) 
         col_pal <- my_palette(k*2)
         layout3d(matrix(1:2, ncol = 2))
         
         next3d()
         plot3d(x[,1], x[,2], x[,3], col = col_pal[y], size = 4)
         title3d("Final membership", font=2,line = 7, cex = 3)
         rgl.spheres(0, col = "transparent", alpha = 0.2)
         
         next3d()
         plot3d(x[,1], x[,2], x[,3], col = col_pal[true_label+k], size = 4)
         title3d("True label", line = 7, cex = 3, font=2)
         rgl.spheres(0 , col = "transparent", alpha = 0.2)
         
      } else {
         
         open3d()
         plot3d(x[,1], x[,2], x[,3], col = y, size = 4)
         rgl.spheres(0, col = "transparent", alpha = 0.2)
         
      }
      
   } else {
      
      pca_result <- PcaLocantore(x)
      pca_data <- data.frame(PC1 = pca_result@scores[,1], 
                             PC2 = pca_result@scores[,2], 
                             PC3 = pca_result@scores[,3])
      pca_data <- pca_data/sqrt(rowSums(pca_data^2))
      pca_data <- data.frame(pca_data, Cluster = y)
      
      
      my_palette <- colorRampPalette(c("orange","blue", "green", "red"))
      col_pal <- my_palette(k*2)
      col_pal <- sample(col_pal, length(col_pal), replace=FALSE)
      
      if(is.null(true_label)){
         open3d()
         plot3d(pca_data$PC1, pca_data$PC2, pca_data$PC3, col = col_pal[y], 
                size = 4, xlab="PC1", ylab="PC2", zlab="PC3")
         title3d("Final membership", line = 7, cex = 3, font=2)
         rgl.spheres(0, col = "transparent", alpha = 0.2)
         
      } else {
         
         layout3d(matrix(1:2, ncol = 2))
         next3d()
         plot3d(pca_data$PC1, pca_data$PC2, pca_data$PC3, col = col_pal[y], 
                size = 4, xlab="PC1", ylab="PC2", zlab="PC3")
         title3d("Final membership", line = 7, cex = 3, font=2)
         rgl.spheres(0, col = "transparent", alpha = 0.2)
         
         next3d()
         plot3d(pca_data$PC1, pca_data$PC2, pca_data$PC3, 
                col = col_pal[true_label+k], size = 4, 
                xlab="PC1", ylab="PC2", zlab="PC3") 
         title3d("True label", line = 7, cex = 3, font=2)
         rgl.spheres(0, col = "transparent", alpha = 0.2)
         
      }
      
      if(pca_res){
         return(pca_result)
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
      labs(title = "Elbow Plot", x = "Number of clusters", 
           y = "Within-cluster sum of squares (WCSS)") +
      theme_minimal()
   
   wcss_values <- data.frame(k = object@input$nClust, wcss = wcss[,2])
   pl_cos <- ggplot(wcss_values, aes(x = k, y = wcss)) +
      geom_line() +
      geom_point() +
      labs(title = "Elbow Plot", x = "Number of clusters", 
           y = "Within-cluster sum of squares (WCSS)") +
      theme_minimal()
   
   fig <- ggarrange(plotlist = list(pl, pl_cos), nrow=1)
   print(fig)

}
#' Cluster spherical observations by mixture of Poisson kernel-based densities
#' 
#' Cluster spherical observations based on mixture of Poisson kernel-based 
#' densities estimated by \code{pkbc}
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
#' @seealso [pkbc()]
#' 
#' @srrstats {G1.a} roxygen2 is used
#' @srrstats {UL3.3} prediction function for pkbc object
#' 
#' @rdname predict.pkbc
#' 
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
#'
#' @param object Object of class \code{pkbc}
#' @param true_label factor or vector of true membership to clusters (if 
#'                   available). It must have the same length of final 
#'                   memberships.
#' @param h Tuning parameter of the k-sample test. (default: 1.5)
#'
#' @details The following evaluation measures are computed: 
#'  In-Group Proportion. If true label are provided, ARI, Average
#'  Silhouette Width, Macro-Precision and Macro-Recall are computed.
#'
#' @return List with the following components:
#' \itemize{
#'    \item \code{metrics} Table of computed evaluation measures.
#'    \item \code{IGP} List of in-group proportions for each value of number of 
#'                     clusters specified.
#' }
#'
#' @references
#' Kapp, A.V., Tibshirani, R. (2007) "Are clusters found in one dataset present 
#' in another dataset?", Biostatistics, 8(1), 9–31, 
#' https://doi.org/10.1093/biostatistics/kxj029
#' 
#' Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis. Journal of Computational and Applied 
#' Mathematics, 20, 53–65.
#'
#' @importFrom mclust adjustedRandIndex
#' @importFrom clusterRepro IGP.clusterRepro
#' @importFrom cluster silhouette
#' @importFrom stats dist
#'
#' @examples
#' #We generate three samples of 100 observations from 3-dimensional
#' #Poisson kernel-based densities with rho=0.8 and different mean directions
#' \donttest{
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
#' pkbc_res<- pkbc(data, 2:4)
#' pkbc_validation(pkbc_res)
#' }
#' 
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G2.0, G2.0a, G2.1, G2.1a,G2.2} input true_label
#' @srrstats {UL3.2} true label can be provided as a separate input 
#' 
#' @export
pkbc_validation <- function(object, true_label=NULL, h=1.5){
   
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
      metrics <- matrix(nrow=2)
   } else {
      metrics <- matrix(nrow=5)
   }
   igp_k <- list()
   
   # For each number of cluster considered
   for(k in object@input$nClust){
      
      # Compute the In-Group Proportion
      if(k>1){
         igp_k[[k]] <- IGP.clusterRepro(as.data.frame(t(x)), 
                              as.data.frame(t(object@res_k[[k]]$params$mu)))$IGP
      }
      # Compute the k-sample test
      #y_k <- as.numeric(object@res_k[[k]]$finalMemb)
      #k_test <- kb.test(x=x, y=y_k,h=h)
      
      # Compute the Average Silhouette Width
      sil <- mean(silhouette(x =object@res_k[[k]]$finalMemb, dist =dist(x))[,3])
      
     # test_res <- cbind(test_res, c(k, k_test@Un[1], k_test@CV_Un[1], 
                               #   k_test@H0_Un[1]))
      
      # If true labels are provided
      if(!is.null(true_label)){
         
         # Compute the Adjusted Rand Index
         ari <- adjustedRandIndex(object@res_k[[k]]$finalMemb, true_label)
         
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
         
         metrics <- cbind(metrics, c(k, sil, ari, macroPrecision, 
                                     macroRecall))
      } else {
         ari <- NA
         metrics <- cbind(metrics, c(k, sil))
      }
      
   }
   metrics <- metrics[,-1]
   #test_res <- test_res[,-1]
   metrics <- as.data.frame(metrics,colnames=NULL)
   #test_res <- as.data.frame(test_res, colnames=NULL)
   #colnames(metrics) <- paste0("k=",object@input$nClust)
   #rownames(test_res) <- c("k","Test statistic", "Critical value", 
                        #  "Reject_H0")
   if(!is.null(true_label)){
      rownames(metrics) <- c("k","ASW", "ARI","Macro_Precision", 
                             "Macro_Recall")
   } else {
      rownames(metrics) <- c("k","ASW")
   }
   

   results <- list(metrics = metrics, IGP = igp_k)
   #results <- list(metrics = metrics, IGP = igp_k, tests = test_res)
   
   return(results)
}
.onLoad <- function(libname, pkgname) {
   if (Sys.getenv("RGL_USE_NULL") == "") {
      Sys.setenv(RGL_USE_NULL = "TRUE")
   }
}