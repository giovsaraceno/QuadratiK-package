#'
#' Poisson kernel-based clustering on the sphere
#'
#' The function \code{pkbc} performs the Poisson kernel-based clustering algorithm 
#' on the sphere based on the Poisson kernel-based densities.
#'
#' @param dat Data matrix of observations to be clustered.
#' @param nClust Number of clusters. It can be a single value or a numeric vector.
#' @param maxIter The maximum number of iterations before a run is terminated.
#' @param stoppingRule String describing the stopping rule to be
#'   used within each run. Currently must be either \code{'max'} (until the change
#'   in the log-likelihood is less than a given threshold $(1e-7)$), \code{'membership'}
#'    (until the membership is unchanged), or \code{'loglik'} (based on a maximum
#'     number of iterations).
#' @param initMethod String describing the initialization method to be used.
#'   Currently must be \code{'sampleData'}.
#' @param numInit Number of initializations.
#'
#' @details The function estimates the parameter of a mixture of Poisson kernel-based
#'  densities. The obtained estimates are used for assigning final memberships, 
#'  identifying the \code{nClust} clusters.
#'
#' @return An S4 object of class \code{pkbc} containing the results of the clustering
#'  procedure based on Poisson kernel-based distributions. The object contains 
#'  the following slots:
#'
#'   \code{res_k}: List of results of the Poisson kernel-based clustering algorithm for each value of number of clusters specified in \code{nClust}. Each object in the list contains:
#'   \itemize{
#'      \item \code{postProbs} Posterior probabilities of each observation for the indicated clusters.
#'      \item \code{LogLik} Maximum value of log-likelihood function
#'      \item \code{wcss} Values of within-cluster sum of squares computed with Euclidean distance and cosine similarity, respectively.
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
#' size <- 100
#' groups<-c(rep(1, size), rep(2, size),rep(3,size))
#' rho=0.8
#' set.seed(081423)
#' data1<-rpkb(size, c(1,0,0),rho,method='rejvmf')
#' data2<-rpkb(size, c(0,1,0),rho,method='rejvmf')
#' data3<-rpkb(size, c(-1,0,0),rho,method='rejvmf')
#' dat<-rbind(data1$x,data2$x, data3$x)
#'
#' #The following code can be used for plotting
#' \donttest{
#' library("rgl")
#' plot3d(dat,col = groups, size=4,xlab = "", ylab = "",zlab = "")
#' legend3d("right", legend = c("Class 1", "Class 2", "Class 3"), col = 1:3, pch = 20)
#' }
#' #Perform the clustering algorithm with number od clusters k=3.
#' pkbd<- pkbc(dat, 3)
#'
#' @export
setGeneric("pkbc",function(dat,
                           nClust,
                           maxIter = 300,
                           stoppingRule = 'loglik',
                           initMethod = 'sampleData',
                           numInit = 10){
   
   standardGeneric("pkbc")
})
#' @rdname pkbc
#' @export
setMethod("pkbc", signature(dat = "ANY"),
          function(dat,
                   nClust,
                   maxIter = 300,
                   stoppingRule = 'loglik',
                   initMethod = 'sampleData',
                   numInit = 10){
             # Constant defining threshold by which log likelihood must change to continue
             # iterations, if applicable.
             LL_TOL <- 1e-7
             
             # validate input
             if(is.null(nClust)){
                stop("Input parameter nClust is require. Provide one specific value or a set of possible values.")
             } else {
                if(any(nClust < 1)){
                   stop('Values in the input parameter nClust must be greater than 0')
                }
             }#else if(length(nClust)==1){
             #if (nClust < 2 ) {
             #  stop('Input parameter nClust must be an integer greater than 1')
             # }
             #}
             if (maxIter < 1) {
                stop('Input parameter maxIter must be greater than 0')
             }
             if ( !(stoppingRule %in% c('max', 'membership', 'loglik')) ) {
                stop(paste('Unrecognized value "', stoppingRule, '" in input parameter
      stoppingRule.', sep=''))
             }
             if ( !(initMethod %in% c('sampleData')) ) {
                stop(paste('Unrecognized value "', initMethod, '" in input parameter
      initMethod.', sep=''))
             }
             if (numInit < 1) {
                stop('Input parameter numInit must be greater than 0')
             }
             
             # set options for stopping rule
             checkMembership <- stoppingRule == 'membership'
             checkLoglik <- stoppingRule == 'loglik'
             
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
                                 'When initMethod = "sampleData", must have more than numClust unique observations.'
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
                      mu_current <- matrix(t(uniqueData[sample(numUniqueObs, size=numClust, replace=FALSE) ,]), ncol=numVar, byrow=TRUE)
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
                      log_probMat_denom <- log(1 + rho_mat_current^2 - 2*rho_mat_current*v_mat)
                      # eq (11) in ICMLD16 paper
                      log_probMat <- log(1 - (rho_mat_current^2)) - log_w_d - (numVar / 2) * log_probMat_denom
                      ######### E step done#############################################
                      ########## M step now#############################################
                      # Denominator of eq (18) in ICMLD16 paper
                      probSum <- matrix(
                         exp(log_probMat) %*% alpha_current,
                         nrow = numData,
                         ncol = numClust
                      )
                      # eq (18) in ICMLD16 paper
                      log_normProbMat_current <- log(alpha_mat_current) + log_probMat - log(probSum)
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
                            (-2*numData*x*alpha_current_h)/(1 - x^2) + numVar*mu_denom_h -numVar*x*sum_h_weightMat
                         }
                         rho_current[h] <- (uniroot(
                            f = root_func,
                            interval = c(0,1),
                            f.lower = root_func(0),
                            f.upper = root_func(1),
                            tol = 0.001
                         ))$root
                      } # for h
                      
                      # Terminate iterations if maxIter is reached. Redundant with while
                      # condition, but this ensures that the stored numIter below is accurate.
                      if (currentIter >= maxIter) {
                         break
                      }
                      
                      # Terminate iterations if membership is unchanged, if applicable
                      if (checkMembership) {
                         # calculate and store group membership
                         membPrevious <- membCurrent
                         membCurrent <- apply(exp(log_normProbMat_current), 1, which.max)
                         if (all(membPrevious == membCurrent)) {
                            break
                         }
                      }
                      # Terminate iterations if log-likelihood is unchanged, if applicable.
                      if (checkLoglik) {
                         loglikPrevious <- loglikCurrent
                         # Calculate log likelihood for the current iteration
                         loglikCurrent <- sum(log( exp(log_probMat) %*% alpha_current ))
                         if (abs(loglikPrevious - loglikCurrent) < LL_TOL) {
                            break
                         }
                      }
                      # Update counter to NEXT iteration number.
                      currentIter <- currentIter + 1
                   } # for currentIter
                   # Compare the current log likelihood to the best stored log likelihood;
                   # if it exceeds the stored, then store it and all associated estimates.
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
                wcss <- sapply(1:numClust, function(cl) {
                   idx <- which(memb_best==cl)
                   
                   if(length(idx)>0){
                      
                      dat_memb <- matrix(t(dat[idx,]),ncol=numVar, byrow=TRUE)
                      mu_memb <- matrix(rep(as.numeric(mu_best[cl,]), times = nrow(dat_memb)), ncol=numVar, byrow=TRUE)
                      return(sum(apply(dat_memb - mu_memb, 1, norm, "2")**2))
                      
                   } else {
                      return(0)
                   }
                })
                wcss <- sum(wcss)
                
                # Calculate the within-cluster sum of squares using cosine similarity
                wcss_cos <- sapply(1:numClust, function(cl) {
                   idx <- which(memb_best==cl)
                   if(length(idx)>0){
                      dat_memb <- matrix(t(dat[idx,]),ncol=numVar, byrow=TRUE)
                      mu_memb <- matrix(rep(mu_best[cl,], nrow(dat_memb)), ncol=numVar,byrow=TRUE)
                      
                      return(sum(dat_memb * mu_memb))
                   } else {
                      return(0)
                   }
                })
                
                wcss_cos <- sum(wcss_cos)
                
                res_k[[numClust]] <- list(postProbs = normprobMat_best,
                                          LogLik = max(logLikVec),
                                          wcss = c(wcss,wcss_cos),
                                          params = list(mu = mu_best, rho = rho_best, alpha = alpha_best),
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
#'
#' Validation of Poisson kernel-based clustering results
#'
#' Method for objects of class \code{pkbc} which displays the elbow plots for the given values of number of clusters and compute evaluation measures for clustering results.
#'
#' @param object Object of class \code{pkbc}
#' @param true_label Vector of true membership to clusters (if available)
#' @param elbow.plot Logical, if TRUE the function returns the elbow plots computed with the Euclidean distance and cosine similarity (default: TRUE).
#' @param h Tuning parameter of the k-sample test. (default: 1.5)
#'
#' @details The function extracts the within-cluster sum of squares and displays the obtained elbow plots. The following evaluation measures are computed: k-sample test, In-Group Proportion. If true label are provided, ARI, Average Silhouette Width, Macro-Precision and Macro-Recall are computed.
#'
#' @return List with the following components:
#' \itemize{
#'    \item \code{metrics} Table of computed evaluation measures.
#'    \item \code{IGP} List of in-group proportions for each value of number of clusters specified.
#'    \item \code{wcss} Table of collected values of within-cluster sum of squares.
#' }
#'
#' @references
#' Kapp, A.V., Tibshirani, R. (2007) "Are clusters found in one dataset present in another dataset?", Biostatistics, 8(1), 9–31, 
#' https://doi.org/10.1093/biostatistics/kxj029
#' 
#' Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. Journal of Computational and Applied Mathematics, 20, 53–65.
#'
#' @importFrom mclust adjustedRandIndex
#' @import clusterRepro
#' @import ggplot2
#' @importFrom cluster silhouette
#'
#' @examples
#' #We generate three sample of 100 observations from 3-dimensional
#' #Poisson kernel-based densities with rho=0.8 and different mean directions
#' \donttest{
#' size <- 20
#' groups<-c(rep(1, size), rep(2, size),rep(3,size))
#' rho=0.8
#' set.seed(081423)
#' data1<-rpkb(size, c(1,0,0),rho,method='rejvmf')
#' data2<-rpkb(size, c(0,1,0),rho,method='rejvmf')
#' data3<-rpkb(size, c(-1,0,0),rho,method='rejvmf')
#' data<-rbind(data1$x,data2$x, data3$x)
#'
#' #Perform the clustering algorithm
#' pkbc_res<- pkbc(data, 2:4)
#' validation(pkbc_res)
#' }
#' 
#' @export
validation <- function(object, true_label=NULL, elbow.plot=TRUE, h=1.5){
   
   x <- object@input$dat
   x <- x/sqrt(rowSums(x^2))
   
   if(is.null(true_label)){
      metrics <- matrix(nrow=5)
   } else {
      metrics <- matrix(nrow=8)
   }
   igp_k <- list()
   wcss <- matrix(ncol=2)
   
   # For each number of cluster considered
   for(k in object@input$nClust){
      
      # Compute the In-Group Proportion
      if(k>1){
         igp_k[[k]] <- IGP.clusterRepro(as.data.frame(t(x)), as.data.frame(t(object@res_k[[k]]$params$mu)))$IGP
      }
      
      # Compute the k-sample test
      k_test <- kb.test(x=x, y=object@res_k[[k]]$finalMemb,h=h,K_threshold=k)
      
      # Compute the Average Silhouette Width
      sil <- mean(silhouette(x = object@res_k[[k]]$finalMemb, dist = dist(x))[,3])
      
      # If true labels are provided
      if(!is.null(true_label)){
         
         # Compute the Adjusted Rand Index
         ari <- adjustedRandIndex(object@res_k[[k]]$finalMemb, true_label)

         m <- max(length(unique(true_label)),length(unique(object@res_k[[k]]$finalMemb)))
         conf_mat <- table(factor(true_label,levels=c(1:m)), factor(object@res_k[[k]]$finalMemb,levels=c(1:m)))
         
         # Compute Macro-Precision and Macro-Recall
         precision <- recall <- numeric(nrow(conf_mat))
         for (i in 1:ncol(conf_mat)) {
            TP <- conf_mat[i, i]
            FP <- sum(conf_mat[, i]) - TP
            FN <- sum(conf_mat[i, ]) - TP
            
            precision[i] <- TP / (TP + FP)
            recall[i] <- TP / (TP + FN)
         }
         
         macroPrecision <- mean(precision,na.rm=T)
         macroRecall <- mean(recall,na.rm=T)
         
         metrics <- cbind(metrics, c(k, k_test@Dn[1], k_test@CV[1], k_test@H0, sil, ari, macroPrecision, macroRecall))
      } else {
         ari <- NA
         metrics <- cbind(metrics, c(k, k_test@Dn[1], k_test@CV[1], k_test@H0, sil))
      }
      
      wcss <- rbind(wcss, object@res_k[[k]]$wcss)
      
   }
   metrics <- metrics[,-1]
   metrics <- as.data.frame(metrics,colnames=NULL)
   colnames(metrics) <- paste0("k=",object@input$nClust)
   if(!is.null(true_label)){
      rownames(metrics) <- c("k","Test statistic", "Critical value", "Reject_H0", "ASW", "ARI","Macro_Precision", "Macro_Recall")
   } else {
      rownames(metrics) <- c("k","Test statistic", "Critical value", "Reject_H0", "ASW")
   }
   
   wcss <- wcss[-1,]
   colnames(wcss) <- c("Euclidean","cosine_similarity")
   rownames(wcss) <- paste0("k=",object@input$nClust)
   
   if(elbow.plot){
      
      wcss_values <- data.frame(k = object@input$nClust, wcss = wcss[,1])
      pl <- ggplot(wcss_values, aes(x = k, y = wcss)) +
         geom_line() +
         geom_point() +
         labs(title = "Elbow Plot", x = "Number of clusters", y = "Within-cluster sum of squares (WCSS)") +
         theme_minimal()
      
      wcss_values <- data.frame(k = object@input$nClust, wcss = wcss[,2])
      pl_cos <- ggplot(wcss_values, aes(x = k, y = wcss)) +
         geom_line() +
         geom_point() +
         labs(title = "Elbow Plot", x = "Number of clusters", y = "Within-cluster sum of squares (WCSS)") +
         theme_minimal()
      
      print(pl)
      print(pl_cos)
   }
   
   results <- list(metrics = metrics, IGP = igp_k, wcss = wcss)
   
   return(results)
}
#'
#' Descriptive statistics for Poisson kernel-based clustering results.
#'
#' Method for objects of class \code{pkbc} which computes and displays some descriptive for each variable with respect to the detected groups. It also generates a plot of data points colored by the final membership.
#'
#' @param object Object of class \code{pkbc}.
#' @param k Number of clusters to be used.
#' @param true_label Vector of true memberships to clusters (default: NULL).
#'
#' @details The function computes mean, standard deviation, median, inter-quantile range, minimum and maximum for each variable in the data set given the final membership assigned by the clustering algorithm. If dimension is equal to 2 or 3, points are displayed on the circle and sphere, respectively. If dimension if greater than 3, the spherical Principal Component procedure proposed by Locantore et al., (1999) is applied for dimensionality reduction and the first three principal components are normalized and displayed on the sphere.
#' 
#' @examples
#' #We generate three samples of 100 observations from 3-dimensional
#' #Poisson kernel-based densities with rho=0.8 and different mean directions
#' size <- 100
#' groups<-c(rep(1, size), rep(2, size),rep(3,size))
#' rho=0.8
#' set.seed(081423)
#' data1<-rpkb(size, c(1,0,0),rho,method='rejvmf')
#' data2<-rpkb(size, c(0,1,0),rho,method='rejvmf')
#' data3<-rpkb(size, c(-1,0,0),rho,method='rejvmf')
#' data<-rbind(data1$x,data2$x, data3$x)
#'
#' #Perform the clustering algorithm
#' pkbc_res<- pkbc(data, 3)
#' summary_stat(pkbc_res, 3)
#' 
#' @references
#' Locantore, N., Marron, J.S., Simpson, D.G. et al. (1999) "Robust principal component
#' analysis for functional data." Test 8, 1–73. https://doi.org/10.1007/BF02595862
#'
#' @return List with computed descriptive statistics for each variable. For d > 3, the complete results from the \code{PcaLocantore} function (package \code{rrcov}) are also returned. 
#'
#' @import rgl
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom rrcov PcaLocantore
#'
#' @export
summary_stat <- function(object, k, true_label=NULL){
   
   if(!(k %in% object@input$nClust)){
      stop("The provided pkbc object does not contain results for the requested number of clusters")
   }
   
   x <- object@input$dat
   y <- as.factor(object@res_k[[k]]$finalMemb)
   
   metrics <- list()
   Results <- list()
   
   for(i in 1:ncol(x)){
      res <- rbind(as.numeric(by(x[,i],y,mean)),
                   as.numeric(by(x[,i],y,sd)),
                   as.numeric(by(x[,i],y,median)),
                   as.numeric(by(x[,i],y,IQR)),
                   as.numeric(by(x[,i],y,min)),
                   as.numeric(by(x[,i],y,max))
      )
      res <- cbind(res, c(mean(x[,i]), sd(x[,i]), median(x[,i]), IQR(x[,i]), min(x[,i]), max(x[,i])))
      res <- as.data.frame(res)
      rownames(res) <- c("mean","sd", "median", "IQR", "min", "max")
      colnames(res) <- c(paste("Group ",seq(1,k),sep=""), "Overall")
      metrics[[length(metrics)+1]] <- res
   }
   
   if (ncol(x) == 2) {
      
      df <- data.frame(V1 = x[,1], V2 = x[,2], clusters = as.factor(y))
      pl <- ggplot(df, aes(x = df$V1, y = df$V2, color = as.factor(df$clusters))) +
         geom_point() +
         theme_minimal() +
         labs(color = "Cluster")
      print(pl)
      
   } else if (ncol(x) == 3) {
      
      if(!is.null(true_label)){
         
         my_palette <- colorRampPalette(c("orange","blue", "green", "red"))  # Gradient from blue to green to red
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
      pca_data <- data.frame(PC1 = pca_result@scores[,1], PC2 = pca_result@scores[,2], PC3 = pca_result@scores[,3])
      #pca_result <- prcomp(x)  # Replace with spherical PCA if available
      #pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], PC3 = pca_result$x[,3])
      pca_data <- pca_data/sqrt(rowSums(pca_data^2))
      pca_data <- data.frame(pca_data, Cluster = y)
      
      Results$pca = pca_result
      
      my_palette <- colorRampPalette(c("orange","blue", "green", "red"))
      col_pal <- my_palette(k*2)
      col_pal <- sample(col_pal, length(col_pal), replace=FALSE)
      
      if(is.null(true_label)){
         open3d()
         plot3d(pca_data$PC1, pca_data$PC2, pca_data$PC3, col = col_pal[y], size = 4, xlab="PC1", ylab="PC2", zlab="PC3")
         title3d("Final membership", line = 7, cex = 3, font=2)
         rgl.spheres(0, col = "transparent", alpha = 0.2)
         
      } else {
         
         layout3d(matrix(1:2, ncol = 2))
         next3d()
         plot3d(pca_data$PC1, pca_data$PC2, pca_data$PC3, col = col_pal[y], size = 4, xlab="PC1", ylab="PC2", zlab="PC3")
         title3d("Final membership", line = 7, cex = 3, font=2)
         rgl.spheres(0, col = "transparent", alpha = 0.2)
         
         next3d()
         plot3d(pca_data$PC1, pca_data$PC2, pca_data$PC3, col = col_pal[true_label+k], size = 4, xlab="PC1", ylab="PC2", zlab="PC3") 
         title3d("True label", line = 7, cex = 3, font=2)
         rgl.spheres(0, col = "transparent", alpha = 0.2)
         
      }
      
   }
   
   Results$metrics = metrics
   return(Results)
}
