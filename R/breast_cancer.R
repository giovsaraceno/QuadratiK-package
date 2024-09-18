#' Breast Cancer Wisconsin (Diagnostic)
#'
#' The \code{breast_cancer} Wisconsin data has 569 rows and 31 columns. The 
#' first 30 variables report the features that are computed from a digitized 
#' image of a fine needle aspirate (FNA) of a breast mass. They describe   
#' characteristics of the cell nuclei present in the image. The last column 
#' indicates the class labels (Benign = 0 or Malignant = 1).
#'
#' @format A data frame of 569 observations and 31 variables. 
#' 
#' @examples
#' data(breast_cancer)
#' summary(breast_cancer)
#'
#' @source
#' Wolberg, W., Mangasarian, O., Street, N., & Street, W. (1993). 
#' Breast Cancer Wisconsin (Diagnostic). 
#' UCI Machine Learning Repository. 
#' https://doi.org/10.24432/C5DW2B.
#'
#' @references
#' Street, W. N., Wolberg, W. H., & Mangasarian, O. L. (1993, July). Nuclear 
#' feature extraction for breast tumor diagnosis. In Biomedical image processing
#' and biomedical visualization (Vol. 1905, pp. 861-870). SPIE.
#' 
#' @srrstats {G1.0} Reference section reports the related literature
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G5.1} the data set is made available
#' 
"breast_cancer"
