#' Wine data set
#'
#' The \code{wine} data frame has 178 rows and 14 columns. The first 13 
#' variables report 13 constituents found in each of the three types of wines. 
#' The last column indicates the class labels (1,2 or 3).
#' 
#' @format A data frame containing the following columns:
#' \itemize{
#'    \item \code{Alcohol} 
#'    \item \code{Malic acid} 
#'    \item \code{Ash} 
#'    \item \code{Alcalinity of ash}
#'    \item \code{Magnesium} 
#'    \item \code{Total phenols}
#'    \item \code{Flavanoids}
#'    \item \code{Nonflavanoid phenols}
#'    \item \code{Proanthocyanins}
#'    \item \code{Color intensity}
#'    \item \code{Hue}
#'    \item \code{OD280/OD315 of diluted wines}
#'    \item \code{Proline}
#'    \item \code{y}: class membership
#' }
#'
#' @details
#' These data are the results of a chemical analysis of wines grown in the same 
#' region in Italy but derived from three different cultivars. The analysis 
#' determined the quantities of 13 constituents found in each of the three types
#' of wines.
#'
#' @examples
#' data(wine)
#' summary(wine)
#'
#' @source
#' Aeberhard, S. and Forina, M. (1991). Wine.
#' UCI Machine Learning Repository. 
#' https://doi.org/10.24432/C5PC7J.
#'
#' @references
#' Aeberhard, S., Coomans, D. and De Vel, O. (1994). Comparative analysis of 
#' statistical pattern recognition methods in high dimensional settings. 
#' Pattern Recognition, 27(8), 1065-1077.
#' 
#' @srrstats {G1.0} Reference section reports the related literature
#' @srrstats {G1.4} roxigen2 is used
#' @srrstats {G5.1} the data set is made available
#' 
"wine"
