#' Wireless Indoor Localization
#'
#' The \code{wireless} data frame has 2000 rows and 8 columns. The first 7 variables
#' report the measurements of the Wi-Fi signal strength received from 7 Wi-Fi routers in an
#' office location in Pittsburgh (USA). The last column indicates the class labels.
#'
#' @format A data frame containing the following columns:
#' \itemize{
#'    \item \code{V1} Signal strength from router 1.
#'    \item \code{V2} Signal strength from router 2.
#'    \item \code{V3} Signal strength from router 3.
#'    \item \code{V4} Signal strength from router 4.
#'    \item \code{V5} Signal strength from router 5.
#'    \item \code{V6} Signal strength from router 6.
#'    \item \code{V7} Signal strength from router 7.
#'    \item \code{V8} Group memberships, from 1 to 4.
#' }
#'
#' @details
#' The Wi-Fi signal strength is measured in dBm, decibel milliwatts, which is expressed
#' as a negative value ranging from -100 to 0.
#' The labels correspond to 4 different rooms.
#' In total, we have 4 groups with 500 observations each.
#'
#' @examples
#' data(wireless)
#' summary(wireless)
#'
#' @source
#' Bhatt,Rajen. (2017). Wireless Indoor Localization. UCI Machine Learning Repository. 
#' 
#' https://doi.org/10.24432/C51880.
#'
#' @references
#' Rohra, J.G., Perumal, B., Narayanan, S.J., Thakur, P., Bhatt, R.B. (2017). "User Localization in an Indoor Environment Using Fuzzy Hybrid of Particle Swarm Optimization & Gravitational Search Algorithm with Neural Networks". In: Deep, K., et al. Proceedings of Sixth International Conference on Soft Computing for Problem Solving. Advances in Intelligent Systems and Computing, vol 546. Springer, Singapore. https://doi.org/10.1007/978-981-10-3322-3_27
#'
"wireless"
