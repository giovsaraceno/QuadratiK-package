#' NA_standards
#'
#'
#'
#' @srrstatsNA {G1.6} there are no alternative implementations
#' @srrstatsNA {G2.4, G2.4a, G2.4c, G2.4d, G2.4e} other conversions are not 
#'             used.
#' @srrstatsNA {G2.5} there are no input of type factor which should be ordered
#' @srrstatsNA {G2.6} we check that each univariate input has the expected type
#' @srrstatsNA {G2.9} there are no conversions needed
#' @srrstatsNA {G2.10} there is no extraction
#' @srrstatsNA {G2.11,G2.12} We consider only data.frame with numeric variables
#' @srrstatsNA {G2.14b,G2.14c} we just report errors in case of missing data
#' @srrstatsNA {G3.1,G3.1a} methods are not based on covariance computation
#' @srrstatsNA {G4.0} no outputs are written in local files
#' @srrstatsNA {G5.0} the functions are tested with simple simulated data sets
#' @srrstatsNA {G5.3} missing and undefined values are not considered
#' @srrstatsNA {G5.4b,G5.4c} no implementation of existing methods 
#' @srrstatsNA {G5.6, G5.6a, G5.6b} There are no parameter recovery methods
#' @srrstatsNA {G5.8d} there are no such data
#' @srrstatsNA {G5.9a,G5.9b} correctness tested on simple cases
#' @srrstatsNA {G5.10, G5.11, G5.11a, G5.12} there are no extended tests
#' 
#' @srrstatsNA {UL1.2,UL1.3,UL1.3a} rownames and colnames are not used
#' @srrstatsNA {UL1.4a,UL1.4b} there are no differences with respect to scales
#' @srrstatsNA {UL2.0} Violations of assumptions are not considered.
#' @srrstatsNA {UL2.2} do not accept NA
#' @srrstatsNA {UL3.0,UL7.2,UL7.3} Labels are not ordered with respect to sizes
#' @srrstatsNA {UL3.1} not a dimensionality reduction software 
#' @srrstatsNA {UL7.5,UL7.5a} there are no batch processes
#' 
#' @noRd
NULL
