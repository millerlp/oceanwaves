# Filename: AlmostZero.R
# 
###############################################################################


#' Define a function to test whether elements in vector x are effectively zero
#' within the square root of machine tolerance for a floating point number.
#' Returns TRUE for all vector elements that are within rounding error of zero.
#' @param x Numeric vector of values to compared against 0
#' @param y Value to be compared against. Default is 0. 
#' @param tolerance The maximum difference between x and y necessary to call a 
#' value non-zero
#' @return TRUE or FALSE for each element of the vector x. 

AlmostZero <- function(x, y = 0, tolerance=sqrt(.Machine$double.eps)) {
	diff <- abs(x - y)
	mag <- pmax(abs(x), abs(y))
	ifelse(mag > tolerance, diff/mag <= tolerance, diff <= tolerance)
}
